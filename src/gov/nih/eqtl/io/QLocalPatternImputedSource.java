/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Bounded-memory nearest local-genotype-pattern imputation.
 *
 * <p>This is deliberately described as a proxy, not phasing or reference-panel imputation.
 * For each missing dosage it finds called donors whose observed flanking dosages are closest
 * to the sample's flanking dosages, and averages tied closest donor dosages. If no comparable
 * donor exists it falls back to the row mean.</p>
 */
public final class QLocalPatternImputedSource implements QMatrixRowSource {
    private static final Pattern LOCUS = Pattern.compile("^([^:]+):(\\d+)(?::.*)?$");
    private static final double TIE_TOLERANCE = 1e-12;

    private record Locus(String chromosome, long position) { }
    private record Row(long index, String id, double[] values, Locus locus) { }

    private final QMatrixRowSource source;
    private final QMissingnessScan scan;
    private final int flanks;
    private final Metadata metadata;

    public QLocalPatternImputedSource(QMatrixRowSource source, QMissingnessScan scan, int flanks) {
        if (source == null || scan == null)
            throw new IllegalArgumentException("source and scan are required");
        if (flanks < 1)
            throw new IllegalArgumentException("flanks must be at least 1");
        if (source.metadata().rowCount() != scan.rowCount())
            throw new IllegalArgumentException("Source and missingness scan row counts differ");
        this.source = source;
        this.scan = scan;
        this.flanks = flanks;
        String tag = "local-pattern-imputation-v1;flanks=" + flanks;
        if (source.metadata().cacheSignatureTag() != null)
            tag = source.metadata().cacheSignatureTag() + ";" + tag;
        metadata = new Metadata(source.metadata().path(), source.metadata().rowCount(),
            source.metadata().columnCount(), source.metadata().sampleIds(), tag);
    }

    @Override
    public Metadata metadata() {
        return metadata;
    }

    @Override
    public BlockReader open(int[] columnOrder) throws IOException {
        return new Reader(source.open(columnOrder));
    }

    private final class Reader implements BlockReader {
        private final QMatrixRowSource.BlockReader inner;
        private final Deque<Row> left = new ArrayDeque<>();
        private final Deque<Row> future = new ArrayDeque<>();
        private final Set<String> completedChromosomes = new HashSet<>();
        private QMatrixRowSource.Block input;
        private int inputRow;
        private String currentChromosome;
        private long previousPosition = -1;
        private long outputOffset;
        private boolean warnedOrderAssumption;
        private boolean warnedMeanFallback;
        private boolean closed;

        Reader(QMatrixRowSource.BlockReader inner) {
            this.inner = inner;
        }

        @Override
        public Block readBlock(int maximumRows) throws IOException {
            if (maximumRows <= 0)
                throw new IllegalArgumentException("maximumRows must be positive");
            if (closed)
                throw new IOException("Matrix reader is closed");
            List<String> ids = new ArrayList<>(maximumRows);
            List<double[]> rows = new ArrayList<>(maximumRows);
            long offset = outputOffset;
            while (rows.size() < maximumRows) {
                fillFuture(flanks + 1);
                Row center = future.pollFirst();
                if (center == null)
                    break;
                List<Row> neighbors = neighbors(center);
                double[] imputed = center.values().clone();
                impute(center, neighbors, imputed);
                ids.add(center.id());
                rows.add(imputed);
                outputOffset++;
                left.addLast(center);
                while (left.size() > flanks)
                    left.removeFirst();
            }
            return rows.isEmpty() ? null
                : new Block(offset, ids.toArray(String[]::new), rows.toArray(double[][]::new));
        }

        private void fillFuture(int requested) throws IOException {
            while (future.size() < requested) {
                Row row = nextRawRow();
                if (row == null)
                   return;
                validateOrder(row);
                future.addLast(row);
            }
        }

        private Row nextRawRow() throws IOException {
            while (input == null || inputRow >= input.rowCount()) {
                input = inner.readBlock(1024);
                inputRow = 0;
                if (input == null)
                    return null;
            }
            long index = input.rowOffset() + inputRow;
            String id = input.rowIds()[inputRow];
            double[] values = input.values()[inputRow].clone();
            inputRow++;
            return new Row(index, id, values, parseLocus(id));
        }

        private void validateOrder(Row row) throws IOException {
            if (row.locus() == null) {
                if (!warnedOrderAssumption) {
                    warnedOrderAssumption = true;
                    System.err.println("WARNING: Predictor row IDs do not all contain chromosome:position; "
                        + "local-pattern imputation assumes the input rows are in genomic order and may cross "
                        + "an unmarked chromosome boundary.");
                }
                return;
            }
            if (!row.locus().chromosome().equals(currentChromosome)) {
                if (currentChromosome != null)
                    completedChromosomes.add(currentChromosome);
                if (completedChromosomes.contains(row.locus().chromosome()))
                    throw new IOException("Chromosome '" + row.locus().chromosome()
                        + "' is not contiguous at predictor row '" + row.id() + "'");
                currentChromosome = row.locus().chromosome();
                previousPosition = -1;
            }
            if (row.locus().position() < previousPosition)
                throw new IOException("Predictor positions are not sorted at row '" + row.id() + "'");
            previousPosition = row.locus().position();
        }

        private List<Row> neighbors(Row center) {
            List<Row> result = new ArrayList<>(2 * flanks);
            for (Row row : left)
                if (sameRegion(center, row))
                    result.add(row);
            int added = 0;
            for (Row row : future) {
                if (added++ >= flanks)
                    break;
                if (sameRegion(center, row))
                    result.add(row);
            }
            return result;
        }

        private boolean sameRegion(Row center, Row neighbor) {
            return center.locus() == null || neighbor.locus() == null
                || center.locus().chromosome().equals(neighbor.locus().chromosome());
        }

        private void impute(Row center, List<Row> neighbors, double[] values) {
            for (int sample = 0; sample < values.length; sample++) {
                if (!QMissingnessScan.isMissing(values[sample]))
                    continue;
                double bestDistance = Double.POSITIVE_INFINITY;
                double donorSum = 0;
                int tiedDonors = 0;
                for (int donor = 0; donor < values.length; donor++) {
                    double target = center.values()[donor];
                    if (QMissingnessScan.isMissing(target))
                        continue;
                    double distance = 0;
                    int compared = 0;
                    for (Row neighbor : neighbors) {
                        double sampleValue = neighbor.values()[sample];
                        double donorValue = neighbor.values()[donor];
                        if (!QMissingnessScan.isMissing(sampleValue)
                            && !QMissingnessScan.isMissing(donorValue)) {
                            double difference = sampleValue - donorValue;
                            distance += difference * difference;
                            compared++;
                        }
                    }
                    if (compared == 0)
                        continue;
                    distance /= compared;
                    if (distance + TIE_TOLERANCE < bestDistance) {
                        bestDistance = distance;
                        donorSum = target;
                        tiedDonors = 1;
                    } else if (Math.abs(distance - bestDistance) <= TIE_TOLERANCE) {
                        donorSum += target;
                        tiedDonors++;
                    }
                }
                if (tiedDonors > 0) {
                    values[sample] = donorSum / tiedDonors;
                } else {
                    values[sample] = scan.rowMean(center.index());
                    if (!Double.isFinite(values[sample]))
                        throw new IllegalArgumentException("Predictor row '" + center.id()
                            + "' has neither a local-pattern donor nor an observed row mean");
                    if (!warnedMeanFallback) {
                        warnedMeanFallback = true;
                        System.err.println("WARNING: At least one local-pattern genotype imputation had no "
                            + "comparable donor and fell back to the predictor-row mean.");
                    }
                }
            }
        }

        @Override
        public void close() throws IOException {
            if (!closed) {
                closed = true;
                inner.close();
            }
        }
    }

    private static Locus parseLocus(String id) {
        Matcher matcher = LOCUS.matcher(id == null ? "" : id.trim());
        if (!matcher.matches())
            return null;
        try {
            return new Locus(matcher.group(1), Long.parseLong(matcher.group(2)));
        } catch (NumberFormatException e) {
            return null;
        }
    }
}
