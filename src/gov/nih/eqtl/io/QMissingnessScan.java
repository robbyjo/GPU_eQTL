/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import static gov.nih.utils.QDataUtils.kUndefinedValue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Exact row/sample missingness metadata collected without materializing a matrix. */
public final class QMissingnessScan {
    public record MissingRow(long rowIndex, String rowId, int missingValues, int patternId) { }

    public record Pattern(int id, BitSet missingSamples, long[] rowIndices) {
        public Pattern {
            missingSamples = (BitSet) missingSamples.clone();
            rowIndices = rowIndices.clone();
        }

        @Override
        public BitSet missingSamples() {
            return (BitSet) missingSamples.clone();
        }

        @Override
        public long[] rowIndices() {
            return rowIndices.clone();
        }
    }

    private static final class MaskKey {
        private final long[] words;
        private final int hash;

        MaskKey(BitSet mask) {
            words = mask.toLongArray();
            hash = Arrays.hashCode(words);
        }

        @Override
        public int hashCode() {
            return hash;
        }

        @Override
        public boolean equals(Object other) {
            return other instanceof MaskKey key && Arrays.equals(words, key.words);
        }
    }

    private static final class PatternBuilder {
        private final int id;
        private final BitSet missing;
        private final List<Long> rows = new ArrayList<>();

        PatternBuilder(int id, BitSet missing) {
            this.id = id;
            this.missing = (BitSet) missing.clone();
        }
    }

    private final String matrixName;
    private final String[] sampleIds;
    private final long rowCount;
    private final long totalMissingValues;
    private final long[] sampleMissingValues;
    private final double[] rowMeans;
    private final List<MissingRow> missingRows;
    private final List<Pattern> patterns;
    private final BitSet rowsWithMissing;

    private QMissingnessScan(String matrixName, String[] sampleIds, long rowCount,
        long totalMissingValues, long[] sampleMissingValues, double[] rowMeans,
        List<MissingRow> missingRows, List<Pattern> patterns, BitSet rowsWithMissing) {
        this.matrixName = matrixName;
        this.sampleIds = sampleIds.clone();
        this.rowCount = rowCount;
        this.totalMissingValues = totalMissingValues;
        this.sampleMissingValues = sampleMissingValues.clone();
        this.rowMeans = rowMeans.clone();
        this.missingRows = List.copyOf(missingRows);
        this.patterns = List.copyOf(patterns);
        this.rowsWithMissing = (BitSet) rowsWithMissing.clone();
    }

    public static QMissingnessScan scan(String matrixName, QMatrixRowSource source,
        int[] columnOrder) throws IOException {
        int selectedColumns = columnOrder == null ? source.metadata().columnCount() : columnOrder.length;
        if (selectedColumns <= 0)
            throw new IllegalArgumentException(matrixName + " has no selected samples");
        if (source.metadata().rowCount() > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Missingness scan supports at most " + Integer.MAX_VALUE + " rows");

        String[] sourceIds = source.metadata().sampleIds();
        String[] selectedIds = new String[selectedColumns];
        if (columnOrder == null) {
            System.arraycopy(sourceIds, 0, selectedIds, 0, selectedColumns);
        } else {
            for (int i = 0; i < selectedColumns; i++)
                selectedIds[i] = sourceIds[columnOrder[i]];
        }

        long[] sampleMissing = new long[selectedColumns];
        double[] means = new double[(int) source.metadata().rowCount()];
        List<MissingRow> missingRows = new ArrayList<>();
        LinkedHashMap<MaskKey, PatternBuilder> grouped = new LinkedHashMap<>();
        BitSet rowsWithMissing = new BitSet(means.length);
        long totalMissing = 0;
        long rowsSeen = 0;

        try (QMatrixRowSource.BlockReader reader = source.open(columnOrder)) {
            QMatrixRowSource.Block block;
            while ((block = reader.readBlock(1024)) != null) {
                for (int row = 0; row < block.rowCount(); row++) {
                    long rowIndex = block.rowOffset() + row;
                    if (rowIndex < 0 || rowIndex >= means.length)
                        throw new IOException("Unexpected row offset while scanning " + matrixName);
                    BitSet missing = new BitSet(selectedColumns);
                    double sum = 0;
                    int observed = 0;
                    for (int sample = 0; sample < selectedColumns; sample++) {
                        double value = block.values()[row][sample];
                        if (isMissing(value)) {
                            missing.set(sample);
                            sampleMissing[sample]++;
                            totalMissing++;
                        } else {
                            sum += value;
                            observed++;
                        }
                    }
                    means[(int) rowIndex] = observed == 0 ? Double.NaN : sum / observed;
                    MaskKey key = new MaskKey(missing);
                    PatternBuilder pattern = grouped.get(key);
                    if (pattern == null) {
                        pattern = new PatternBuilder(grouped.size(), missing);
                        grouped.put(key, pattern);
                    }
                    pattern.rows.add(rowIndex);
                    if (!missing.isEmpty()) {
                        rowsWithMissing.set((int) rowIndex);
                        missingRows.add(new MissingRow(rowIndex, block.rowIds()[row],
                            missing.cardinality(), pattern.id));
                    }
                    rowsSeen++;
                }
            }
        }
        if (rowsSeen != source.metadata().rowCount())
            throw new IOException(matrixName + " row count changed during missingness scan: expected "
                + source.metadata().rowCount() + ", found " + rowsSeen);

        List<Pattern> patterns = new ArrayList<>(grouped.size());
        for (PatternBuilder builder : grouped.values()) {
            long[] rowIndices = new long[builder.rows.size()];
            for (int i = 0; i < rowIndices.length; i++)
                rowIndices[i] = builder.rows.get(i);
            patterns.add(new Pattern(builder.id, builder.missing, rowIndices));
        }
        return new QMissingnessScan(matrixName, selectedIds, rowsSeen, totalMissing,
            sampleMissing, means, missingRows, patterns, rowsWithMissing);
    }

    public static boolean isMissing(double value) {
        return value == kUndefinedValue || !Double.isFinite(value);
    }

    public String matrixName() { return matrixName; }
    public String[] sampleIds() { return sampleIds.clone(); }
    public long rowCount() { return rowCount; }
    public int sampleCount() { return sampleIds.length; }
    public long totalMissingValues() { return totalMissingValues; }
    public long totalValues() { return Math.multiplyExact(rowCount, sampleIds.length); }
    public boolean hasMissingValues() { return totalMissingValues > 0; }
    public long[] sampleMissingValues() { return sampleMissingValues.clone(); }
    public List<MissingRow> missingRows() { return missingRows; }
    public List<Pattern> patterns() { return patterns; }
    public BitSet rowsWithMissing() { return (BitSet) rowsWithMissing.clone(); }

    public double rowMean(long rowIndex) {
        if (rowIndex < 0 || rowIndex >= rowMeans.length)
            throw new IllegalArgumentException("Row index is outside the missingness scan: " + rowIndex);
        return rowMeans[(int) rowIndex];
    }

}
