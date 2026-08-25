/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HexFormat;
import java.util.List;
import java.util.zip.CheckedInputStream;
import java.util.zip.CheckedOutputStream;
import java.util.zip.CRC32;

import gov.nih.eqtl.QMissingValuePolicy;

/**
 * Exact pattern-specific genotype summaries and missing-value filling over selected samples.
 * The durable statistics avoid repeating frequency scans when the same trait mask is reused.
 */
public final class QPatternVariantSource implements QMatrixRowSource {
    private static final int MAGIC = 0x51505653; // QPVS
    private static final int VERSION = 2;
    private static final int READ_BLOCK_ROWS = 1024;
    private static final double TOLERANCE = 1e-12;

    public record Summary(long inputVariants, long includedVariants, long monomorphicVariants,
        long belowMinimumMaf, long belowMinimumMac, long noCallVariants,
        long missingGenotypes, Path cachePath, boolean reused) { }

    private record Statistics(BitSet selected, double[] means, Summary summary) { }

    private final QMatrixRowSource source;
    private final int[] patternColumns;
    private final QMissingValuePolicy missingPolicy;
    private final BitSet selected;
    private final double[] means;
    private final Metadata metadata;
    private final Summary summary;

    public static QPatternVariantSource openOrBuild(Path cacheDirectory, String label,
        QMatrixRowSource source, int[] patternColumns, QMissingValuePolicy missingPolicy,
        boolean applyFrequencyFilters, double minimumMaf, double minimumMac,
        boolean rebuild) throws IOException {
        if (cacheDirectory == null || source == null || patternColumns == null
            || missingPolicy == null)
            throw new IllegalArgumentException("Pattern statistics require a cache, source, columns, and policy");
        if (patternColumns.length < 2)
            throw new IllegalArgumentException("Pattern-specific analysis requires at least two samples");
        if (!Double.isFinite(minimumMaf) || minimumMaf < 0 || minimumMaf > 0.5
            || !Double.isFinite(minimumMac) || minimumMac < 0)
            throw new IllegalArgumentException("Invalid pattern MAF/MAC threshold");
        String signature = signature(source, patternColumns, missingPolicy,
            applyFrequencyFilters, minimumMaf, minimumMac);
        Files.createDirectories(cacheDirectory);
        String safeLabel = label == null ? "pattern" : label.replaceAll("[^A-Za-z0-9_.-]", "_");
        Path cachePath = cacheDirectory.resolve(safeLabel + "-" + signature + ".qpvs")
            .toAbsolutePath().normalize();
        Statistics statistics = null;
        if (!rebuild && Files.isRegularFile(cachePath)) {
            try {
                statistics = read(cachePath, signature, source.metadata().rowCount(),
                    patternColumns.length);
                System.out.println("Reusing pattern-specific variant statistics: " + cachePath);
            } catch (IOException e) {
                System.err.println("Pattern-specific statistics cache is unusable and will be rebuilt: "
                    + e.getMessage());
            }
        }
        if (statistics == null) {
            System.out.println("Building pattern-specific variant statistics: " + cachePath);
            statistics = build(cachePath, signature, source, patternColumns, missingPolicy,
                applyFrequencyFilters, minimumMaf, minimumMac);
        }
        return new QPatternVariantSource(source, patternColumns, missingPolicy,
            statistics, signature);
    }

    private QPatternVariantSource(QMatrixRowSource source, int[] patternColumns,
        QMissingValuePolicy missingPolicy, Statistics statistics, String signature) {
        this.source = source;
        this.patternColumns = patternColumns.clone();
        this.missingPolicy = missingPolicy;
        selected = (BitSet) statistics.selected().clone();
        means = statistics.means().clone();
        summary = statistics.summary();
        String tag = source.metadata().cacheSignatureTag();
        tag = (tag == null ? "" : tag + ";") + "pattern-variant-v1;signature=" + signature
            + ";selected=" + selected.cardinality();
        metadata = new Metadata(source.metadata().path(), selected.cardinality(),
            source.metadata().columnCount(), source.metadata().sampleIds(), tag);
    }

    public Summary summary() {
        return summary;
    }

    @Override
    public Metadata metadata() {
        return metadata;
    }

    @Override
    public BlockReader open(int[] columnOrder) throws IOException {
        int[] requested = columnOrder == null ? identity(source.metadata().columnCount()) : columnOrder;
        if (!Arrays.equals(requested, patternColumns))
            throw new IOException("Pattern-specific variant statistics were built for a different sample order");
        return new Reader(source.open(requested));
    }

    private final class Reader implements BlockReader {
        private final BlockReader inner;
        private Block current;
        private int currentRow;
        private long outputOffset;
        private boolean closed;

        Reader(BlockReader inner) {
            this.inner = inner;
        }

        @Override
        public Block readBlock(int maximumRows) throws IOException {
            if (maximumRows <= 0)
                throw new IllegalArgumentException("maximumRows must be positive");
            if (closed)
                throw new IOException("Pattern variant reader is closed");
            List<String> ids = new ArrayList<>(maximumRows);
            List<double[]> values = new ArrayList<>(maximumRows);
            long offset = outputOffset;
            while (values.size() < maximumRows) {
                if (current == null || currentRow >= current.rowCount()) {
                    current = inner.readBlock(READ_BLOCK_ROWS);
                    currentRow = 0;
                    if (current == null)
                        break;
                }
                int row = Math.toIntExact(current.rowOffset() + currentRow);
                String id = current.rowIds()[currentRow];
                double[] rowValues = current.values()[currentRow++];
                if (!selected.get(row))
                    continue;
                applyMissingPolicy(rowValues, means[row], id);
                ids.add(id);
                values.add(rowValues);
                outputOffset++;
            }
            if (values.isEmpty())
                return null;
            return new Block(offset, ids.toArray(String[]::new), values.toArray(double[][]::new));
        }

        @Override
        public void close() throws IOException {
            if (!closed) {
                closed = true;
                inner.close();
            }
        }
    }

    private void applyMissingPolicy(double[] values, double mean, String rowId) {
        for (int sample = 0; sample < values.length; sample++) {
            if (!QMissingnessScan.isMissing(values[sample]))
                continue;
            if (missingPolicy == QMissingValuePolicy.ZERO) {
                values[sample] = 0;
            } else if (missingPolicy == QMissingValuePolicy.MEAN) {
                if (!Double.isFinite(mean))
                    throw new IllegalArgumentException("Variant '" + rowId
                        + "' has no pattern-observed dosage for mean imputation");
                values[sample] = mean;
            } else {
                throw new IllegalArgumentException("Variant '" + rowId
                    + "' contains a missing dosage under policy " + missingPolicy.optionName());
            }
        }
    }

    private static Statistics build(Path cachePath, String signature, QMatrixRowSource source,
        int[] columns, QMissingValuePolicy policy, boolean applyFrequencyFilters,
        double minimumMaf, double minimumMac) throws IOException {
        long rowCountLong = source.metadata().rowCount();
        if (rowCountLong <= 0 || rowCountLong > Integer.MAX_VALUE)
            throw new IOException("Unsupported variant row count " + rowCountLong);
        int rowCount = (int) rowCountLong;
        BitSet selected = new BitSet(rowCount);
        double[] means = new double[rowCount];
        Arrays.fill(means, Double.NaN);
        long input = 0;
        long monomorphic = 0;
        long belowMaf = 0;
        long belowMac = 0;
        long noCalls = 0;
        long missingValues = 0;
        Path parent = cachePath.getParent();
        Path temporary = Files.createTempFile(parent, cachePath.getFileName().toString(), ".partial");
        boolean complete = false;
        CheckedOutputStream checked = new CheckedOutputStream(
            new BufferedOutputStream(Files.newOutputStream(temporary)), new CRC32());
        try (DataOutputStream output = new DataOutputStream(checked);
             BlockReader reader = source.open(columns)) {
            output.writeInt(MAGIC);
            output.writeInt(VERSION);
            output.writeUTF(signature);
            output.writeInt(rowCount);
            output.writeInt(columns.length);
            Block block;
            while ((block = reader.readBlock(READ_BLOCK_ROWS)) != null) {
                for (int blockRow = 0; blockRow < block.rowCount(); blockRow++) {
                    int row = Math.toIntExact(block.rowOffset() + blockRow);
                    if (row != input)
                        throw new IOException("Variant row order changed while building pattern statistics");
                    double sum = 0;
                    double sumSquares = 0;
                    int called = 0;
                    int missing = 0;
                    for (double value : block.values()[blockRow]) {
                        if (QMissingnessScan.isMissing(value)) {
                            missing++;
                        } else {
                            sum += value;
                            sumSquares += value * value;
                            called++;
                        }
                    }
                    missingValues += missing;
                    double mean = called == 0 ? Double.NaN : sum / called;
                    means[row] = mean;
                    double alleleNumber = 2.0 * called;
                    double eaf = called == 0 ? Double.NaN : sum / alleleNumber;
                    double maf = called == 0 ? Double.NaN : Math.min(eaf, 1 - eaf);
                    double mac = called == 0 ? Double.NaN : Math.min(sum, alleleNumber - sum);
                    boolean include = true;
                    int reason = 0;
                    if (called == 0) {
                        include = false;
                        noCalls++;
                        reason |= 1;
                    } else {
                        int varianceSamples = policy == QMissingValuePolicy.ZERO
                            ? columns.length : called;
                        double centered = sumSquares - sum * sum / varianceSamples;
                        if (!(centered > TOLERANCE) || !Double.isFinite(centered)) {
                            include = false;
                            monomorphic++;
                            reason |= 2;
                        }
                        if (maf + TOLERANCE < minimumMaf) {
                            belowMaf++;
                            reason |= 4;
                            if (applyFrequencyFilters)
                                include = false;
                        }
                        if (mac + TOLERANCE < minimumMac) {
                            belowMac++;
                            reason |= 8;
                            if (applyFrequencyFilters)
                                include = false;
                        }
                    }
                    if (missing > 0 && policy != QMissingValuePolicy.MEAN
                        && policy != QMissingValuePolicy.ZERO) {
                        throw new IOException("Variant '" + block.rowIds()[blockRow]
                            + "' contains a missing dosage under policy " + policy.optionName());
                    }
                    if (include)
                        selected.set(row);
                    output.writeBoolean(include);
                    output.writeByte(reason);
                    output.writeInt(called);
                    output.writeInt(missing);
                    output.writeDouble(sum);
                    output.writeDouble(sumSquares);
                    output.writeDouble(mean);
                    output.writeDouble(eaf);
                    output.writeDouble(maf);
                    output.writeDouble(mac);
                    input++;
                }
            }
            if (input != rowCount)
                throw new IOException("Variant row count changed while building pattern statistics: expected "
                    + rowCount + ", found " + input);
            output.writeLong(checked.getChecksum().getValue());
            output.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(temporary, cachePath);
            else
                Files.deleteIfExists(temporary);
        }
        Summary summary = new Summary(input, selected.cardinality(), monomorphic,
            belowMaf, belowMac, noCalls, missingValues, cachePath, false);
        return new Statistics(selected, means, summary);
    }

    private static Statistics read(Path cachePath, String signature, long expectedRows,
        int expectedSamples) throws IOException {
        BitSet selected = new BitSet(Math.toIntExact(expectedRows));
        double[] means = new double[Math.toIntExact(expectedRows)];
        long monomorphic = 0;
        long belowMaf = 0;
        long belowMac = 0;
        long noCalls = 0;
        long missingValues = 0;
        CheckedInputStream checked = new CheckedInputStream(
            new BufferedInputStream(Files.newInputStream(cachePath)), new CRC32());
        try (DataInputStream input = new DataInputStream(checked)) {
            if (input.readInt() != MAGIC || input.readInt() != VERSION)
                throw new IOException("Unsupported pattern statistics cache " + cachePath);
            if (!input.readUTF().equals(signature))
                throw new IOException("Pattern statistics signature mismatch in " + cachePath);
            int rows = input.readInt();
            int samples = input.readInt();
            if (rows != expectedRows || samples != expectedSamples)
                throw new IOException("Pattern statistics dimensions differ in " + cachePath);
            for (int row = 0; row < rows; row++) {
                boolean include = input.readBoolean();
                int reason = input.readUnsignedByte();
                input.readInt(); // called
                int missing = input.readInt();
                input.readDouble(); // sum
                input.readDouble(); // sum of squares
                means[row] = input.readDouble();
                input.readDouble(); // EAF
                input.readDouble(); // MAF
                input.readDouble(); // MAC
                if (include) selected.set(row);
                if ((reason & 1) != 0) noCalls++;
                if ((reason & 2) != 0) monomorphic++;
                if ((reason & 4) != 0) belowMaf++;
                if ((reason & 8) != 0) belowMac++;
                missingValues += missing;
            }
            long observedChecksum = checked.getChecksum().getValue();
            long expectedChecksum = input.readLong();
            if (observedChecksum != expectedChecksum)
                throw new IOException("Pattern statistics checksum failure in " + cachePath);
            if (input.read() != -1)
                throw new IOException("Pattern statistics cache has trailing data: " + cachePath);
        } catch (EOFException e) {
            throw new IOException("Truncated pattern statistics cache " + cachePath, e);
        }
        Summary summary = new Summary(expectedRows, selected.cardinality(), monomorphic,
            belowMaf, belowMac, noCalls, missingValues, cachePath, true);
        return new Statistics(selected, means, summary);
    }

    private static String signature(QMatrixRowSource source, int[] columns,
        QMissingValuePolicy policy, boolean applyFrequencyFilters,
        double minimumMaf, double minimumMac) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            update(digest, "gpu-eqtl-pattern-variant-statistics-v1");
            Path path = source.metadata().path().toAbsolutePath().normalize();
            update(digest, path.toString());
            if (Files.isRegularFile(path)) {
                update(digest, Long.toString(Files.size(path)));
                update(digest, Long.toString(Files.getLastModifiedTime(path).toMillis()));
            }
            update(digest, Long.toString(source.metadata().rowCount()));
            update(digest, Integer.toString(source.metadata().columnCount()));
            update(digest, String.valueOf(source.metadata().cacheSignatureTag()));
            update(digest, policy.name());
            update(digest, Boolean.toString(applyFrequencyFilters));
            update(digest, Double.toHexString(minimumMaf));
            update(digest, Double.toHexString(minimumMac));
            for (int column : columns) {
                update(digest, Integer.toString(column));
                update(digest, source.metadata().sampleIds()[column]);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
    }

    private static void update(MessageDigest digest, String value) {
        digest.update((byte) 0);
        digest.update(value.getBytes(StandardCharsets.UTF_8));
    }

    private static int[] identity(int count) {
        int[] values = new int[count];
        for (int i = 0; i < count; i++)
            values[i] = i;
        return values;
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
