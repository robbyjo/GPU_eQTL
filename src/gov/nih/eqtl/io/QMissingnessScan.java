/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import static gov.nih.utils.QDataUtils.kUndefinedValue;

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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.zip.CheckedInputStream;
import java.util.zip.CheckedOutputStream;
import java.util.zip.CRC32;

/** Exact row/sample missingness metadata collected without materializing a matrix. */
public final class QMissingnessScan {
    private static final int CACHE_MAGIC = 0x514d4953; // QMIS
    private static final int CACHE_VERSION = 1;

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

    /** Load a signed exact scan when possible, otherwise scan once and commit it atomically. */
    public static QMissingnessScan scanOrLoad(String matrixName, QMatrixRowSource source,
        int[] columnOrder, Path cacheDirectory, boolean rebuild) throws IOException {
        if (cacheDirectory == null)
            return scan(matrixName, source, columnOrder);
        String signature = cacheSignature(source, columnOrder);
        Path path = cachePath(cacheDirectory, source, signature);
        if (!rebuild && Files.isRegularFile(path)) {
            try {
                QMissingnessScan cached = readCache(matrixName, path, signature, source, columnOrder);
                System.out.println("Reusing " + matrixName + " missingness scan: " + path);
                return cached;
            } catch (IOException | RuntimeException e) {
                System.err.println("Missingness scan cache is unusable and will be rebuilt: "
                    + e.getMessage());
            }
        }
        QMissingnessScan result = scan(matrixName, source, columnOrder);
        Files.createDirectories(path.getParent());
        writeCache(path, signature, result);
        System.out.println("Saved " + matrixName + " missingness scan: " + path);
        return result;
    }

    static String cacheSignature(QMatrixRowSource source, int[] columnOrder) throws IOException {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
        update(digest, "gpu-eqtl-missingness-cache-v" + CACHE_VERSION);
        Path sourcePath = source.metadata().path().toAbsolutePath().normalize();
        update(digest, sourcePath.toString());
        if (Files.isRegularFile(sourcePath)) {
            update(digest, Files.size(sourcePath));
            update(digest, Files.getLastModifiedTime(sourcePath).toMillis());
        }
        update(digest, source.metadata().rowCount());
        update(digest, source.metadata().columnCount());
        update(digest, source.metadata().cacheSignatureTag());
        int[] selected = columnOrder == null ? identity(source.metadata().columnCount()) : columnOrder;
        String[] ids = source.metadata().sampleIds();
        for (int column : selected) {
            if (column < 0 || column >= ids.length)
                throw new IllegalArgumentException("Invalid missingness-cache column index " + column);
            update(digest, column);
            update(digest, ids[column]);
        }
        return HexFormat.of().formatHex(digest.digest());
    }

    private static void writeCache(Path path, String signature, QMissingnessScan scan)
        throws IOException {
        Path temporary = Files.createTempFile(path.getParent(), path.getFileName().toString(), ".partial");
        boolean complete = false;
        CRC32 checksum = new CRC32();
        try (DataOutputStream output = new DataOutputStream(new CheckedOutputStream(
            new BufferedOutputStream(Files.newOutputStream(temporary), 1 << 20), checksum))) {
            output.writeInt(CACHE_MAGIC);
            output.writeInt(CACHE_VERSION);
            writeString(output, signature);
            output.writeLong(scan.rowCount);
            output.writeInt(scan.sampleIds.length);
            output.writeLong(scan.totalMissingValues);
            for (String id : scan.sampleIds)
                writeString(output, id);
            for (long missing : scan.sampleMissingValues)
                output.writeLong(missing);
            for (double mean : scan.rowMeans)
                output.writeDouble(mean);
            output.writeInt(scan.missingRows.size());
            for (MissingRow row : scan.missingRows) {
                output.writeLong(row.rowIndex());
                writeString(output, row.rowId());
                output.writeInt(row.missingValues());
                output.writeInt(row.patternId());
            }
            output.writeInt(scan.patterns.size());
            for (Pattern pattern : scan.patterns) {
                output.writeInt(pattern.id());
                long[] words = pattern.missingSamples().toLongArray();
                output.writeInt(words.length);
                for (long word : words)
                    output.writeLong(word);
                long[] rows = pattern.rowIndices();
                output.writeInt(rows.length);
                for (long row : rows)
                    output.writeLong(row);
            }
            long expectedChecksum = checksum.getValue();
            output.writeLong(expectedChecksum);
            output.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(temporary, path);
            else
                Files.deleteIfExists(temporary);
        }
    }

    private static QMissingnessScan readCache(String matrixName, Path path, String signature,
        QMatrixRowSource source, int[] columnOrder) throws IOException {
        CRC32 checksum = new CRC32();
        try (DataInputStream input = new DataInputStream(new CheckedInputStream(
            new BufferedInputStream(Files.newInputStream(path), 1 << 20), checksum))) {
            if (input.readInt() != CACHE_MAGIC || input.readInt() != CACHE_VERSION)
                throw new IOException("Unsupported missingness scan cache " + path);
            if (!signature.equals(readString(input)))
                throw new IOException("Missingness scan cache signature mismatch in " + path);
            long rows = input.readLong();
            int samples = input.readInt();
            long totalMissing = input.readLong();
            if (rows != source.metadata().rowCount() || rows <= 0 || rows > Integer.MAX_VALUE
                || samples <= 0 || totalMissing < 0
                || totalMissing > Math.multiplyExact(rows, samples))
                throw new IOException("Invalid missingness scan dimensions in " + path);
            String[] expectedIds = selectedSampleIds(source, columnOrder);
            if (samples != expectedIds.length)
                throw new IOException("Missingness scan sample count changed in " + path);
            String[] sampleIds = new String[samples];
            for (int sample = 0; sample < samples; sample++) {
                sampleIds[sample] = readString(input);
                if (!sampleIds[sample].equals(expectedIds[sample]))
                    throw new IOException("Missingness scan sample order changed in " + path);
            }
            long[] sampleMissing = new long[samples];
            long sampleMissingTotal = 0;
            for (int sample = 0; sample < samples; sample++) {
                sampleMissing[sample] = input.readLong();
                if (sampleMissing[sample] < 0 || sampleMissing[sample] > rows)
                    throw new IOException("Invalid sample missingness count in " + path);
                sampleMissingTotal = Math.addExact(sampleMissingTotal, sampleMissing[sample]);
            }
            if (sampleMissingTotal != totalMissing)
                throw new IOException("Missingness scan totals disagree in " + path);
            double[] rowMeans = new double[(int) rows];
            for (int row = 0; row < rowMeans.length; row++)
                rowMeans[row] = input.readDouble();

            int missingRowCount = checkedCount(input.readInt(), rows, "missing row", path);
            List<MissingRow> missingRows = new ArrayList<>(missingRowCount);
            BitSet rowsWithMissing = new BitSet((int) rows);
            long rowMissingTotal = 0;
            for (int row = 0; row < missingRowCount; row++) {
                long rowIndex = input.readLong();
                String rowId = readString(input);
                int missingValues = input.readInt();
                int patternId = input.readInt();
                if (rowIndex < 0 || rowIndex >= rows || rowsWithMissing.get((int) rowIndex)
                    || missingValues <= 0 || missingValues > samples || patternId < 0)
                    throw new IOException("Invalid missing row in " + path);
                rowsWithMissing.set((int) rowIndex);
                rowMissingTotal = Math.addExact(rowMissingTotal, missingValues);
                missingRows.add(new MissingRow(rowIndex, rowId, missingValues, patternId));
            }
            if (rowMissingTotal != totalMissing)
                throw new IOException("Missingness row totals disagree in " + path);

            int patternCount = checkedCount(input.readInt(), rows, "pattern", path);
            if (patternCount == 0)
                throw new IOException("Missingness scan cache has no patterns: " + path);
            List<Pattern> patterns = new ArrayList<>(patternCount);
            BitSet assignedRows = new BitSet((int) rows);
            int[] rowPatterns = new int[(int) rows];
            Arrays.fill(rowPatterns, -1);
            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                int id = input.readInt();
                if (id != patternIndex)
                    throw new IOException("Missingness pattern order changed in " + path);
                int wordCount = checkedCount(input.readInt(), (samples + 63L) / 64L,
                    "pattern-mask word", path);
                long[] words = new long[wordCount];
                for (int word = 0; word < wordCount; word++)
                    words[word] = input.readLong();
                BitSet mask = BitSet.valueOf(words);
                if (mask.length() > samples)
                    throw new IOException("Missingness mask exceeds the sample count in " + path);
                int patternRows = checkedCount(input.readInt(), rows, "pattern row", path);
                if (patternRows == 0)
                    throw new IOException("Missingness pattern has no rows in " + path);
                long[] rowIndices = new long[patternRows];
                for (int row = 0; row < patternRows; row++) {
                    long rowIndex = input.readLong();
                    if (rowIndex < 0 || rowIndex >= rows || assignedRows.get((int) rowIndex))
                        throw new IOException("Invalid or duplicate pattern row in " + path);
                    assignedRows.set((int) rowIndex);
                    rowPatterns[(int) rowIndex] = id;
                    rowIndices[row] = rowIndex;
                }
                patterns.add(new Pattern(id, mask, rowIndices));
            }
            if (assignedRows.cardinality() != rows)
                throw new IOException("Missingness patterns do not cover every row in " + path);
            for (MissingRow row : missingRows) {
                if (row.patternId() >= patterns.size()
                    || rowPatterns[(int) row.rowIndex()] != row.patternId()
                    || patterns.get(row.patternId()).missingSamples().cardinality() != row.missingValues())
                    throw new IOException("Missing row/pattern metadata disagree in " + path);
            }
            long calculatedChecksum = checksum.getValue();
            long expectedChecksum = input.readLong();
            if (calculatedChecksum != expectedChecksum || input.read() != -1)
                throw new IOException("Missingness scan checksum failure in " + path);
            return new QMissingnessScan(matrixName, sampleIds, rows, totalMissing,
                sampleMissing, rowMeans, missingRows, patterns, rowsWithMissing);
        } catch (EOFException e) {
            throw new IOException("Truncated missingness scan cache " + path, e);
        } catch (ArithmeticException e) {
            throw new IOException("Invalid missingness scan counts in " + path, e);
        }
    }

    private static int checkedCount(int count, long maximum, String label, Path path)
        throws IOException {
        if (count < 0 || count > maximum)
            throw new IOException("Invalid " + label + " count in " + path);
        return count;
    }

    private static String[] selectedSampleIds(QMatrixRowSource source, int[] columnOrder) {
        String[] ids = source.metadata().sampleIds();
        int[] selected = columnOrder == null ? identity(ids.length) : columnOrder;
        String[] result = new String[selected.length];
        for (int i = 0; i < selected.length; i++)
            result[i] = ids[selected[i]];
        return result;
    }

    private static Path cachePath(Path directory, QMatrixRowSource source, String signature) {
        String sourceName = source.metadata().path().getFileName().toString()
            .replaceAll("[^A-Za-z0-9_.-]", "_");
        return directory.toAbsolutePath().normalize().resolve(
            sourceName + "." + signature.substring(0, 20) + ".qmiss");
    }

    private static void writeString(DataOutputStream output, String value) throws IOException {
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        output.writeInt(bytes.length);
        output.write(bytes);
    }

    private static String readString(DataInputStream input) throws IOException {
        int length = input.readInt();
        if (length < 0 || length > 16 * 1024 * 1024)
            throw new IOException("Invalid missingness-cache string length " + length);
        byte[] bytes = input.readNBytes(length);
        if (bytes.length != length)
            throw new EOFException("Truncated missingness-cache string");
        return new String(bytes, StandardCharsets.UTF_8);
    }

    private static int[] identity(int count) {
        int[] values = new int[count];
        for (int i = 0; i < count; i++)
            values[i] = i;
        return values;
    }

    private static void update(MessageDigest digest, String value) {
        if (value == null) {
            update(digest, -1);
            return;
        }
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        update(digest, bytes.length);
        digest.update(bytes);
    }

    private static void update(MessageDigest digest, long value) {
        for (int shift = 56; shift >= 0; shift -= 8)
            digest.update((byte) (value >>> shift));
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
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
