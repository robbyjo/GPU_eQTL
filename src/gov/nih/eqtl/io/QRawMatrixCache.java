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
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;
import java.util.zip.CRC32;

/** Durable aligned FP64 rows before residualization/standardization. */
public final class QRawMatrixCache implements QMatrixRowSource, AutoCloseable {
    private static final int MAGIC = 0x51524157; // QRAW
    private static final int VERSION = 1;
    private static final int INDEX_MAGIC = 0x51524958; // QRIX
    private static final int INDEX_VERSION = 2;

    private final Path path;
    private final String signature;
    private final Metadata metadata;
    private final AtomicLong indexedSelectionCalls = new AtomicLong();
    private final AtomicLong indexedRowsRead = new AtomicLong();
    private final AtomicLong indexedNumericBytesRead = new AtomicLong();
    private volatile RowIndex rowIndex;
    private volatile boolean persistentIndexReused;

    public record IndexedReadStatistics(long selectionCalls, long selectedRows,
        long numericBytesRead, long indexedRows, boolean persistentIndexReused) { }

    private record RowIndex(String[] ids, long[] offsets) { }

    public static QRawMatrixCache openOrBuild(Path cacheDirectory, String signature,
        QMatrixRowSource source, int[] columnOrder, int rowsPerBlock, boolean rebuild)
        throws IOException {
        if (rowsPerBlock <= 0)
            throw new IllegalArgumentException("Raw cache build block must be positive");
        Files.createDirectories(cacheDirectory);
        Path cachePath = cachePath(cacheDirectory, signature, source);
        if (!rebuild && Files.isRegularFile(cachePath)) {
            try {
                QRawMatrixCache cache = open(cachePath, signature);
                System.out.println("Reusing aligned raw predictor cache: " + cachePath);
                return cache;
            } catch (IOException e) {
                System.err.println("Aligned raw predictor cache is unusable and will be rebuilt: "
                    + e.getMessage());
            }
        }
        System.out.println("Building aligned raw predictor cache: " + cachePath);
        build(cachePath, signature, source, columnOrder, rowsPerBlock);
        return open(cachePath, signature);
    }

    public static QRawMatrixCache openIfPresent(Path cacheDirectory, String signature,
        QMatrixRowSource source) throws IOException {
        Path cachePath = cachePath(cacheDirectory, signature, source);
        if (!Files.isRegularFile(cachePath))
            return null;
        QRawMatrixCache cache = open(cachePath, signature);
        System.out.println("Reusing preprocessed aligned variant cache: " + cachePath);
        return cache;
    }

    public static QRawMatrixCache open(Path path, String expectedSignature) throws IOException {
        return new QRawMatrixCache(path, expectedSignature);
    }

    public static String signature(QMatrixRowSource source, int[] columnOrder) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            update(digest, "gpu-eqtl-aligned-raw-cache-v1");
            Path sourcePath = source.metadata().path().toAbsolutePath().normalize();
            update(digest, sourcePath.toString());
            if (Files.isRegularFile(sourcePath)) {
                update(digest, Long.toString(Files.size(sourcePath)));
                update(digest, Long.toString(Files.getLastModifiedTime(sourcePath).toMillis()));
            }
            update(digest, Long.toString(source.metadata().rowCount()));
            update(digest, String.valueOf(source.metadata().cacheSignatureTag()));
            int[] selected = columnOrder == null ? identity(source.metadata().columnCount()) : columnOrder;
            String[] ids = source.metadata().sampleIds();
            for (int column : selected) {
                if (column < 0 || column >= ids.length)
                    throw new IllegalArgumentException("Invalid raw-cache column index " + column);
                update(digest, Integer.toString(column));
                update(digest, ids[column]);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
    }

    private QRawMatrixCache(Path path, String expectedSignature) throws IOException {
        this.path = path.toAbsolutePath().normalize();
        try (DataInputStream input = openInput(this.path)) {
            if (input.readInt() != MAGIC || input.readInt() != VERSION)
                throw new IOException("Unsupported aligned raw matrix cache " + path);
            signature = readString(input);
            if (!signature.equals(expectedSignature))
                throw new IOException("Aligned raw matrix cache signature mismatch in " + path);
            long rows = input.readLong();
            int columns = input.readInt();
            if (rows <= 0 || rows > Integer.MAX_VALUE || columns <= 0)
                throw new IOException("Invalid aligned raw matrix cache dimensions in " + path);
            String[] sampleIds = new String[columns];
            for (int column = 0; column < columns; column++)
                sampleIds[column] = readString(input);
            metadata = new Metadata(this.path, rows, columns, sampleIds,
                "aligned-raw-v1;signature=" + signature);
        }
    }

    public Path path() {
        return path;
    }

    @Override
    public Metadata metadata() {
        return metadata;
    }

    @Override
    public BlockReader open(int[] columnOrder) throws IOException {
        int[] selected = validateColumns(columnOrder);
        return new Reader(selected);
    }

    /**
     * Read only requested rows, preserving cache/source order and validating every
     * returned row checksum. The lazy row-offset index does not alter the version-1
     * on-disk format.
     */
    public Block readSelected(Set<String> requested, int[] columnOrder) throws IOException {
        if (requested == null)
            throw new IllegalArgumentException("Requested raw-cache row IDs are required");
        int[] selected = validateColumns(columnOrder);
        int[] outputIndex = new int[metadata.columnCount()];
        java.util.Arrays.fill(outputIndex, -1);
        for (int output = 0; output < selected.length; output++)
            outputIndex[selected[output]] = output;
        RowIndex index = rowIndex();
        ArrayList<String> ids = new ArrayList<>();
        ArrayList<double[]> values = new ArrayList<>();
        Set<String> seen = new HashSet<>();
        byte[] numericBytes = new byte[Math.multiplyExact(metadata.columnCount(), Double.BYTES)];
        try (RandomAccessFile input = new RandomAccessFile(path.toFile(), "r")) {
            long nextOffset = -1;
            for (int row = 0; row < index.ids().length; row++) {
                String indexedId = index.ids()[row];
                if (!requested.contains(indexedId)) continue;
                if (!seen.add(indexedId))
                    throw new IOException("Duplicate aligned raw-cache row ID '" + indexedId
                        + "' in " + path);
                if (nextOffset != index.offsets()[row])
                    input.seek(index.offsets()[row]);
                byte[] idBytes = readBytes(input);
                String id = new String(idBytes, StandardCharsets.UTF_8);
                if (!id.equals(indexedId))
                    throw new IOException("Aligned raw-cache row index changed at row " + row
                        + " in " + path);
                CRC32 checksum = new CRC32();
                checksum.update(idBytes);
                input.readFully(numericBytes);
                checksum.update(numericBytes);
                ByteBuffer numeric = ByteBuffer.wrap(numericBytes);
                double[] selectedValues = new double[selected.length];
                for (int column = 0; column < metadata.columnCount(); column++) {
                    int output = outputIndex[column];
                    if (output >= 0)
                        selectedValues[output] = Double.longBitsToDouble(
                            numeric.getLong(column * Double.BYTES));
                }
                int expected = input.readInt();
                if ((int) checksum.getValue() != expected)
                    throw new IOException("Aligned raw cache checksum failure at row " + row
                        + " in " + path);
                ids.add(id);
                values.add(selectedValues);
                nextOffset = input.getFilePointer();
            }
        } catch (EOFException e) {
            throw new IOException("Truncated aligned raw cache " + path, e);
        }
        indexedSelectionCalls.incrementAndGet();
        indexedRowsRead.addAndGet(ids.size());
        indexedNumericBytesRead.addAndGet(Math.multiplyExact((long) ids.size(),
            Math.multiplyExact((long) metadata.columnCount(), Double.BYTES)));
        return new Block(0, ids.toArray(String[]::new), values.toArray(double[][]::new));
    }

    public IndexedReadStatistics indexedReadStatistics() {
        RowIndex index = rowIndex;
        return new IndexedReadStatistics(indexedSelectionCalls.get(), indexedRowsRead.get(),
            indexedNumericBytesRead.get(), index == null ? 0 : index.ids().length,
            persistentIndexReused);
    }

    private int[] validateColumns(int[] columnOrder) {
        int[] selected = columnOrder == null ? identity(metadata.columnCount()) : columnOrder.clone();
        boolean[] seen = new boolean[metadata.columnCount()];
        for (int column : selected) {
            if (column < 0 || column >= seen.length || seen[column])
                throw new IllegalArgumentException("Invalid or duplicate aligned raw-cache column index");
            seen[column] = true;
        }
        return selected;
    }

    @Override public void close() { }

    private RowIndex rowIndex() throws IOException {
        RowIndex cached = rowIndex;
        if (cached != null) return cached;
        synchronized (this) {
            if (rowIndex != null) return rowIndex;
            Path indexPath = indexPath();
            if (Files.isRegularFile(indexPath)) {
                try {
                    rowIndex = loadIndex(indexPath);
                    persistentIndexReused = true;
                    return rowIndex;
                } catch (IOException e) {
                    System.err.println("Aligned raw row index is unusable and will be rebuilt: "
                        + e.getMessage());
                }
            }
            String[] ids = new String[Math.toIntExact(metadata.rowCount())];
            long[] offsets = new long[ids.length];
            Set<String> seen = new HashSet<>();
            try (RandomAccessFile input = new RandomAccessFile(path.toFile(), "r")) {
                validateHeader(input);
                long numericBytes = Math.multiplyExact((long) metadata.columnCount(),
                    Double.BYTES);
                for (int row = 0; row < ids.length; row++) {
                    offsets[row] = input.getFilePointer();
                    ids[row] = readString(input);
                    if (!seen.add(ids[row]))
                        throw new IOException("Duplicate aligned raw-cache row ID '" + ids[row]
                            + "' at row " + row + " in " + path);
                    long checksumOffset = Math.addExact(input.getFilePointer(), numericBytes);
                    if (checksumOffset + Integer.BYTES > input.length())
                        throw new IOException("Truncated aligned raw cache " + path);
                    input.seek(checksumOffset + Integer.BYTES);
                }
                if (input.getFilePointer() != input.length())
                    throw new IOException("Aligned raw cache has trailing data: " + path);
            }
            rowIndex = new RowIndex(ids, offsets);
            writeIndex(indexPath, rowIndex);
            return rowIndex;
        }
    }

    private RowIndex loadIndex(Path indexPath) throws IOException {
        try (DataInputStream input = openInput(indexPath)) {
            if (input.readInt() != INDEX_MAGIC || input.readInt() != INDEX_VERSION)
                throw new IOException("Unsupported aligned raw row index " + indexPath);
            if (!readString(input).equals(signature)
                || input.readLong() != Files.size(path)
                || input.readLong() != Files.getLastModifiedTime(path).toMillis())
                throw new IOException("Aligned raw row index source identity mismatch in "
                    + indexPath);
            int rows = input.readInt();
            if (rows != metadata.rowCount())
                throw new IOException("Aligned raw row index dimension mismatch in " + indexPath);
            String[] ids = new String[rows];
            long[] offsets = new long[rows];
            Set<String> seen = new HashSet<>();
            long previous = -1;
            long sourceLength = Files.size(path);
            for (int row = 0; row < rows; row++) {
                byte[] idBytes = readBytes(input);
                ids[row] = new String(idBytes, StandardCharsets.UTF_8);
                offsets[row] = input.readLong();
                CRC32 checksum = new CRC32();
                checksum.update(idBytes);
                update(checksum, offsets[row]);
                if ((int) checksum.getValue() != input.readInt())
                    throw new IOException("Aligned raw row index checksum failure at row "
                        + row + " in " + indexPath);
                if (!seen.add(ids[row]) || offsets[row] <= previous || offsets[row] >= sourceLength)
                    throw new IOException("Invalid aligned raw row index entry at row " + row
                        + " in " + indexPath);
                previous = offsets[row];
            }
            if (input.read() != -1)
                throw new IOException("Aligned raw row index has trailing data: " + indexPath);
            return new RowIndex(ids, offsets);
        } catch (EOFException e) {
            throw new IOException("Truncated aligned raw row index " + indexPath, e);
        }
    }

    private void writeIndex(Path indexPath, RowIndex index) throws IOException {
        Path temporary = Files.createTempFile(indexPath.getParent(),
            indexPath.getFileName().toString(), ".partial");
        boolean complete = false;
        try (DataOutputStream output = new DataOutputStream(new BufferedOutputStream(
                Files.newOutputStream(temporary)))) {
            output.writeInt(INDEX_MAGIC);
            output.writeInt(INDEX_VERSION);
            writeString(output, signature);
            output.writeLong(Files.size(path));
            output.writeLong(Files.getLastModifiedTime(path).toMillis());
            output.writeInt(index.ids().length);
            for (int row = 0; row < index.ids().length; row++) {
                byte[] idBytes = index.ids()[row].getBytes(StandardCharsets.UTF_8);
                writeBytes(output, idBytes);
                output.writeLong(index.offsets()[row]);
                CRC32 checksum = new CRC32();
                checksum.update(idBytes);
                update(checksum, index.offsets()[row]);
                output.writeInt((int) checksum.getValue());
            }
            output.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(temporary, indexPath);
            else
                Files.deleteIfExists(temporary);
        }
    }

    private Path indexPath() {
        return path.resolveSibling(path.getFileName().toString() + ".idx");
    }

    private void validateHeader(RandomAccessFile input) throws IOException {
        if (input.readInt() != MAGIC || input.readInt() != VERSION
            || !readString(input).equals(signature))
            throw new IOException("Aligned raw cache header changed: " + path);
        long rows = input.readLong();
        int columns = input.readInt();
        if (rows != metadata.rowCount() || columns != metadata.columnCount())
            throw new IOException("Aligned raw cache dimensions changed: " + path);
        String[] sampleIds = metadata.sampleIds();
        for (int column = 0; column < columns; column++)
            if (!readString(input).equals(sampleIds[column]))
                throw new IOException("Aligned raw cache sample IDs changed: " + path);
    }

    private final class Reader implements BlockReader {
        private final DataInputStream input;
        private final int[] selectedColumns;
        private final int[] outputIndex;
        private final byte[] numericBytes;
        private long rowsRead;
        private boolean closed;

        Reader(int[] selectedColumns) throws IOException {
            this.selectedColumns = selectedColumns;
            outputIndex = new int[metadata.columnCount()];
            java.util.Arrays.fill(outputIndex, -1);
            for (int output = 0; output < selectedColumns.length; output++)
                outputIndex[selectedColumns[output]] = output;
            numericBytes = new byte[Math.multiplyExact(metadata.columnCount(), Double.BYTES)];
            input = openInput(path);
            if (input.readInt() != MAGIC || input.readInt() != VERSION
                || !readString(input).equals(signature))
                throw new IOException("Aligned raw cache header changed: " + path);
            long rows = input.readLong();
            int columns = input.readInt();
            if (rows != metadata.rowCount() || columns != metadata.columnCount())
                throw new IOException("Aligned raw cache dimensions changed: " + path);
            for (int column = 0; column < columns; column++)
                if (!readString(input).equals(metadata.sampleIds()[column]))
                    throw new IOException("Aligned raw cache sample IDs changed: " + path);
        }

        @Override
        public Block readBlock(int maximumRows) throws IOException {
            if (maximumRows <= 0)
                throw new IllegalArgumentException("maximumRows must be positive");
            if (closed)
                throw new IOException("Aligned raw cache reader is closed");
            if (rowsRead == metadata.rowCount()) {
                if (input.read() != -1)
                    throw new IOException("Aligned raw cache has trailing data: " + path);
                return null;
            }
            int count = (int) Math.min(maximumRows, metadata.rowCount() - rowsRead);
            String[] ids = new String[count];
            double[][] values = new double[count][selectedColumns.length];
            long offset = rowsRead;
            try {
                for (int row = 0; row < count; row++) {
                    byte[] idBytes = readBytes(input);
                    ids[row] = new String(idBytes, StandardCharsets.UTF_8);
                    CRC32 checksum = new CRC32();
                    checksum.update(idBytes);
                    input.readFully(numericBytes);
                    checksum.update(numericBytes);
                    ByteBuffer numeric = ByteBuffer.wrap(numericBytes);
                    for (int column = 0; column < metadata.columnCount(); column++) {
                        int output = outputIndex[column];
                        if (output >= 0)
                            values[row][output] = Double.longBitsToDouble(
                                numeric.getLong(column * Double.BYTES));
                    }
                    int expected = input.readInt();
                    if ((int) checksum.getValue() != expected)
                        throw new IOException("Aligned raw cache checksum failure at row "
                            + (rowsRead + row) + " in " + path);
                }
            } catch (EOFException e) {
                throw new IOException("Truncated aligned raw cache " + path, e);
            }
            rowsRead += count;
            return new Block(offset, ids, values);
        }

        @Override
        public void close() throws IOException {
            if (!closed) {
                closed = true;
                input.close();
            }
        }
    }

    private static void build(Path cachePath, String signature, QMatrixRowSource source,
        int[] columnOrder, int rowsPerBlock) throws IOException {
        int[] selected = columnOrder == null ? identity(source.metadata().columnCount()) : columnOrder;
        String[] sourceIds = source.metadata().sampleIds();
        String[] selectedIds = new String[selected.length];
        for (int i = 0; i < selected.length; i++)
            selectedIds[i] = sourceIds[selected[i]];
        Path temporary = Files.createTempFile(cachePath.getParent(),
            cachePath.getFileName().toString(), ".partial");
        boolean complete = false;
        long written = 0;
        try (DataOutputStream output = new DataOutputStream(new BufferedOutputStream(
                Files.newOutputStream(temporary)));
             BlockReader reader = source.open(selected)) {
            output.writeInt(MAGIC);
            output.writeInt(VERSION);
            writeString(output, signature);
            output.writeLong(source.metadata().rowCount());
            output.writeInt(selected.length);
            for (String sampleId : selectedIds)
                writeString(output, sampleId);
            Block block;
            while ((block = reader.readBlock(rowsPerBlock)) != null) {
                if (block.rowOffset() != written)
                    throw new IOException("Source row order changed while building aligned raw cache");
                for (int row = 0; row < block.rowCount(); row++) {
                    byte[] id = block.rowIds()[row].getBytes(StandardCharsets.UTF_8);
                    writeBytes(output, id);
                    CRC32 checksum = new CRC32();
                    checksum.update(id);
                    if (block.values()[row].length != selected.length)
                        throw new IOException("Source sample count changed while building aligned raw cache");
                    for (double value : block.values()[row]) {
                        long bits = Double.doubleToLongBits(value);
                        output.writeLong(bits);
                        update(checksum, bits);
                    }
                    output.writeInt((int) checksum.getValue());
                    written++;
                }
            }
            if (written != source.metadata().rowCount())
                throw new IOException("Source row count changed while building aligned raw cache: expected "
                    + source.metadata().rowCount() + ", found " + written);
            output.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(temporary, cachePath);
            else
                Files.deleteIfExists(temporary);
        }
    }

    private static DataInputStream openInput(Path path) throws IOException {
        return new DataInputStream(new BufferedInputStream(Files.newInputStream(path), 1 << 20));
    }

    private static void writeString(DataOutputStream output, String value) throws IOException {
        writeBytes(output, value.getBytes(StandardCharsets.UTF_8));
    }

    private static String readString(DataInputStream input) throws IOException {
        return new String(readBytes(input), StandardCharsets.UTF_8);
    }

    private static void writeBytes(DataOutputStream output, byte[] value) throws IOException {
        output.writeInt(value.length);
        output.write(value);
    }

    private static byte[] readBytes(DataInputStream input) throws IOException {
        int length = input.readInt();
        if (length < 0 || length > 16 * 1024 * 1024)
            throw new IOException("Invalid aligned raw cache string length " + length);
        byte[] value = input.readNBytes(length);
        if (value.length != length)
            throw new EOFException("Truncated aligned raw cache string");
        return value;
    }

    private static String readString(RandomAccessFile input) throws IOException {
        return new String(readBytes(input), StandardCharsets.UTF_8);
    }

    private static byte[] readBytes(RandomAccessFile input) throws IOException {
        int length = input.readInt();
        if (length < 0 || length > 16 * 1024 * 1024)
            throw new IOException("Invalid aligned raw cache string length " + length);
        byte[] value = new byte[length];
        input.readFully(value);
        return value;
    }

    private static void update(CRC32 checksum, long value) {
        for (int shift = 56; shift >= 0; shift -= 8)
            checksum.update((int) (value >>> shift) & 0xff);
    }

    private static void update(MessageDigest digest, String value) {
        digest.update((byte) 0);
        digest.update(value.getBytes(StandardCharsets.UTF_8));
    }

    private static Path cachePath(Path cacheDirectory, String signature,
        QMatrixRowSource source) {
        String sourceName = source.metadata().path().getFileName().toString()
            .replaceAll("[^A-Za-z0-9_.-]", "_");
        return cacheDirectory.resolve(sourceName + "-" + signature + ".qraw")
            .toAbsolutePath().normalize();
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
