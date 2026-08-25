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
import java.util.HexFormat;
import java.util.zip.CRC32;

/** Durable aligned FP64 rows before residualization/standardization. */
public final class QRawMatrixCache implements QMatrixRowSource, AutoCloseable {
    private static final int MAGIC = 0x51524157; // QRAW
    private static final int VERSION = 1;

    private final Path path;
    private final String signature;
    private final Metadata metadata;

    public static QRawMatrixCache openOrBuild(Path cacheDirectory, String signature,
        QMatrixRowSource source, int[] columnOrder, int rowsPerBlock, boolean rebuild)
        throws IOException {
        if (rowsPerBlock <= 0)
            throw new IllegalArgumentException("Raw cache build block must be positive");
        Files.createDirectories(cacheDirectory);
        String sourceName = source.metadata().path().getFileName().toString()
            .replaceAll("[^A-Za-z0-9_.-]", "_");
        Path cachePath = cacheDirectory.resolve(sourceName + "-" + signature + ".qraw")
            .toAbsolutePath().normalize();
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
        int[] selected = columnOrder == null ? identity(metadata.columnCount()) : columnOrder.clone();
        boolean[] seen = new boolean[metadata.columnCount()];
        for (int column : selected) {
            if (column < 0 || column >= seen.length || seen[column])
                throw new IllegalArgumentException("Invalid or duplicate aligned raw-cache column index");
            seen[column] = true;
        }
        return new Reader(selected);
    }

    @Override public void close() { }

    private final class Reader implements BlockReader {
        private final DataInputStream input;
        private final int[] selectedColumns;
        private final int[] outputIndex;
        private long rowsRead;
        private boolean closed;

        Reader(int[] selectedColumns) throws IOException {
            this.selectedColumns = selectedColumns;
            outputIndex = new int[metadata.columnCount()];
            java.util.Arrays.fill(outputIndex, -1);
            for (int output = 0; output < selectedColumns.length; output++)
                outputIndex[selectedColumns[output]] = output;
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
                    for (int column = 0; column < metadata.columnCount(); column++) {
                        long bits = input.readLong();
                        int output = outputIndex[column];
                        if (output >= 0)
                            values[row][output] = Double.longBitsToDouble(bits);
                        update(checksum, bits);
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

    private static void update(CRC32 checksum, long value) {
        for (int shift = 56; shift >= 0; shift -= 8)
            checksum.update((int) (value >>> shift) & 0xff);
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
