/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
package gov.nih.eqtl.io;

import java.io.EOFException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.FileTime;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HexFormat;
import java.util.zip.CRC32;

import gov.nih.eqtl.QeQTLPreprocessor;
import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;

/** Indexed, reusable FP64 rows after alignment, residualization, and standardization. */
public final class QBinaryMatrixCache implements AutoCloseable {
    private static final int MAGIC = 0x514d4348; // QMCH
    private static final int INDEX_MAGIC = 0x514d4958; // QMIX
    private static final int VERSION = 1;

    private final Path path;
    private final RandomAccessFile file;
    private final String kind;
    private final String signature;
    private final int sampleCount;
    private final long rowCount;
    private final long[] rowOffsets;

    private QBinaryMatrixCache(Path path, String expectedKind, String expectedSignature) throws IOException {
        this.path = path.toAbsolutePath().normalize();
        file = new RandomAccessFile(this.path.toFile(), "r");
        try {
            if (file.readInt() != MAGIC)
                throw new IOException("Not a GPU eQTL matrix cache: " + path);
            int version = file.readInt();
            if (version != VERSION)
                throw new IOException("Unsupported matrix cache version " + version + " in " + path);
            sampleCount = file.readInt();
            rowCount = file.readLong();
            long indexOffset = file.readLong();
            kind = readString(file);
            signature = readString(file);
            if (!kind.equals(expectedKind))
                throw new IOException("Matrix cache kind mismatch in " + path);
            if (!signature.equals(expectedSignature))
                throw new IOException("Matrix cache signature mismatch in " + path);
            if (sampleCount <= 0 || rowCount <= 0 || rowCount > Integer.MAX_VALUE
                || indexOffset <= file.getFilePointer() || indexOffset >= file.length())
                throw new IOException("Invalid matrix cache header in " + path);

            file.seek(indexOffset);
            if (file.readInt() != INDEX_MAGIC || file.readLong() != rowCount)
                throw new IOException("Invalid matrix cache index in " + path);
            rowOffsets = new long[(int) rowCount];
            long previous = -1;
            for (int i = 0; i < rowOffsets.length; i++) {
                rowOffsets[i] = file.readLong();
                if (rowOffsets[i] <= previous || rowOffsets[i] >= indexOffset)
                    throw new IOException("Invalid row offset in matrix cache " + path);
                previous = rowOffsets[i];
            }
        } catch (IOException | RuntimeException e) {
            file.close();
            throw e;
        }
    }

    public static QBinaryMatrixCache open(Path path, String kind, String signature) throws IOException {
        return new QBinaryMatrixCache(path, kind, signature);
    }

    public static QBinaryMatrixCache openOrBuild(Path cacheDirectory, String kind, String signature,
        QDelimitedMatrixSource source, int[] columnOrder, double[][] covariateQ,
        int rowsPerBuildBlock, boolean rebuild) throws IOException {
        Files.createDirectories(cacheDirectory);
        Path cachePath = cachePath(cacheDirectory, source.metadata().path(), kind, signature);
        if (!rebuild && Files.isRegularFile(cachePath)) {
            try {
                QBinaryMatrixCache cache = open(cachePath, kind, signature);
                System.out.println("Reusing " + kind + " cache: " + cachePath);
                return cache;
            } catch (IOException e) {
                System.err.println("Matrix cache is unusable and will be rebuilt: " + e.getMessage());
            }
        }
        System.out.println("Building " + kind + " cache: " + cachePath);
        build(cachePath, kind, signature, source, columnOrder, covariateQ, rowsPerBuildBlock);
        return open(cachePath, kind, signature);
    }

    private static void build(Path cachePath, String kind, String signature,
        QDelimitedMatrixSource source, int[] columnOrder, double[][] covariateQ,
        int rowsPerBuildBlock) throws IOException {
        Path parent = cachePath.toAbsolutePath().getParent();
        Path temporary = Files.createTempFile(parent, cachePath.getFileName().toString(), ".partial");
        boolean complete = false;
        try (RandomAccessFile output = new RandomAccessFile(temporary.toFile(), "rw");
             QDelimitedMatrixSource.BlockReader reader = source.open(columnOrder)) {
            output.writeInt(MAGIC);
            output.writeInt(VERSION);
            output.writeInt(source.metadata().columnCount());
            long rowCountPosition = output.getFilePointer();
            output.writeLong(0);
            long indexOffsetPosition = output.getFilePointer();
            output.writeLong(0);
            writeString(output, kind);
            writeString(output, signature);

            long expectedRows = source.metadata().rowCount();
            if (expectedRows <= 0 || expectedRows > Integer.MAX_VALUE)
                throw new IOException("Unsupported source row count " + expectedRows);
            long[] offsets = new long[(int) expectedRows];
            int written = 0;
            QDelimitedMatrixSource.Block raw;
            while ((raw = reader.readBlock(rowsPerBuildBlock)) != null) {
                PreparedBlock prepared = QeQTLPreprocessor.prepare(raw, covariateQ, kind);
                for (int row = 0; row < prepared.values().length; row++) {
                    if (written >= offsets.length)
                        throw new IOException("Source grew while building cache: " + source.metadata().path());
                    offsets[written++] = output.getFilePointer();
                    writeRow(output, prepared.rowIds()[row], prepared.standardDeviations()[row], prepared.values()[row]);
                }
            }
            if (written != offsets.length)
                throw new IOException("Source row count changed while building cache: expected "
                    + offsets.length + ", found " + written);

            long indexOffset = output.getFilePointer();
            output.writeInt(INDEX_MAGIC);
            output.writeLong(written);
            for (long offset : offsets)
                output.writeLong(offset);
            output.seek(rowCountPosition);
            output.writeLong(written);
            output.seek(indexOffsetPosition);
            output.writeLong(indexOffset);
            output.getFD().sync();
            complete = true;
        } finally {
            if (complete) {
                moveAtomically(temporary, cachePath);
            } else {
                Files.deleteIfExists(temporary);
            }
        }
    }

    public synchronized PreparedBlock readBlock(long rowOffset, int maximumRows) throws IOException {
        if (rowOffset < 0 || rowOffset >= rowCount)
            throw new IllegalArgumentException("Row offset is outside the cache: " + rowOffset);
        if (maximumRows <= 0)
            throw new IllegalArgumentException("maximumRows must be positive");
        int count = (int) Math.min(maximumRows, rowCount - rowOffset);
        String[] rowIds = new String[count];
        double[][] values = new double[count][sampleCount];
        double[] standardDeviations = new double[count];
        for (int row = 0; row < count; row++) {
            file.seek(rowOffsets[(int) rowOffset + row]);
            CRC32 checksum = new CRC32();
            rowIds[row] = readChecksummedString(file, checksum);
            standardDeviations[row] = readChecksummedDouble(file, checksum);
            for (int sample = 0; sample < sampleCount; sample++)
                values[row][sample] = readChecksummedDouble(file, checksum);
            long expectedChecksum = file.readLong();
            if (checksum.getValue() != expectedChecksum)
                throw new IOException("Checksum failure at cache row " + (rowOffset + row) + " in " + path);
        }
        return new PreparedBlock(rowOffset, rowIds, values, standardDeviations);
    }

    public Path path() {
        return path;
    }

    public String signature() {
        return signature;
    }

    public int sampleCount() {
        return sampleCount;
    }

    public long rowCount() {
        return rowCount;
    }

    @Override
    public void close() throws IOException {
        file.close();
    }

    public static String signature(String kind, QDelimitedMatrixSource source,
        int[] columnOrder, double[][] covariateQ) throws IOException {
        MessageDigest digest = sha256();
        update(digest, "gpu-eqtl-cache-v" + VERSION);
        update(digest, kind);
        Path path = source.metadata().path().toAbsolutePath().normalize();
        update(digest, path.toString());
        update(digest, Files.size(path));
        FileTime modified = Files.getLastModifiedTime(path);
        update(digest, modified.toMillis());
        update(digest, source.metadata().rowCount());
        update(digest, source.metadata().columnCount());
        for (String sampleId : source.metadata().sampleIds())
            update(digest, sampleId);
        for (int column : columnOrder)
            update(digest, column);
        if (covariateQ == null) {
            update(digest, -1);
        } else {
            update(digest, covariateQ.length);
            update(digest, covariateQ[0].length);
            for (double[] row : covariateQ)
                for (double value : row)
                    update(digest, Double.doubleToLongBits(value));
        }
        return HexFormat.of().formatHex(digest.digest());
    }

    public static String analysisSignature(QBinaryMatrixCache genotype, QBinaryMatrixCache expression,
        int genotypeRows, int expressionRows, int degreesOfFreedomOffset,
        int errorDegreesOfFreedom, double rSquaredThreshold, boolean simplify, boolean rSquaredOnly) {
        MessageDigest digest = sha256();
        update(digest, "gpu-eqtl-checkpoint-v1");
        update(digest, genotype.signature());
        update(digest, expression.signature());
        update(digest, genotypeRows);
        update(digest, expressionRows);
        update(digest, degreesOfFreedomOffset);
        update(digest, errorDegreesOfFreedom);
        update(digest, Double.doubleToLongBits(rSquaredThreshold));
        update(digest, simplify ? 1 : 0);
        update(digest, rSquaredOnly ? 1 : 0);
        return HexFormat.of().formatHex(digest.digest());
    }

    private static Path cachePath(Path directory, Path source, String kind, String signature) {
        String name = source.getFileName().toString().replaceAll("[^A-Za-z0-9._-]", "_");
        return directory.resolve(name + "." + kind.toLowerCase() + "." + signature.substring(0, 20) + ".qcache");
    }

    private static void writeRow(RandomAccessFile output, String rowId, double standardDeviation,
        double[] values) throws IOException {
        CRC32 checksum = new CRC32();
        writeChecksummedString(output, checksum, rowId);
        writeChecksummedDouble(output, checksum, standardDeviation);
        for (double value : values)
            writeChecksummedDouble(output, checksum, value);
        output.writeLong(checksum.getValue());
    }

    private static void writeString(RandomAccessFile output, String value) throws IOException {
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        output.writeInt(bytes.length);
        output.write(bytes);
    }

    private static String readString(RandomAccessFile input) throws IOException {
        int length = input.readInt();
        if (length < 0 || length > 16 * 1024 * 1024)
            throw new IOException("Invalid string length in matrix cache");
        byte[] bytes = new byte[length];
        input.readFully(bytes);
        return new String(bytes, StandardCharsets.UTF_8);
    }

    private static void writeChecksummedString(RandomAccessFile output, CRC32 checksum,
        String value) throws IOException {
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        output.writeInt(bytes.length);
        update(checksum, bytes.length);
        output.write(bytes);
        checksum.update(bytes);
    }

    private static String readChecksummedString(RandomAccessFile input, CRC32 checksum) throws IOException {
        int length = input.readInt();
        if (length < 0 || length > 16 * 1024 * 1024)
            throw new EOFException("Invalid row identifier length in matrix cache");
        update(checksum, length);
        byte[] bytes = new byte[length];
        input.readFully(bytes);
        checksum.update(bytes);
        return new String(bytes, StandardCharsets.UTF_8);
    }

    private static void writeChecksummedDouble(RandomAccessFile output, CRC32 checksum,
        double value) throws IOException {
        long bits = Double.doubleToRawLongBits(value);
        output.writeLong(bits);
        update(checksum, bits);
    }

    private static double readChecksummedDouble(RandomAccessFile input, CRC32 checksum) throws IOException {
        long bits = input.readLong();
        update(checksum, bits);
        return Double.longBitsToDouble(bits);
    }

    private static void update(CRC32 checksum, int value) {
        checksum.update((value >>> 24) & 0xff);
        checksum.update((value >>> 16) & 0xff);
        checksum.update((value >>> 8) & 0xff);
        checksum.update(value & 0xff);
    }

    private static void update(CRC32 checksum, long value) {
        for (int shift = 56; shift >= 0; shift -= 8)
            checksum.update((int) (value >>> shift) & 0xff);
    }

    private static MessageDigest sha256() {
        try {
            return MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
    }

    private static void update(MessageDigest digest, String value) {
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        update(digest, bytes.length);
        digest.update(bytes);
    }

    private static void update(MessageDigest digest, int value) {
        digest.update((byte) (value >>> 24));
        digest.update((byte) (value >>> 16));
        digest.update((byte) (value >>> 8));
        digest.update((byte) value);
    }

    private static void update(MessageDigest digest, long value) {
        for (int shift = 56; shift >= 0; shift -= 8)
            digest.update((byte) (value >>> shift));
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
