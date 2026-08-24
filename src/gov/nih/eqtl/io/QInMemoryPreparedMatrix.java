/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.IOException;
import java.util.Arrays;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;

/** Immutable, concurrently readable in-memory view of one verified prepared cache. */
public final class QInMemoryPreparedMatrix implements QPreparedMatrix, AutoCloseable {
    private final String signature;
    private final int sampleCount;
    private String[] rowIds;
    private double[][] values;
    private double[] standardDeviations;

    private QInMemoryPreparedMatrix(String signature, int sampleCount, int rowCount) {
        this.signature = signature;
        this.sampleCount = sampleCount;
        rowIds = new String[rowCount];
        values = new double[rowCount][];
        standardDeviations = new double[rowCount];
    }

    public static QInMemoryPreparedMatrix load(QBinaryMatrixCache cache, int rowsPerBlock)
        throws IOException {
        if (rowsPerBlock <= 0)
            throw new IllegalArgumentException("rowsPerBlock must be positive");
        int rowCount = Math.toIntExact(cache.rowCount());
        QInMemoryPreparedMatrix result = new QInMemoryPreparedMatrix(
            cache.signature(), cache.sampleCount(), rowCount);
        int destination = 0;
        for (long offset = 0; offset < cache.rowCount(); offset += rowsPerBlock) {
            PreparedBlock block = cache.readBlock(offset, rowsPerBlock);
            for (int row = 0; row < block.values().length; row++) {
                result.rowIds[destination] = block.rowIds()[row];
                result.values[destination] = block.values()[row];
                result.standardDeviations[destination] = block.standardDeviations()[row];
                destination++;
            }
        }
        if (destination != rowCount)
            throw new IOException("Prepared trait cache row count changed while loading memory");
        return result;
    }

    public static long estimateResidentBytes(long rowCount, int sampleCount) {
        try {
            long numeric = Math.multiplyExact(Math.multiplyExact(rowCount, sampleCount), Double.BYTES);
            long rowOverhead = Math.multiplyExact(rowCount, 64L);
            return Math.addExact(Math.addExact(numeric, numeric / 8),
                Math.addExact(rowOverhead, 1024L * 1024L));
        } catch (ArithmeticException e) {
            return Long.MAX_VALUE;
        }
    }

    @Override
    public String signature() {
        return signature;
    }

    @Override
    public int sampleCount() {
        return sampleCount;
    }

    @Override
    public long rowCount() {
        ensureOpen();
        return values.length;
    }

    @Override
    public PreparedBlock readBlock(long rowOffset, int maximumRows) {
        ensureOpen();
        if (rowOffset < 0 || rowOffset >= values.length)
            throw new IllegalArgumentException("Row offset is outside the prepared matrix: " + rowOffset);
        if (maximumRows <= 0)
            throw new IllegalArgumentException("maximumRows must be positive");
        int start = Math.toIntExact(rowOffset);
        int end = (int) Math.min(values.length, rowOffset + (long) maximumRows);
        return new PreparedBlock(rowOffset,
            Arrays.copyOfRange(rowIds, start, end),
            Arrays.copyOfRange(values, start, end),
            Arrays.copyOfRange(standardDeviations, start, end));
    }

    @Override
    public void close() {
        rowIds = null;
        values = null;
        standardDeviations = null;
    }

    private void ensureOpen() {
        if (values == null)
            throw new IllegalStateException("Prepared in-memory matrix is closed");
    }
}
