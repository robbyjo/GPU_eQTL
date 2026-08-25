/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import gov.nih.eqtl.QMissingValuePolicy;

/** Applies a declared numeric missing-value policy and optional row selection while streaming. */
public final class QPolicyMatrixSource implements QMatrixRowSource {
    private static final int INNER_BLOCK_ROWS = 1024;

    private final QMatrixRowSource source;
    private final QMissingnessScan scan;
    private final QMissingValuePolicy policy;
    private final BitSet selectedRows;
    private final Metadata metadata;

    public QPolicyMatrixSource(QMatrixRowSource source, QMissingnessScan scan,
        QMissingValuePolicy policy) {
        this(source, scan, policy, null);
    }

    public QPolicyMatrixSource(QMatrixRowSource source, QMissingnessScan scan,
        QMissingValuePolicy policy, BitSet requestedRows) {
        if (source == null || scan == null || policy == null)
            throw new IllegalArgumentException("source, scan, and policy are required");
        if (source.metadata().rowCount() != scan.rowCount())
            throw new IllegalArgumentException("Source and missingness scan row counts differ");
        this.source = source;
        this.scan = scan;
        this.policy = policy;
        selectedRows = requestedRows == null ? allRows(scan.rowCount()) : (BitSet) requestedRows.clone();
        if (selectedRows.length() > scan.rowCount())
            throw new IllegalArgumentException("Selected row index exceeds the source row count");
        if (policy == QMissingValuePolicy.EXCLUDE_ROW)
            selectedRows.andNot(scan.rowsWithMissing());
        if (selectedRows.isEmpty())
            throw new IllegalArgumentException("Missing-value policy excludes every row of " + scan.matrixName());
        String tag = "missing-policy-v2;policy=" + policy.optionName() + ";rows="
            + selectedRows.cardinality() + ";selection-hash=" + java.util.Arrays.hashCode(selectedRows.toLongArray());
        if (source.metadata().cacheSignatureTag() != null)
            tag = source.metadata().cacheSignatureTag() + ";" + tag;
        metadata = new Metadata(source.metadata().path(), selectedRows.cardinality(),
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
        private QMatrixRowSource.Block current;
        private int currentRow;
        private long outputOffset;
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
            List<String> rowIds = new ArrayList<>(maximumRows);
            List<double[]> values = new ArrayList<>(maximumRows);
            long blockOffset = outputOffset;
            while (values.size() < maximumRows) {
                if (current == null || currentRow >= current.rowCount()) {
                    current = inner.readBlock(INNER_BLOCK_ROWS);
                    currentRow = 0;
                    if (current == null)
                        break;
                }
                long sourceRow = current.rowOffset() + currentRow;
                String rowId = current.rowIds()[currentRow];
                double[] row = current.values()[currentRow];
                currentRow++;
                if (!selectedRows.get((int) sourceRow))
                    continue;
                applyPolicy(row, sourceRow, rowId);
                rowIds.add(rowId);
                values.add(row);
                outputOffset++;
            }
            if (values.isEmpty())
                return null;
            return new Block(blockOffset, rowIds.toArray(String[]::new), values.toArray(double[][]::new));
        }

        @Override
        public void close() throws IOException {
            if (!closed) {
                closed = true;
                inner.close();
            }
        }
    }

    private void applyPolicy(double[] values, long sourceRow, String rowId) {
        boolean missing = false;
        for (double value : values)
            missing |= QMissingnessScan.isMissing(value);
        if (!missing)
            return;
        if (policy == QMissingValuePolicy.ERROR || policy == QMissingValuePolicy.PATTERN
            || policy == QMissingValuePolicy.EXCLUDE_ROW)
            throw new IllegalArgumentException(scan.matrixName() + " row '" + rowId
                + "' still contains a missing value under policy " + policy.optionName());
        double replacement = policy == QMissingValuePolicy.ZERO ? 0.0 : observedMean(values);
        if (!Double.isFinite(replacement))
            throw new IllegalArgumentException(scan.matrixName() + " row '" + rowId
                + "' has no observed value available for " + policy.optionName() + " imputation");
        for (int i = 0; i < values.length; i++)
            if (QMissingnessScan.isMissing(values[i]))
                values[i] = replacement;
    }

    private static double observedMean(double[] values) {
        double sum = 0;
        int observed = 0;
        for (double value : values) {
            if (!QMissingnessScan.isMissing(value)) {
                sum += value;
                observed++;
            }
        }
        return observed == 0 ? Double.NaN : sum / observed;
    }

    private static BitSet allRows(long rowCount) {
        if (rowCount <= 0 || rowCount > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Unsupported row count " + rowCount);
        BitSet rows = new BitSet((int) rowCount);
        rows.set(0, (int) rowCount);
        return rows;
    }
}
