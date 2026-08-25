/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.List;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.eqtl.io.QPreparedMatrix;

/** Sequential builder source for globally padded, pattern-residualized trait rows. */
final class QPatternPreparedTraitSource implements QPreparedMatrix, AutoCloseable {
    private static final int INNER_ROWS = 256;

    private final QMatrixRowSource source;
    private final int[] columnOrder;
    private final QTraitPatternModelSet models;
    private final String signature;
    private QMatrixRowSource.BlockReader reader;
    private QMatrixRowSource.Block current;
    private int currentRow;
    private long sourceRowsRead;
    private long outputRowsRead;
    private boolean closed;

    QPatternPreparedTraitSource(QMatrixRowSource source, int[] columnOrder,
        QTraitPatternModelSet models) throws IOException {
        this.source = source;
        this.columnOrder = columnOrder.clone();
        this.models = models;
        signature = signature(source, columnOrder, models.signature());
    }

    @Override public String signature() { return signature; }
    @Override public int sampleCount() { return models.sampleCount(); }
    @Override public long rowCount() { return models.estimableTraitRows(); }

    @Override
    public PreparedBlock readBlock(long rowOffset, int maximumRows) throws IOException {
        if (rowOffset != outputRowsRead)
            throw new IllegalArgumentException("Pattern trait builder requires sequential row access");
        if (maximumRows <= 0)
            throw new IllegalArgumentException("maximumRows must be positive");
        if (closed)
            throw new IOException("Pattern trait builder is closed");
        if (reader == null)
            reader = source.open(columnOrder);
        List<String> ids = new ArrayList<>(maximumRows);
        List<double[]> values = new ArrayList<>(maximumRows);
        List<Double> standardDeviations = new ArrayList<>(maximumRows);
        while (values.size() < maximumRows) {
            if (current == null || currentRow >= current.rowCount()) {
                current = reader.readBlock(INNER_ROWS);
                currentRow = 0;
                if (current == null) break;
                if (current.rowOffset() != sourceRowsRead)
                    throw new IOException("Trait source row order changed during pattern preparation");
            }
            long sourceRow = current.rowOffset() + currentRow;
            String id = current.rowIds()[currentRow];
            double[] raw = current.values()[currentRow];
            currentRow++;
            sourceRowsRead++;
            QTraitPatternModelSet.Model model = models.model(models.patternForSourceRow(sourceRow));
            if (!model.estimable())
                continue;
            PreparedBlock prepared = models.prepareTrait(outputRowsRead, sourceRow, id, raw);
            ids.add(id);
            values.add(prepared.values()[0]);
            standardDeviations.add(prepared.standardDeviations()[0]);
            outputRowsRead++;
        }
        if (values.isEmpty()) {
            if (sourceRowsRead != source.metadata().rowCount()
                || outputRowsRead != models.estimableTraitRows())
                throw new IOException("Trait source dimensions changed during pattern preparation");
            close();
            return null;
        }
        double[] sd = new double[standardDeviations.size()];
        for (int i = 0; i < sd.length; i++) sd[i] = standardDeviations.get(i);
        return new PreparedBlock(rowOffset, ids.toArray(String[]::new),
            values.toArray(double[][]::new), sd);
    }

    @Override
    public void close() throws IOException {
        if (!closed) {
            closed = true;
            if (reader != null) reader.close();
        }
    }

    private static String signature(QMatrixRowSource source, int[] columns, String modelSignature)
        throws IOException {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
        update(digest, "gpu-eqtl-pattern-trait-prepared-v2");
        update(digest, gov.nih.eqtl.io.QBinaryMatrixCache.signature(
            "PatternTraitSource", source, columns, null));
        update(digest, modelSignature);
        return HexFormat.of().formatHex(digest.digest());
    }

    private static void update(MessageDigest digest, String value) {
        byte[] bytes = value.getBytes(java.nio.charset.StandardCharsets.UTF_8);
        digest.update((byte) (bytes.length >>> 24));
        digest.update((byte) (bytes.length >>> 16));
        digest.update((byte) (bytes.length >>> 8));
        digest.update((byte) bytes.length);
        digest.update(bytes);
    }
}
