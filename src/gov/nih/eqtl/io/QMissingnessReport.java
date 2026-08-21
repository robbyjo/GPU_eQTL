/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Map;

import gov.nih.eqtl.QDataType;
import gov.nih.eqtl.QMissingValuePolicy;

/** Atomic, tab-separated audit of matrix and selected-covariate missingness. */
public final class QMissingnessReport {
    private QMissingnessReport() { }

    public static void write(Path output, QMissingnessScan predictor, QDataType predictorType,
        QMissingValuePolicy predictorPolicy, QMissingnessScan trait, QDataType traitType,
        QMissingValuePolicy traitPolicy, QCovariateTable.Missingness covariates,
        String covariatePolicy, String[] canonicalSampleIds) throws IOException {
        Path normalized = output.toAbsolutePath().normalize();
        Path parent = normalized.getParent();
        if (parent != null)
            Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, normalized.getFileName().toString(), ".partial");
        boolean complete = false;
        try (BufferedWriter writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {
            writer.write("record_type\tmatrix\tentity_index\tentity_id\ttotal_values\tmissing_values\tmissing_rate\tpattern_id\tpattern_rows\tpolicy\tdetail\n");
            writeScan(writer, predictor, predictorType.optionName(), predictorPolicy.optionName());
            writeScan(writer, trait, traitType.optionName(), traitPolicy.optionName());
            if (covariates != null) {
                long total = (long) covariates.rowMissingCounts().length
                    * covariates.columnMissingCounts().size();
                summary(writer, "covariates", "selected", total, covariates.totalMissing(),
                    covariatePolicy, "selected_columns=" + covariates.columnMissingCounts().size()
                        + ";complete_samples=" + covariates.completeSampleCount());
                int columnIndex = 0;
                for (Map.Entry<String, Integer> entry : covariates.columnMissingCounts().entrySet()) {
                    if (entry.getValue() > 0)
                        record(writer, "COVARIATE_COLUMN", "covariates", columnIndex,
                            entry.getKey(), covariates.rowMissingCounts().length, entry.getValue(),
                            -1, -1, covariatePolicy, "");
                    columnIndex++;
                }
                int[] missing = covariates.rowMissingCounts();
                for (int sample = 0; sample < missing.length; sample++)
                    if (missing[sample] > 0)
                        record(writer, "COVARIATE_SAMPLE", "covariates", sample,
                            canonicalSampleIds == null ? Integer.toString(sample) : canonicalSampleIds[sample],
                            covariates.columnMissingCounts().size(), missing[sample], -1, -1,
                            covariatePolicy, "");
            }
            writer.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(temporary, normalized);
            else
                Files.deleteIfExists(temporary);
        }
    }

    private static void writeScan(BufferedWriter writer, QMissingnessScan scan,
        String dataType, String policy) throws IOException {
        summary(writer, scan.matrixName(), dataType, scan.totalValues(), scan.totalMissingValues(),
            policy, "rows=" + scan.rowCount() + ";samples=" + scan.sampleCount()
                + ";patterns=" + scan.patterns().size());
        for (QMissingnessScan.Pattern pattern : scan.patterns())
            record(writer, "PATTERN", scan.matrixName(), pattern.id(), "pattern-" + pattern.id(),
                scan.sampleCount(), pattern.missingSamples().cardinality(), pattern.id(),
                pattern.rowIndices().length, policy,
                "observed_samples=" + (scan.sampleCount() - pattern.missingSamples().cardinality()));
        for (QMissingnessScan.MissingRow row : scan.missingRows())
            record(writer, "ROW", scan.matrixName(), row.rowIndex(), row.rowId(),
                scan.sampleCount(), row.missingValues(), row.patternId(), -1, policy, "");
        long[] sampleMissing = scan.sampleMissingValues();
        String[] sampleIds = scan.sampleIds();
        for (int sample = 0; sample < sampleMissing.length; sample++)
            if (sampleMissing[sample] > 0)
                record(writer, "SAMPLE", scan.matrixName(), sample, sampleIds[sample],
                    scan.rowCount(), sampleMissing[sample], -1, -1, policy, "");
    }

    private static void summary(BufferedWriter writer, String matrix, String id,
        long total, long missing, String policy, String detail) throws IOException {
        record(writer, "SUMMARY", matrix, -1, id, total, missing, -1, -1, policy, detail);
    }

    private static void record(BufferedWriter writer, String type, String matrix, long index,
        String id, long total, long missing, int pattern, long patternRows,
        String policy, String detail) throws IOException {
        double rate = total == 0 ? Double.NaN : missing / (double) total;
        writer.write(tsv(type)); writer.write('\t');
        writer.write(tsv(matrix)); writer.write('\t');
        writer.write(index < 0 ? "" : Long.toString(index)); writer.write('\t');
        writer.write(tsv(id)); writer.write('\t');
        writer.write(Long.toString(total)); writer.write('\t');
        writer.write(Long.toString(missing)); writer.write('\t');
        writer.write(Double.isFinite(rate) ? Double.toString(rate) : "NA"); writer.write('\t');
        writer.write(pattern < 0 ? "" : Integer.toString(pattern)); writer.write('\t');
        writer.write(patternRows < 0 ? "" : Long.toString(patternRows)); writer.write('\t');
        writer.write(tsv(policy)); writer.write('\t');
        writer.write(tsv(detail)); writer.write('\n');
    }

    private static String tsv(String value) {
        if (value == null)
            return "";
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
