/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.jama.CholeskyDecomposition;
import gov.nih.jama.Matrix;
import gov.nih.utils.matrix.EMultiplicationMode;
import gov.nih.utils.matrix.QMatrixUtils;

/** Cohort-specific fixed effects and optional within-cohort covariance whitening. */
final class QCohortModel {
    private static final String FORMAT = "gpu-eqtl-cohort-model-v1";

    record Definition(String cohort, String[] fixedCovariates, String[] factorCovariates,
        String subjectColumn, String repeatedPolicy, double withinSubjectCorrelation,
        Path covarianceFile) { }

    static final class Definitions {
        private final Path path;
        private final LinkedHashMap<String, Definition> values;

        private Definitions(Path path, LinkedHashMap<String, Definition> values) {
            this.path = path;
            this.values = values;
        }

        static Definitions load(Path path) throws IOException {
            Path normalized = path.toAbsolutePath().normalize();
            List<String> lines = Files.readAllLines(normalized, StandardCharsets.UTF_8);
            if (lines.isEmpty()) throw new IOException("Cohort model file is empty: " + normalized);
            String[] header = lines.get(0).split("\t", -1);
            Map<String, Integer> columns = new LinkedHashMap<>();
            for (int i = 0; i < header.length; i++) {
                String name = header[i].trim().toUpperCase(Locale.ROOT);
                if (name.isEmpty() || columns.put(name, i) != null)
                    throw new IOException("Blank or duplicate cohort-model column in " + normalized);
            }
            require(columns, "COHORT", normalized);
            require(columns, "FIXED_COVARIATES", normalized);
            LinkedHashMap<String, Definition> definitions = new LinkedHashMap<>();
            for (int line = 1; line < lines.size(); line++) {
                if (lines.get(line).isBlank() || lines.get(line).stripLeading().startsWith("#"))
                    continue;
                String[] fields = lines.get(line).split("\t", -1);
                if (fields.length != header.length)
                    throw new IOException("Expected " + header.length + " cohort-model fields at line "
                        + (line + 1) + " in " + normalized);
                String cohort = field(fields, columns, "COHORT");
                if (cohort.isEmpty() || definitions.containsKey(cohort))
                    throw new IOException("Blank or duplicate cohort at line " + (line + 1));
                String subject = field(fields, columns, "SUBJECT_COLUMN");
                String policy = field(fields, columns, "REPEATED_POLICY");
                if (policy.isEmpty()) policy = "independent";
                policy = policy.toLowerCase(Locale.ROOT);
                if (!(policy.equals("independent") || policy.equals("compound-symmetry")
                    || policy.equals("one-per-subject")))
                    throw new IOException("Unsupported repeated policy '" + policy + "' for " + cohort);
                String correlationText = field(fields, columns, "WITHIN_SUBJECT_CORRELATION");
                double correlation = correlationText.isEmpty() ? 0 : Double.parseDouble(correlationText);
                if (!Double.isFinite(correlation) || correlation < 0 || correlation >= 1)
                    throw new IOException("Within-subject correlation must be in [0,1) for " + cohort);
                if ((policy.equals("compound-symmetry") || policy.equals("one-per-subject"))
                    && subject.isEmpty())
                    throw new IOException("Repeated policy " + policy
                        + " requires SUBJECT_COLUMN for " + cohort);
                String covariance = field(fields, columns, "COVARIANCE_FILE");
                Path covariancePath = covariance.isEmpty() ? null
                    : normalized.getParent().resolve(covariance).normalize();
                if (covariancePath != null && policy.equals("compound-symmetry"))
                    throw new IOException("Specify COVARIANCE_FILE or compound-symmetry, not both, for " + cohort);
                definitions.put(cohort, new Definition(cohort,
                    list(field(fields, columns, "FIXED_COVARIATES")),
                    list(field(fields, columns, "FACTOR_COVARIATES")),
                    subject.isEmpty() ? null : subject, policy, correlation, covariancePath));
            }
            if (definitions.isEmpty()) throw new IOException("Cohort model contains no definitions");
            return new Definitions(normalized, definitions);
        }

        String[] missingnessTerms(String cohortColumn, String[] commonTerms) {
            LinkedHashMap<String, Boolean> terms = new LinkedHashMap<>();
            terms.put(cohortColumn, true);
            addTerms(terms, commonTerms);
            for (Definition definition : values.values()) {
                addTerms(terms, definition.fixedCovariates());
                if (definition.subjectColumn() != null) terms.put(definition.subjectColumn(), true);
            }
            return terms.keySet().toArray(String[]::new);
        }

        private static void addTerms(Map<String, Boolean> terms, String[] values) {
            if (values == null) return;
            for (String term : values)
                for (String component : term.split("\\*"))
                    if (!component.isBlank()) terms.put(component.trim(), true);
        }
    }

    private record CohortBlock(String name, int[] samples, double[][] lower,
        int designColumns, String repeatedPolicy, Path covarianceFile) { }

    private final double[][] design;
    private final String[] columnNames;
    private final CohortBlock[] blocks;
    private final String signature;

    static QCohortModel create(Definitions definitions, QCovariateTable covariates,
        String cohortColumn, boolean[] retainedRows, String[] canonicalSampleIds,
        String[] commonFixed, String[] commonFactors) throws IOException {
        if (definitions == null || covariates == null || cohortColumn == null
            || cohortColumn.isBlank() || retainedRows == null || canonicalSampleIds == null)
            throw new IllegalArgumentException("Cohort definitions, covariates, selection, and IDs are required");
        String[] cohorts = covariates.columnValues(cohortColumn, retainedRows);
        if (cohorts.length != canonicalSampleIds.length)
            throw new IllegalArgumentException("Cohort labels do not match aligned sample count");
        LinkedHashMap<String, List<Integer>> groups = new LinkedHashMap<>();
        for (int sample = 0; sample < cohorts.length; sample++)
            groups.computeIfAbsent(cohorts[sample], ignored -> new ArrayList<>()).add(sample);
        if (!groups.keySet().equals(definitions.values.keySet()))
            throw new IllegalArgumentException("Cohort-model definitions " + definitions.values.keySet()
                + " do not exactly match aligned cohort levels " + groups.keySet());

        List<QCovariateTable.ModelMatrix> localModels = new ArrayList<>();
        List<int[]> localSamples = new ArrayList<>();
        List<CohortBlock> cohortBlocks = new ArrayList<>();
        int totalColumns = 0;
        for (Definition definition : definitions.values.values()) {
            boolean[] selected = new boolean[retainedRows.length];
            int retainedIndex = 0;
            for (int row = 0; row < retainedRows.length; row++) {
                if (!retainedRows[row]) continue;
                if (cohorts[retainedIndex].equals(definition.cohort())) selected[row] = true;
                retainedIndex++;
            }
            String[] fixed = combine(commonFixed, definition.fixedCovariates());
            String[] factors = combine(commonFactors, definition.factorCovariates());
            QCovariateTable.ModelMatrix model = covariates.buildModelMatrix(fixed, factors, selected);
            int[] samples = groups.get(definition.cohort()).stream().mapToInt(Integer::intValue).toArray();
            if (model.values().length != samples.length)
                throw new IllegalArgumentException("Cohort model row count mismatch for " + definition.cohort());
            String[] subjects = definition.subjectColumn() == null ? null
                : selectedValues(covariates, definition.subjectColumn(), selected);
            if (definition.repeatedPolicy().equals("one-per-subject")) requireUniqueSubjects(
                definition.cohort(), subjects);
            double[][] covariance = null;
            if (definition.covarianceFile() != null)
                covariance = readCovariance(definition.covarianceFile(), select(canonicalSampleIds, samples));
            else if (definition.repeatedPolicy().equals("compound-symmetry"))
                covariance = compoundSymmetry(subjects, definition.withinSubjectCorrelation());
            double[][] lower = covariance == null ? null : cholesky(definition.cohort(), covariance);
            localModels.add(model);
            localSamples.add(samples);
            cohortBlocks.add(new CohortBlock(definition.cohort(), samples, lower,
                model.columnNames().length, definition.repeatedPolicy(), definition.covarianceFile()));
            totalColumns = Math.addExact(totalColumns, model.columnNames().length);
        }

        double[][] rawDesign = new double[cohorts.length][totalColumns];
        String[] names = new String[totalColumns];
        int columnOffset = 0;
        int cohortIndex = 0;
        for (Definition definition : definitions.values.values()) {
            QCovariateTable.ModelMatrix model = localModels.get(cohortIndex);
            int[] samples = localSamples.get(cohortIndex);
            for (int column = 0; column < model.columnNames().length; column++) {
                names[columnOffset + column] = definition.cohort() + ":" + model.columnNames()[column];
                for (int row = 0; row < samples.length; row++)
                    rawDesign[samples[row]][columnOffset + column] = model.values()[row][column];
            }
            columnOffset += model.columnNames().length;
            cohortIndex++;
        }
        CohortBlock[] blocks = cohortBlocks.toArray(CohortBlock[]::new);
        double[][] whitenedDesign = transformDesign(rawDesign, blocks);
        String signature = signature(definitions.path, cohortColumn, canonicalSampleIds,
            whitenedDesign, blocks);
        return new QCohortModel(whitenedDesign, names, blocks, signature);
    }

    private QCohortModel(double[][] design, String[] columnNames, CohortBlock[] blocks,
        String signature) {
        this.design = design;
        this.columnNames = columnNames;
        this.blocks = blocks;
        this.signature = signature;
    }

    double[][] design() { return design; }
    String[] columnNames() { return columnNames.clone(); }
    String signature() { return signature; }

    boolean hasWhitening() {
        for (CohortBlock block : blocks)
            if (block.lower() != null) return true;
        return false;
    }

    QeQTLPreprocessor.Residualizer wrap(QeQTLPreprocessor.Residualizer delegate) {
        return new QeQTLPreprocessor.Residualizer() {
            @Override
            public double[][] residualize(double[][] values, double[][] covariateQ,
                String matrixName) {
                double[][] transformed = transformRows(values, blocks);
                return delegate == null
                    ? QMatrixUtils.parallelMatrixMultiplication(transformed, covariateQ, null, 1,
                        transformed.length, transformed[0].length, EMultiplicationMode.XMinusXYYt)
                    : delegate.residualize(transformed, covariateQ, matrixName);
            }

            @Override public int concurrency() { return delegate == null ? 1 : delegate.concurrency(); }
            @Override public int estimatedHostBytesPerValue() {
                return Math.addExact(Double.BYTES,
                    delegate == null ? 0 : delegate.estimatedHostBytesPerValue());
            }
            @Override public String cacheSignatureTag() {
                return "cohort-whitened-fp64-v1;" + signature + ";delegate="
                    + (delegate == null ? "cpu-fp64-v1" : delegate.cacheSignatureTag());
            }
        };
    }

    void writeAudit(Path output) throws IOException {
        Path normalized = output.toAbsolutePath().normalize();
        if (normalized.getParent() != null) Files.createDirectories(normalized.getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(normalized, StandardCharsets.UTF_8)) {
            writer.write("cohort\tsamples\tdesign_columns\trepeated_policy\tcovariance\tsignature");
            writer.newLine();
            for (CohortBlock block : blocks) {
                writer.write(block.name()); writer.write('\t');
                writer.write(Integer.toString(block.samples().length)); writer.write('\t');
                writer.write(Integer.toString(block.designColumns())); writer.write('\t');
                writer.write(block.repeatedPolicy()); writer.write('\t');
                writer.write(block.covarianceFile() == null
                    ? (block.lower() == null ? "identity" : "compound-symmetry")
                    : block.covarianceFile().toString());
                writer.write('\t'); writer.write(signature); writer.newLine();
            }
        }
    }

    static double stackedDot(double[] left, double[] right) {
        if (left.length != right.length) throw new IllegalArgumentException("Sample counts differ");
        double value = 0;
        for (int i = 0; i < left.length; i++) value += left[i] * right[i];
        return value;
    }

    static double blockwiseDot(double[] left, double[] right, int[][] cohortSamples) {
        double value = 0;
        for (int[] cohort : cohortSamples)
            for (int sample : cohort) value += left[sample] * right[sample];
        return value;
    }

    private static double[][] transformRows(double[][] values, CohortBlock[] blocks) {
        double[][] result = new double[values.length][];
        for (int row = 0; row < values.length; row++) {
            result[row] = values[row].clone();
            transformVector(result[row], blocks);
        }
        return result;
    }

    private static double[][] transformDesign(double[][] design, CohortBlock[] blocks) {
        double[][] result = new double[design.length][design[0].length];
        for (int column = 0; column < design[0].length; column++) {
            double[] values = new double[design.length];
            for (int row = 0; row < design.length; row++) values[row] = design[row][column];
            transformVector(values, blocks);
            for (int row = 0; row < design.length; row++) result[row][column] = values[row];
        }
        return result;
    }

    private static void transformVector(double[] values, CohortBlock[] blocks) {
        for (CohortBlock block : blocks) {
            if (block.lower() == null) continue;
            int[] samples = block.samples();
            double[][] lower = block.lower();
            double[] solved = new double[samples.length];
            for (int row = 0; row < samples.length; row++) {
                double value = values[samples[row]];
                for (int column = 0; column < row; column++)
                    value -= lower[row][column] * solved[column];
                solved[row] = value / lower[row][row];
            }
            for (int row = 0; row < samples.length; row++) values[samples[row]] = solved[row];
        }
    }

    private static double[][] cholesky(String cohort, double[][] covariance) {
        CholeskyDecomposition decomposition = new CholeskyDecomposition(new Matrix(covariance));
        if (!decomposition.isSPD())
            throw new IllegalArgumentException("Covariance is not symmetric positive definite for " + cohort);
        return decomposition.getL().getArray();
    }

    private static double[][] compoundSymmetry(String[] subjects, double correlation) {
        if (subjects == null) throw new IllegalArgumentException("Subject IDs are required");
        double[][] covariance = new double[subjects.length][subjects.length];
        for (int row = 0; row < subjects.length; row++)
            for (int column = 0; column < subjects.length; column++)
                covariance[row][column] = row == column ? 1
                    : subjects[row].equals(subjects[column]) ? correlation : 0;
        return covariance;
    }

    private static void requireUniqueSubjects(String cohort, String[] subjects) {
        Set<String> seen = new HashSet<>();
        for (String subject : subjects)
            if (!seen.add(subject))
                throw new IllegalArgumentException("Cohort " + cohort
                    + " uses one-per-subject but subject '" + subject + "' is repeated");
    }

    private static double[][] readCovariance(Path path, String[] expectedIds) throws IOException {
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        if (lines.size() < 2) throw new IOException("Covariance file is empty: " + path);
        String delimiter = lines.get(0).contains("\t") ? "\t" : ",";
        String[] header = lines.get(0).split(delimiter, -1);
        if (header.length != expectedIds.length + 1)
            throw new IOException("Covariance header size does not match cohort samples: " + path);
        Map<String, Integer> columns = new HashMap<>();
        for (int i = 1; i < header.length; i++)
            if (header[i].isBlank() || columns.put(header[i].trim(), i - 1) != null)
                throw new IOException("Blank or duplicate covariance column ID in " + path);
        Map<String, double[]> rows = new HashMap<>();
        for (int line = 1; line < lines.size(); line++) {
            String[] fields = lines.get(line).split(delimiter, -1);
            if (fields.length != header.length) throw new IOException("Ragged covariance row in " + path);
            double[] values = new double[expectedIds.length];
            for (int column = 1; column < fields.length; column++)
                values[column - 1] = Double.parseDouble(fields[column]);
            if (rows.put(fields[0].trim(), values) != null)
                throw new IOException("Duplicate covariance row ID in " + path);
        }
        double[][] covariance = new double[expectedIds.length][expectedIds.length];
        for (int row = 0; row < expectedIds.length; row++) {
            double[] source = rows.get(expectedIds[row]);
            if (source == null) throw new IOException("Covariance omits sample " + expectedIds[row]);
            for (int column = 0; column < expectedIds.length; column++) {
                Integer sourceColumn = columns.get(expectedIds[column]);
                if (sourceColumn == null) throw new IOException("Covariance omits sample " + expectedIds[column]);
                covariance[row][column] = source[sourceColumn];
            }
        }
        if (rows.size() != expectedIds.length)
            throw new IOException("Covariance contains rows outside the cohort sample set: " + path);
        for (int row = 0; row < covariance.length; row++) {
            if (!(covariance[row][row] > 0)) throw new IOException("Non-positive covariance diagonal");
            for (int column = 0; column < row; column++) {
                if (Math.abs(covariance[row][column] - covariance[column][row]) > 1e-12)
                    throw new IOException("Covariance is not symmetric within tolerance: " + path);
                double average = (covariance[row][column] + covariance[column][row]) / 2;
                covariance[row][column] = average;
                covariance[column][row] = average;
            }
        }
        return covariance;
    }

    private static String signature(Path definitions, String cohortColumn, String[] ids,
        double[][] design, CohortBlock[] blocks) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            update(digest, FORMAT); update(digest, definitions.toString());
            update(digest, cohortColumn);
            for (String id : ids) update(digest, id);
            for (double[] row : design)
                for (double value : row) update(digest, Double.toHexString(value));
            for (CohortBlock block : blocks) {
                update(digest, block.name()); update(digest, block.repeatedPolicy());
                if (block.lower() != null)
                    for (double[] row : block.lower())
                        for (double value : row) update(digest, Double.toHexString(value));
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException(e);
        }
    }

    private static void update(MessageDigest digest, String value) {
        digest.update(value.getBytes(StandardCharsets.UTF_8)); digest.update((byte) 0);
    }

    private static String[] selectedValues(QCovariateTable table, String column,
        boolean[] selected) {
        return table.columnValues(column, selected);
    }

    private static String[] select(String[] values, int[] indices) {
        String[] selected = new String[indices.length];
        for (int i = 0; i < indices.length; i++) selected[i] = values[indices[i]];
        return selected;
    }

    private static String[] combine(String[] first, String[] second) {
        List<String> values = new ArrayList<>();
        if (first != null) values.addAll(Arrays.asList(first));
        if (second != null) values.addAll(Arrays.asList(second));
        return values.toArray(String[]::new);
    }

    private static String[] list(String value) {
        if (value == null || value.isBlank()) return new String[0];
        return Arrays.stream(value.split(",")).map(String::trim)
            .filter(item -> !item.isEmpty()).toArray(String[]::new);
    }

    private static String field(String[] fields, Map<String, Integer> columns, String name) {
        Integer column = columns.get(name);
        return column == null ? "" : fields[column].trim();
    }

    private static void require(Map<String, Integer> columns, String name, Path path)
        throws IOException {
        if (!columns.containsKey(name)) throw new IOException("Cohort model " + path
            + " is missing required column " + name);
    }
}
