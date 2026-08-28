/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

import com.csvreader.CsvReader;

/** A small mixed-type covariate table with explicit sample alignment. */
public final class QCovariateTable {
    public record ModelMatrix(double[][] values, String[] columnNames, String[] automaticFactors) {
        public ModelMatrix {
            columnNames = columnNames.clone();
            automaticFactors = automaticFactors.clone();
        }

        @Override
        public String[] columnNames() {
            return columnNames.clone();
        }

        @Override
        public String[] automaticFactors() {
            return automaticFactors.clone();
        }
    }

    public record Missingness(boolean[] completeRows, int[] rowMissingCounts,
        Map<String, Integer> columnMissingCounts, int totalMissing) {
        public Missingness {
            completeRows = completeRows.clone();
            rowMissingCounts = rowMissingCounts.clone();
            columnMissingCounts = Collections.unmodifiableMap(new LinkedHashMap<>(columnMissingCounts));
        }

        @Override
        public boolean[] completeRows() { return completeRows.clone(); }

        @Override
        public int[] rowMissingCounts() { return rowMissingCounts.clone(); }

        public int completeSampleCount() {
            int count = 0;
            for (boolean complete : completeRows)
                if (complete) count++;
            return count;
        }
    }

    private record EncodedColumns(double[][] columns, String[] names, boolean automaticFactor) { }
    private record NormalizedIds(String[] values, int transformedCount) { }

    private final Path path;
    private final String[] columnNames;
    private final String[][] rows;
    private final Map<String, Integer> columnIndices;

    private QCovariateTable(Path path, String[] columnNames, String[][] rows) {
        this.path = path;
        this.columnNames = columnNames;
        this.rows = rows;
        columnIndices = new LinkedHashMap<>();
        for (int i = 0; i < columnNames.length; i++)
            columnIndices.put(columnNames[i], i);
    }

    public static QCovariateTable load(Path path, char[] delimiters, String commentMarkers) throws IOException {
        CsvReader reader = null;
        try {
            reader = new CsvReader(QDelimitedMatrixSource.openReader(path));
            reader.setTrimWhitespace(true);
            reader.setDelimiter(delimiters);
            if (commentMarkers != null && !commentMarkers.isEmpty()) {
                reader.setUseComments(true);
                reader.setComment(commentMarkers.charAt(0));
            }
            if (!reader.readRecord())
                throw new IOException("Covariate file is empty: " + path);
            String[] header = reader.getValues();
            validateHeader(header, path);
            List<String[]> rows = new ArrayList<>();
            while (reader.readRecord()) {
                String[] row = reader.getValues();
                if (row.length != header.length)
                    throw new IOException("Expected " + header.length + " covariate fields but found "
                        + row.length + " at data row " + (rows.size() + 1) + " in " + path);
                rows.add(row);
            }
            if (rows.isEmpty())
                throw new IOException("Covariate file contains no data rows: " + path);
            return new QCovariateTable(path.toAbsolutePath().normalize(), header.clone(), rows.toArray(String[][]::new));
        } finally {
            if (reader != null)
                reader.close();
        }
    }

    public int rowCount() {
        return rows.length;
    }

    public String[] columnNames() {
        return columnNames.clone();
    }

    public String[] columnValues(String name, boolean[] includedRows) {
        requireColumn(name);
        if (includedRows != null && includedRows.length != rows.length)
            throw new IllegalArgumentException("Covariate row selection has the wrong length");
        int column = columnIndices.get(name);
        List<String> values = new ArrayList<>();
        for (int row = 0; row < rows.length; row++)
            if (includedRows == null || includedRows[row]) {
                String value = rows[row][column].trim();
                requirePresent(value, name, row);
                values.add(value);
            }
        return values.toArray(String[]::new);
    }

    public QSampleAlignment align(String[] genotypeIds, String[] expressionIds,
        String requestedGenotypeIdColumn, String requestedExpressionIdColumn) {
        return align(genotypeIds, expressionIds, requestedGenotypeIdColumn,
            requestedExpressionIdColumn, QSampleAlignmentPolicy.STRICT, null, null);
    }

    public QSampleAlignment align(String[] genotypeIds, String[] expressionIds,
        String requestedGenotypeIdColumn, String requestedExpressionIdColumn,
        QSampleAlignmentPolicy policy, String genotypeIdStripPrefix,
        String expressionIdStripPrefix) {
        if (policy == null)
            throw new IllegalArgumentException("Sample alignment policy is required");
        NormalizedIds genotype = normalizeIds(genotypeIds, genotypeIdStripPrefix, "genotype sample");
        NormalizedIds expression = normalizeIds(expressionIds, expressionIdStripPrefix, "expression sample");
        if (policy == QSampleAlignmentPolicy.STRICT
            && (genotypeIds.length != rows.length || expressionIds.length != rows.length))
            throw new IllegalArgumentException("Sample counts do not match under strict alignment: genotype="
                + genotypeIds.length + ", expression=" + expressionIds.length + ", covariates=" + rows.length
                + "; use --sample-alignment covariate-subset only when matrix-only samples should be excluded");
        if (genotypeIds.length < rows.length || expressionIds.length < rows.length)
            throw new IllegalArgumentException("A matrix has fewer samples than the covariate table: genotype="
                + genotypeIds.length + ", expression=" + expressionIds.length + ", covariates=" + rows.length);

        int genotypeIdIndex = resolveIdColumn(genotype.values(), requestedGenotypeIdColumn,
            "genotype", policy);
        int expressionIdIndex = resolveIdColumn(expression.values(), requestedExpressionIdColumn,
            "expression", policy);
        int[] genotypeOrder = orderFor(genotype.values(), genotypeIdIndex, "genotype");
        int[] expressionOrder = orderFor(expression.values(), expressionIdIndex, "expression");
        return new QSampleAlignment(genotypeOrder, expressionOrder,
            columnNames[genotypeIdIndex], columnNames[expressionIdIndex],
            reorderedCount(genotypeOrder), reorderedCount(expressionOrder), policy,
            normalizedPrefix(genotypeIdStripPrefix), normalizedPrefix(expressionIdStripPrefix),
            genotype.transformedCount(), expression.transformedCount(),
            genotypeIds.length - rows.length, expressionIds.length - rows.length);
    }

    public static QSampleAlignment alignDirectly(String[] genotypeIds, String[] expressionIds) {
        return alignDirectly(genotypeIds, expressionIds, null, null);
    }

    public static QSampleAlignment alignDirectly(String[] genotypeIds, String[] expressionIds,
        String genotypeIdStripPrefix, String expressionIdStripPrefix) {
        NormalizedIds genotype = normalizeIds(genotypeIds, genotypeIdStripPrefix, "genotype sample");
        NormalizedIds expression = normalizeIds(expressionIds, expressionIdStripPrefix, "expression sample");
        if (genotypeIds.length != expressionIds.length)
            throw new IllegalArgumentException("Genotype and expression sample counts differ");
        Map<String, Integer> expressionIndex = index(expression.values(), "expression sample");
        int[] identity = new int[genotypeIds.length];
        int[] expressionOrder = new int[genotypeIds.length];
        for (int i = 0; i < genotypeIds.length; i++) {
            identity[i] = i;
            Integer source = expressionIndex.get(genotype.values()[i]);
            if (source == null)
                throw new IllegalArgumentException("Expression data are missing genotype sample '"
                    + genotype.values()[i] + "'");
            expressionOrder[i] = source;
        }
        return new QSampleAlignment(identity, expressionOrder, "matrix header", "matrix header",
            0, reorderedCount(expressionOrder), QSampleAlignmentPolicy.STRICT,
            normalizedPrefix(genotypeIdStripPrefix), normalizedPrefix(expressionIdStripPrefix),
            genotype.transformedCount(), expression.transformedCount(), 0, 0);
    }

    public ModelMatrix buildModelMatrix(String[] terms, String[] forcedFactors) {
        Set<String> forced = forcedFactors == null ? Set.of() : new HashSet<>(Arrays.asList(forcedFactors));
        if (forcedFactors != null)
            for (String factor : forcedFactors)
                requireColumn(factor);

        Map<String, EncodedColumns> encoded = new HashMap<>();
        List<double[]> modelColumns = new ArrayList<>();
        List<String> modelNames = new ArrayList<>();
        List<String> automaticFactors = new ArrayList<>();
        double[] intercept = new double[rows.length];
        Arrays.fill(intercept, 1.0);
        modelColumns.add(intercept);
        modelNames.add("(Intercept)");

        if (terms != null) {
            for (String term : terms) {
                if (term == null || term.isBlank())
                    continue;
                String[] factors = term.split("\\*");
                EncodedColumns product = null;
                for (String factor : factors) {
                    String name = factor.trim();
                    requireColumn(name);
                    EncodedColumns current = encoded.get(name);
                    if (current == null) {
                        current = encode(name, forced.contains(name));
                        encoded.put(name, current);
                        if (current.automaticFactor())
                            automaticFactors.add(name);
                    }
                    product = product == null ? current : interact(product, current);
                }
                if (product == null || product.columns().length == 0)
                    throw new IllegalArgumentException("Covariate term '" + term + "' has no estimable columns");
                modelColumns.addAll(Arrays.asList(product.columns()));
                modelNames.addAll(Arrays.asList(product.names()));
            }
        }

        double[][] values = new double[rows.length][modelColumns.size()];
        for (int column = 0; column < modelColumns.size(); column++)
            for (int row = 0; row < rows.length; row++)
                values[row][column] = modelColumns.get(column)[row];
        return new ModelMatrix(values, modelNames.toArray(String[]::new), automaticFactors.toArray(String[]::new));
    }

    public ModelMatrix buildModelMatrix(String[] terms, String[] forcedFactors, boolean[] includedRows) {
        if (includedRows == null)
            return buildModelMatrix(terms, forcedFactors);
        if (includedRows.length != rows.length)
            throw new IllegalArgumentException("Covariate row selection has the wrong length");
        List<String[]> selected = new ArrayList<>();
        for (int row = 0; row < rows.length; row++)
            if (includedRows[row])
                selected.add(rows[row]);
        if (selected.isEmpty())
            throw new IllegalArgumentException("Covariate missingness policy removed every sample");
        return new QCovariateTable(path, columnNames, selected.toArray(String[][]::new))
            .buildModelMatrix(terms, forcedFactors);
    }

    public Missingness inspectMissingness(String[] terms) {
        LinkedHashMap<String, Integer> selectedColumns = new LinkedHashMap<>();
        if (terms != null) {
            for (String term : terms) {
                if (term == null || term.isBlank())
                    continue;
                for (String raw : term.split("\\*")) {
                    String name = raw.trim();
                    requireColumn(name);
                    selectedColumns.putIfAbsent(name, 0);
                }
            }
        }
        boolean[] complete = new boolean[rows.length];
        Arrays.fill(complete, true);
        int[] rowMissing = new int[rows.length];
        int total = 0;
        for (String name : new ArrayList<>(selectedColumns.keySet())) {
            int column = columnIndices.get(name);
            int columnMissing = 0;
            for (int row = 0; row < rows.length; row++) {
                if (isMissingToken(rows[row][column])) {
                    complete[row] = false;
                    rowMissing[row]++;
                    columnMissing++;
                    total++;
                }
            }
            selectedColumns.put(name, columnMissing);
        }
        return new Missingness(complete, rowMissing, selectedColumns, total);
    }

    private EncodedColumns encode(String name, boolean forceFactor) {
        int index = columnIndices.get(name);
        boolean numeric = !forceFactor;
        double[] numericValues = new double[rows.length];
        if (numeric) {
            for (int row = 0; row < rows.length; row++) {
                String raw = rows[row][index].trim();
                requirePresent(raw, name, row);
                try {
                    numericValues[row] = Double.parseDouble(raw);
                    if (!Double.isFinite(numericValues[row]))
                        throw new NumberFormatException("not finite");
                } catch (NumberFormatException e) {
                    numeric = false;
                    break;
                }
            }
        }
        if (numeric)
            return new EncodedColumns(new double[][] { numericValues }, new String[] { name }, false);

        TreeSet<String> levels = new TreeSet<>();
        for (int row = 0; row < rows.length; row++) {
            String raw = rows[row][index].trim();
            requirePresent(raw, name, row);
            levels.add(raw);
        }
        if (levels.size() < 2)
            throw new IllegalArgumentException("Categorical covariate '" + name + "' has fewer than two levels");
        String reference = levels.first();
        List<String> nonReference = new ArrayList<>(levels);
        nonReference.remove(0);
        double[][] columns = new double[nonReference.size()][rows.length];
        String[] names = new String[nonReference.size()];
        for (int level = 0; level < nonReference.size(); level++) {
            String value = nonReference.get(level);
            names[level] = name + "[" + value + "]";
            for (int row = 0; row < rows.length; row++)
                columns[level][row] = value.equals(rows[row][index].trim()) ? 1.0 : 0.0;
        }
        if (reference.isEmpty())
            throw new IllegalArgumentException("Blank reference level for covariate '" + name + "'");
        return new EncodedColumns(columns, names, !forceFactor);
    }

    private static EncodedColumns interact(EncodedColumns left, EncodedColumns right) {
        int count = left.columns().length * right.columns().length;
        double[][] columns = new double[count][];
        String[] names = new String[count];
        int output = 0;
        for (int i = 0; i < left.columns().length; i++) {
            for (int j = 0; j < right.columns().length; j++) {
                double[] values = new double[left.columns()[i].length];
                for (int row = 0; row < values.length; row++)
                    values[row] = left.columns()[i][row] * right.columns()[j][row];
                columns[output] = values;
                names[output] = left.names()[i] + "*" + right.names()[j];
                output++;
            }
        }
        return new EncodedColumns(columns, names, left.automaticFactor() || right.automaticFactor());
    }

    private int resolveIdColumn(String[] matrixIds, String requested, String matrixName,
        QSampleAlignmentPolicy policy) {
        if (requested != null) {
            requireColumn(requested);
            if (!columnMatches(requested, matrixIds, policy))
                throw new IllegalArgumentException("Covariate column '" + requested + "' does not match "
                    + matrixName + " sample IDs under " + policy.optionName() + " alignment");
            return columnIndices.get(requested);
        }
        List<Integer> matches = new ArrayList<>();
        for (int column = 0; column < columnNames.length; column++)
            if (columnMatches(columnNames[column], matrixIds, policy))
                matches.add(column);
        if (matches.isEmpty())
            throw new IllegalArgumentException("No covariate column matches " + matrixName
                + " sample IDs; specify --" + matrixName + "-id-column");
        if (matches.size() > 1)
            throw new IllegalArgumentException("More than one covariate column matches " + matrixName
                + " sample IDs (" + matchingNames(matches) + "); specify --" + matrixName + "-id-column");
        return matches.get(0);
    }

    private boolean columnMatches(String columnName, String[] identifiers,
        QSampleAlignmentPolicy policy) {
        int column = columnIndices.get(columnName);
        Set<String> expected = new HashSet<>(Arrays.asList(identifiers));
        Set<String> actual = new HashSet<>();
        for (String[] row : rows) {
            String value = row[column].trim();
            if (value.isEmpty() || !actual.add(value))
                return false;
        }
        return policy == QSampleAlignmentPolicy.STRICT
            ? actual.equals(expected) : expected.containsAll(actual);
    }

    private int[] orderFor(String[] matrixIds, int covariateIdColumn, String matrixName) {
        Map<String, Integer> matrixIndex = index(matrixIds, matrixName + " sample");
        int[] order = new int[rows.length];
        Set<String> covariateIds = new HashSet<>();
        for (int row = 0; row < rows.length; row++) {
            String id = rows[row][covariateIdColumn].trim();
            if (id.isEmpty() || !covariateIds.add(id))
                throw new IllegalArgumentException("Blank or duplicate ID in covariate column '"
                    + columnNames[covariateIdColumn] + "'");
            Integer source = matrixIndex.get(id);
            if (source == null)
                throw new IllegalArgumentException("Covariate " + matrixName + " ID '" + id
                    + "' is absent from the matrix header");
            order[row] = source;
        }
        return order;
    }

    private void requireColumn(String name) {
        if (!columnIndices.containsKey(name))
            throw new IllegalArgumentException("Covariate '" + name + "' is not present in " + path);
    }

    private static void requirePresent(String value, String name, int row) {
        if (isMissingToken(value))
            throw new IllegalArgumentException("Missing value for covariate '" + name + "' at data row " + (row + 1));
    }

    public static boolean isMissingToken(String value) {
        String normalized = value == null ? "" : value.trim().toLowerCase(Locale.ROOT);
        return normalized.isEmpty() || normalized.equals("na") || normalized.equals("n/a")
            || normalized.equals("null") || normalized.equals("nan") || normalized.equals(".");
    }

    private static void validateHeader(String[] header, Path path) throws IOException {
        Set<String> seen = new HashSet<>();
        for (String raw : header) {
            String value = raw.trim();
            if (value.isEmpty())
                throw new IOException("Blank covariate column name in " + path);
            if (!seen.add(value))
                throw new IOException("Duplicate covariate column name '" + value + "' in " + path);
        }
    }

    private static void validateUnique(String[] identifiers, String kind) {
        if (identifiers == null)
            throw new IllegalArgumentException("Missing " + kind + " identifiers");
        index(identifiers, kind);
    }

    private static NormalizedIds normalizeIds(String[] identifiers, String prefix, String kind) {
        if (identifiers == null)
            throw new IllegalArgumentException("Missing " + kind + " identifiers");
        String literalPrefix = normalizedPrefix(prefix);
        String[] normalized = new String[identifiers.length];
        int transformed = 0;
        for (int i = 0; i < identifiers.length; i++) {
            String id = identifiers[i] == null ? "" : identifiers[i].trim();
            if (!literalPrefix.isEmpty() && id.startsWith(literalPrefix)) {
                id = id.substring(literalPrefix.length());
                transformed++;
            }
            normalized[i] = id;
        }
        index(normalized, kind + " after prefix normalization");
        return new NormalizedIds(normalized, transformed);
    }

    private static String normalizedPrefix(String prefix) {
        return prefix == null ? "" : prefix;
    }

    private static Map<String, Integer> index(String[] identifiers, String kind) {
        Map<String, Integer> result = new HashMap<>();
        for (int i = 0; i < identifiers.length; i++) {
            String id = identifiers[i] == null ? "" : identifiers[i].trim();
            if (id.isEmpty())
                throw new IllegalArgumentException("Blank " + kind + " identifier");
            if (result.put(id, i) != null)
                throw new IllegalArgumentException("Duplicate " + kind + " identifier '" + id + "'");
        }
        return result;
    }

    private String matchingNames(List<Integer> matches) {
        List<String> names = new ArrayList<>();
        for (int match : matches)
            names.add(columnNames[match]);
        return String.join(", ", names);
    }

    private static int reorderedCount(int[] order) {
        int count = 0;
        for (int i = 0; i < order.length; i++)
            if (order[i] != i)
                count++;
        return count;
    }
}
