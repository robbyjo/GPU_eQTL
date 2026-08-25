/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import static gov.nih.utils.QStatsUtils.calcStdDevAndStandardize;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HexFormat;
import java.util.Locale;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QMissingnessScan;
import gov.nih.jama.QRDecomposition;

/** Covariate/rank/DF state shared by exact trait masks in genotype-outer scheduling. */
final class QTraitPatternModelSet {
    private record BuiltModel(Model model) { }

    static final class Model {
        final int id;
        final long traitRows;
        final int[] observed;
        final double[] designSums;
        final double[][] upperR;
        final int rank;
        final int errorDegreesOfFreedom;
        final double rSquaredThreshold;
        final String exclusionReason;

        Model(int id, long traitRows, int[] observed, double[] designSums,
            double[][] upperR, int rank,
            int errorDegreesOfFreedom, double rSquaredThreshold, String exclusionReason) {
            this.id = id;
            this.traitRows = traitRows;
            this.observed = observed;
            this.designSums = designSums;
            this.upperR = upperR;
            this.rank = rank;
            this.errorDegreesOfFreedom = errorDegreesOfFreedom;
            this.rSquaredThreshold = rSquaredThreshold;
            this.exclusionReason = exclusionReason;
        }

        boolean estimable() { return exclusionReason == null; }
        int reportedDegreesOfFreedom(int offset) { return errorDegreesOfFreedom - offset; }
    }

    private final QMissingnessScan scan;
    private final double[][] design;
    private final Model[] models;
    private final int[] sourceRowPattern;
    private final int[] preparedRowPatterns;
    private final long estimableTraitRows;
    private final long excludedTraitRows;
    private final String firstExclusion;
    private final String signature;

    static QTraitPatternModelSet create(QMissingnessScan scan, double[][] covariateModel,
        int degreesOfFreedomOffset, String thresholdType, double threshold,
        QUnestimableTraitPolicy policy) {
        int samples = scan.sampleCount();
        double[][] design = covariateModel;
        if (design == null) {
            design = new double[samples][1];
            for (double[] row : design) row[0] = 1;
        }
        if (design.length != samples || design[0].length == 0)
            throw new IllegalArgumentException("Trait-pattern covariate design has the wrong dimensions");
        int columns = design[0].length;
        for (double[] row : design)
            if (row.length != columns)
                throw new IllegalArgumentException("Trait-pattern covariate design is ragged");

        Model[] models = new Model[scan.patterns().size()];
        int[] rowPattern = new int[Math.toIntExact(scan.rowCount())];
        Arrays.fill(rowPattern, -1);
        long started = System.nanoTime();
        long lastReport = started;
        int workers = Math.min(models.length,
            Math.max(1, Runtime.getRuntime().availableProcessors() - 1));
        System.out.println("Trait-pattern model preflight started: " + models.length
            + " exact missingness pattern(s), " + workers
            + " CPU worker(s); progress approximately every 15 seconds.");
        ExecutorService executor = Executors.newFixedThreadPool(workers);
        CompletionService<BuiltModel> completed = new ExecutorCompletionService<>(executor);
        for (QMissingnessScan.Pattern pattern : scan.patterns()) {
            double[][] immutableDesign = design;
            completed.submit(() -> new BuiltModel(buildModel(pattern, immutableDesign, samples,
                columns, degreesOfFreedomOffset, thresholdType, threshold)));
        }
        try {
            for (int completedPatterns = 1; completedPatterns <= models.length; completedPatterns++) {
                BuiltModel built = completed.take().get();
                models[built.model().id] = built.model();
                long now = System.nanoTime();
                if (completedPatterns == models.length || now - lastReport >= 15_000_000_000L) {
                    System.out.println("Trait-pattern model preflight "
                        + (completedPatterns == models.length ? "complete: " : "progress: ")
                        + String.format(Locale.ROOT, "%,d/%,d (%.1f%%); elapsed=%.1fs",
                            completedPatterns, models.length,
                            100.0 * completedPatterns / models.length,
                            (now - started) / 1_000_000_000.0));
                    lastReport = now;
                }
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IllegalStateException("Interrupted during trait-pattern model preflight", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof RuntimeException runtime) throw runtime;
            if (cause instanceof Error error) throw error;
            throw new IllegalStateException("Trait-pattern model preflight failed", cause);
        } finally {
            executor.shutdownNow();
            try {
                if (!executor.awaitTermination(1, TimeUnit.MINUTES))
                    System.err.println("WARNING: Trait-pattern model workers did not stop within one minute.");
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }

        long estimableRows = 0;
        long excludedRows = 0;
        String firstFailure = null;
        for (QMissingnessScan.Pattern pattern : scan.patterns()) {
            Model model = models[pattern.id()];
            String reason = model.exclusionReason;
            if (reason == null)
                estimableRows += pattern.rowIndices().length;
            else {
                excludedRows += pattern.rowIndices().length;
                if (firstFailure == null)
                    firstFailure = "pattern " + pattern.id() + " (traits="
                        + pattern.rowIndices().length + ", N=" + model.observed.length
                        + ", rank=" + model.rank + ", reason=" + reason + ")";
            }
            for (long row : pattern.rowIndices())
                rowPattern[Math.toIntExact(row)] = pattern.id();
        }
        for (int row = 0; row < rowPattern.length; row++)
            if (rowPattern[row] < 0)
                throw new IllegalArgumentException("Trait pattern metadata omits row " + row);
        int[] preparedPatterns = new int[Math.toIntExact(estimableRows)];
        int prepared = 0;
        for (int pattern : rowPattern)
            if (models[pattern].estimable())
                preparedPatterns[prepared++] = pattern;
        String signature = signature(scan, design, models, degreesOfFreedomOffset,
            thresholdType, threshold, policy);
        return new QTraitPatternModelSet(scan, design, models, rowPattern,
            preparedPatterns, estimableRows, excludedRows, firstFailure, signature);
    }

    private static Model buildModel(QMissingnessScan.Pattern pattern, double[][] design,
        int samples, int columns, int degreesOfFreedomOffset, String thresholdType,
        double threshold) {
        int[] observed = complement(pattern.missingSamples(), samples);
        double[][] selected = selectRows(design, observed);
        double[] designSums = new double[columns];
        for (double[] row : selected)
            for (int column = 0; column < columns; column++)
                designSums[column] += row[column];
        QRDecomposition qr = new QRDecomposition(selected);
        int rank = qr.getRank();
        int errorDf = observed.length - rank - 1;
        String reason = null;
        if (!qr.isFullRank())
            reason = "rank-deficient-design";
        else if (errorDf - degreesOfFreedomOffset <= 0)
            reason = "non-positive-residual-df";
        double rsqThreshold = reason == null
            ? QeQTLAnalysis.rSquaredThreshold(thresholdType, threshold,
                errorDf, degreesOfFreedomOffset) : Double.NaN;
        double[][] upperR = reason == null ? qr.getR().getArray() : null;
        return new Model(pattern.id(), pattern.rowIndices().length,
            observed, designSums, upperR, rank, errorDf, rsqThreshold, reason);
    }

    private QTraitPatternModelSet(QMissingnessScan scan, double[][] design, Model[] models,
        int[] sourceRowPattern, int[] preparedRowPatterns, long estimableTraitRows,
        long excludedTraitRows, String firstExclusion, String signature) {
        this.scan = scan;
        this.design = design;
        this.models = models;
        this.sourceRowPattern = sourceRowPattern;
        this.preparedRowPatterns = preparedRowPatterns;
        this.estimableTraitRows = estimableTraitRows;
        this.excludedTraitRows = excludedTraitRows;
        this.firstExclusion = firstExclusion;
        this.signature = signature;
    }

    Model model(int id) { return models[id]; }
    Model[] models() { return models; }
    int patternForSourceRow(long row) { return sourceRowPattern[Math.toIntExact(row)]; }
    int patternForPreparedRow(long row) { return preparedRowPatterns[Math.toIntExact(row)]; }
    int sampleCount() { return scan.sampleCount(); }
    int designColumns() { return design[0].length; }
    long estimableTraitRows() { return estimableTraitRows; }
    long excludedTraitRows() { return excludedTraitRows; }
    String firstExclusion() { return firstExclusion; }
    String signature() { return signature; }
    double designValue(int sample, int column) { return design[sample][column]; }

    PreparedBlock prepareTrait(long outputRowOffset, long sourceRow, String rowId,
        double[] rawValues) {
        Model model = model(patternForSourceRow(sourceRow));
        if (rawValues.length != sampleCount())
            throw new IllegalArgumentException("Trait-pattern preparation dimensions differ");
        double[] selected = new double[model.observed.length];
        for (int i = 0; i < model.observed.length; i++) {
            double value = rawValues[model.observed[i]];
            if (QMissingnessScan.isMissing(value))
                throw new IllegalArgumentException("Trait row '" + rowId
                    + "' is missing at an observed pattern sample");
            selected[i] = value;
        }
        double[] xTranspose = new double[model.rank];
        for (int sample = 0; sample < selected.length; sample++)
            for (int column = 0; column < xTranspose.length; column++)
                xTranspose[column] += design[model.observed[sample]][column] * selected[sample];
        double[] qTranspose = new double[model.rank];
        double[] coefficients = new double[model.rank];
        QPatternSufficientStatistics.leastSquaresCoefficients(xTranspose, model.upperR,
            qTranspose, coefficients);
        for (int sample = 0; sample < selected.length; sample++)
            for (int column = 0; column < coefficients.length; column++)
                selected[sample] -= design[model.observed[sample]][column] * coefficients[column];
        double standardDeviation = calcStdDevAndStandardize(selected);
        if (!(standardDeviation > 0) || !Double.isFinite(standardDeviation))
            throw new IllegalArgumentException("Trait row '" + rowId
                + "' has zero or invalid variance after pattern covariate adjustment");
        double[] padded = new double[sampleCount()];
        for (int i = 0; i < model.observed.length; i++)
            padded[model.observed[i]] = selected[i];
        return new PreparedBlock(outputRowOffset, new String[] {rowId},
            new double[][] {padded}, new double[] {standardDeviation});
    }

    void writeAudit(Path output, int degreesOfFreedomOffset, String scheduler) throws IOException {
        Path normalized = output.toAbsolutePath().normalize();
        Path parent = normalized.getParent();
        if (parent != null) Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, normalized.getFileName().toString(), ".partial");
        boolean complete = false;
        try (BufferedWriter writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {
            writer.write("pattern_id\ttrait_rows\tobserved_samples\tcovariate_rank\tresidual_df\tstatus\treason\tscheduler\n");
            for (Model model : models) {
                writer.write(Integer.toString(model.id)); writer.write('\t');
                writer.write(Long.toString(model.traitRows)); writer.write('\t');
                writer.write(Integer.toString(model.observed.length)); writer.write('\t');
                writer.write(Integer.toString(model.rank)); writer.write('\t');
                writer.write(Integer.toString(model.reportedDegreesOfFreedom(degreesOfFreedomOffset)));
                writer.write('\t'); writer.write(model.estimable() ? "included" : "excluded");
                writer.write('\t'); writer.write(model.exclusionReason == null ? "" : model.exclusionReason);
                writer.write('\t'); writer.write(scheduler); writer.newLine();
            }
            writer.flush();
            complete = true;
        } finally {
            if (complete) moveAtomically(temporary, normalized);
            else Files.deleteIfExists(temporary);
        }
    }

    private static String signature(QMissingnessScan scan, double[][] design, Model[] models,
        int dfOffset, String thresholdType, double threshold, QUnestimableTraitPolicy policy) {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
        update(digest, "gpu-eqtl-trait-pattern-model-v1");
        update(digest, scan.rowCount()); update(digest, scan.sampleCount());
        for (double[] row : design)
            for (double value : row) update(digest, Double.doubleToLongBits(value));
        for (Model model : models) {
            update(digest, model.id); update(digest, model.traitRows);
            for (int sample : model.observed) update(digest, sample);
            update(digest, model.rank); update(digest, model.errorDegreesOfFreedom);
            update(digest, model.exclusionReason);
        }
        update(digest, dfOffset); update(digest, thresholdType);
        update(digest, Double.doubleToLongBits(threshold)); update(digest, policy.optionName());
        return HexFormat.of().formatHex(digest.digest());
    }

    private static int[] complement(BitSet missing, int sampleCount) {
        int[] result = new int[sampleCount - missing.cardinality()];
        int output = 0;
        for (int sample = missing.nextClearBit(0); sample >= 0 && sample < sampleCount;
             sample = missing.nextClearBit(sample + 1))
            result[output++] = sample;
        return result;
    }

    private static double[][] selectRows(double[][] values, int[] rows) {
        double[][] selected = new double[rows.length][];
        for (int row = 0; row < rows.length; row++)
            selected[row] = values[rows[row]].clone();
        return selected;
    }

    private static void update(MessageDigest digest, String value) {
        if (value == null) { update(digest, -1); return; }
        byte[] bytes = value.getBytes(StandardCharsets.UTF_8);
        update(digest, bytes.length); digest.update(bytes);
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
}
