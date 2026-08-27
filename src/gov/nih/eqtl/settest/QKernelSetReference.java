/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import static gov.nih.utils.QDataUtils.kUndefinedValue;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import gov.nih.eqtl.QMissingValuePolicy;
import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.settest.QBurdenReference.SetAudit;
import gov.nih.eqtl.settest.QBurdenReference.Variant;
import gov.nih.eqtl.settest.QSetTestPolicy.FailurePolicy;
import gov.nih.eqtl.settest.QVariantSetTable.Entry;
import gov.nih.eqtl.settest.QVariantSetTable.SetDefinition;
import gov.nih.utils.matrix.EMultiplicationMode;
import gov.nih.utils.matrix.QMatrixUtils;

/** Scalar FP64 set preparation and kernel-test reference used by SKAT and SKAT-O. */
public final class QKernelSetReference {
    private static final double DOSAGE_TOLERANCE = 1e-12;

    public record Options(QSetTestMethod method, double[] rhoGrid,
        int simulations, long seed) {
        public Options {
            if (method != QSetTestMethod.SKAT && method != QSetTestMethod.SKAT_O)
                throw new IllegalArgumentException("Kernel reference requires skat or skat-o");
            rhoGrid = rhoGrid == null ? null : rhoGrid.clone();
            if (method == QSetTestMethod.SKAT_O
                && (rhoGrid == null || rhoGrid.length == 0 || simulations < 1))
                throw new IllegalArgumentException("SKAT-O requires rho values and simulations");
            if (method == QSetTestMethod.SKAT_O) {
                if (rhoGrid[0] != 0 || rhoGrid[rhoGrid.length - 1] != 1)
                    throw new IllegalArgumentException("SKAT-O rho grid must start at 0 and end at 1");
                for (int i = 0; i < rhoGrid.length; i++)
                    if (!Double.isFinite(rhoGrid[i]) || rhoGrid[i] < 0 || rhoGrid[i] > 1
                        || (i > 0 && rhoGrid[i] <= rhoGrid[i - 1]))
                        throw new IllegalArgumentException(
                            "SKAT-O rho grid must be strictly increasing within [0,1]");
            }
        }
        @Override public double[] rhoGrid() { return rhoGrid == null ? null : rhoGrid.clone(); }
    }

    public record Result(String setId, String traitId, int variantCount, int sampleCount,
        int residualDegreesOfFreedom, double statistic, double pValue, double log10P,
        String pValueMethod, double minimumComponentP, String rhoComponentPValues) { }

    public record Analysis(List<SetAudit> audits, List<Result> results) {
        public Analysis {
            audits = List.copyOf(audits);
            results = List.copyOf(results);
        }
    }

    private QKernelSetReference() { }

    public static Analysis analyze(QVariantSetTable table, List<Variant> variants,
        QSetTestNullModel nullModel, QSetTestPolicy policy, Options options) {
        if (table == null || variants == null || nullModel == null || policy == null
            || options == null)
            throw new IllegalArgumentException("Kernel set-test inputs are required");
        Map<String, Variant> source = indexVariants(variants, nullModel.sampleCount());
        List<SetAudit> audits = new ArrayList<>();
        List<Result> results = new ArrayList<>();
        for (SetDefinition set : table.sets())
            analyzeSet(set, source, nullModel, policy, options, audits, results);
        return new Analysis(audits, results);
    }

    private static void analyzeSet(SetDefinition set, Map<String, Variant> source,
        QSetTestNullModel nullModel, QSetTestPolicy policy, Options options,
        List<SetAudit> audits, List<Result> results) {
        List<double[]> dosages = new ArrayList<>();
        List<Double> weights = new ArrayList<>();
        List<String> includedIds = new ArrayList<>();
        int absent = 0;
        int frequencyExcluded = 0;
        for (Entry entry : set.entries()) {
            Variant variant = source.get(entry.variantId());
            if (variant == null) {
                absent++;
                if (policy.absentVariant() == FailurePolicy.ERROR)
                    throw new IllegalArgumentException("Variant '" + entry.variantId()
                        + "' from set '" + set.id() + "' is absent from the aligned source");
                continue;
            }
            entry.requireExactAlleles(variant.ref(), variant.alt());
            Frequency frequency = frequency(variant);
            if (!policy.includes(frequency.maf(), frequency.mac())) {
                frequencyExcluded++;
                continue;
            }
            dosages.add(effectDosages(variant, entry, frequency, policy.predictorMissing()));
            weights.add(entry.weight());
            includedIds.add(entry.variantId());
        }
        if (dosages.isEmpty()) {
            handleDegenerate(set, absent, frequencyExcluded, includedIds, policy,
                audits, "skipped-empty", "has no variants after source and frequency filtering");
            return;
        }
        double[][] residualVariants = residualize(dosages.toArray(double[][]::new),
            nullModel.covariateQ());
        double[] numericWeights = weights.stream().mapToDouble(Double::doubleValue).toArray();
        List<Result> setResults = new ArrayList<>();
        PreparedBlock traits = nullModel.traits();
        try {
            for (int trait = 0; trait < traits.values().length; trait++) {
                double[] residualTrait = unstandardize(traits.values()[trait],
                    traits.standardDeviations()[trait]);
                if (options.method() == QSetTestMethod.SKAT) {
                    QSkatReference.Result result = QSkatReference.calculate(residualVariants,
                        numericWeights, residualTrait, nullModel.nullResidualDegreesOfFreedom());
                    setResults.add(new Result(set.id(), traits.rowIds()[trait], includedIds.size(),
                        nullModel.sampleCount(), nullModel.nullResidualDegreesOfFreedom(),
                        result.statistic(), result.pValue(), result.log10P(), result.pValueMethod(),
                        Double.NaN, ""));
                } else {
                    long traitSeed = options.seed() ^ stableHash(set.id())
                        ^ Long.rotateLeft(stableHash(traits.rowIds()[trait]), 23);
                    QSkatReference.OmnibusResult result = QSkatReference.calculateOmnibus(
                        residualVariants, numericWeights, residualTrait,
                        nullModel.nullResidualDegreesOfFreedom(), options.rhoGrid(),
                        options.simulations(), traitSeed);
                    setResults.add(new Result(set.id(), traits.rowIds()[trait], includedIds.size(),
                        nullModel.sampleCount(), nullModel.nullResidualDegreesOfFreedom(),
                        Double.NaN, result.adjustedP(), Math.log10(result.adjustedP()),
                        result.adjustmentMethod() + ";simulations=" + result.simulations()
                            + ";seed=" + traitSeed,
                        result.minimumComponentP(), componentValues(result.components())));
                }
            }
        } catch (IllegalArgumentException e) {
            if (policy.degenerateSet() == FailurePolicy.ERROR)
                throw new IllegalArgumentException("Set '" + set.id()
                    + "' has a zero or invalid projected kernel", e);
            audits.add(new SetAudit(set.id(), set.entries().size(), absent,
                frequencyExcluded, includedIds.size(), includedIds, "skipped-degenerate-kernel"));
            return;
        }
        audits.add(new SetAudit(set.id(), set.entries().size(), absent,
            frequencyExcluded, includedIds.size(), includedIds, "analyzed"));
        results.addAll(setResults);
    }

    private static void handleDegenerate(SetDefinition set, int absent, int excluded,
        List<String> includedIds, QSetTestPolicy policy, List<SetAudit> audits,
        String status, String message) {
        if (policy.degenerateSet() == FailurePolicy.ERROR)
            throw new IllegalArgumentException("Set '" + set.id() + "' " + message);
        audits.add(new SetAudit(set.id(), set.entries().size(), absent, excluded,
            includedIds.size(), includedIds, status));
    }

    private static Map<String, Variant> indexVariants(List<Variant> variants, int samples) {
        Map<String, Variant> indexed = new LinkedHashMap<>();
        for (Variant variant : variants) {
            if (variant == null || variant.altDosagesInternal().length != samples)
                throw new IllegalArgumentException("Kernel-reference variant has wrong sample count");
            if (indexed.putIfAbsent(variant.id(), variant) != null)
                throw new IllegalArgumentException("Duplicate aligned variant ID '" + variant.id() + "'");
        }
        return indexed;
    }

    private static Frequency frequency(Variant variant) {
        int called = 0;
        double sum = 0;
        for (double dosage : variant.altDosagesInternal()) {
            if (isMissing(dosage)) continue;
            sum += validateDosage(variant.id(), dosage);
            called++;
        }
        if (called == 0) return new Frequency(Double.NaN, Double.NaN, Double.NaN);
        double eaf = sum / (2.0 * called);
        double maf = Math.min(eaf, 1 - eaf);
        return new Frequency(eaf, maf, 2.0 * called * maf);
    }

    private static double[] effectDosages(Variant variant, Entry entry, Frequency frequency,
        QMissingValuePolicy missingPolicy) {
        double[] source = variant.altDosagesInternal();
        double[] result = new double[source.length];
        double mean = 2 * (entry.effectIsAlt() ? frequency.eaf() : 1 - frequency.eaf());
        for (int sample = 0; sample < source.length; sample++) {
            if (isMissing(source[sample])) {
                if (missingPolicy == QMissingValuePolicy.ERROR)
                    throw new IllegalArgumentException("Variant '" + variant.id()
                        + "' contains a missing dosage under set-test policy error");
                result[sample] = missingPolicy == QMissingValuePolicy.MEAN ? mean : 0;
            } else {
                double alt = validateDosage(variant.id(), source[sample]);
                result[sample] = entry.effectIsAlt() ? alt : 2 - alt;
            }
        }
        return result;
    }

    private static double[][] residualize(double[][] rows, double[][] q) {
        return QMatrixUtils.parallelMatrixMultiplication(rows, q, null, 1,
            rows.length, rows[0].length, EMultiplicationMode.XMinusXYYt);
    }

    private static double[] unstandardize(double[] values, double standardDeviation) {
        double[] result = values.clone();
        for (int i = 0; i < result.length; i++) result[i] *= standardDeviation;
        return result;
    }

    private static String componentValues(List<QSkatReference.Component> components) {
        StringBuilder result = new StringBuilder();
        for (QSkatReference.Component component : components) {
            if (result.length() > 0) result.append(';');
            result.append(component.rho()).append(':').append(component.result().pValue());
        }
        return result.toString();
    }

    private static long stableHash(String value) {
        long hash = 0xcbf29ce484222325L;
        for (int i = 0; i < value.length(); i++) {
            hash ^= value.charAt(i);
            hash *= 0x100000001b3L;
        }
        return hash;
    }

    private static double validateDosage(String id, double dosage) {
        if (!Double.isFinite(dosage) || dosage < -DOSAGE_TOLERANCE
            || dosage > 2 + DOSAGE_TOLERANCE)
            throw new IllegalArgumentException("Variant '" + id
                + "' has a non-finite or out-of-range additive dosage " + dosage);
        return Math.max(0, Math.min(2, dosage));
    }

    private static boolean isMissing(double dosage) {
        return Double.isNaN(dosage) || dosage == kUndefinedValue;
    }

    private record Frequency(double eaf, double maf, double mac) { }
}
