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
import gov.nih.eqtl.QeQTLPreprocessor;
import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.QeQTLStatistics;
import gov.nih.eqtl.io.QMatrixRowSource.Block;
import gov.nih.eqtl.settest.QSetTestPolicy.FailurePolicy;
import gov.nih.eqtl.settest.QVariantSetTable.Entry;
import gov.nih.eqtl.settest.QVariantSetTable.SetDefinition;
import gov.nih.utils.matrix.EMultiplicationMode;
import gov.nih.utils.matrix.QMatrixUtils;

/** Deterministic scalar/CPU weighted-burden reference for continuous traits. */
public final class QBurdenReference {
    private static final double DOSAGE_TOLERANCE = 1e-12;

    public record Variant(String id, String ref, String alt, double[] altDosages) {
        public Variant {
            if (id == null || id.isBlank() || ref == null || ref.isBlank()
                || alt == null || alt.isBlank() || altDosages == null || altDosages.length == 0)
                throw new IllegalArgumentException("A burden-reference variant is incomplete");
            id = id.trim();
            ref = ref.trim().toUpperCase(java.util.Locale.ROOT);
            alt = alt.trim().toUpperCase(java.util.Locale.ROOT);
            altDosages = altDosages.clone();
        }

        @Override
        public double[] altDosages() { return altDosages.clone(); }
    }

    public record SetAudit(String setId, int requestedVariants, int absentVariants,
        int frequencyExcludedVariants, int includedVariants, List<String> includedVariantIds,
        String status) {
        public SetAudit {
            includedVariantIds = List.copyOf(includedVariantIds);
        }
    }

    public record Result(String setId, String traitId, int variantCount, int sampleCount,
        int residualDegreesOfFreedom, double rSquared, double effect,
        double tStatistic, double log10P) { }

    public record Analysis(List<SetAudit> audits, List<Result> results) {
        public Analysis {
            audits = List.copyOf(audits);
            results = List.copyOf(results);
        }
    }

    private QBurdenReference() { }

    public static Analysis analyze(QVariantSetTable table, List<Variant> variants,
        QSetTestNullModel nullModel, QSetTestPolicy policy) {
        return analyze(table, variants, nullModel, policy, false);
    }

    /** Production FP64 path: prepares each burden once and batches all set-by-trait products. */
    public static Analysis analyzeBatched(QVariantSetTable table, List<Variant> variants,
        QSetTestNullModel nullModel, QSetTestPolicy policy) {
        return analyze(table, variants, nullModel, policy, true);
    }

    private static Analysis analyze(QVariantSetTable table, List<Variant> variants,
        QSetTestNullModel nullModel, QSetTestPolicy policy, boolean batched) {
        if (table == null || variants == null || nullModel == null || policy == null)
            throw new IllegalArgumentException("Set definitions, variants, null model, and policy are required");
        Map<String, Variant> source = indexVariants(variants, nullModel.sampleCount());
        List<SetAudit> audits = new ArrayList<>();
        List<Result> results = new ArrayList<>();
        List<PreparedSet> preparedSets = new ArrayList<>();
        for (SetDefinition set : table.sets()) {
            PreparedSet prepared = prepareSet(set, source, nullModel, policy, audits);
            if (prepared != null) preparedSets.add(prepared);
        }
        if (batched) calculateBatched(preparedSets, nullModel, results);
        else calculateScalar(preparedSets, nullModel, results);
        return new Analysis(audits, results);
    }

    private static PreparedSet prepareSet(SetDefinition set, Map<String, Variant> source,
        QSetTestNullModel nullModel, QSetTestPolicy policy, List<SetAudit> audits) {
        double[] burden = new double[nullModel.sampleCount()];
        int absent = 0;
        int frequencyExcluded = 0;
        List<String> includedIds = new ArrayList<>();
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
            addVariant(burden, variant, entry, frequency, policy.predictorMissing());
            includedIds.add(entry.variantId());
        }
        if (includedIds.isEmpty()) {
            if (policy.degenerateSet() == FailurePolicy.ERROR)
                throw new IllegalArgumentException("Set '" + set.id()
                    + "' has no variants after source and frequency filtering");
            audits.add(new SetAudit(set.id(), set.entries().size(), absent,
                frequencyExcluded, 0, List.of(), "skipped-empty"));
            return null;
        }
        PreparedBlock prepared;
        try {
            prepared = QeQTLPreprocessor.prepare(new Block(0, new String[] {set.id()},
                new double[][] {burden}), nullModel.covariateQ(), "Burden");
        } catch (IllegalArgumentException e) {
            if (policy.degenerateSet() == FailurePolicy.ERROR)
                throw new IllegalArgumentException("Set '" + set.id()
                    + "' has zero or invalid burden variance after covariate adjustment", e);
            audits.add(new SetAudit(set.id(), set.entries().size(), absent,
                frequencyExcluded, includedIds.size(), includedIds, "skipped-monomorphic"));
            return null;
        }
        audits.add(new SetAudit(set.id(), set.entries().size(), absent,
            frequencyExcluded, includedIds.size(), includedIds, "analyzed"));
        return new PreparedSet(set.id(), includedIds.size(), prepared.values()[0],
            prepared.standardDeviations()[0]);
    }

    private static void calculateScalar(List<PreparedSet> sets,
        QSetTestNullModel nullModel, List<Result> results) {
        PreparedBlock traits = nullModel.traits();
        for (PreparedSet set : sets)
            for (int trait = 0; trait < traits.values().length; trait++)
                addResult(results, set, traits, trait,
                    QeQTLPreprocessor.correlation(set.standardizedBurden(),
                        traits.values()[trait]), nullModel);
    }

    private static void calculateBatched(List<PreparedSet> sets,
        QSetTestNullModel nullModel, List<Result> results) {
        if (sets.isEmpty()) return;
        PreparedBlock traits = nullModel.traits();
        double[][] burdens = new double[sets.size()][];
        for (int set = 0; set < sets.size(); set++) burdens[set] = sets.get(set).standardizedBurden();
        double[][] products = QMatrixUtils.parallelMatrixMultiplication(burdens,
            QMatrixUtils.transpose(traits.values()), null, 1, burdens.length,
            traits.values().length, EMultiplicationMode.XY);
        double denominator = nullModel.sampleCount() - 1.0;
        for (int set = 0; set < sets.size(); set++)
            for (int trait = 0; trait < traits.values().length; trait++)
                addResult(results, sets.get(set), traits, trait,
                    products[set][trait] / denominator, nullModel);
    }

    private static void addResult(List<Result> results, PreparedSet set,
        PreparedBlock traits, int trait, double correlation, QSetTestNullModel nullModel) {
        QeQTLStatistics.Result statistic = QeQTLStatistics.calculate(correlation,
            traits.standardDeviations()[trait], set.standardDeviation(),
            nullModel.residualDegreesOfFreedom(), 0);
        results.add(new Result(set.id(), traits.rowIds()[trait], set.variantCount(),
            nullModel.sampleCount(), nullModel.residualDegreesOfFreedom(),
            statistic.rSquared(), statistic.effect(), statistic.tStatistic(),
            statistic.log10P()));
    }

    private static Map<String, Variant> indexVariants(List<Variant> variants, int samples) {
        LinkedHashMap<String, Variant> indexed = new LinkedHashMap<>();
        for (Variant variant : variants) {
            if (variant == null || variant.altDosages.length != samples)
                throw new IllegalArgumentException("Burden-reference variant has the wrong sample count");
            if (indexed.putIfAbsent(variant.id(), variant) != null)
                throw new IllegalArgumentException("Duplicate aligned variant ID '" + variant.id() + "'");
        }
        return indexed;
    }

    private static Frequency frequency(Variant variant) {
        int called = 0;
        double sum = 0;
        for (double dosage : variant.altDosages) {
            if (isMissing(dosage))
                continue;
            double valid = validateDosage(variant.id(), dosage);
            called++;
            sum += valid;
        }
        if (called == 0)
            return new Frequency(0, Double.NaN, Double.NaN, Double.NaN);
        double eaf = sum / (2.0 * called);
        double maf = Math.min(eaf, 1 - eaf);
        return new Frequency(called, eaf, maf, 2.0 * called * maf);
    }

    private static void addVariant(double[] burden, Variant variant, Entry entry,
        Frequency frequency, QMissingValuePolicy missingPolicy) {
        double meanEffectDosage = entry.effectIsAlt()
            ? frequency.eaf() * 2 : (1 - frequency.eaf()) * 2;
        for (int sample = 0; sample < burden.length; sample++) {
            double dosage = variant.altDosages[sample];
            double effectDosage;
            if (isMissing(dosage)) {
                if (missingPolicy == QMissingValuePolicy.ERROR)
                    throw new IllegalArgumentException("Variant '" + variant.id()
                        + "' contains a missing dosage under set-test policy error");
                effectDosage = missingPolicy == QMissingValuePolicy.MEAN
                    ? meanEffectDosage : 0.0;
            } else {
                double altDosage = validateDosage(variant.id(), dosage);
                effectDosage = entry.effectIsAlt() ? altDosage : 2.0 - altDosage;
            }
            burden[sample] += entry.weight() * effectDosage;
        }
    }

    private static double validateDosage(String variantId, double dosage) {
        if (!Double.isFinite(dosage) || dosage < -DOSAGE_TOLERANCE
            || dosage > 2 + DOSAGE_TOLERANCE)
            throw new IllegalArgumentException("Variant '" + variantId
                + "' has a non-finite or out-of-range additive dosage " + dosage);
        return Math.max(0, Math.min(2, dosage));
    }

    private static boolean isMissing(double dosage) {
        return Double.isNaN(dosage) || dosage == kUndefinedValue;
    }

    private record Frequency(int called, double eaf, double maf, double mac) { }
    private record PreparedSet(String id, int variantCount, double[] standardizedBurden,
        double standardDeviation) { }
}
