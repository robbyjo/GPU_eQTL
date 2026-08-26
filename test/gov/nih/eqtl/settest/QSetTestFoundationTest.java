/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QMissingValuePolicy;
import gov.nih.eqtl.settest.QBurdenReference.Analysis;
import gov.nih.eqtl.settest.QBurdenReference.Result;
import gov.nih.eqtl.settest.QBurdenReference.SetAudit;
import gov.nih.eqtl.settest.QBurdenReference.Variant;
import gov.nih.eqtl.settest.QSetTestPolicy.FailurePolicy;

class QSetTestFoundationTest {
    private static final Path FIXTURE = Path.of("test/resources/set-test-reference");

    @Test
    void parsesDeterministicWeightedOverlappingMembership() throws Exception {
        QVariantSetTable first = QVariantSetTable.load(FIXTURE.resolve("sets.tsv"));
        QVariantSetTable second = QVariantSetTable.load(FIXTURE.resolve("sets.tsv"));

        assertArrayEquals(new String[] {"geneA", "geneB", "geneEmpty", "mono"},
            first.sets().stream().map(QVariantSetTable.SetDefinition::id).toArray(String[]::new));
        assertArrayEquals(new String[] {"v1", "v2"}, first.sets().get(0).entries().stream()
            .map(QVariantSetTable.Entry::variantId).toArray(String[]::new));
        assertEquals(0.5, first.sets().get(0).entries().get(1).weight());
        assertEquals("C", first.sets().get(0).entries().get(1).effectAllele());
        assertEquals("v2", first.sets().get(1).entries().get(0).variantId());
        assertEquals(first.signature(), second.signature());
        assertEquals(64, first.signature().length());
    }

    @Test
    void weightedBurdenMatchesIndependentFp64FixtureAndAuditsFiltering() throws Exception {
        QVariantSetTable sets = QVariantSetTable.load(FIXTURE.resolve("sets.tsv"));
        QSetTestNullModel nullModel = nullModel();
        QSetTestPolicy policy = new QSetTestPolicy(0, 0.2, 0, Double.POSITIVE_INFINITY,
            QMissingValuePolicy.MEAN, FailurePolicy.SKIP, FailurePolicy.SKIP);

        Analysis analysis = QBurdenReference.analyze(sets, variants(), nullModel, policy);

        assertEquals(2, nullModel.covariateRank());
        assertEquals(5, nullModel.residualDegreesOfFreedom());
        assertEquals(4, analysis.audits().size());
        assertAudit(analysis.audits().get(0), "geneA", 2, 0, 0, 2,
            List.of("v1", "v2"), "analyzed");
        assertAudit(analysis.audits().get(1), "geneB", 2, 0, 1, 1,
            List.of("v2"), "analyzed");
        assertAudit(analysis.audits().get(2), "geneEmpty", 1, 1, 0, 0,
            List.of(), "skipped-empty");
        assertAudit(analysis.audits().get(3), "mono", 1, 0, 0, 1,
            List.of("v4"), "skipped-monomorphic");

        assertArrayEquals(new String[] {"geneA/trait1", "geneA/trait2",
            "geneB/trait1", "geneB/trait2"}, analysis.results().stream()
                .map(result -> result.setId() + "/" + result.traitId()).toArray(String[]::new));
        for (Result result : analysis.results()) {
            assertEquals(8, result.sampleCount());
            assertEquals(5, result.residualDegreesOfFreedom());
        }
        assertArrayEquals(new double[] {
            0.6747780738826216, 0.5395555683137198,
            0.01083266499989602, 0.009789919608688065
        }, analysis.results().stream().mapToDouble(Result::rSquared).toArray(), 1e-12);
        assertArrayEquals(new double[] {
            0.406515205649314, -0.10457115235072852,
            0.14854046495368042, -0.04062226883412121
        }, analysis.results().stream().mapToDouble(Result::effect).toArray(), 1e-12);
        assertArrayEquals(new double[] {
            3.220887637744619, -2.4205523460484324,
            0.23400103342370288, -0.22233656583953693
        }, analysis.results().stream().mapToDouble(Result::tStatistic).toArray(), 1e-12);
        assertArrayEquals(new double[] {
            -1.6299970085493116, -1.2212956736106053,
            -0.08393386934985025, -0.07943403355250683
        }, analysis.results().stream().mapToDouble(Result::log10P).toArray(), 1e-12);
    }

    @Test
    void rejectsAlleleMismatchMissingDosageDuplicateMembershipAndRankLoss(
        @TempDir Path directory) throws Exception {
        QVariantSetTable sets = QVariantSetTable.load(FIXTURE.resolve("sets.tsv"));
        QSetTestNullModel nullModel = nullModel();
        QSetTestPolicy missingError = new QSetTestPolicy(0, 0.2, 0,
            Double.POSITIVE_INFINITY, QMissingValuePolicy.ERROR,
            FailurePolicy.SKIP, FailurePolicy.SKIP);
        IllegalArgumentException missing = assertThrows(IllegalArgumentException.class,
            () -> QBurdenReference.analyze(sets, variants(), nullModel, missingError));
        assertTrue(missing.getMessage().contains("missing dosage"));

        List<Variant> mismatched = List.of(
            new Variant("v1", "A", "T", new double[] {0, 0, 0, 1, 0, 0, 1, 0}),
            variants().get(1), variants().get(2), variants().get(3));
        IllegalArgumentException alleles = assertThrows(IllegalArgumentException.class,
            () -> QBurdenReference.analyze(sets, mismatched, nullModel,
                new QSetTestPolicy(0, 0.2, 0, Double.POSITIVE_INFINITY,
                    QMissingValuePolicy.MEAN, FailurePolicy.SKIP, FailurePolicy.SKIP)));
        assertTrue(alleles.getMessage().contains("Allele mismatch"));

        Path duplicate = directory.resolve("duplicate.tsv");
        Files.writeString(duplicate, "SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\n"
            + "gene\tv1\tA\tG\tG\n"
            + "gene\tv1\tA\tG\tG\n");
        IOException duplicateFailure = assertThrows(IOException.class,
            () -> QVariantSetTable.load(duplicate));
        assertTrue(duplicateFailure.getMessage().contains("Duplicate variant"));

        double[][] rankDeficient = new double[8][2];
        for (int sample = 0; sample < rankDeficient.length; sample++) {
            rankDeficient[sample][0] = 1;
            rankDeficient[sample][1] = 1;
        }
        IllegalArgumentException rank = assertThrows(IllegalArgumentException.class,
            () -> QSetTestNullModel.create(new String[] {"trait"},
                new double[][] {{1, 2, 3, 4, 5, 6, 7, 8}}, rankDeficient));
        assertTrue(rank.getMessage().contains("rank deficient"));
    }

    private static QSetTestNullModel nullModel() throws IOException {
        List<String> covariateLines = Files.readAllLines(FIXTURE.resolve("covariates.csv"));
        double[][] design = new double[covariateLines.size() - 1][2];
        for (int sample = 0; sample < design.length; sample++) {
            String[] fields = covariateLines.get(sample + 1).split(",", -1);
            design[sample][0] = 1;
            design[sample][1] = Double.parseDouble(fields[1]);
        }
        List<String> traitLines = Files.readAllLines(FIXTURE.resolve("traits.csv"));
        String[] traitIds = new String[traitLines.size() - 1];
        double[][] traitValues = new double[traitIds.length][design.length];
        for (int trait = 0; trait < traitIds.length; trait++) {
            String[] fields = traitLines.get(trait + 1).split(",", -1);
            traitIds[trait] = fields[0];
            for (int sample = 0; sample < design.length; sample++)
                traitValues[trait][sample] = Double.parseDouble(fields[sample + 1]);
        }
        return QSetTestNullModel.create(traitIds, traitValues, design);
    }

    private static List<Variant> variants() throws IOException {
        List<String> lines = Files.readAllLines(FIXTURE.resolve("variants.csv"));
        int samples = lines.get(0).split(",", -1).length - 3;
        List<Variant> variants = new ArrayList<>();
        for (int row = 1; row < lines.size(); row++) {
            String[] fields = lines.get(row).split(",", -1);
            double[] dosages = new double[samples];
            for (int sample = 0; sample < samples; sample++)
                dosages[sample] = fields[sample + 3].isBlank()
                    ? Double.NaN : Double.parseDouble(fields[sample + 3]);
            variants.add(new Variant(fields[0], fields[1], fields[2], dosages));
        }
        return List.copyOf(variants);
    }

    private static void assertAudit(SetAudit actual, String setId, int requested,
        int absent, int frequencyExcluded, int included, List<String> ids, String status) {
        assertEquals(setId, actual.setId());
        assertEquals(requested, actual.requestedVariants());
        assertEquals(absent, actual.absentVariants());
        assertEquals(frequencyExcluded, actual.frequencyExcludedVariants());
        assertEquals(included, actual.includedVariants());
        assertEquals(ids, actual.includedVariantIds());
        assertEquals(status, actual.status());
    }
}
