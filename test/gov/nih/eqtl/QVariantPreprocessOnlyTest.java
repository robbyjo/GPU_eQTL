/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.FileTime;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QVariantPreprocessOnlyTest {
    @Test
    void defaultMacTwentyIsAppliedToTheAlignedVcfCohort(@TempDir Path directory)
        throws Exception {
        Path genotype = Path.of("test/resources/variant-reference/genotype.vcf")
            .toAbsolutePath().normalize();
        Path expression = Path.of("test/resources/variant-reference/expression.csv")
            .toAbsolutePath().normalize();
        Map<String, String> values = new LinkedHashMap<>();
        values.put("ini.path", directory.toString() + java.io.File.separator);
        values.put("genotype_file", genotype.toString());
        values.put("expression_file", expression.toString());
        values.put("genotype_format", "vcf");
        values.put("preprocess_only", "true");
        values.put("trait_missing", "error");
        values.put("variant_qc_output", directory.resolve("variants.tsv").toString());
        QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
        QeQTLAnalysis.profiler = new QeQTLProfiler(false);

        IOException error = assertThrows(IOException.class,
            () -> QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
                expression.toString(), null, null, null, "none", Double.NaN, 0,
                true, 0, 1, "vcf"));
        assertTrue(error.getMessage().contains("No variants remain"));
        assertTrue(error.getMessage().contains("--min-mac/--min-maf"));
        assertTrue(error.getMessage().contains("min_mac=20"));
        assertEquals(20, QeQTLAnalysis.config.getMinimumMac());
    }

    @Test
    void alignedVcfPreprocessingCreatesAndReusesDurableRawCache(@TempDir Path directory)
        throws Exception {
        Path genotype = Path.of("test/resources/variant-reference/genotype.vcf")
            .toAbsolutePath().normalize();
        Path expression = Path.of("test/resources/variant-reference/expression.csv")
            .toAbsolutePath().normalize();
        Path cache = directory.resolve("cache");
        Path qc = directory.resolve("variants.tsv");
        Path missingness = directory.resolve("missingness.tsv");
        Map<String, String> values = new LinkedHashMap<>();
        values.put("ini.path", directory.toString() + java.io.File.separator);
        values.put("genotype_file", genotype.toString());
        values.put("expression_file", expression.toString());
        values.put("genotype_format", "vcf");
        values.put("preprocess_only", "true");
        values.put("trait_missing", "error");
        values.put("min_mac", "0");
        values.put("cache_dir", cache.toString());
        values.put("variant_qc_output", qc.toString());
        values.put("missingness_qc_output", missingness.toString());
        values.put("genotype_block_rows", "2");
        QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
        QeQTLAnalysis.profiler = new QeQTLProfiler(false);

        QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
            expression.toString(), null, null, null, "none", Double.NaN, 0,
            true, 0, 1, "vcf");

        assertTrue(Files.isRegularFile(qc));
        assertTrue(Files.isRegularFile(missingness));
        List<Path> rawCaches;
        try (var paths = Files.list(cache.resolve("aligned-raw"))) {
            rawCaches = paths.filter(path -> path.getFileName().toString().endsWith(".qraw")).toList();
        }
        assertEquals(1, rawCaches.size());
        assertFalse(Files.exists(directory.resolve("results.csv")));
        byte[] first = Files.readAllBytes(rawCaches.get(0));
        FileTime modified = Files.getLastModifiedTime(rawCaches.get(0));

        QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
            expression.toString(), null, null, null, "none", Double.NaN, 0,
            true, 0, 1, "vcf");
        assertEquals(modified, Files.getLastModifiedTime(rawCaches.get(0)));
        assertTrue(java.util.Arrays.equals(first, Files.readAllBytes(rawCaches.get(0))));
    }

    @Test
    void alignedRawCacheDoesNotReapplySourceColumnIndicesToCachedSampleIds(
        @TempDir Path directory) throws Exception {
        Path genotype = Path.of("test/resources/variant-reference/genotype.vcf")
            .toAbsolutePath().normalize();
        Path expression = Path.of("test/resources/variant-reference/expression.csv")
            .toAbsolutePath().normalize();
        Path covariates = directory.resolve("covariates.csv");
        Files.writeString(covariates,
            "genotype_id,expression_id\nS4,S4\nS1,S1\nS2,S2\n");
        Path cache = directory.resolve("cache");
        Path missingness = directory.resolve("missingness.tsv");
        Map<String, String> values = new LinkedHashMap<>();
        values.put("ini.path", directory.toString() + java.io.File.separator);
        values.put("genotype_file", genotype.toString());
        values.put("expression_file", expression.toString());
        values.put("covariate_file", covariates.toString());
        values.put("genotype_format", "vcf");
        values.put("genotype_id_column", "genotype_id");
        values.put("expression_id_column", "expression_id");
        values.put("sample_alignment", "covariate-subset");
        values.put("preprocess_only", "true");
        values.put("trait_missing", "error");
        values.put("min_mac", "0");
        values.put("cache_dir", cache.toString());
        values.put("variant_qc_output", directory.resolve("variants.tsv").toString());
        values.put("missingness_qc_output", missingness.toString());
        values.put("genotype_block_rows", "2");
        QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
        QeQTLAnalysis.profiler = new QeQTLProfiler(false);

        QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
            expression.toString(), covariates.toString(), null, null,
            "none", Double.NaN, 0, true, 0, 1, "vcf");
        QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
            expression.toString(), covariates.toString(), null, null,
            "none", Double.NaN, 0, true, 0, 1, "vcf");

        String report = Files.readString(missingness);
        assertTrue(report.contains("selected_samples=3;matrix_only_samples_excluded=1"));
        assertTrue(report.contains("\nSAMPLE\tpredictor\t1\tS1\t"));
    }
}
