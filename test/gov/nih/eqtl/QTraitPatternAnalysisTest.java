/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuPrecision;
import gov.nih.gpu.cpu.CpuBackend;

class QTraitPatternAnalysisTest {
    private static final String FAIL_PROPERTY = "eqtl.test.trait.pattern.fail.after";
    private static final String GENOTYPE_FAIL_PROPERTY =
        "eqtl.test.genotype.outer.fail.before.assembly";

    @AfterEach
    void restoreState() {
        System.clearProperty(FAIL_PROPERTY);
        System.clearProperty(GENOTYPE_FAIL_PROPERTY);
        System.clearProperty(CpuBackend.BLAS_PROPERTY);
        if (QeQTLAnalysis.mContexts != null)
            for (GpuContext context : QeQTLAnalysis.mContexts)
                if (context != null) context.close();
        QeQTLAnalysis.mContexts = null;
        QeQTLAnalysis.simplifyResult = false;
        QeQTLAnalysis.rsqOnly = false;
        QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        QeQTLAnalysis.residualizationMode = QResidualizationMode.AUTO;
    }

    @Test
    void exactPatternOutputMatchesCpuReferenceAndResumeIsByteIdentical(
        @TempDir Path directory) throws Exception {
        Path baselineOutput = directory.resolve("baseline.csv");
        runAnalysis(directory.resolve("baseline"), baselineOutput, false, false);
        List<String> rows = Files.readAllLines(baselineOutput);
        assertEquals("Rs_ID,ProbesetID,RSq,Fx,T,log10P,N,DF", rows.get(0));
        assertEquals(7, rows.size());
        String[] expectedIds = {"rs1/iso1", "rs2/iso1", "rs3/iso1",
            "rs1/iso2", "rs2/iso2", "rs3/iso2"};
        double[] expectedRsq = {0.33429286679254133, 0.18924585599861957,
            0.2735057531379482, 0.22047136051758234, 0.03308959550638096,
            0.35070209723441076};
        double[] expectedEffects = {0.2594575510480451, -0.24749432530757393,
            0.38742926434923264, -0.11012349969527942, 0.05316436390087009,
            -0.23150143426852185};
        double[] expectedT = {1.4172678817212558, -0.9662702228081002,
            1.2271485511712632, -0.921129477588173, 0.32041529561311494,
            -1.272939634032138};
        double[] expectedLog10P = {-0.6394636593437162, -0.410466712182518,
            -0.5420393432793345, -0.3716960473068875, -0.11369289994852227,
            -0.533533124937666};
        for (int i = 0; i < expectedIds.length; i++) {
            String[] fields = rows.get(i + 1).split(",", -1);
            assertEquals(expectedIds[i], fields[0] + "/" + fields[1]);
            assertEquals(expectedRsq[i], Double.parseDouble(fields[2]), 1e-12);
            assertEquals(expectedEffects[i], Double.parseDouble(fields[3]), 1e-12);
            assertEquals(expectedT[i], Double.parseDouble(fields[4]), 1e-12);
            assertEquals(expectedLog10P[i], Double.parseDouble(fields[5]), 1e-12);
            assertEquals(i < 3 ? 8 : 7, Integer.parseInt(fields[6]));
            assertEquals(i < 3 ? 4 : 3, Integer.parseInt(fields[7]));
        }

        Path resumedRoot = directory.resolve("resumed");
        Path resumedOutput = directory.resolve("resumed.csv");
        System.setProperty(FAIL_PROPERTY, "1");
        Exception interruption = assertThrows(Exception.class,
            () -> runAnalysis(resumedRoot, resumedOutput, false, true));
        assertTrue(messageChain(interruption).contains("Injected test failure after trait pattern 1"));
        System.clearProperty(FAIL_PROPERTY);
        assertTrue(Files.isRegularFile(resumedRoot.resolve("checkpoint")
            .resolve("pattern-00000000.results.part")));

        runAnalysis(resumedRoot, resumedOutput, true, true);
        assertArrayEquals(Files.readAllBytes(baselineOutput), Files.readAllBytes(resumedOutput));
        assertTrue(Files.isRegularFile(resumedRoot.resolve("checkpoint")
            .resolve("pattern-00000001.results.part")));
    }

	@Test
    void patternSafetyLimitFailsBeforePatternCachesAreCreated(@TempDir Path directory)
		throws Exception {
		Path root = directory.resolve("limited");
		Exception error = assertThrows(IllegalArgumentException.class,
			() -> runAnalysis(root, directory.resolve("limited.csv"), false, false, 1,
				"pattern"));
		assertTrue(messageChain(error).contains("2 patterns exceed --max-trait-patterns 1"));
		assertTrue(!Files.exists(root.resolve("cache").resolve("trait-patterns")));
		assertTrue(!Files.exists(root.resolve("checkpoint")));
	}

	@Test
	void genotypeOuterSchedulerMatchesPatternOuterStatisticsWithoutPredictorPatternCaches(
		@TempDir Path directory) throws Exception {
		Path patternOutput = directory.resolve("pattern.csv");
		runAnalysis(directory.resolve("pattern"), patternOutput, false, false);
		Path genotypeOutput = directory.resolve("genotype.csv");
		Path genotypeRoot = directory.resolve("genotype");
		runAnalysis(genotypeRoot, genotypeOutput, false, false, null, "genotype");
		Map<String, double[]> expected = numericRows(patternOutput);
		Map<String, double[]> actual = numericRows(genotypeOutput);
		assertEquals(expected.keySet(), actual.keySet());
		for (String id : expected.keySet())
			assertArrayEquals(expected.get(id), actual.get(id), 2e-12, id);
		assertTrue(Files.isRegularFile(Path.of(genotypeOutput + ".trait-patterns.tsv")));
		assertTrue(!Files.exists(genotypeRoot.resolve("cache").resolve("trait-patterns")));
		assertTrue(Files.isDirectory(genotypeRoot.resolve("cache").resolve("trait-pattern-global")));

		Path resumedRoot = directory.resolve("genotype-resumed");
		Path resumedOutput = directory.resolve("genotype-resumed.csv");
		System.setProperty(GENOTYPE_FAIL_PROPERTY, "true");
		Exception interruption = assertThrows(Exception.class,
			() -> runAnalysis(resumedRoot, resumedOutput, false, true, null, "genotype"));
		assertTrue(messageChain(interruption).contains(
			"Injected test failure before genotype-outer assembly"));
		System.clearProperty(GENOTYPE_FAIL_PROPERTY);
		assertTrue(!Files.exists(resumedOutput));
		assertTrue(Files.isRegularFile(resumedRoot.resolve("checkpoint").resolve("block-00000000.part")));
		assertTrue(Files.isRegularFile(resumedRoot.resolve("checkpoint.pattern-qc")
			.resolve("block-00000000.qc")));
		runAnalysis(resumedRoot, resumedOutput, true, true, null, "genotype");
		assertArrayEquals(Files.readAllBytes(genotypeOutput), Files.readAllBytes(resumedOutput));
		assertEquals(qcRows(Path.of(genotypeOutput + ".pattern-variant-qc.tsv")),
			qcRows(Path.of(resumedOutput + ".pattern-variant-qc.tsv")));
	}

	@Test
	void unestimablePatternsAreAuditedBeforeErrorOrExplicitlySkipped(@TempDir Path directory)
		throws Exception {
		Path expression = directory.resolve("unestimable-expression.csv");
		Files.writeString(expression, ",T2,T8,T1,T5,T4,T7,T3,T6\n"
			+ "good,2.2,4.9,1.0,3.6,2.9,4.3,2.4,3.8\n"
			+ "bad,NA,NA,NA,NA,NA,1.9,3.7,2.3\n");
		Path errorOutput = directory.resolve("error.csv");
		Exception error = assertThrows(IllegalArgumentException.class,
			() -> runAnalysisWithExpression(directory.resolve("error"), errorOutput,
				expression, "error"));
		assertTrue(messageChain(error).contains("Unestimable trait missingness patterns"));
		Path errorAudit = Path.of(errorOutput + ".trait-patterns.tsv");
		assertTrue(Files.isRegularFile(errorAudit));
		assertTrue(Files.readString(errorAudit).contains("excluded\tnon-positive-residual-df"));

		Path skipOutput = directory.resolve("skip.csv");
		runAnalysisWithExpression(directory.resolve("skip"), skipOutput, expression, "skip");
		List<String> rows = Files.readAllLines(skipOutput);
		assertEquals(4, rows.size());
		assertTrue(rows.stream().skip(1).allMatch(row -> row.contains(",good,")));
		assertTrue(Files.readString(Path.of(skipOutput + ".trait-patterns.tsv"))
			.contains("excluded\tnon-positive-residual-df"));
	}

	@Test
	void genotypeOuterPatternMeanAndZeroImputationMatchPatternOuterVcfStatistics(
		@TempDir Path directory) throws Exception {
		for (String policy : List.of("mean", "zero")) {
			Path patternOutput = directory.resolve("vcf-pattern-" + policy + ".csv");
			runVariantPatternAnalysis(directory.resolve("vcf-pattern-" + policy),
				patternOutput, "pattern", policy);
			Path genotypeOutput = directory.resolve("vcf-genotype-" + policy + ".csv");
			runVariantPatternAnalysis(directory.resolve("vcf-genotype-" + policy),
				genotypeOutput, "genotype", policy);
			Map<String, double[]> expected = numericRows(patternOutput);
			Map<String, double[]> actual = numericRows(genotypeOutput);
			assertEquals(expected.keySet(), actual.keySet());
			for (String id : expected.keySet())
				assertArrayEquals(expected.get(id), actual.get(id), 2e-12,
					policy + ": " + id);
			assertEquals(qcRows(Path.of(patternOutput + ".pattern-variant-qc.tsv")),
				qcRows(Path.of(genotypeOutput + ".pattern-variant-qc.tsv")), policy);
		}
	}

    private static void runAnalysis(Path root, Path output, boolean resume,
        boolean keepCheckpoints) throws Exception {
		runAnalysis(root, output, resume, keepCheckpoints, null);
	}

	private static void runAnalysis(Path root, Path output, boolean resume,
		boolean keepCheckpoints, Integer maximumTraitPatterns) throws Exception {
		runAnalysis(root, output, resume, keepCheckpoints, maximumTraitPatterns, null);
	}

	private static void runAnalysis(Path root, Path output, boolean resume,
		boolean keepCheckpoints, Integer maximumTraitPatterns, String scheduler) throws Exception {
		runAnalysis(root, output, resume, keepCheckpoints, maximumTraitPatterns, scheduler,
			null, null);
	}

	private static void runAnalysisWithExpression(Path root, Path output, Path expression,
		String unestimablePolicy) throws Exception {
		runAnalysis(root, output, false, false, null, "genotype", expression,
			unestimablePolicy);
	}

	private static void runAnalysis(Path root, Path output, boolean resume,
		boolean keepCheckpoints, Integer maximumTraitPatterns, String scheduler,
		Path expressionOverride, String unestimablePolicy) throws Exception {
        Files.createDirectories(root);
        Path genotype = Path.of("test/resources/eqtl-reference/genotype.csv")
            .toAbsolutePath().normalize();
		Path expression = expressionOverride == null
			? Path.of("test/resources/missing-reference/expression.csv").toAbsolutePath().normalize()
			: expressionOverride.toAbsolutePath().normalize();
        Path covariates = Path.of("test/resources/eqtl-reference/covariates.csv")
            .toAbsolutePath().normalize();
        Map<String, String> values = new LinkedHashMap<>();
        values.put("ini.path", root.toString() + java.io.File.separator);
        values.put("genotype_file", genotype.toString());
        values.put("expression_file", expression.toString());
        values.put("covariate_file", covariates.toString());
        values.put("output_file", output.toString());
        values.put("genotype_format", "csv");
        values.put("trait_missing", "pattern");
        values.put("predictor_missing", "error");
        values.put("covariate_missing", "error");
        values.put("genotype_id_column", "genotype_id");
        values.put("expression_id_column", "expression_id");
        values.put("sample_alignment", "strict");
        values.put("genotype_block_rows", "16");
        values.put("expression_block_rows", "16");
        values.put("trait_cache", "disk");
        values.put("cache_dir", root.resolve("cache").toString());
        values.put("checkpoint_dir", root.resolve("checkpoint").toString());
        values.put("num_threads", "1");
        values.put("precision", "fp64");
        values.put("residualization", "cpu");
        if (resume) values.put("resume", "true");
        if (keepCheckpoints) values.put("keep_checkpoints", "true");
		if (maximumTraitPatterns != null)
			values.put("max_trait_patterns", maximumTraitPatterns.toString());
		if (scheduler != null)
			values.put("trait_pattern_scheduler", scheduler);
		if (unestimablePolicy != null)
			values.put("unestimable_trait_patterns", unestimablePolicy);
        QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
        QeQTLAnalysis.profiler = new QeQTLProfiler(false);
        QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        QeQTLAnalysis.residualizationMode = QResidualizationMode.CPU;
        QeQTLAnalysis.simplifyResult = false;
        QeQTLAnalysis.rsqOnly = false;
        System.setProperty(CpuBackend.BLAS_PROPERTY, "java");
        QeQTLAnalysis.mContexts = new GpuContext[] {
            new CpuBackend().discoverGpuDevices().get(0).openContext()
        };
		QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
            expression.toString(), covariates.toString(), new String[] {"Age", "Batch"},
            null, "none", Double.NaN, 0, true, 16, 1, "csv");
    }

    private static String messageChain(Throwable error) {
        StringBuilder result = new StringBuilder();
        for (Throwable current = error; current != null; current = current.getCause())
            result.append(current.getMessage()).append('\n');
        return result.toString();
    }

	private static Map<String, double[]> numericRows(Path output) throws Exception {
		Map<String, double[]> result = new TreeMap<>();
		for (String row : Files.readAllLines(output).subList(1,
			Files.readAllLines(output).size())) {
			String[] fields = row.split(",", -1);
			double[] numeric = new double[fields.length - 2];
			for (int i = 2; i < fields.length; i++) numeric[i - 2] = Double.parseDouble(fields[i]);
			result.put(fields[0] + "/" + fields[1], numeric);
		}
		return result;
	}

	private static List<String> qcRows(Path output) throws Exception {
		return Files.readAllLines(output).stream().skip(1).map(row -> {
			String[] fields = row.split("\t", -1);
			return String.join("\t", java.util.Arrays.copyOf(fields, 11));
		}).toList();
	}

	private static void runVariantPatternAnalysis(Path root, Path output, String scheduler,
		String predictorMissing)
		throws Exception {
		Files.createDirectories(root);
		Path genotype = Path.of("test/resources/variant-reference/genotype.vcf")
			.toAbsolutePath().normalize();
		Path expression = Path.of("test/resources/variant-pattern-reference/expression.csv")
			.toAbsolutePath().normalize();
		Map<String, String> values = new LinkedHashMap<>();
		values.put("ini.path", root.toString() + java.io.File.separator);
		values.put("genotype_file", genotype.toString());
		values.put("expression_file", expression.toString());
		values.put("output_file", output.toString());
		values.put("genotype_format", "vcf");
		values.put("predictor_missing", predictorMissing);
		values.put("trait_missing", "pattern");
		values.put("sample_alignment", "strict");
		values.put("min_mac", "0");
		values.put("genotype_block_rows", "16");
		values.put("expression_block_rows", "32");
		values.put("trait_cache", "disk");
		values.put("cache_dir", root.resolve("cache").toString());
		values.put("checkpoint_dir", root.resolve("checkpoint").toString());
		values.put("variant_qc_output", root.resolve("variants.tsv").toString());
		values.put("num_threads", "1");
		values.put("precision", "fp64");
		values.put("residualization", "cpu");
		values.put("trait_pattern_scheduler", scheduler);
		QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
		QeQTLAnalysis.profiler = new QeQTLProfiler(false);
		QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
		QeQTLAnalysis.residualizationMode = QResidualizationMode.CPU;
		QeQTLAnalysis.simplifyResult = false;
		QeQTLAnalysis.rsqOnly = false;
		System.setProperty(CpuBackend.BLAS_PROPERTY, "java");
		QeQTLAnalysis.mContexts = new GpuContext[] {
			new CpuBackend().discoverGpuDevices().get(0).openContext()
		};
		QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
			expression.toString(), null, null, null, "none", Double.NaN, 0, true,
			0, 1, "vcf");
	}
}
