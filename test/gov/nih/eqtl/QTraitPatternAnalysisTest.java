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

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuPrecision;
import gov.nih.gpu.cpu.CpuBackend;

class QTraitPatternAnalysisTest {
    private static final String FAIL_PROPERTY = "eqtl.test.trait.pattern.fail.after";

    @AfterEach
    void restoreState() {
        System.clearProperty(FAIL_PROPERTY);
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
			() -> runAnalysis(root, directory.resolve("limited.csv"), false, false, 1));
		assertTrue(messageChain(error).contains("2 patterns exceed --max-trait-patterns 1"));
		assertTrue(!Files.exists(root.resolve("cache").resolve("trait-patterns")));
		assertTrue(!Files.exists(root.resolve("checkpoint")));
	}

    private static void runAnalysis(Path root, Path output, boolean resume,
        boolean keepCheckpoints) throws Exception {
		runAnalysis(root, output, resume, keepCheckpoints, null);
	}

	private static void runAnalysis(Path root, Path output, boolean resume,
		boolean keepCheckpoints, Integer maximumTraitPatterns) throws Exception {
        Files.createDirectories(root);
        Path genotype = Path.of("test/resources/eqtl-reference/genotype.csv")
            .toAbsolutePath().normalize();
        Path expression = Path.of("test/resources/missing-reference/expression.csv")
            .toAbsolutePath().normalize();
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
}
