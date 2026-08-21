/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QeQTLCommandLineTest {
    @Test
    void legacyIniCanBeOverriddenByArguments(@TempDir Path directory) throws Exception {
        Path ini = directory.resolve("analysis.ini");
        Files.writeString(ini, "genotype_file = genotype.csv\n"
            + "expression_file = expression.csv\noutput_file = result.csv\nnum_threads = 2\n");

        QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
            "--config", ini.toString(), "--threads", "7", "--genotype-block-rows", "64",
            "--fixed-covariates", "Age,Batch", "--gpu-backend", "cuda", "--precision", "fp32",
			"--residualization", "gpu",
            "--cache-dir", directory.resolve("cache").toString(),
            "--checkpoint-dir", directory.resolve("checkpoint").toString(),
            "--profile-output", directory.resolve("profile.csv").toString(),
            "--rebuild-cache", "--resume", "--keep-checkpoints", "--profile" });
        assertEquals(7, result.config().getNumThreads());
        assertEquals(64, result.config().getGenotypeBlockRows());
        assertArrayEquals(new String[] { "Age", "Batch" }, result.config().getFixedCovariates());
        assertEquals("cuda", result.gpuBackend());
        assertEquals(gov.nih.gpu.GpuPrecision.FP32, result.config().getGpuPrecision());
		assertEquals(QResidualizationMode.GPU, result.config().getResidualizationMode());
        assertTrue(result.config().getGenotypeFilename().startsWith(directory.toString()));
        assertTrue(result.config().getCacheDirectory().startsWith(directory.toString()));
        assertTrue(result.config().getCheckpointDirectory().startsWith(directory.toString()));
        assertTrue(result.config().getRebuildCache());
        assertTrue(result.config().getResume());
        assertTrue(result.config().getKeepCheckpoints());
        assertTrue(result.config().getProfile());
        assertTrue(result.config().getProfileOutputFilename().startsWith(directory.toString()));
    }

    @Test
    void directArgumentsDoNotRequireAnIniFile(@TempDir Path directory) throws Exception {
        QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
            "--genotype", directory.resolve("g.csv").toString(),
            "--expression", directory.resolve("e.csv").toString(),
            "--output", directory.resolve("o.csv").toString(),
            "--threshold", "pval", "1e-4", "--validate-only" });
        assertEquals("pval", result.config().getThresholdType());
        assertEquals(1e-4, result.config().getThresholdValue());
        assertEquals("auto", result.config().getGenotypeFileFormat());
        assertTrue(result.config().getValidateOnly());
        assertEquals(0, result.config().getNumThreads());
		assertEquals(QResidualizationMode.AUTO, result.config().getResidualizationMode());
    }

    @Test
    void variantQcArgumentsAreParsed(@TempDir Path directory) throws Exception {
        Path qc = directory.resolve("variants.tsv");
        QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
            "--genotype", directory.resolve("g.vcf.gz").toString(),
            "--expression", directory.resolve("e.csv").toString(),
            "--output", directory.resolve("o.csv").toString(),
            "--genotype-format", "vcf.gz", "--genotype-field", "GT",
            "--genotype-missing", "mean", "--multiallelic", "error",
            "--min-maf", "0.01", "--min-mac", "5", "--variant-qc-output", qc.toString() });
        assertEquals("vcf.gz", result.config().getGenotypeFileFormat());
        assertEquals("GT", result.config().getGenotypeField());
        assertEquals("mean", result.config().getGenotypeMissingPolicy());
        assertEquals("error", result.config().getMultiallelicPolicy());
        assertEquals(0.01, result.config().getMinimumMaf());
        assertEquals(5, result.config().getMinimumMac());
        assertEquals(qc.toAbsolutePath().normalize().toString(), result.config().getVariantQcOutputFilename());
    }

	@Test
	void explicitMatrixTypesAndMissingnessPoliciesAreParsed(@TempDir Path directory) throws Exception {
		Path qc = directory.resolve("missingness.tsv");
		QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
			"--predictor", directory.resolve("methylation.csv").toString(),
			"--traits", directory.resolve("proteomics.csv").toString(),
			"--output", directory.resolve("result.csv").toString(),
			"--predictor-type", "methylation", "--trait-type", "proteomics",
			"--predictor-missing", "mean", "--trait-missing", "pattern",
			"--predictor-flanks", "2", "--covariate-missing", "complete-samples",
			"--missingness-qc-output", qc.toString(), "--inspect-missingness" });
		assertEquals(QDataType.METHYLATION, result.config().getPredictorDataType());
		assertEquals(QDataType.PROTEOMICS, result.config().getTraitDataType());
		assertEquals(QMissingValuePolicy.MEAN, result.config().getPredictorMissingPolicy());
		assertEquals(QMissingValuePolicy.PATTERN, result.config().getTraitMissingPolicy());
		assertEquals(2, result.config().getPredictorFlankCount());
		assertTrue(result.config().getInspectMissingness());
		assertEquals(qc.toAbsolutePath().normalize().toString(),
			result.config().getMissingnessQcOutputFilename());
	}

	@Test
	void genericMatrixNamesRequireExplicitBiologicalTypes(@TempDir Path directory) {
		assertThrows(IllegalArgumentException.class, () -> QeQTLCommandLine.parse(new String[] {
			"--predictor", directory.resolve("p.csv").toString(),
			"--expression", directory.resolve("e.csv").toString(),
			"--output", directory.resolve("o.csv").toString() }));
		assertThrows(IllegalArgumentException.class, () -> QeQTLCommandLine.parse(new String[] {
			"--genotype", directory.resolve("g.csv").toString(),
			"--traits", directory.resolve("t.csv").toString(),
			"--output", directory.resolve("o.csv").toString() }));
	}
}
