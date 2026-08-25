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
			"--max-trait-patterns", "17",
			"--trait-pattern-scheduler", "genotype",
			"--unestimable-trait-patterns", "skip",
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
		assertEquals(17, result.config().getMaximumTraitPatterns());
		assertEquals(QTraitPatternScheduler.GENOTYPE,
			result.config().getTraitPatternScheduler());
		assertEquals(QUnestimableTraitPolicy.SKIP,
			result.config().getUnestimableTraitPolicy());
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
        assertEquals(QeQTLAnalysisConfig.DEFAULT_MINIMUM_MAC,
            result.config().getMinimumMac());
		assertEquals(QeQTLAnalysisConfig.DEFAULT_MAXIMUM_TRAIT_PATTERNS,
			result.config().getMaximumTraitPatterns());
		assertEquals(QTraitPatternScheduler.AUTO,
			result.config().getTraitPatternScheduler());
		assertEquals(QUnestimableTraitPolicy.ERROR,
			result.config().getUnestimableTraitPolicy());
    }

    @Test
    void preprocessOnlyAndMacDisableAreParsedWithoutAssociationOutput(@TempDir Path directory)
        throws Exception {
        QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
            "--genotype", directory.resolve("g.vcf.gz").toString(),
            "--expression", directory.resolve("e.csv").toString(),
            "--genotype-format", "vcf", "--preprocess-only", "--min-mac", "0",
            "--cache-dir", directory.resolve("cache").toString() });
        assertTrue(result.config().getPreprocessOnly());
        assertEquals(0, result.config().getMinimumMac());
        assertEquals(null, result.config().getOutputFilename());
    }

	@Test
	void backendAliasSelectsCpu(@TempDir Path directory) throws Exception {
		QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
			"--genotype", directory.resolve("g.csv").toString(),
			"--expression", directory.resolve("e.csv").toString(),
			"--output", directory.resolve("o.csv").toString(),
			"--backend", "cpu", "--printbackendinfo" });
		assertEquals("cpu", result.gpuBackend());
		assertTrue(result.printGpuInfo());
	}

    @Test
    void variantQcArgumentsAreParsed(@TempDir Path directory) throws Exception {
        Path qc = directory.resolve("variants.tsv");
        Path checkpoint = directory.resolve("variant-checkpoint");
        Path index = directory.resolve("g.vcf.gz.tbi");
        Path regions = directory.resolve("regions.tsv");
        QeQTLCommandLine.Result result = QeQTLCommandLine.parse(new String[] {
            "--genotype", directory.resolve("g.vcf.gz").toString(),
            "--expression", directory.resolve("e.csv").toString(),
            "--output", directory.resolve("o.csv").toString(),
            "--genotype-format", "vcf.gz", "--genotype-field", "GT",
            "--genotype-missing", "mean", "--multiallelic", "error",
            "--min-maf", "0.01", "--min-mac", "5", "--variant-qc-output", qc.toString(),
            "--variant-qc-threads", "3", "--variant-qc-checkpoint", checkpoint.toString(),
            "--variant-index", index.toString(), "--region", "geneA=1:100-200",
            "--region", "geneB=2:300-400", "--regions-file", regions.toString(),
            "--region-coordinates", "bed", "--frequency-scope", "pattern" });
        assertEquals("vcf.gz", result.config().getGenotypeFileFormat());
        assertEquals("GT", result.config().getGenotypeField());
        assertEquals("mean", result.config().getGenotypeMissingPolicy());
        assertEquals("error", result.config().getMultiallelicPolicy());
        assertEquals(0.01, result.config().getMinimumMaf());
        assertEquals(5, result.config().getMinimumMac());
        assertEquals(3, result.config().getVariantQcThreads());
        assertEquals(qc.toAbsolutePath().normalize().toString(), result.config().getVariantQcOutputFilename());
        assertEquals(checkpoint.toAbsolutePath().normalize().toString(),
            result.config().getVariantQcCheckpointDirectory());
        assertEquals(index.toAbsolutePath().normalize().toString(),
            result.config().getVariantIndexFilename());
        assertEquals("geneA=1:100-200;geneB=2:300-400", result.config().getRegions());
        assertEquals(regions.toAbsolutePath().normalize().toString(),
            result.config().getRegionsFilename());
        assertEquals("bed", result.config().getRegionCoordinates());
        assertEquals("pattern", result.config().getFrequencyScope());
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
			"--sample-alignment", "covariate-subset",
			"--predictor-id-strip-prefix", "P", "--trait-id-strip-prefix", "X",
			"--trait-cache", "memory",
			"--missingness-qc-output", qc.toString(), "--inspect-missingness" });
		assertEquals(QDataType.METHYLATION, result.config().getPredictorDataType());
		assertEquals(QDataType.PROTEOMICS, result.config().getTraitDataType());
		assertEquals(QMissingValuePolicy.MEAN, result.config().getPredictorMissingPolicy());
		assertEquals(QMissingValuePolicy.PATTERN, result.config().getTraitMissingPolicy());
		assertEquals(2, result.config().getPredictorFlankCount());
		assertEquals(gov.nih.eqtl.io.QSampleAlignmentPolicy.COVARIATE_SUBSET,
			result.config().getSampleAlignmentPolicy());
		assertEquals("P", result.config().getGenotypeIdStripPrefix());
		assertEquals("X", result.config().getExpressionIdStripPrefix());
		assertEquals(QTraitCacheMode.MEMORY, result.config().getTraitCacheMode());
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
