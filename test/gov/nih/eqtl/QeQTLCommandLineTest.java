/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

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
            "--fixed-covariates", "Age,Batch", "--gpu-backend", "cuda" });
        assertEquals(7, result.config().getNumThreads());
        assertEquals(64, result.config().getGenotypeBlockRows());
        assertArrayEquals(new String[] { "Age", "Batch" }, result.config().getFixedCovariates());
        assertEquals("cuda", result.gpuBackend());
        assertTrue(result.config().getGenotypeFilename().startsWith(directory.toString()));
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
        assertEquals("csv", result.config().getGenotypeFileFormat());
        assertTrue(result.config().getValidateOnly());
    }
}
