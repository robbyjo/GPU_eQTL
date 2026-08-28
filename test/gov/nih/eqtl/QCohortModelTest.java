/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuPrecision;
import gov.nih.gpu.cpu.CpuBackend;

class QCohortModelTest {
    @Test
    void buildsBlockDesignAndAppliesTheSameWhiteningToBothSides(@TempDir Path root)
        throws Exception {
        Path covariates = writeCovariates(root);
        Path definitions = writeDefinitions(root);
        QCovariateTable table = QCovariateTable.load(covariates, new char[] {','}, "#");
        boolean[] selected = new boolean[8];
        java.util.Arrays.fill(selected, true);
        String[] ids = {"G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"};
        QCohortModel model = QCohortModel.create(QCohortModel.Definitions.load(definitions),
            table, "cohort", selected, ids, null, null);

        assertEquals(4, model.design()[0].length);
        assertTrue(model.hasWhitening());
        assertArrayEquals(new String[] {"A:(Intercept)", "A:Age",
            "B:(Intercept)", "B:Age"}, model.columnNames());
        double[] predictor = {0, 1, 2, 3, 0, 1, 2, 3};
        double[] trait = {3, 1, 4, 2, 2, 5, 1, 6};
        assertEquals(QCohortModel.stackedDot(predictor, trait),
            QCohortModel.blockwiseDot(predictor, trait,
                new int[][] {{0, 1, 2, 3}, {4, 5, 6, 7}}), 0);

        QeQTLPreprocessor.Residualizer whitening = model.wrap(null);
        double[][] left = whitening.residualize(new double[][] {predictor},
            orthogonal(model.design()), "predictor");
        double[][] right = whitening.residualize(new double[][] {trait},
            orthogonal(model.design()), "trait");
        assertEquals(1, left.length);
        assertEquals(8, right[0].length);
        assertTrue(whitening.cacheSignatureTag().contains(model.signature()));
    }

    @Test
    void streamedAndResidentCohortAnalysesAreByteIdentical(@TempDir Path root)
        throws Exception {
        Path covariates = writeCovariates(root);
        Path definitions = writeDefinitions(root);
        Path resident = root.resolve("resident.csv");
        Path streamed = root.resolve("streamed.csv");
        run(root.resolve("resident-run"), covariates, definitions, resident, false);
        run(root.resolve("streamed-run"), covariates, definitions, streamed, true);
        assertArrayEquals(Files.readAllBytes(resident), Files.readAllBytes(streamed));
        assertTrue(Files.readString(Path.of(streamed + ".cohorts.tsv"))
            .contains("B\t4\t2\tcompound-symmetry"));
    }

    private static double[][] orthogonal(double[][] design) {
        return new gov.nih.jama.QRDecomposition(design).getQ().getArray();
    }

    private static Path writeCovariates(Path root) throws Exception {
        Path path = root.resolve("covariates.csv");
        Files.writeString(path, "subject_id,genotype_id,expression_id,Age,cohort,repeat_id\n"
            + "S1,G1,T1,20,A,A1\nS2,G2,T2,22,A,A2\n"
            + "S3,G3,T3,25,A,A3\nS4,G4,T4,27,A,A4\n"
            + "S5,G5,T5,30,B,B1\nS6,G6,T6,33,B,B1\n"
            + "S7,G7,T7,36,B,B2\nS8,G8,T8,40,B,B2\n");
        return path;
    }

    private static Path writeDefinitions(Path root) throws Exception {
        Path path = root.resolve("cohorts.tsv");
        Files.writeString(path, "COHORT\tFIXED_COVARIATES\tSUBJECT_COLUMN\tREPEATED_POLICY"
            + "\tWITHIN_SUBJECT_CORRELATION\n"
            + "A\tAge\t\tindependent\t\n"
            + "B\tAge\trepeat_id\tcompound-symmetry\t0.25\n");
        return path;
    }

    private static void run(Path runRoot, Path covariates, Path definitions,
        Path output, boolean stream) throws Exception {
        Files.createDirectories(runRoot);
        Path genotype = Path.of("test/resources/eqtl-reference/genotype.csv")
            .toAbsolutePath().normalize();
        Path expression = Path.of("test/resources/eqtl-reference/expression.csv")
            .toAbsolutePath().normalize();
        Map<String, String> values = new LinkedHashMap<>();
        values.put("ini.path", runRoot.toString() + java.io.File.separator);
        values.put("genotype_file", genotype.toString());
        values.put("expression_file", expression.toString());
        values.put("covariate_file", covariates.toString());
        values.put("cohort_model", definitions.toString());
        values.put("cohort_column", "cohort");
        values.put("output_file", output.toString());
        values.put("genotype_format", "csv");
        values.put("trait_missing", "error");
        values.put("predictor_missing", "error");
        values.put("covariate_missing", "error");
        values.put("genotype_id_column", "genotype_id");
        values.put("expression_id_column", "expression_id");
        values.put("sample_alignment", "strict");
        values.put("cache_dir", runRoot.resolve("cache").toString());
		values.put("block_size", "16");
        values.put("num_threads", "1");
        values.put("precision", "fp64");
        values.put("residualization", "cpu");
        if (stream) {
            values.put("genotype_block_rows", "16");
            values.put("expression_block_rows", "16");
        }
        QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
        QeQTLAnalysis.profiler = new QeQTLProfiler(false);
        QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        QeQTLAnalysis.residualizationMode = QResidualizationMode.CPU;
        QeQTLAnalysis.simplifyResult = false;
        QeQTLAnalysis.rsqOnly = false;
        System.setProperty(CpuBackend.BLAS_PROPERTY, "java");
        GpuContext context = new CpuBackend().discoverGpuDevices().get(0).openContext();
        QeQTLAnalysis.mContexts = new GpuContext[] {context};
        try {
            QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
                expression.toString(), covariates.toString(), null, null, "none",
                Double.NaN, 0, true, 16, 1, "csv");
        } finally {
            context.close();
        }
    }
}
