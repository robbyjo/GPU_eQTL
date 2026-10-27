/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Proxy;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QDelimitedMatrixSource;
import gov.nih.eqtl.io.QMissingnessScan;
import gov.nih.eqtl.io.QPreparedMatrix;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;
import gov.nih.gpu.GpuPatternStatisticsResult;
import gov.nih.gpu.cpu.CpuBackend;

class QGenotypeOuterPatternJobTest {
    private static final double NOISE = 0.0003;
    private static final double REPORTED_ABSOLUTE = 1.000000061191563;

    @Test
    void failedBlockResumesWithScalarVerifiedStatisticsAndDeterministicOutput(
        @TempDir Path directory) throws Exception {
        Path resumed = directory.resolve("resumed");
        RuntimeException failure = assertThrows(RuntimeException.class,
            () -> runBlock(resumed, 1.01, false, 0));
        assertTrue(failure.getCause().getMessage().contains("predictor 'rs1', trait 'positive'"));
        assertTrue(failure.getCause().getMessage().contains("pattern 0, N=8, DF=5"));
        assertFalse(Files.exists(resumed.resolve("checkpoint/block-00000000.part")));
        assertFalse(Files.exists(resumed.resolve("qc/block-00000000.qc")));

        runBlock(resumed, REPORTED_ABSOLUTE, true, 0);
        Path clean = directory.resolve("clean");
        runBlock(clean, REPORTED_ABSOLUTE, false, 0);
        byte[] rows = Files.readAllBytes(resumed.resolve("checkpoint/block-00000000.part"));
        assertArrayEquals(Files.readAllBytes(clean.resolve("checkpoint/block-00000000.part")), rows);
        String[] lines = new String(rows, java.nio.charset.StandardCharsets.UTF_8).strip().split("\\R");
        assertEquals(2, lines.length);
        for (int i = 0; i < lines.length; i++) {
            double sign = i == 0 ? 1 : -1;
            String[] fields = lines[i].split(",");
            assertEquals("rs1", fields[0]);
            assertEquals(i == 0 ? "positive" : "negative", fields[1]);
            double correlation = sign / Math.sqrt(1 + NOISE * NOISE);
            var expected = QeQTLStatistics.calculate(correlation,
                Math.sqrt(6 * (1 + NOISE * NOISE) / 7), Math.sqrt(6.0 / 7), 5, 0);
            assertEquals(expected.rSquared(), Double.parseDouble(fields[2]), 3e-15);
            assertEquals(expected.effect(), Double.parseDouble(fields[3]), 3e-13);
            assertEquals(expected.tStatistic(), Double.parseDouble(fields[4]),
                Math.abs(expected.tStatistic()) * 1e-8);
            assertEquals(expected.log10P(), Double.parseDouble(fields[5]), 1e-7);
            assertEquals("8", fields[6]);
            assertEquals("5", fields[7]);
        }
    }

    @Test
    void thresholdUsesVerifiedCorrelationRatherThanAnArtificialPerfectFit(
        @TempDir Path directory) throws Exception {
        runBlock(directory, REPORTED_ABSOLUTE, false, 0.99999999);
        assertEquals("", Files.readString(directory.resolve("checkpoint/block-00000000.part")));
    }

    private static void runBlock(Path root, double targetCorrelation, boolean resume,
        double threshold) throws Exception {
        Files.createDirectories(root);
        Path traitFile = root.resolve("traits.csv");
        if (!Files.exists(traitFile))
            Files.writeString(traitFile, ",S0,S1,S2,S3,S4,S5,S6,S7\n"
                + "positive,1,-2,1,0,0,0.0003,-0.0006,0.0003\n"
                + "negative,-1,2,-1,0,0,-0.0003,0.0006,-0.0003\n");
        double[][] design = new double[8][2];
        double[] genotype = new double[8];
        double[] u = {1, -2, 1, 0, 0, 0, 0, 0};
        for (int i = 0; i < 8; i++) {
            design[i] = new double[] {1, i - 3};
            genotype[i] = 10000 + 2 * (i - 3) + u[i];
        }
        QTraitPatternModelSet models;
        PreparedBlock block;
        {
            var source = new QDelimitedMatrixSource(traitFile, new char[] {','}, "#");
            models = QTraitPatternModelSet.create(QMissingnessScan.scan("trait", source, null),
                design, 0, "rsq", threshold, QUnestimableTraitPolicy.ERROR);
            try (var reader = source.open(null)) {
                var raw = reader.readBlock(16);
                PreparedBlock positive = models.prepareTrait(0, 0, "positive", raw.values()[0]);
                PreparedBlock negative = models.prepareTrait(1, 1, "negative", raw.values()[1]);
                block = new PreparedBlock(0, new String[] {"positive", "negative"},
                    new double[][] {positive.values()[0], negative.values()[0]},
                    new double[] {positive.standardDeviations()[0], negative.standardDeviations()[0]});
            }
        }
        QPreparedMatrix traits = new QPreparedMatrix() {
            public String signature() { return "fixture"; }
            public int sampleCount() { return 8; }
            public long rowCount() { return 2; }
            public PreparedBlock readBlock(long offset, int maximumRows) { return block; }
        };
        var checkpoint = QAnalysisCheckpoint.open(root.resolve("checkpoint"), "fixture", 1, resume, true);
        var qc = QGenotypeOuterPatternQcCheckpoint.open(root.resolve("qc"), "fixture", 1, 1, resume, true);
        String previous = System.getProperty(CpuBackend.BLAS_PROPERTY);
        boolean previousRsq = QeQTLAnalysis.rsqOnly;
        boolean previousSimplify = QeQTLAnalysis.simplifyResult;
        System.setProperty(CpuBackend.BLAS_PROPERTY, "java");
        QeQTLAnalysis.rsqOnly = false;
        QeQTLAnalysis.simplifyResult = false;
        try (GpuContext delegate = new CpuBackend().discoverGpuDevices().get(0).openContext()) {
            delegate.compileKernel("#define DATATYPE double\n#define N_MIN_1 1\n"
                + QeQTLAnalysis.eqtlReal, "eqtlReal");
            // Simulate cancellation in compact variance state while retaining real
            // production Java matrix products and the actual checkpoint/job path.
            GpuContext context = (GpuContext) Proxy.newProxyInstance(GpuContext.class.getClassLoader(),
                new Class<?>[] {GpuContext.class}, (proxy, method, args) -> {
                    try {
                        Object result = method.invoke(delegate, args);
                        if (method.getName().equals("computePatternStatisticsDouble")) {
                            var statistics = (GpuPatternStatisticsResult) result;
                            double[] values = statistics.values().clone();
                            double trueCorrelation = 1 / Math.sqrt(1 + NOISE * NOISE);
                            values[GpuPatternStatisticsResult.RESIDUAL_SUM_SQUARES] =
                                6 * Math.pow(trueCorrelation / targetCorrelation, 2);
                            return new GpuPatternStatisticsResult(1, 1, values);
                        }
                        return result;
                    } catch (InvocationTargetException e) {
                        throw e.getCause();
                    }
                });
            try (GpuContextPool pool = GpuContextPool.borrowed(new GpuContext[] {context});
                QAnalysisProgress progress = QAnalysisProgress.dynamic("fixture", 0, 2)) {
                new QGenotypeOuterPatternJob(new QGenotypeOuterPatternJob.QMatrixBlock(
                    0, new String[] {"rs1"}, new double[][] {genotype}), traits, models,
                    checkpoint, 0, qc, pool, 16, 16, QMissingValuePolicy.MEAN, 0,
                    new QeQTLProfiler(false), progress, 0, 0, false).run();
            }
        } finally {
            if (previous == null) System.clearProperty(CpuBackend.BLAS_PROPERTY);
            else System.setProperty(CpuBackend.BLAS_PROPERTY, previous);
            QeQTLAnalysis.rsqOnly = previousRsq;
            QeQTLAnalysis.simplifyResult = previousSimplify;
        }
    }
}
