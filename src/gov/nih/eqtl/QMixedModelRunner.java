/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.*;
import java.util.*;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.gpu.*;

/** One fitted trait at a time, bounded predictor tiles, atomic trait checkpoints. */
public final class QMixedModelRunner {
    private QMixedModelRunner() { }
    private static final String HEADER = "SNP,Probe,N,DF,Effect,SE,Statistic,log10P,Test,Status,VarianceComponents,REMLLogLik,Iterations,Singular";

    public static void run(QeQTLAnalysisConfig config, QMixedModel model,
        QMatrixRowSource genotypes, int[] genotypeColumns, QMatrixRowSource traits,
        int[] traitColumns, GpuContext[] contexts) throws Exception {
        String sourceStamp = sourceStamp(config);
        String execution = execution(contexts);
        String signature = signature(config, genotypeColumns, traitColumns, execution);
        requireUnchanged(config,sourceStamp);
        Path output = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
        Path directory = config.getCheckpointDirectory() == null
            ? Path.of(output + ".mixed-checkpoint") : Path.of(config.getCheckpointDirectory());
        int traitCount = Math.toIntExact(traits.metadata().rowCount());
        boolean score = config.getMixedModel().equals("score");
        int blockRows = config.getGenotypeBlockRows() > 0 ? config.getGenotypeBlockRows() : 256;
        QAnalysisCheckpoint checkpoint = QAnalysisCheckpoint.open(directory, signature, traitCount,
            config.getResume(), config.getKeepCheckpoints());
        System.out.println("Gaussian mixed model: " + (score ? "REML null + exact score variance (GENESIS-style)"
            : "per-predictor REML refit + " + config.getRemlDf() + " Wald inference")
            + "; traits=" + traitCount + "; predictor tile=" + blockRows);
        try (GpuContextPool pool = new GpuContextPool(contexts);
             QMatrixRowSource.BlockReader reader = traits.open(traitColumns)) {
            if (score && pool.getNumContexts() > 0) {
                String kernel = "#define BLOCK_SIZE 16\n#define DATATYPE double\n#define N_MIN_1 1\n"
                    + "#if defined(cl_khr_fp64)\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
                    + "#elif defined(cl_amd_fp64)\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n"
                    + QeQTLAnalysis.eqtlReal;
                for (GpuContext context : pool.getAllContexts()) context.compileKernel(kernel, "eqtlReal", GpuPrecision.FP64);
            }
            long started = System.nanoTime(), comparisons = 0;
            for (int trait = 0; trait < traitCount; trait++) {
                var row = reader.readBlock(1);
                if (row == null || row.rowCount() != 1 || row.rowOffset() != trait)
                    throw new IOException("Mixed-model trait source changed during analysis");
                if (checkpoint.isComplete(trait)) continue;
                requireUnchanged(config,sourceStamp);
                double[] y = row.values()[0];
                QMixedModel.NullState state = model.fitNull(y);
                System.out.println("REML null " + row.rowIds()[0] + ": variances="
                    + Arrays.toString(state.fit().varianceComponents()) + "; logLik="
                    + state.fit().restrictedLogLikelihood() + "; iterations=" + state.fit().iterations());
                double[] packedProjection = score && pool.getNumContexts() > 0 ? packProjection(state, y.length) : null;
                checkpoint.writeBlock(trait, writer -> {
                    try (var genotypeReader = genotypes.open(genotypeColumns)) {
                        QMatrixRowSource.Block block;
                        long processed = 0, lastReport = System.nanoTime();
                        while ((block = genotypeReader.readBlock(blockRows)) != null) {
                            double[][] centered = new double[block.rowCount()][];
                            for (int g = 0; g < centered.length; g++) centered[g] = state.center(block.values()[g]);
                            double[][] projected = score ? project(state, centered, packedProjection, pool) : null;
                            for (int g = 0; g < block.rowCount(); g++) {
                                var result = score ? state.score(centered[g], projected[g]) : model.refit(y, block.values()[g], state);
                                if (config.getThresholdType().equals("pval")
                                    && result.log10P() > Math.log10(config.getThresholdValue())) continue;
                                writer.write(csv(block.rowIds()[g]) + "," + csv(row.rowIds()[0]) + "," + y.length
                                    + "," + result.df() + "," + result.beta() + "," + result.standardError()
                                    + "," + result.statistic() + "," + result.log10P() + "," + result.test()
                                    + "," + result.status() + "," + csv(variances(result.fit()))
                                    + "," + result.fit().restrictedLogLikelihood() + "," + result.fit().iterations()
                                    + "," + result.fit().isSingular(1e-6) + "\n");
                            }
                            processed += block.rowCount();
                            if (System.nanoTime()-lastReport > 5_000_000_000L) {
                                System.out.println("Mixed-model trait " + row.rowIds()[0] + ": " + processed + " comparisons");
                                lastReport = System.nanoTime();
                            }
                        }
                    }
                    requireUnchanged(config,sourceStamp);
                });
                comparisons += genotypes.metadata().rowCount();
                System.out.println("Mixed-model progress " + (trait+1) + "/" + traitCount + " traits; "
                    + comparisons + " current-run comparisons in " + ((System.nanoTime()-started)/1e9) + " seconds");
            }
            if (reader.readBlock(1) != null) throw new IOException("Mixed-model trait source grew during analysis");
        }
        if (!signature.equals(signature(config, genotypeColumns, traitColumns, execution)))
            throw new IOException("Mixed-model input changed during analysis; output was not assembled");
        checkpoint.assemble(output, HEADER);
    }

    static double[] packProjection(QMixedModel.NullState state, int n) {
        int width = round(n,64), rows = round(n,16);
        double[] result = new double[Math.multiplyExact(rows,width)], p = state.projection();
        for (int i = 0; i < n; i++) System.arraycopy(p,i*n,result,i*width,n);
        return result;
    }
    static double[][] project(QMixedModel.NullState state, double[][] variants,
        double[] projection, GpuContextPool pool) {
        double[][] result = new double[variants.length][];
        if (projection == null) {
            for (int g = 0; g < variants.length; g++) result[g] = state.project(variants[g]);
            return result;
        }
        int n = variants[0].length, width = round(n,64), columns = round(variants.length,16), rows = round(n,16);
        double[] packed = new double[Math.multiplyExact(width,columns)];
        for (int g = 0; g < variants.length; g++) for (int i = 0; i < n; i++) packed[i*columns+g] = variants[g][i];
        GpuContext context = pool.reserveContext();
        double[] product;
        try {
            product = context.executeDoubleKernel(projection, packed, Math.multiplyExact(rows,columns),
                17L*64*8, width, columns, new long[]{columns,rows}, new long[]{16,16});
        } finally { pool.releaseContext(context); }
        for (int g = 0; g < variants.length; g++) {
            result[g] = new double[n];
            for (int i = 0; i < n; i++) {
                double value = product[i*columns+g];
                if (!Double.isFinite(value)) throw new IllegalArgumentException("Backend returned non-finite mixed-model projection");
                result[g][i] = value;
            }
        }
        return result;
    }
    private static int round(int n, int multiple) { return Math.multiplyExact(Math.addExact(n,multiple-1)/multiple,multiple); }
    private static String csv(String value) { return "\"" + value.replace("\"", "\"\"") + "\""; }
    private static String variances(org.jlinalg.reml.RemlResult fit) {
        StringJoiner joiner = new StringJoiner(";"); double[] values = fit.varianceComponents();
        for (int i = 0; i < values.length; i++) joiner.add(fit.componentNames().get(i)+"="+values[i]);
        return joiner.toString();
    }
    private static String signature(QeQTLAnalysisConfig config, int[] genotypeColumns, int[] traitColumns, String execution) throws Exception {
        MessageDigest digest = MessageDigest.getInstance("SHA-256");
        digest.update("mixed-v1-jlinalg-0.3.0-jdistlib-0.10.2-fp64".getBytes(StandardCharsets.UTF_8));
        digest.update(execution.getBytes(StandardCharsets.UTF_8));
        Map<String,String> settings = new TreeMap<>(config.asMap());
        for (String key : List.of("resume", "keep_checkpoints", "rebuild_cache", "profile", "profile_output")) settings.remove(key);
        for (var entry : settings.entrySet()) digest.update((entry.getKey()+"="+entry.getValue()+"\n").getBytes(StandardCharsets.UTF_8));
        digest.update(Arrays.toString(genotypeColumns).getBytes(StandardCharsets.UTF_8));
        digest.update(Arrays.toString(traitColumns).getBytes(StandardCharsets.UTF_8));
        for (String filename : new String[]{config.getGenotypeFilename(), config.getExpressionFilename(),
                config.getCovariateFilename(), config.getPedigreeFilename(), config.getVariantIndexFilename(), config.getRegionsFilename()}) {
            if (filename == null) continue;
            Path path = Path.of(filename).toAbsolutePath().normalize();
            digest.update(path.toString().getBytes(StandardCharsets.UTF_8));
            try (var input = new DigestInputStream(Files.newInputStream(path),digest)) { input.transferTo(OutputStream.nullOutputStream()); }
        }
        return HexFormat.of().formatHex(digest.digest());
    }
    private static String execution(GpuContext[] contexts) {
        StringBuilder value = new StringBuilder("cpu-blas=").append(System.getProperty("eqtl.cpu.blas","auto"))
            .append(";threads=").append(System.getProperty("eqtl.cpu.threads","auto"));
        if (contexts != null) for (GpuContext context : contexts) {
            GpuDevice device = context.getDevice();
            value.append(';').append(device.getBackendName()).append(';').append(device.getName())
                .append(';').append(device.getDriverVersion()).append(';').append(device.getComputeApiVersion());
        }
        return value.toString();
    }
    private static void requireUnchanged(QeQTLAnalysisConfig config, String stamp) throws IOException {
        if (!stamp.equals(sourceStamp(config)))
            throw new IOException("Mixed-model input changed before checkpoint completion");
    }
    private static String sourceStamp(QeQTLAnalysisConfig config) throws IOException {
        StringBuilder stamp = new StringBuilder();
        for (String filename : new String[]{config.getGenotypeFilename(),config.getExpressionFilename(),
                config.getCovariateFilename(),config.getPedigreeFilename(),config.getVariantIndexFilename(),config.getRegionsFilename()}) {
            if (filename == null) continue;
            var attrs = Files.readAttributes(Path.of(filename),java.nio.file.attribute.BasicFileAttributes.class);
            stamp.append(filename).append(':').append(attrs.fileKey()).append(':').append(attrs.size())
                .append(':').append(attrs.lastModifiedTime()).append('\n');
        }
        return stamp.toString();
    }
}
