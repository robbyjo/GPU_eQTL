/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import java.io.IOException;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.List;
import java.util.Locale;
import java.util.Set;

import gov.nih.eqtl.QAnalysisCheckpoint;
import gov.nih.eqtl.QMissingValuePolicy;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.eqtl.io.QMatrixRowSource.Block;
import gov.nih.eqtl.io.QRawMatrixCache;
import gov.nih.eqtl.settest.QBurdenReference.Analysis;
import gov.nih.eqtl.settest.QBurdenReference.Result;
import gov.nih.eqtl.settest.QBurdenReference.SetAudit;
import gov.nih.eqtl.settest.QBurdenReference.Variant;

/** Bounded deterministic production runner shared by command-line set tests. */
public final class QSetTestRunner {
    private static final int MAXIMUM_AUTOMATIC_SET_BLOCK_SIZE = 256;
    private static final long MINIMUM_AUTOMATIC_BUDGET = 64L * 1024 * 1024;
    private static final long AUTOMATIC_FIXED_RESERVE = 256L * 1024 * 1024;

    public record Options(QSetTestMethod method, Path definitions, Path output, Path auditOutput,
        Path checkpointDirectory, int traitBlockRows, QSetTestPolicy policy,
        boolean resume, boolean keepCheckpoints, String thresholdType, double threshold,
        int setBlockSize, double[] skatORhoGrid, int skatOSimulations, long skatOSeed) {
        public Options {
            if (method == null || !method.isSetTest() || output == null
                || auditOutput == null || checkpointDirectory == null || traitBlockRows <= 0
                || setBlockSize < 0 || policy == null || thresholdType == null)
                throw new IllegalArgumentException("Incomplete set-test runner options");
            skatORhoGrid = skatORhoGrid == null ? null : skatORhoGrid.clone();
            if (method == QSetTestMethod.SKAT_O
                && (skatORhoGrid == null || skatORhoGrid.length == 0 || skatOSimulations < 1))
                throw new IllegalArgumentException("SKAT-O requires rho values and simulations");
            if (method == QSetTestMethod.SKAT_O) validateRhoGrid(skatORhoGrid);
        }
        @Override public double[] skatORhoGrid() {
            return skatORhoGrid == null ? null : skatORhoGrid.clone();
        }
    }

    record SetBlockRecommendation(int blockSize, long heapLimitBytes,
        long tileBudgetBytes, long estimatedTileBytes) { }

    private static void validateRhoGrid(double[] values) {
        if (values[0] != 0 || values[values.length - 1] != 1)
            throw new IllegalArgumentException("SKAT-O rho grid must start at 0 and end at 1");
        for (int i = 0; i < values.length; i++) {
            if (!Double.isFinite(values[i]) || values[i] < 0 || values[i] > 1
                || (i > 0 && values[i] <= values[i - 1]))
                throw new IllegalArgumentException(
                    "SKAT-O rho grid must be strictly increasing within [0,1]");
        }
    }

    private QSetTestRunner() { }

    public static void run(QMatrixRowSource variants, int[] variantColumns,
        QMatrixRowSource traits, int[] traitColumns, double[][] covariateDesign,
        Options options) throws Exception {
        QVariantSetTable definitions = QVariantSetTable.load(options.definitions());
        run(variants, variantColumns, traits, traitColumns, covariateDesign,
            definitions, options);
    }

    public static void run(QMatrixRowSource variants, int[] variantColumns,
        QMatrixRowSource traits, int[] traitColumns, double[][] covariateDesign,
        QVariantSetTable definitions, Options options) throws Exception {
        if (definitions == null) throw new IllegalArgumentException("Variant-set definitions are required");
        int effectiveSetBlockSize = options.setBlockSize();
        if (effectiveSetBlockSize == 0) {
            SetBlockRecommendation recommendation = recommendSetBlockSize(definitions,
                variantColumns.length, options.method(),
                Runtime.getRuntime().maxMemory());
            effectiveSetBlockSize = recommendation.blockSize();
            System.out.println("Automatic set block size = " + effectiveSetBlockSize
                + "; JVM heap limit=" + recommendation.heapLimitBytes()
                + " bytes; set-tile budget=" + recommendation.tileBudgetBytes()
                + " bytes; worst estimated tile=" + recommendation.estimatedTileBytes()
                + " bytes");
        }
        long traitRows = traits.metadata().rowCount();
        int traitBlocks = Math.toIntExact((traitRows + options.traitBlockRows() - 1)
            / options.traitBlockRows());
        int setBlocks = (definitions.sets().size() + effectiveSetBlockSize - 1)
            / effectiveSetBlockSize;
        int totalBlocks = Math.multiplyExact(setBlocks, traitBlocks);
        System.out.println("Set-test schedule: sets=" + definitions.sets().size()
            + ", set_blocks=" + setBlocks + ", set_block_size=" + effectiveSetBlockSize
            + ", traits=" + traitRows + ", trait_blocks=" + traitBlocks
            + ", trait_block_rows=" + options.traitBlockRows());
        long observedHeapPeak = usedHeapBytes();
        String signature = signature(variants, variantColumns, traits, traitColumns,
            covariateDesign, definitions, options, effectiveSetBlockSize);
        QAnalysisCheckpoint checkpoint = QAnalysisCheckpoint.open(
            options.checkpointDirectory(), signature, totalBlocks,
            options.resume(), options.keepCheckpoints());
        List<SetAudit> audits = new ArrayList<>();
        for (int setBlock = 0; setBlock < setBlocks; setBlock++) {
            int setFrom = setBlock * effectiveSetBlockSize;
            QVariantSetTable setTile = definitions.subset(setFrom,
                Math.min(definitions.sets().size(), setFrom + effectiveSetBlockSize));
            List<Variant> selectedVariants = loadSelectedVariants(variants, variantColumns,
                setTile, options.traitBlockRows());
            observedHeapPeak = Math.max(observedHeapPeak, usedHeapBytes());
            List<SetAudit> tileAudit = null;
            try (QMatrixRowSource.BlockReader reader = traits.open(traitColumns)) {
                for (int traitBlock = 0; traitBlock < traitBlocks; traitBlock++) {
                    Block block = reader.readBlock(options.traitBlockRows());
                    if (block == null || block.rowOffset()
                        != (long) traitBlock * options.traitBlockRows())
                        throw new IOException("Trait source ended before set-test trait block "
                            + traitBlock);
                    int blockNumber = setBlock * traitBlocks + traitBlock;
                    if (checkpoint.isComplete(blockNumber)) continue;
                    QSetTestNullModel nullModel = QSetTestNullModel.create(
                        block.rowIds(), block.values(), covariateDesign);
                    if (options.method() == QSetTestMethod.BURDEN) {
                        Analysis analysis = QBurdenReference.analyzeBatched(setTile,
                            selectedVariants, nullModel, options.policy());
                        if (tileAudit == null) tileAudit = analysis.audits();
                        checkpoint.writeBlock(blockNumber,
                            writer -> writeBurdenResults(writer, analysis.results(), options));
                    } else {
                        QKernelSetReference.Analysis analysis = QKernelSetReference.analyze(
                            setTile, selectedVariants, nullModel, options.policy(),
                            kernelOptions(options));
                        if (tileAudit == null) tileAudit = analysis.audits();
                        checkpoint.writeBlock(blockNumber,
                            writer -> writeKernelResults(writer, analysis.results(), options));
                    }
                }
            }
            if (tileAudit == null)
                tileAudit = regenerateAudit(traits, traitColumns, covariateDesign,
                    setTile, selectedVariants, options);
            audits.addAll(tileAudit);
            observedHeapPeak = Math.max(observedHeapPeak, usedHeapBytes());
        }
        writeAudit(options.auditOutput(), audits, definitions.signature(), options);
        checkpoint.assemble(options.output(),
            "Set_ID,Trait_ID,Method,Variants,N,DF,Statistic,R_squared,Effect,T,P_value,"
                + "log10P,P_value_method,Minimum_component_P,Rho_component_P_values");
        if (variants instanceof QRawMatrixCache raw) {
            QRawMatrixCache.IndexedReadStatistics statistics = raw.indexedReadStatistics();
            System.out.println("Indexed set-test variant reads: selections="
                + statistics.selectionCalls() + ", indexed_rows=" + statistics.indexedRows()
                + ", selected_rows=" + statistics.selectedRows() + ", numeric_bytes="
                + statistics.numericBytesRead() + ", persistent_index_reused="
                + statistics.persistentIndexReused());
        }
        System.out.println("Set-test observed JVM heap peak at tile boundaries: "
            + observedHeapPeak + " bytes");
    }

    static SetBlockRecommendation recommendSetBlockSize(QVariantSetTable definitions,
        int sampleCount, QSetTestMethod method, long heapLimitBytes) {
        if (definitions == null || definitions.sets().isEmpty() || sampleCount < 1
            || method == null || !method.isSetTest() || heapLimitBytes < 1)
            throw new IllegalArgumentException("Invalid automatic set-block tuning inputs");
        long halfHeap = heapLimitBytes / 2;
        long desired = Math.max(MINIMUM_AUTOMATIC_BUDGET,
            Math.max(0, halfHeap - AUTOMATIC_FIXED_RESERVE));
        long budget = Math.max(1, Math.min(halfHeap, desired));
        int candidate = Math.min(MAXIMUM_AUTOMATIC_SET_BLOCK_SIZE,
            definitions.sets().size());
        long estimate;
        while (true) {
            estimate = worstTileBytes(definitions, candidate, sampleCount, method);
            if (estimate <= budget || candidate == 1)
                return new SetBlockRecommendation(candidate, heapLimitBytes, budget, estimate);
            double ratio = budget / (double) estimate;
            int reduced = Math.max(1, (int) Math.floor(candidate * ratio * 0.9));
            candidate = reduced < candidate ? reduced : candidate - 1;
        }
    }

    private static long worstTileBytes(QVariantSetTable definitions, int blockSize,
        int sampleCount, QSetTestMethod method) {
        long worst = 0;
        for (int from = 0; from < definitions.sets().size(); from += blockSize) {
            int to = Math.min(definitions.sets().size(), from + blockSize);
            Set<String> variants = new HashSet<>();
            int largestSet = 0;
            for (int set = from; set < to; set++) {
                List<QVariantSetTable.Entry> entries = definitions.sets().get(set).entries();
                largestSet = Math.max(largestSet, entries.size());
                for (QVariantSetTable.Entry entry : entries)
                    variants.add(entry.variantId());
            }
            long rowBytes = saturatedAdd(saturatedMultiply(sampleCount, Double.BYTES), 176);
            long selectedRows = saturatedMultiply(variants.size(), rowBytes);
            long methodWorkspace;
            if (method == QSetTestMethod.BURDEN) {
                methodWorkspace = saturatedMultiply(to - from,
                    saturatedMultiply(sampleCount, 2L * Double.BYTES));
            } else {
                long dosageWorkspace = saturatedMultiply(largestSet,
                    saturatedMultiply(sampleCount, 2L * Double.BYTES));
                long covarianceWorkspace = saturatedMultiply(largestSet,
                    saturatedMultiply(largestSet, 3L * Double.BYTES));
                methodWorkspace = saturatedAdd(dosageWorkspace, covarianceWorkspace);
            }
            long rawEstimate = saturatedAdd(selectedRows, methodWorkspace);
            long safetyAdjusted = saturatedAdd(rawEstimate, rawEstimate / 4);
            worst = Math.max(worst, safetyAdjusted);
        }
        return worst;
    }

    private static long saturatedMultiply(long left, long right) {
        if (left == 0 || right == 0) return 0;
        if (left > Long.MAX_VALUE / right) return Long.MAX_VALUE;
        return left * right;
    }

    private static long saturatedAdd(long left, long right) {
        if (left > Long.MAX_VALUE - right) return Long.MAX_VALUE;
        return left + right;
    }

    private static List<SetAudit> regenerateAudit(QMatrixRowSource traits, int[] traitColumns,
        double[][] covariateDesign, QVariantSetTable definitions, List<Variant> selectedVariants,
        Options options) throws Exception {
        try (QMatrixRowSource.BlockReader reader = traits.open(traitColumns)) {
            Block block = reader.readBlock(options.traitBlockRows());
            if (block == null) throw new IOException("Trait source is empty");
            QSetTestNullModel nullModel = QSetTestNullModel.create(block.rowIds(),
                block.values(), covariateDesign);
            return options.method() == QSetTestMethod.BURDEN
                ? QBurdenReference.analyze(definitions, selectedVariants,
                    nullModel, options.policy()).audits()
                : QKernelSetReference.analyze(definitions, selectedVariants, nullModel,
                    options.policy(), kernelOptions(options)).audits();
        }
    }

    private static QKernelSetReference.Options kernelOptions(Options options) {
        return new QKernelSetReference.Options(options.method(), options.skatORhoGrid(),
            options.skatOSimulations(), options.skatOSeed());
    }

    private static List<Variant> loadSelectedVariants(QMatrixRowSource source,
        int[] columnOrder, QVariantSetTable definitions, int blockRows) throws IOException {
        Set<String> requested = new HashSet<>();
        for (QVariantSetTable.SetDefinition set : definitions.sets())
            for (QVariantSetTable.Entry entry : set.entries())
                requested.add(entry.variantId());
        if (source instanceof QRawMatrixCache raw) {
            Block block = raw.readSelected(requested, columnOrder);
            List<Variant> selected = new ArrayList<>(block.rowCount());
            for (int row = 0; row < block.rowCount(); row++) {
                String id = block.rowIds()[row];
                String[] alleles = allelesFromIdentifier(id);
                selected.add(Variant.takeOwnership(id, alleles[0], alleles[1],
                    block.values()[row]));
            }
            return List.copyOf(selected);
        }
        List<Variant> selected = new ArrayList<>();
        Set<String> seen = new HashSet<>();
        try (QMatrixRowSource.BlockReader reader = source.open(columnOrder)) {
            Block block;
            while ((block = reader.readBlock(Math.max(16, blockRows))) != null) {
                for (int row = 0; row < block.rowCount(); row++) {
                    String id = block.rowIds()[row];
                    if (!requested.contains(id)) continue;
                    if (!seen.add(id))
                        throw new IOException("Duplicate aligned variant ID '" + id + "'");
                    String[] alleles = allelesFromIdentifier(id);
                    selected.add(new Variant(id, alleles[0], alleles[1], block.values()[row]));
                }
            }
        }
        return List.copyOf(selected);
    }

    private static String[] allelesFromIdentifier(String id) {
        String[] fields = id.split(":", -1);
        if (fields.length < 4 || fields[2].isBlank() || fields[3].isBlank())
            throw new IllegalArgumentException("Set-test source variant ID '" + id
                + "' must encode CHROM:POS:REF:ALT so alleles can be verified");
        return new String[] {fields[2], fields[3]};
    }

    private static void writeBurdenResults(Writer writer, List<Result> results, Options options)
        throws IOException {
        for (Result result : results) {
            if (!retained(result.log10P(), result.rSquared(), options)) continue;
            writer.write(csv(result.setId())); writer.write(',');
            writer.write(csv(result.traitId())); writer.write(',');
            writer.write(options.method().optionName()); writer.write(',');
            writer.write(Integer.toString(result.variantCount())); writer.write(',');
            writer.write(Integer.toString(result.sampleCount())); writer.write(',');
            writer.write(Integer.toString(result.residualDegreesOfFreedom())); writer.write(',');
            writer.write(Double.toString(result.tStatistic() * result.tStatistic())); writer.write(',');
            writer.write(Double.toString(result.rSquared())); writer.write(',');
            writer.write(Double.toString(result.effect())); writer.write(',');
            writer.write(Double.toString(result.tStatistic())); writer.write(',');
            writer.write(Double.toString(Math.pow(10, result.log10P()))); writer.write(',');
            writer.write(Double.toString(result.log10P())); writer.write(',');
            writer.write("student-t"); writer.write(',');
            writer.write("NaN,");
            writer.write(System.lineSeparator());
        }
    }

    private static void writeKernelResults(Writer writer,
        List<QKernelSetReference.Result> results, Options options) throws IOException {
        for (QKernelSetReference.Result result : results) {
            if (!retained(result.log10P(), Double.NaN, options)) continue;
            writer.write(csv(result.setId())); writer.write(',');
            writer.write(csv(result.traitId())); writer.write(',');
            writer.write(options.method().optionName()); writer.write(',');
            writer.write(Integer.toString(result.variantCount())); writer.write(',');
            writer.write(Integer.toString(result.sampleCount())); writer.write(',');
            writer.write(Integer.toString(result.residualDegreesOfFreedom())); writer.write(',');
            writer.write(Double.toString(result.statistic())); writer.write(',');
            writer.write("NaN,NaN,NaN,");
            writer.write(Double.toString(result.pValue())); writer.write(',');
            writer.write(Double.toString(result.log10P())); writer.write(',');
            writer.write(csv(result.pValueMethod())); writer.write(',');
            writer.write(Double.toString(result.minimumComponentP())); writer.write(',');
            writer.write(csv(result.rhoComponentPValues()));
            writer.write(System.lineSeparator());
        }
    }

    private static boolean retained(double log10P, double rSquared, Options options) {
        return switch (options.thresholdType().toLowerCase(Locale.ROOT)) {
            case "none" -> true;
            case "rsq" -> {
                if (options.method() != QSetTestMethod.BURDEN)
                    throw new IllegalArgumentException("The rsq threshold is only valid for burden tests");
                yield rSquared >= options.threshold();
            }
            case "pval" -> options.threshold() > 0
                && log10P <= Math.log10(options.threshold());
            default -> throw new IllegalArgumentException("Unsupported set-test threshold "
                + options.thresholdType());
        };
    }

    private static void writeAudit(Path path, List<SetAudit> audits, String definitionSignature,
        Options options) throws IOException {
        Path output = path.toAbsolutePath().normalize();
        Path parent = output.getParent();
        if (parent != null) Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, output.getFileName().toString(), ".partial");
        List<String> lines = new ArrayList<>();
        lines.add("Set_ID\tRequested\tAbsent\tFrequency_excluded\tIncluded\tStatus"
            + "\tIncluded_variant_IDs\tDefinition_signature\tMethod");
        for (SetAudit audit : audits)
            lines.add(tsv(audit.setId()) + "\t" + audit.requestedVariants() + "\t"
                + audit.absentVariants() + "\t" + audit.frequencyExcludedVariants() + "\t"
                + audit.includedVariants() + "\t" + audit.status() + "\t"
                + tsv(String.join(";", audit.includedVariantIds())) + "\t"
                + definitionSignature + "\t" + options.method().optionName());
        Files.write(temporary, lines, StandardCharsets.UTF_8);
        try {
            Files.move(temporary, output, StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(temporary, output, StandardCopyOption.REPLACE_EXISTING);
        }
    }

    private static String signature(QMatrixRowSource variants, int[] variantColumns,
        QMatrixRowSource traits, int[] traitColumns, double[][] design,
        QVariantSetTable definitions, Options options, int effectiveSetBlockSize) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            update(digest, "gpu-eqtl-set-runner-v1");
            update(digest, sourceSignature(variants)); update(digest, sourceSignature(traits));
            update(digest, definitions.signature()); update(digest, options.method().optionName());
            update(digest, options.policy().toString()); update(digest, options.thresholdType());
            update(digest, Long.toString(Double.doubleToLongBits(options.threshold())));
            update(digest, Integer.toString(options.traitBlockRows()));
            update(digest, Integer.toString(effectiveSetBlockSize));
            if (options.skatORhoGrid() != null)
                for (double rho : options.skatORhoGrid())
                    update(digest, Long.toString(Double.doubleToLongBits(rho)));
            update(digest, Integer.toString(options.skatOSimulations()));
            update(digest, Long.toString(options.skatOSeed()));
            for (int value : variantColumns) update(digest, Integer.toString(value));
            for (int value : traitColumns) update(digest, Integer.toString(value));
            for (double[] row : design)
                for (double value : row)
                    update(digest, Long.toString(Double.doubleToLongBits(value)));
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
    }

    private static String sourceSignature(QMatrixRowSource source) {
        QMatrixRowSource.Metadata metadata = source.metadata();
        Path path = metadata.path().toAbsolutePath().normalize();
        try {
            return path + "|" + Files.size(path) + "|" + Files.getLastModifiedTime(path).toMillis()
                + "|" + metadata.rowCount() + "|" + metadata.columnCount() + "|"
                + metadata.cacheSignatureTag();
        } catch (IOException e) {
            throw new IllegalArgumentException("Could not sign set-test source " + path, e);
        }
    }

    private static String csv(String value) {
        if (value.indexOf(',') < 0 && value.indexOf('"') < 0 && value.indexOf('\n') < 0
            && value.indexOf('\r') < 0) return value;
        return '"' + value.replace("\"", "\"\"") + '"';
    }

    private static String tsv(String value) {
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    private static void update(MessageDigest digest, String value) {
        digest.update((byte) 0);
        digest.update(value.getBytes(StandardCharsets.UTF_8));
    }

    private static long usedHeapBytes() {
        Runtime runtime = Runtime.getRuntime();
        return runtime.totalMemory() - runtime.freeMemory();
    }
}
