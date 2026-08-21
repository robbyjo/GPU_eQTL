/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import static gov.nih.gpu.GpuRuntime.DEFAULT_ALIGNMENT;
import static gov.nih.utils.QStringUtils.sLn;

import java.io.IOException;
import java.io.Writer;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QBinaryMatrixCache;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;
import gov.nih.gpu.GpuExecutionMetrics;
import gov.nih.gpu.GpuPrecision;

/** One bounded-RAM genotype block, paired with sequential expression blocks. */
final class QeQTLStreamedJobReal implements Runnable {
    private static final int OUTPUT_BUFFER_CHARACTERS = 1024 * 1024;

    private final PreparedBlock genotype;
    private final QBinaryMatrixCache expressionCache;
    private final QAnalysisCheckpoint checkpoint;
    private final int genotypeBlockNumber;
    private final GpuContextPool contextPool;
    private final int genotypeCapacity;
    private final int expressionCapacity;
    private final int localBlockSize;
    private final int errorDegreesOfFreedom;
    private final int degreesOfFreedomOffset;
    private final double rSquaredThreshold;
    private final GpuPrecision precision;
    private final QeQTLProfiler profiler;
    private final boolean includeSampleStatistics;
    private final int sampleCount;

    QeQTLStreamedJobReal(PreparedBlock genotype, QBinaryMatrixCache expressionCache,
        QAnalysisCheckpoint checkpoint, int genotypeBlockNumber, GpuContextPool contextPool,
        int genotypeCapacity, int expressionCapacity, int localBlockSize,
        int degreesOfFreedomOffset, int errorDegreesOfFreedom, double rSquaredThreshold,
        GpuPrecision precision, QeQTLProfiler profiler, boolean includeSampleStatistics) {
        this.genotype = genotype;
        this.expressionCache = expressionCache;
        this.checkpoint = checkpoint;
        this.genotypeBlockNumber = genotypeBlockNumber;
        this.contextPool = contextPool;
        this.genotypeCapacity = genotypeCapacity;
        this.expressionCapacity = expressionCapacity;
        this.localBlockSize = localBlockSize;
        this.degreesOfFreedomOffset = degreesOfFreedomOffset;
        this.errorDegreesOfFreedom = errorDegreesOfFreedom;
        this.rSquaredThreshold = rSquaredThreshold;
        this.precision = precision;
        this.profiler = profiler;
        this.includeSampleStatistics = includeSampleStatistics;
        sampleCount = genotype.values()[0].length;
    }

    @Override
    public void run() {
        try {
            checkpoint.writeBlock(genotypeBlockNumber, this::runWithWriter);
        } catch (Exception e) {
            throw new RuntimeException("Failed genotype checkpoint block " + genotypeBlockNumber, e);
        }
    }

    private void runWithWriter(Writer writer) throws IOException {
        long packingStarted = profiler.start();
        int sampleCount = genotype.values()[0].length;
        int paddedSampleCount = roundUp(sampleCount, DEFAULT_ALIGNMENT);
        double[] flattenedGenotypes = precision == GpuPrecision.FP64
            ? new double[genotypeCapacity * paddedSampleCount] : null;
        float[] flattenedGenotypes32 = precision == GpuPrecision.FP32
            ? new float[genotypeCapacity * paddedSampleCount] : null;
        for (int snp = 0; snp < genotype.values().length; snp++) {
            for (int sample = 0; sample < sampleCount; sample++) {
                int destination = sample * genotypeCapacity + snp;
                double value = genotype.values()[snp][sample] - 1;
                if (precision == GpuPrecision.FP32)
                    flattenedGenotypes32[destination] = (float) value;
                else
                    flattenedGenotypes[destination] = value;
            }
        }
        profiler.record(QeQTLProfiler.Phase.GENOTYPE_PACK, packingStarted,
            genotype.values().length, (long) genotype.values().length * sampleCount * precision.bytes());

        for (long offset = 0; offset < expressionCache.rowCount(); offset += expressionCapacity) {
            long cacheReadStarted = profiler.start();
            PreparedBlock expression = expressionCache.readBlock(offset, expressionCapacity);
            profiler.record(QeQTLProfiler.Phase.EXPRESSION_CACHE_READ, cacheReadStarted,
                expression.values().length,
                (long) expression.values().length * sampleCount * Double.BYTES);
            executePair(genotype, expression, flattenedGenotypes, flattenedGenotypes32,
                paddedSampleCount, writer);
        }
    }

    private void executePair(PreparedBlock snps, PreparedBlock traits, double[] flattenedGenotypes,
        float[] flattenedGenotypes32, int paddedSampleCount, Writer writer) {
        if (QeQTLAnalysis.DEBUG)
            System.out.println(snps.rowOffset() + "," + traits.rowOffset());
        long packingStarted = profiler.start();
        double[] flattenedTraits = precision == GpuPrecision.FP64
            ? new double[expressionCapacity * paddedSampleCount] : null;
        float[] flattenedTraits32 = precision == GpuPrecision.FP32
            ? new float[expressionCapacity * paddedSampleCount] : null;
        for (int trait = 0; trait < traits.values().length; trait++) {
            int destination = paddedSampleCount * trait;
            if (precision == GpuPrecision.FP32) {
                for (int sample = 0; sample < traits.values()[trait].length; sample++)
                    flattenedTraits32[destination + sample] = (float) traits.values()[trait][sample];
            } else {
                System.arraycopy(traits.values()[trait], 0, flattenedTraits,
                    destination, traits.values()[trait].length);
            }
        }
        profiler.record(QeQTLProfiler.Phase.EXPRESSION_PACK, packingStarted,
            traits.values().length, (long) traits.values().length * traits.values()[0].length * precision.bytes());

        double[] products = null;
        float[] products32 = null;
        long waitStarted = profiler.start();
        GpuContext context = contextPool.reserveContext();
        profiler.record(QeQTLProfiler.Phase.GPU_CONTEXT_WAIT, waitStarted, 1, 0);
        try {
            long localMemoryBytes = (long) (localBlockSize + 1) * (4L * localBlockSize) * precision.bytes();
            int activeSnps = roundUp(snps.values().length, localBlockSize);
            int activeTraits = roundUp(traits.values().length, localBlockSize);
            if (precision == GpuPrecision.FP32)
                products32 = context.executeFloatKernel(flattenedTraits32, flattenedGenotypes32,
                    expressionCapacity * genotypeCapacity, localMemoryBytes, paddedSampleCount,
                    genotypeCapacity, new long[] { activeSnps, activeTraits },
                    new long[] { localBlockSize, localBlockSize });
            else
                products = context.executeDoubleKernel(flattenedTraits, flattenedGenotypes,
                    expressionCapacity * genotypeCapacity, localMemoryBytes, paddedSampleCount,
                    genotypeCapacity, new long[] { activeSnps, activeTraits },
                    new long[] { localBlockSize, localBlockSize });
            recordGpuMetrics(context.getLastExecutionMetrics());
        } finally {
            contextPool.releaseContext(context);
        }

        long resultsStarted = profiler.start();
        StringBuilder output = new StringBuilder(Math.min(OUTPUT_BUFFER_CHARACTERS, 64 * 1024));
        for (int snp = 0; snp < snps.values().length; snp++) {
            for (int trait = 0; trait < traits.values().length; trait++) {
                int resultOffset = trait * genotypeCapacity + snp;
                double correlation = QeQTLStatistics.validateCorrelation(
                    precision == GpuPrecision.FP32 ? products32[resultOffset] : products[resultOffset], precision);
                double rSquared = correlation * correlation;
                if (rSquared < rSquaredThreshold)
                    continue;
                appendResult(output, snps, traits, snp, trait, correlation);
                if (output.length() >= OUTPUT_BUFFER_CHARACTERS)
                    flush(output, writer);
            }
        }
        flush(output, writer);
        profiler.record(QeQTLProfiler.Phase.CPU_RESULTS_AND_WRITE, resultsStarted,
            (long) snps.values().length * traits.values().length, 0);
    }

    private void recordGpuMetrics(GpuExecutionMetrics metrics) {
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_BUFFER_SETUP,
            metrics.bufferSetupNanoseconds(), 1, 0);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_UPLOAD,
            metrics.uploadNanoseconds(), 1, metrics.uploadedBytes());
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_COMPUTE,
            metrics.computeNanoseconds(), 1, 0);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_DOWNLOAD,
            metrics.downloadNanoseconds(), 1, metrics.downloadedBytes());
    }

    private void appendResult(StringBuilder output, PreparedBlock snps, PreparedBlock traits,
        int snp, int trait, double correlation) {
        QeQTLStatistics.Result result = QeQTLStatistics.calculate(correlation,
            traits.standardDeviations()[trait], snps.standardDeviations()[snp],
            errorDegreesOfFreedom, degreesOfFreedomOffset);
        double rSquared = result.rSquared();
        if (!QeQTLAnalysis.rsqOnly) {
            double effect = result.effect();
            double t = result.tStatistic();
            double p = result.log10P();
            if (QeQTLAnalysis.simplifyResult) {
                rSquared = Math.round(rSquared * 10000) / 10000.0;
                effect = Math.round(effect * 10000) / 10000.0;
                t = Math.round(Math.abs(t) * 10000) / 10000.0;
                p = Math.round(p * 10000) / 10000.0;
                t = Math.copySign(t, correlation);
            }
            output.append(snps.rowIds()[snp]).append(',').append(traits.rowIds()[trait]).append(',')
                .append(rSquared).append(',').append(effect).append(',').append(t).append(',').append(p);
        } else {
            if (QeQTLAnalysis.simplifyResult)
                rSquared = Math.round(rSquared * 10000) / 10000.0;
            output.append(snps.rowIds()[snp]).append(',').append(traits.rowIds()[trait]).append(',')
                .append(rSquared).append(correlation < 0 ? ",-" : ",+");
        }
        if (includeSampleStatistics)
            output.append(',').append(sampleCount).append(',')
                .append(errorDegreesOfFreedom - degreesOfFreedomOffset);
        output.append(sLn);
    }

    private void flush(StringBuilder output, Writer writer) {
        if (output.length() == 0)
            return;
        try {
            writer.write(output.toString());
        } catch (IOException e) {
            throw new RuntimeException("Cannot write eQTL results", e);
        }
        output.setLength(0);
    }

    private static int roundUp(int number, int multiple) {
        int remainder = number % multiple;
        return remainder == 0 ? number : number + multiple - remainder;
    }
}
