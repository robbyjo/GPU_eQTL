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
import java.util.Arrays;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QDelimitedMatrixSource;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;

/** One bounded-RAM genotype block, paired with sequential expression blocks. */
final class QeQTLStreamedJobReal implements Runnable {
    private static final int OUTPUT_BUFFER_CHARACTERS = 1024 * 1024;

    private final PreparedBlock genotype;
    private final QDelimitedMatrixSource expressionSource;
    private final int[] expressionColumnOrder;
    private final double[][] covariateQ;
    private final GpuContextPool contextPool;
    private final int genotypeCapacity;
    private final int expressionCapacity;
    private final int localBlockSize;
    private final int errorDegreesOfFreedom;
    private final int degreesOfFreedomOffset;
    private final double rSquaredThreshold;
    private final Writer writer;

    QeQTLStreamedJobReal(PreparedBlock genotype, QDelimitedMatrixSource expressionSource,
        int[] expressionColumnOrder, double[][] covariateQ, GpuContextPool contextPool,
        int genotypeCapacity, int expressionCapacity, int localBlockSize,
        int degreesOfFreedomOffset, int errorDegreesOfFreedom, double rSquaredThreshold, Writer writer) {
        this.genotype = genotype;
        this.expressionSource = expressionSource;
        this.expressionColumnOrder = expressionColumnOrder.clone();
        this.covariateQ = covariateQ;
        this.contextPool = contextPool;
        this.genotypeCapacity = genotypeCapacity;
        this.expressionCapacity = expressionCapacity;
        this.localBlockSize = localBlockSize;
        this.degreesOfFreedomOffset = degreesOfFreedomOffset;
        this.errorDegreesOfFreedom = errorDegreesOfFreedom;
        this.rSquaredThreshold = rSquaredThreshold;
        this.writer = writer;
    }

    @Override
    public void run() {
        int sampleCount = genotype.values()[0].length;
        int paddedSampleCount = roundUp(sampleCount, DEFAULT_ALIGNMENT);
        double[] flattenedGenotypes = new double[genotypeCapacity * paddedSampleCount];
        for (int snp = 0; snp < genotype.values().length; snp++)
            for (int sample = 0; sample < sampleCount; sample++)
                flattenedGenotypes[sample * genotypeCapacity + snp] = genotype.values()[snp][sample] - 1;

        try (QDelimitedMatrixSource.BlockReader expressionReader = expressionSource.open(expressionColumnOrder)) {
            QDelimitedMatrixSource.Block rawExpression;
            while ((rawExpression = expressionReader.readBlock(expressionCapacity)) != null) {
                PreparedBlock expression = QeQTLPreprocessor.prepare(rawExpression, covariateQ, "Expression");
                executePair(genotype, expression, flattenedGenotypes, paddedSampleCount);
            }
        } catch (IOException e) {
            throw new RuntimeException("Cannot stream expression matrix", e);
        }
    }

    private void executePair(PreparedBlock snps, PreparedBlock traits, double[] flattenedGenotypes,
        int paddedSampleCount) {
        System.out.println(snps.rowOffset() + "," + traits.rowOffset());
        double[] flattenedTraits = new double[expressionCapacity * paddedSampleCount];
        for (int trait = 0; trait < traits.values().length; trait++)
            System.arraycopy(traits.values()[trait], 0, flattenedTraits,
                paddedSampleCount * trait, traits.values()[trait].length);

        double[] products;
        GpuContext context = contextPool.reserveContext();
        try {
            long localMemoryBytes = (long) (localBlockSize + 1) * (4L * localBlockSize) * Double.BYTES;
            int activeSnps = roundUp(snps.values().length, localBlockSize);
            int activeTraits = roundUp(traits.values().length, localBlockSize);
            products = context.executeDoubleKernel(flattenedTraits, flattenedGenotypes,
                expressionCapacity * genotypeCapacity, localMemoryBytes, paddedSampleCount,
                genotypeCapacity, new long[] { activeSnps, activeTraits },
                new long[] { localBlockSize, localBlockSize });
        } finally {
            contextPool.releaseContext(context);
        }

        StringBuilder output = new StringBuilder(Math.min(OUTPUT_BUFFER_CHARACTERS, 64 * 1024));
        for (int snp = 0; snp < snps.values().length; snp++) {
            for (int trait = 0; trait < traits.values().length; trait++) {
                double correlation = products[trait * genotypeCapacity + snp];
                double rSquared = correlation * correlation;
                if (rSquared < rSquaredThreshold)
                    continue;
                appendResult(output, snps, traits, snp, trait, correlation);
                if (output.length() >= OUTPUT_BUFFER_CHARACTERS)
                    flush(output);
            }
        }
        flush(output);
        Arrays.fill(products, 0);
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
                .append(rSquared).append(',').append(effect).append(',').append(t).append(',').append(p).append(sLn);
        } else {
            if (QeQTLAnalysis.simplifyResult)
                rSquared = Math.round(rSquared * 10000) / 10000.0;
            output.append(snps.rowIds()[snp]).append(',').append(traits.rowIds()[trait]).append(',')
                .append(rSquared).append(correlation < 0 ? ",-" : ",+").append(sLn);
        }
    }

    private void flush(StringBuilder output) {
        if (output.length() == 0)
            return;
        synchronized (writer) {
            try {
                writer.write(output.toString());
                writer.flush();
            } catch (IOException e) {
                throw new RuntimeException("Cannot write eQTL results", e);
            }
        }
        output.setLength(0);
    }

    private static int roundUp(int number, int multiple) {
        int remainder = number % multiple;
        return remainder == 0 ? number : number + multiple - remainder;
    }
}
