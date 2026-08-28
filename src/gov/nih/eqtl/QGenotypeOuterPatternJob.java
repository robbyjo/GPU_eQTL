/*
 * Copyright 2026 Roby Joehanes
 *
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
import gov.nih.eqtl.io.QPreparedMatrix;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;
import gov.nih.gpu.GpuExecutionMetrics;
import gov.nih.gpu.GpuPatternStatisticsPlan;
import gov.nih.gpu.GpuPatternStatisticsResult;
import gov.nih.gpu.GpuPrecision;

/** One raw genotype block evaluated against all exact trait masks. */
final class QGenotypeOuterPatternJob implements Runnable {
    private static final int LOCAL_BLOCK_SIZE = 16;
    private static final int OUTPUT_BUFFER_CHARACTERS = 1024 * 1024;
    private static final long DESIGN_OUTPUT_TARGET_BYTES = 64L * 1024 * 1024;

    private final QMatrixBlock genotype;
    private final QPreparedMatrix traits;
    private final QTraitPatternModelSet models;
    private final QAnalysisCheckpoint checkpoint;
    private final QGenotypeOuterPatternQcCheckpoint qcCheckpoint;
    private final int blockNumber;
    private final GpuContextPool contextPool;
    private final int genotypeCapacity;
    private final int traitCapacity;
    private final QMissingValuePolicy predictorMissingPolicy;
    private final int degreesOfFreedomOffset;
    private final QeQTLProfiler profiler;
    private final QAnalysisProgress progress;
    private final double minimumMaf;
    private final double minimumMac;
    private final boolean applyPatternFrequencyFilters;

    record QMatrixBlock(long rowOffset, String[] rowIds, double[][] values) { }

    QGenotypeOuterPatternJob(QMatrixBlock genotype, QPreparedMatrix traits,
        QTraitPatternModelSet models, QAnalysisCheckpoint checkpoint, int blockNumber,
        QGenotypeOuterPatternQcCheckpoint qcCheckpoint,
        GpuContextPool contextPool, int genotypeCapacity, int traitCapacity,
        QMissingValuePolicy predictorMissingPolicy, int degreesOfFreedomOffset,
        QeQTLProfiler profiler, QAnalysisProgress progress, double minimumMaf,
        double minimumMac, boolean applyPatternFrequencyFilters) {
        this.genotype = genotype;
        this.traits = traits;
        this.models = models;
        this.checkpoint = checkpoint;
        this.qcCheckpoint = qcCheckpoint;
        this.blockNumber = blockNumber;
        this.contextPool = contextPool;
        this.genotypeCapacity = genotypeCapacity;
        this.traitCapacity = traitCapacity;
        this.predictorMissingPolicy = predictorMissingPolicy;
        this.degreesOfFreedomOffset = degreesOfFreedomOffset;
        this.profiler = profiler;
        this.progress = progress;
        this.minimumMaf = minimumMaf;
        this.minimumMac = minimumMac;
        this.applyPatternFrequencyFilters = applyPatternFrequencyFilters;
    }

    @Override
    public void run() {
        try {
            checkpoint.writeBlock(blockNumber, this::runWithWriter);
        } catch (Exception e) {
            throw new RuntimeException("Failed genotype-outer checkpoint block " + blockNumber, e);
        }
    }

    private void runWithWriter(Writer writer) throws IOException {
        int samples = models.sampleCount();
        int paddedSamples = roundUp(samples, DEFAULT_ALIGNMENT);
        int activeVariants = genotype.values().length;
        int doubledCapacity = Math.multiplyExact(genotypeCapacity, 2);
        int tripledCapacity = Math.multiplyExact(genotypeCapacity, 3);
        long packingStarted = profiler.start();
        double[] traitInputs = new double[Math.multiplyExact(paddedSamples, doubledCapacity)];
        double[] aggregateInputs = new double[Math.multiplyExact(paddedSamples, tripledCapacity)];
        for (int sample = 0; sample < samples; sample++) {
            for (int variant = 0; variant < activeVariants; variant++) {
                double value = genotype.values()[variant][sample];
                boolean missing = gov.nih.eqtl.io.QMissingnessScan.isMissing(value);
                if (missing && predictorMissingPolicy != QMissingValuePolicy.MEAN
                    && predictorMissingPolicy != QMissingValuePolicy.ZERO)
                    throw new IllegalArgumentException("Predictor row '" + genotype.rowIds()[variant]
                        + "' contains missing/non-finite data after applying policy "
                        + predictorMissingPolicy.optionName());
                int traitBase = sample * doubledCapacity + 2 * variant;
                int aggregateBase = sample * tripledCapacity + 3 * variant;
                if (missing) {
                    traitInputs[traitBase + 1] = 1.0;
                } else {
                    traitInputs[traitBase] = value;
                    aggregateInputs[aggregateBase] = 1.0;
                    aggregateInputs[aggregateBase + 1] = value;
                    aggregateInputs[aggregateBase + 2] = value * value;
                }
            }
        }
        profiler.record(QeQTLProfiler.Phase.GENOTYPE_PACK, packingStarted,
            activeVariants, (long) activeVariants * samples * Double.BYTES * 5);

        GpuContext context = contextPool.reserveContext();
        try {
            PatternVariantStatistics patternStatistics = calculatePatternStatistics(
                context, aggregateInputs, paddedSamples, activeVariants, tripledCapacity);
            progress.registerComparisons(patternStatistics.activeComparisons());
            for (long offset = 0; offset < traits.rowCount(); offset += traitCapacity) {
                long readStarted = profiler.start();
                PreparedBlock traitBlock = traits.readBlock(offset, traitCapacity);
                profiler.record(QeQTLProfiler.Phase.EXPRESSION_CACHE_READ, readStarted,
                    traitBlock.values().length,
                    (long) traitBlock.values().length * samples * Double.BYTES);
                calculateAndWriteAssociations(context, traitInputs, paddedSamples,
                    activeVariants, doubledCapacity, traitBlock, patternStatistics, writer);
            }
            qcCheckpoint.writeBlock(blockNumber, patternStatistics.qcCounts());
        } finally {
            contextPool.releaseContext(context);
        }
    }

    private record PatternVariantStatistics(double[][] standardDeviations,
        double[][] replacements, boolean[][] included, int[] includedVariantCounts,
        long activeComparisons,
        QGenotypeOuterPatternQcCheckpoint.BlockCounts qcCounts) { }

	private PatternVariantStatistics calculatePatternStatistics(GpuContext context,
		double[] aggregateInputs, int paddedSamples, int activeVariants, int tripledCapacity) {
		QTraitPatternModelSet.Model[] allModels = models.models();
		double[][] standardDeviations = new double[allModels.length][];
		double[][] replacements = new double[allModels.length][];
		boolean[][] included = new boolean[allModels.length][];
		int[] includedVariantCounts = new int[allModels.length];
		QGenotypeOuterPatternQcCheckpoint.BlockCounts qcCounts =
			new QGenotypeOuterPatternQcCheckpoint.BlockCounts(allModels.length);
		int rowsPerPattern = models.designColumns() + 1;
		long bytesPerPattern = (long) rowsPerPattern * tripledCapacity * Double.BYTES;
		int patternsPerBatch = (int) Math.max(1,
			Math.min(64, DESIGN_OUTPUT_TARGET_BYTES / Math.max(1, bytesPerPattern)));
		GpuPatternStatisticsPlan plan = models.patternStatisticsPlan();
		GpuPatternStatisticsResult compact = context.computePatternStatisticsDouble(
			aggregateInputs, paddedSamples, activeVariants, genotypeCapacity, plan,
			predictorMissingPolicy == QMissingValuePolicy.MEAN, patternsPerBatch,
			LOCAL_BLOCK_SIZE);
		recordGpuMetrics(context.getLastExecutionMetrics());
		for (int pattern = 0; pattern < plan.patternCount(); pattern++) {
			QTraitPatternModelSet.Model model = models.model(plan.patternId(pattern));
			double[] sd = new double[activeVariants];
			double[] replacementValues = new double[activeVariants];
			boolean[] includedVariants = new boolean[activeVariants];
			for (int variant = 0; variant < activeVariants; variant++) {
				double calledValue = compact.value(pattern, variant,
					GpuPatternStatisticsResult.CALLED_COUNT);
				long called = Math.round(calledValue);
				if (called < 0 || called > model.observed.length
					|| Math.abs(calledValue - called) > 1e-8)
					throw new IllegalStateException("Invalid called genotype count " + calledValue
						+ " for pattern " + model.id + ", variant " + variant);
				double dosageSum = compact.value(pattern, variant,
					GpuPatternStatisticsResult.DOSAGE_SUM);
				long missing = model.observed.length - called;
				qcCounts.inputVariants()[model.id]++;
				qcCounts.missingGenotypes()[model.id] = Math.addExact(
					qcCounts.missingGenotypes()[model.id], missing);
				double replacement = compact.value(pattern, variant,
					GpuPatternStatisticsResult.REPLACEMENT);
				replacementValues[variant] = replacement;
				if (called == 0) {
					qcCounts.noCallVariants()[model.id]++;
				} else {
					double filledSum = dosageSum + missing * replacement;
					double centered = compact.value(pattern, variant,
						GpuPatternStatisticsResult.FILLED_SUM_SQUARES)
						- filledSum * filledSum / model.observed.length;
					boolean monomorphic = !(centered > 1e-12) || !Double.isFinite(centered);
					if (monomorphic) qcCounts.monomorphicVariants()[model.id]++;
					double alleleNumber = 2.0 * called;
					double eaf = dosageSum / alleleNumber;
					double maf = Math.min(eaf, 1.0 - eaf);
					double mac = Math.min(dosageSum, alleleNumber - dosageSum);
					boolean belowMaf = maf + 1e-12 < minimumMaf;
					boolean belowMac = mac + 1e-12 < minimumMac;
					if (belowMaf)
						qcCounts.belowMinimumMaf()[model.id]++;
					if (belowMac)
						qcCounts.belowMinimumMac()[model.id]++;
					boolean include = !monomorphic
						&& (!applyPatternFrequencyFilters || (!belowMaf && !belowMac));
					if (include) {
						includedVariants[variant] = true;
						includedVariantCounts[model.id]++;
						qcCounts.includedVariants()[model.id]++;
					}
				}
				double residual = QPatternSufficientStatistics.validateResidualSumSquares(
					compact.value(pattern, variant,
						GpuPatternStatisticsResult.FILLED_SUM_SQUARES),
					compact.value(pattern, variant,
						GpuPatternStatisticsResult.RESIDUAL_SUM_SQUARES));
				sd[variant] = residual > 0
					? QPatternSufficientStatistics.standardDeviation(residual, model.observed.length)
					: Double.NaN;
            }
			standardDeviations[model.id] = sd;
			replacements[model.id] = replacementValues;
			included[model.id] = includedVariants;
		}
		long activeComparisons = 0;
		for (QTraitPatternModelSet.Model model : allModels)
			if (model.estimable())
				activeComparisons = Math.addExact(activeComparisons,
					Math.multiplyExact(model.traitRows,
						(long) includedVariantCounts[model.id]));
		return new PatternVariantStatistics(standardDeviations, replacements, included,
			includedVariantCounts, activeComparisons, qcCounts);
	}

    private void calculateAndWriteAssociations(GpuContext context, double[] traitInputs,
        int paddedSamples, int activeVariants, int doubledCapacity, PreparedBlock traitBlock,
        PatternVariantStatistics patternStatistics, Writer writer) throws IOException {
        long packingStarted = profiler.start();
        double[] flattenedTraits = new double[Math.multiplyExact(traitCapacity, paddedSamples)];
        for (int trait = 0; trait < traitBlock.values().length; trait++)
            System.arraycopy(traitBlock.values()[trait], 0, flattenedTraits,
                trait * paddedSamples, traitBlock.values()[trait].length);
        profiler.record(QeQTLProfiler.Phase.EXPRESSION_PACK, packingStarted,
            traitBlock.values().length,
            (long) traitBlock.values().length * models.sampleCount() * Double.BYTES);
        long localBytes = (long) (LOCAL_BLOCK_SIZE + 1) * (4L * LOCAL_BLOCK_SIZE)
            * Double.BYTES;
        double[] products = context.executeDoubleKernel(flattenedTraits, traitInputs,
            Math.multiplyExact(traitCapacity, doubledCapacity), localBytes, paddedSamples,
            doubledCapacity, new long[] {roundUp(activeVariants * 2, LOCAL_BLOCK_SIZE),
                roundUp(traitBlock.values().length, LOCAL_BLOCK_SIZE)},
            new long[] {LOCAL_BLOCK_SIZE, LOCAL_BLOCK_SIZE});
        recordGpuMetrics(context.getLastExecutionMetrics());

        long resultsStarted = profiler.start();
        StringBuilder output = new StringBuilder(64 * 1024);
        for (int variant = 0; variant < activeVariants; variant++) {
            for (int trait = 0; trait < traitBlock.values().length; trait++) {
                long preparedRow = traitBlock.rowOffset() + trait;
                QTraitPatternModelSet.Model model = models.model(
                    models.patternForPreparedRow(preparedRow));
                if (!patternStatistics.included()[model.id][variant])
                    continue;
                double genotypeSd = patternStatistics.standardDeviations()[model.id][variant];
                if (!(genotypeSd > 0) || !Double.isFinite(genotypeSd))
                    continue;
                int productBase = trait * doubledCapacity + 2 * variant;
                double rawDot = products[productBase]
                    + patternStatistics.replacements()[model.id][variant]
                        * products[productBase + 1];
                double correlation = QeQTLStatistics.validateCorrelation(
                    QPatternSufficientStatistics.correlation(rawDot, genotypeSd,
                        model.observed.length), GpuPrecision.FP64);
                if (correlation * correlation < model.rSquaredThreshold)
                    continue;
                appendResult(output, variant, trait, traitBlock, model, genotypeSd, correlation);
                if (output.length() >= OUTPUT_BUFFER_CHARACTERS) {
                    writer.write(output.toString());
                    output.setLength(0);
                }
            }
        }
        if (!output.isEmpty()) writer.write(output.toString());
        long comparisons = 0;
        for (int trait = 0; trait < traitBlock.values().length; trait++) {
            QTraitPatternModelSet.Model model = models.model(
                models.patternForPreparedRow(traitBlock.rowOffset() + trait));
            comparisons = Math.addExact(comparisons,
                patternStatistics.includedVariantCounts()[model.id]);
        }
        profiler.record(QeQTLProfiler.Phase.CPU_RESULTS_AND_WRITE,
            resultsStarted, comparisons, 0);
        if (comparisons > 0) progress.addComparisons(comparisons);
    }

    private void appendResult(StringBuilder output, int variant, int trait,
        PreparedBlock traitBlock, QTraitPatternModelSet.Model model,
        double genotypeSd, double correlation) {
        QeQTLStatistics.Result result = QeQTLStatistics.calculate(correlation,
            traitBlock.standardDeviations()[trait], genotypeSd,
            model.errorDegreesOfFreedom, degreesOfFreedomOffset);
        double rSquared = result.rSquared();
        if (!QeQTLAnalysis.rsqOnly) {
            double effect = result.effect();
            double t = result.tStatistic();
            double p = result.log10P();
            if (QeQTLAnalysis.simplifyResult) {
                rSquared = Math.round(rSquared * 10000) / 10000.0;
                effect = Math.round(effect * 10000) / 10000.0;
                t = Math.copySign(Math.round(Math.abs(t) * 10000) / 10000.0, correlation);
                p = Math.round(p * 10000) / 10000.0;
            }
            output.append(genotype.rowIds()[variant]).append(',').append(traitBlock.rowIds()[trait])
                .append(',').append(rSquared).append(',').append(effect).append(',').append(t)
                .append(',').append(p);
        } else {
            if (QeQTLAnalysis.simplifyResult)
                rSquared = Math.round(rSquared * 10000) / 10000.0;
            output.append(genotype.rowIds()[variant]).append(',').append(traitBlock.rowIds()[trait])
                .append(',').append(rSquared).append(correlation < 0 ? ",-" : ",+");
        }
        output.append(',').append(model.observed.length).append(',')
            .append(model.reportedDegreesOfFreedom(degreesOfFreedomOffset)).append(sLn);
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

    private static int roundUp(int number, int multiple) {
        int remainder = number % multiple;
        return remainder == 0 ? number : number + multiple - remainder;
    }
}
