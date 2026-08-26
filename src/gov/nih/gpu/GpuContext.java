/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

/**
 * An exclusive execution context for a GPU.
 *
 * <p>The small API deliberately keeps vendor handles out of the analysis code.
 * A CUDA, HIP, Vulkan, or Level Zero backend can implement the same contract
 * without changing the eQTL scheduling and result code.</p>
 */
public interface GpuContext extends AutoCloseable {
	GpuDevice getDevice();

	/** Enable synchronization and phase timing for subsequent calls. */
	default void setProfilingEnabled(boolean enabled) { }

	/** Metrics from the most recent call while this exclusive context was reserved. */
	default GpuExecutionMetrics getLastExecutionMetrics() {
		return GpuExecutionMetrics.EMPTY;
	}

	void compileKernel(String source, String kernelName);

	default void compileKernel(String source, String kernelName, GpuPrecision precision) {
		if (precision != GpuPrecision.FP64)
			throw new UnsupportedOperationException("This GPU context does not implement FP32");
		compileKernel(source, kernelName);
	}

	double[] executeDoubleKernel(
		double[] inputA,
		double[] inputB,
		int outputElements,
		long localMemoryBytes,
		int widthA,
		int widthB,
		long[] globalWorkSize,
		long[] localWorkSize);

	/**
	 * Compute compact FP64 pattern-specific predictor replacement and variance state.
	 * Implementations may retain device allocations while this exclusive context is open.
	 */
	default GpuPatternStatisticsResult computePatternStatisticsDouble(
		double[] aggregateInputs, int paddedSamples, int activeVariants, int variantCapacity,
		GpuPatternStatisticsPlan plan, boolean meanFill, int patternsPerBatch,
		int workGroupSize) {
		return GpuPatternStatisticsSupport.calculate(this, aggregateInputs, paddedSamples,
			activeVariants, variantCapacity, plan, meanFill, patternsPerBatch, workGroupSize);
	}

	default float[] executeFloatKernel(
		float[] inputA,
		float[] inputB,
		int outputElements,
		long localMemoryBytes,
		int widthA,
		int widthB,
		long[] globalWorkSize,
		long[] localWorkSize) {
		throw new UnsupportedOperationException("This GPU context does not implement FP32");
	}

	/**
	 * Residualize row-major observations as {@code Y - (Y Q) Q^T} in FP64.
	 * The supplied Q matrix is row-major with {@code sampleCount} rows and
	 * {@code covariateRank} orthonormal columns.
	 */
	default double[] residualizeDoubleRows(double[] rowMajorValues, double[] rowMajorQ,
		int rowCount, int sampleCount, int covariateRank) {
		throw new UnsupportedOperationException("This GPU context does not implement FP64 residualization");
	}

	/** FP32 counterpart of {@link #residualizeDoubleRows(double[], double[], int, int, int)}. */
	default float[] residualizeFloatRows(float[] rowMajorValues, float[] rowMajorQ,
		int rowCount, int sampleCount, int covariateRank) {
		throw new UnsupportedOperationException("This GPU context does not implement FP32 residualization");
	}

	/** Release projection-specific buffers without closing this reusable context. */
	default void releaseResidualizationResources() { }

	@Override
	void close();
}
