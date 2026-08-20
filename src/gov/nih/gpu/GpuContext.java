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

	@Override
	void close();
}
