/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cpu;

import java.util.Arrays;

/** Dependency-free deterministic fallback when a native BLAS cannot load. */
final class JavaMatrixEngine implements CpuMatrixEngine {
	private static final int COLUMN_BLOCK = 64;

	@Override
	public String description() {
		return "portable Java FP64/FP32 matrix engine (native OpenBLAS unavailable or disabled)";
	}

	@Override
	public int threadCount() {
		return 1;
	}

	@Override
	public void multiplyDouble(double[] left, double[] right, double[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, double alpha, double beta) {
		scale(output, rows, columns, outputStride, beta);
		for (int row = 0; row < rows; row++) {
			int leftBase = row * leftStride;
			int outputBase = row * outputStride;
			for (int columnBase = 0; columnBase < columns; columnBase += COLUMN_BLOCK) {
				int columnEnd = Math.min(columns, columnBase + COLUMN_BLOCK);
				for (int k = 0; k < shared; k++) {
					double scaledLeft = alpha * left[leftBase + k];
					int rightBase = k * rightStride;
					for (int column = columnBase; column < columnEnd; column++)
						output[outputBase + column] += scaledLeft * right[rightBase + column];
				}
			}
		}
	}

	@Override
	public void multiplyFloat(float[] left, float[] right, float[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, float alpha, float beta) {
		scale(output, rows, columns, outputStride, beta);
		for (int row = 0; row < rows; row++) {
			int leftBase = row * leftStride;
			int outputBase = row * outputStride;
			for (int columnBase = 0; columnBase < columns; columnBase += COLUMN_BLOCK) {
				int columnEnd = Math.min(columns, columnBase + COLUMN_BLOCK);
				for (int k = 0; k < shared; k++) {
					float scaledLeft = alpha * left[leftBase + k];
					int rightBase = k * rightStride;
					for (int column = columnBase; column < columnEnd; column++)
						output[outputBase + column] += scaledLeft * right[rightBase + column];
				}
			}
		}
	}

	@Override
	public void multiplyDoubleRightTransposed(double[] left, double[] right, double[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, double alpha, double beta) {
		scale(output, rows, columns, outputStride, beta);
		for (int row = 0; row < rows; row++) {
			int leftBase = row * leftStride;
			int outputBase = row * outputStride;
			for (int column = 0; column < columns; column++) {
				int rightBase = column * rightStride;
				double sum = 0;
				for (int k = 0; k < shared; k++)
					sum += left[leftBase + k] * right[rightBase + k];
				output[outputBase + column] += alpha * sum;
			}
		}
	}

	@Override
	public void multiplyFloatRightTransposed(float[] left, float[] right, float[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, float alpha, float beta) {
		scale(output, rows, columns, outputStride, beta);
		for (int row = 0; row < rows; row++) {
			int leftBase = row * leftStride;
			int outputBase = row * outputStride;
			for (int column = 0; column < columns; column++) {
				int rightBase = column * rightStride;
				float sum = 0;
				for (int k = 0; k < shared; k++)
					sum += left[leftBase + k] * right[rightBase + k];
				output[outputBase + column] += alpha * sum;
			}
		}
	}

	private static void scale(double[] values, int rows, int columns, int stride, double factor) {
		if (factor == 0) {
			if (columns == stride) Arrays.fill(values, 0);
			else for (int row = 0; row < rows; row++)
				Arrays.fill(values, row * stride, row * stride + columns, 0);
		} else if (factor != 1) {
			for (int row = 0; row < rows; row++)
				for (int column = 0; column < columns; column++)
					values[row * stride + column] *= factor;
		}
	}

	private static void scale(float[] values, int rows, int columns, int stride, float factor) {
		if (factor == 0) {
			if (columns == stride) Arrays.fill(values, 0);
			else for (int row = 0; row < rows; row++)
				Arrays.fill(values, row * stride, row * stride + columns, 0);
		} else if (factor != 1) {
			for (int row = 0; row < rows; row++)
				for (int column = 0; column < columns; column++)
					values[row * stride + column] *= factor;
		}
	}
}
