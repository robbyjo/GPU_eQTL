/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cpu;

/** Host matrix operations used by the CPU compute context. */
interface CpuMatrixEngine {
	String description();

	int threadCount();

	void multiplyDouble(double[] left, double[] right, double[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, double alpha, double beta);

	void multiplyFloat(float[] left, float[] right, float[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, float alpha, float beta);

	void multiplyDoubleRightTransposed(double[] left, double[] right, double[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, double alpha, double beta);

	void multiplyFloatRightTransposed(float[] left, float[] right, float[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, float alpha, float beta);
}
