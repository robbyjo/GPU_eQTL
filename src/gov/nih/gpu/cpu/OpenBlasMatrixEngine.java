/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cpu;

import org.bytedeco.javacpp.Loader;

import static org.bytedeco.openblas.global.openblas.CblasNoTrans;
import static org.bytedeco.openblas.global.openblas.CblasRowMajor;
import static org.bytedeco.openblas.global.openblas.CblasTrans;
import static org.bytedeco.openblas.global.openblas.cblas_dgemm;
import static org.bytedeco.openblas.global.openblas.cblas_sgemm;
import static org.bytedeco.openblas.presets.openblas_nolapack.blas_get_num_threads;
import static org.bytedeco.openblas.presets.openblas_nolapack.blas_set_num_threads;

/** Native OpenBLAS implementation supplied by the JavaCPP preset. */
final class OpenBlasMatrixEngine implements CpuMatrixEngine {
	private final int threads;
	private final String nativePath;

	OpenBlasMatrixEngine(int requestedThreads) {
		if (requestedThreads <= 0)
			throw new IllegalArgumentException("OpenBLAS thread count must be positive");
		nativePath = Loader.load(org.bytedeco.openblas.global.openblas.class);
		blas_set_num_threads(requestedThreads);
		threads = Math.max(1, blas_get_num_threads());
	}

	@Override
	public String description() {
		return "OpenBLAS 0.3.34 via JavaCPP 1.5.14, " + threads + " BLAS thread"
			+ (threads == 1 ? "" : "s") + (nativePath == null ? "" : " (bundled native library)");
	}

	@Override
	public int threadCount() {
		return threads;
	}

	@Override
	public void multiplyDouble(double[] left, double[] right, double[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, double alpha, double beta) {
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			rows, columns, shared, alpha, left, leftStride,
			right, rightStride, beta, output, outputStride);
	}

	@Override
	public void multiplyFloat(float[] left, float[] right, float[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, float alpha, float beta) {
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			rows, columns, shared, alpha, left, leftStride,
			right, rightStride, beta, output, outputStride);
	}

	@Override
	public void multiplyDoubleRightTransposed(double[] left, double[] right, double[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, double alpha, double beta) {
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			rows, columns, shared, alpha, left, leftStride,
			right, rightStride, beta, output, outputStride);
	}

	@Override
	public void multiplyFloatRightTransposed(float[] left, float[] right, float[] output,
		int rows, int columns, int shared, int leftStride, int rightStride,
		int outputStride, float alpha, float beta) {
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			rows, columns, shared, alpha, left, leftStride,
			right, rightStride, beta, output, outputStride);
	}
}
