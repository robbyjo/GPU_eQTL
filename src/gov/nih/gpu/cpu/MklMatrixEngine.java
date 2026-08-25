/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 * The oneMKL linking permission in LICENSE_EXCEPTION applies to this file.
 */
package gov.nih.gpu.cpu;

import org.bytedeco.javacpp.Loader;

import java.nio.charset.StandardCharsets;

import static org.bytedeco.mkl.global.mkl_rt.CblasNoTrans;
import static org.bytedeco.mkl.global.mkl_rt.CblasRowMajor;
import static org.bytedeco.mkl.global.mkl_rt.CblasTrans;
import static org.bytedeco.mkl.global.mkl_rt.MKL_Get_Max_Threads;
import static org.bytedeco.mkl.global.mkl_rt.MKL_Get_Version_String;
import static org.bytedeco.mkl.global.mkl_rt.MKL_Set_Num_Threads;
import static org.bytedeco.mkl.global.mkl_rt.cblas_dgemm;
import static org.bytedeco.mkl.global.mkl_rt.cblas_sgemm;

/** Native Intel oneMKL implementation supplied by the JavaCPP preset. */
final class MklMatrixEngine implements CpuMatrixEngine {
	private final int threads;
	private final String version;
	private final String nativePath;

	MklMatrixEngine(int requestedThreads) {
		if (requestedThreads <= 0)
			throw new IllegalArgumentException("oneMKL thread count must be positive");
		nativePath = Loader.load(org.bytedeco.mkl.global.mkl_rt.class);
		MKL_Set_Num_Threads(requestedThreads);
		threads = Math.max(1, MKL_Get_Max_Threads());
		byte[] buffer = new byte[256];
		MKL_Get_Version_String(buffer, buffer.length);
		int length = 0;
		while (length < buffer.length && buffer[length] != 0) length++;
		String detected = new String(buffer, 0, length, StandardCharsets.UTF_8).trim();
		version = detected.isEmpty() ? "Intel oneMKL 2026.1" : detected;
	}

	@Override
	public String description() {
		return version + ", " + threads + " BLAS thread" + (threads == 1 ? "" : "s")
			+ (nativePath == null ? "" : " (bundled native library)");
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
