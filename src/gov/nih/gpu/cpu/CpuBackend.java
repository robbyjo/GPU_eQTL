/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cpu;

import gov.nih.gpu.GpuBackend;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuException;

import java.util.List;
import java.util.Locale;

/** Always-available host backend with native BLAS acceleration and a Java fallback. */
public final class CpuBackend implements GpuBackend {
	public static final String BLAS_PROPERTY = "eqtl.cpu.blas";
	public static final String THREADS_PROPERTY = "eqtl.cpu.threads";
	private final String requestedEngine;
	private final CpuDevice device;

	public CpuBackend() {
		requestedEngine = parseEngine(System.getProperty(BLAS_PROPERTY, "auto"));
		device = new CpuDevice(this);
	}

	@Override
	public String getName() {
		return "cpu";
	}

	@Override
	public String getRuntimeDescription() {
		return "host FP64/FP32 backend; engine=" + requestedEngine
			+ "; optional oneMKL 2026.1 on Windows/Linux x64; bundled OpenBLAS 0.3.34"
			+ " via JavaCPP 1.5.14; portable Java fallback";
	}

	@Override
	public List<GpuDevice> discoverGpuDevices() {
		return List.of(device);
	}

	CpuMatrixEngine openEngine() {
		int requestedThreads = requestedThreadCount();
		if ("java".equals(requestedEngine))
			return new JavaMatrixEngine();
		if ("mkl".equals(requestedEngine)) {
			try {
				return new MklMatrixEngine(requestedThreads);
			} catch (RuntimeException | LinkageError e) {
				throw new GpuException("The explicitly selected oneMKL CPU engine could not be loaded", e);
			}
		}
		if ("openblas".equals(requestedEngine)) {
			try {
				return new OpenBlasMatrixEngine(requestedThreads);
			} catch (RuntimeException | LinkageError e) {
				throw new GpuException("The explicitly selected OpenBLAS CPU engine could not be loaded", e);
			}
		}
		Throwable mklFailure;
		try {
			return new MklMatrixEngine(requestedThreads);
		} catch (RuntimeException | LinkageError e) {
			mklFailure = e;
		}
		try {
			return new OpenBlasMatrixEngine(requestedThreads);
		} catch (RuntimeException | LinkageError openBlasFailure) {
			System.err.println("WARNING: Neither oneMKL nor bundled OpenBLAS could be loaded; using the slower "
				+ "portable Java CPU engine. oneMKL: " + conciseMessage(mklFailure)
				+ "; OpenBLAS: " + conciseMessage(openBlasFailure));
			return new JavaMatrixEngine();
		}
	}

	private static String parseEngine(String value) {
		String normalized = value == null ? "auto" : value.trim().toLowerCase(Locale.ROOT);
		if (normalized.isEmpty()) normalized = "auto";
		if (!("auto".equals(normalized) || "mkl".equals(normalized)
			|| "openblas".equals(normalized) || "java".equals(normalized)))
			throw new GpuException(BLAS_PROPERTY + " must be auto, mkl, openblas, or java, not '" + value + "'");
		return normalized;
	}

	private static int requestedThreadCount() {
		String configured = System.getProperty(THREADS_PROPERTY);
		if (configured == null || configured.isBlank())
			return Math.max(1, Runtime.getRuntime().availableProcessors() - 1);
		try {
			int value = Integer.parseInt(configured.trim());
			if (value <= 0)
				throw new NumberFormatException();
			return value;
		} catch (NumberFormatException e) {
			throw new GpuException(THREADS_PROPERTY + " must be a positive integer, not '" + configured + "'");
		}
	}

	private static String conciseMessage(Throwable error) {
		String message = error.getMessage();
		return message == null || message.isBlank() ? error.getClass().getSimpleName() : message;
	}
}
