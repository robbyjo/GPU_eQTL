/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cpu;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuExecutionMetrics;
import gov.nih.gpu.GpuException;
import gov.nih.gpu.GpuPrecision;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/** CPU implementation of the production association and residualization operations. */
final class CpuContext implements GpuContext {
	private static final Pattern N_MINUS_ONE = Pattern.compile(
		"(?m)^\\s*#define\\s+N_MIN_1\\s+(\\d+)\\s*$");
	private static final Pattern DATA_TYPE = Pattern.compile(
		"(?m)^\\s*#define\\s+DATATYPE\\s+(float|double)\\s*$");

	private final CpuDevice device;
	private final CpuMatrixEngine engine;
	private GpuPrecision precision = GpuPrecision.FP64;
	private double normalization;
	private boolean kernelReady;
	private boolean profilingEnabled;
	private boolean closed;
	private GpuExecutionMetrics lastExecutionMetrics = GpuExecutionMetrics.EMPTY;

	CpuContext(CpuDevice device, CpuMatrixEngine engine) {
		this.device = device;
		this.engine = engine;
	}

	@Override
	public GpuDevice getDevice() {
		return device;
	}

	@Override
	public synchronized void setProfilingEnabled(boolean enabled) {
		profilingEnabled = enabled;
		lastExecutionMetrics = GpuExecutionMetrics.EMPTY;
	}

	@Override
	public synchronized GpuExecutionMetrics getLastExecutionMetrics() {
		return lastExecutionMetrics;
	}

	@Override
	public synchronized void compileKernel(String source, String kernelName) {
		compileKernel(source, kernelName, GpuPrecision.FP64);
	}

	@Override
	public synchronized void compileKernel(String source, String kernelName, GpuPrecision precision) {
		ensureOpen();
		if (!"eqtlReal".equals(kernelName))
			throw new GpuException("The CPU backend currently supports the production eqtlReal operation only");
		if (source == null || precision == null)
			throw new IllegalArgumentException("Kernel source and precision must not be null");
		Matcher dataType = DATA_TYPE.matcher(source);
		String expectedType = precision == GpuPrecision.FP32 ? "float" : "double";
		if (!dataType.find() || !expectedType.equals(dataType.group(1)))
			throw new GpuException("Kernel DATATYPE does not match " + precision.optionName());
		Matcher nMinusOne = N_MINUS_ONE.matcher(source);
		if (!nMinusOne.find())
			throw new GpuException("The eqtlReal operation is missing a positive N_MIN_1 definition");
		try {
			long value = Long.parseLong(nMinusOne.group(1));
			if (value <= 0) throw new NumberFormatException();
			normalization = 1.0 / value;
		} catch (NumberFormatException e) {
			throw new GpuException("Invalid N_MIN_1 value in the eqtlReal operation", e);
		}
		this.precision = precision;
		kernelReady = true;
	}

	@Override
	public synchronized double[] executeDoubleKernel(double[] inputA, double[] inputB,
		int outputElements, long localMemoryBytes, int widthA, int widthB,
		long[] globalWorkSize, long[] localWorkSize) {
		ensureReady(GpuPrecision.FP64);
		validateArguments(inputA, inputB, outputElements, widthA, widthB,
			globalWorkSize, localWorkSize);
		int columns = Math.toIntExact(globalWorkSize[0]);
		int rows = Math.toIntExact(globalWorkSize[1]);
		long started = profilingEnabled ? System.nanoTime() : 0;
		double[] output = new double[outputElements];
		long setup = elapsed(started);
		started = profilingEnabled ? System.nanoTime() : 0;
		try {
			engine.multiplyDouble(inputA, inputB, output, rows, columns, widthA,
				widthA, widthB, widthB, normalization, 0);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CPU association execution failed using " + engine.description(), e);
		}
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setup, 0, elapsed(started), 0, 0, 0)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	@Override
	public synchronized float[] executeFloatKernel(float[] inputA, float[] inputB,
		int outputElements, long localMemoryBytes, int widthA, int widthB,
		long[] globalWorkSize, long[] localWorkSize) {
		ensureReady(GpuPrecision.FP32);
		validateArguments(inputA, inputB, outputElements, widthA, widthB,
			globalWorkSize, localWorkSize);
		int columns = Math.toIntExact(globalWorkSize[0]);
		int rows = Math.toIntExact(globalWorkSize[1]);
		long started = profilingEnabled ? System.nanoTime() : 0;
		float[] output = new float[outputElements];
		long setup = elapsed(started);
		started = profilingEnabled ? System.nanoTime() : 0;
		try {
			engine.multiplyFloat(inputA, inputB, output, rows, columns, widthA,
				widthA, widthB, widthB, (float) normalization, 0);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CPU FP32 association execution failed using " + engine.description(), e);
		}
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setup, 0, elapsed(started), 0, 0, 0)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	@Override
	public synchronized double[] residualizeDoubleRows(double[] values, double[] q,
		int rowCount, int sampleCount, int covariateRank) {
		ensureOpen();
		validateResidualization(values, q, rowCount, sampleCount, covariateRank);
		long started = profilingEnabled ? System.nanoTime() : 0;
		double[] output = values.clone();
		double[] coefficients = new double[Math.multiplyExact(rowCount, covariateRank)];
		long setup = elapsed(started);
		started = profilingEnabled ? System.nanoTime() : 0;
		try {
			engine.multiplyDouble(values, q, coefficients, rowCount, covariateRank,
				sampleCount, sampleCount, covariateRank, covariateRank, 1, 0);
			engine.multiplyDoubleRightTransposed(coefficients, q, output, rowCount,
				sampleCount, covariateRank, covariateRank, covariateRank, sampleCount, -1, 1);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CPU FP64 residualization failed using " + engine.description(), e);
		}
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setup, 0, elapsed(started), 0, 0, 0)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	@Override
	public synchronized float[] residualizeFloatRows(float[] values, float[] q,
		int rowCount, int sampleCount, int covariateRank) {
		ensureOpen();
		validateResidualization(values, q, rowCount, sampleCount, covariateRank);
		long started = profilingEnabled ? System.nanoTime() : 0;
		float[] output = values.clone();
		float[] coefficients = new float[Math.multiplyExact(rowCount, covariateRank)];
		long setup = elapsed(started);
		started = profilingEnabled ? System.nanoTime() : 0;
		try {
			engine.multiplyFloat(values, q, coefficients, rowCount, covariateRank,
				sampleCount, sampleCount, covariateRank, covariateRank, 1, 0);
			engine.multiplyFloatRightTransposed(coefficients, q, output, rowCount,
				sampleCount, covariateRank, covariateRank, covariateRank, sampleCount, -1, 1);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CPU FP32 residualization failed using " + engine.description(), e);
		}
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setup, 0, elapsed(started), 0, 0, 0)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	@Override
	public synchronized void close() {
		closed = true;
		kernelReady = false;
		lastExecutionMetrics = GpuExecutionMetrics.EMPTY;
	}

	private void ensureReady(GpuPrecision expected) {
		ensureOpen();
		if (!kernelReady)
			throw new GpuException("No CPU eQTL operation has been prepared");
		if (precision != expected)
			throw new GpuException("The CPU operation was prepared for " + precision.optionName());
	}

	private void ensureOpen() {
		if (closed) throw new GpuException("CPU compute context is closed");
	}

	private long elapsed(long started) {
		return profilingEnabled ? System.nanoTime() - started : 0;
	}

	private static void validateArguments(Object inputA, Object inputB, int outputElements,
		int widthA, int widthB, long[] globalWorkSize, long[] localWorkSize) {
		if (inputA == null || inputB == null || outputElements <= 0 || widthA <= 0 || widthB <= 0)
			throw new IllegalArgumentException("Matrix buffers and widths must be non-empty");
		if (globalWorkSize == null || localWorkSize == null
			|| globalWorkSize.length != 2 || localWorkSize.length != 2)
			throw new IllegalArgumentException("The eQTL operation requires two-dimensional work sizes");
		if (globalWorkSize[0] <= 0 || globalWorkSize[1] <= 0
			|| globalWorkSize[0] > Integer.MAX_VALUE || globalWorkSize[1] > Integer.MAX_VALUE)
			throw new IllegalArgumentException("Active work dimensions must be positive 32-bit values");
		int inputALength = inputA instanceof double[] values ? values.length : ((float[]) inputA).length;
		int inputBLength = inputB instanceof double[] values ? values.length : ((float[]) inputB).length;
		long rows = globalWorkSize[1];
		if (globalWorkSize[0] > widthB
			|| Math.multiplyExact(rows, widthA) > inputALength
			|| Math.multiplyExact((long) widthA, widthB) > inputBLength
			|| Math.multiplyExact(rows, widthB) > outputElements)
			throw new IllegalArgumentException("Work dimensions exceed the supplied eQTL matrix buffers");
	}

	private static void validateResidualization(Object values, Object q,
		int rows, int samples, int rank) {
		if (values == null || q == null || rows <= 0 || samples <= 0 || rank <= 0 || rank > samples)
			throw new IllegalArgumentException("Invalid residualization matrix dimensions");
		int valueLength = values instanceof double[] array ? array.length : ((float[]) values).length;
		int qLength = q instanceof double[] array ? array.length : ((float[]) q).length;
		if (valueLength != Math.multiplyExact(rows, samples)
			|| qLength != Math.multiplyExact(samples, rank))
			throw new IllegalArgumentException("Residualization buffers do not match their matrix dimensions");
	}
}
