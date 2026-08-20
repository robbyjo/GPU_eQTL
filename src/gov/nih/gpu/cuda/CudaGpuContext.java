/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cuda;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuExecutionMetrics;
import gov.nih.gpu.GpuException;
import gov.nih.gpu.GpuPrecision;

import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas2;
import jcuda.jcublas.cublasHandle;
import jcuda.runtime.JCuda;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static jcuda.jcublas.cublasOperation.CUBLAS_OP_N;
import static jcuda.jcublas.cublasOperation.CUBLAS_OP_T;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;

/** CUDA execution context backed by cuBLAS SGEMM/DGEMM. */
final class CudaGpuContext implements GpuContext {
	private static final Pattern N_MINUS_ONE = Pattern.compile(
		"(?m)^\\s*#define\\s+N_MIN_1\\s+(\\d+)\\s*$");
	private static final Pattern DATA_TYPE = Pattern.compile(
		"(?m)^\\s*#define\\s+DATATYPE\\s+(float|double)\\s*$");

	private final CudaGpuDevice device;
	private cublasHandle handle;
	private Pointer inputABuffer;
	private Pointer inputBBuffer;
	private Pointer outputBuffer;
	private long inputACapacity;
	private long inputBCapacity;
	private long outputCapacity;
	private Pointer residualValuesBuffer;
	private Pointer residualQBuffer;
	private Pointer residualCoefficientsBuffer;
	private long residualValuesCapacity;
	private long residualQCapacity;
	private long residualCoefficientsCapacity;
	private Object residualProjectionSource;
	private GpuPrecision residualProjectionPrecision;
	private int residualProjectionSamples;
	private int residualProjectionRank;
	private double normalization;
	private GpuPrecision precision = GpuPrecision.FP64;
	private boolean profilingEnabled;
	private GpuExecutionMetrics lastExecutionMetrics = GpuExecutionMetrics.EMPTY;
	private boolean kernelReady;
	private boolean closed;

	CudaGpuContext(CudaGpuDevice device) {
		this.device = device;
		try {
			selectDevice();
			handle = new cublasHandle();
			JCublas2.cublasCreate(handle);
		} catch (RuntimeException | LinkageError e) {
			handle = null;
			throw new GpuException("Could not create a cuBLAS context for " + device.getName(), e);
		}
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
		if (precision == null)
			throw new IllegalArgumentException("precision must not be null");
		if (!"eqtlReal".equals(kernelName)) {
			throw new GpuException("The CUDA backend currently supports the production eqtlReal operation only; "
				+ "select -Deqtl.gpu.backend=opencl for the legacy categorical-SNP kernel");
		}
		if (source == null) {
			throw new IllegalArgumentException("Kernel source must not be null");
		}
		Matcher dataType = DATA_TYPE.matcher(source);
		String expectedType = precision == GpuPrecision.FP32 ? "float" : "double";
		if (!dataType.find() || !expectedType.equals(dataType.group(1)))
			throw new GpuException("Kernel DATATYPE does not match " + precision.optionName());
		Matcher matcher = N_MINUS_ONE.matcher(source);
		if (!matcher.find()) {
			throw new GpuException("The eqtlReal operation is missing a positive N_MIN_1 definition");
		}
		long nMinusOne;
		try {
			nMinusOne = Long.parseLong(matcher.group(1));
		} catch (NumberFormatException e) {
			throw new GpuException("Invalid N_MIN_1 value in the eqtlReal operation", e);
		}
		if (nMinusOne <= 0) {
			throw new GpuException("N_MIN_1 must be positive for the eqtlReal operation");
		}
		normalization = 1.0 / nMinusOne;
		this.precision = precision;
		kernelReady = true;
	}

	@Override
	public synchronized double[] executeDoubleKernel(
		double[] inputA,
		double[] inputB,
		int outputElements,
		long localMemoryBytes,
		int widthA,
		int widthB,
		long[] globalWorkSize,
		long[] localWorkSize) {

		ensureOpen();
		if (!kernelReady) {
			throw new GpuException("No CUDA eQTL operation has been prepared for " + device.getName());
		}
		if (precision != GpuPrecision.FP64)
			throw new GpuException("The CUDA operation was prepared for " + precision.optionName());
		if (inputA == null || inputB == null)
			throw new IllegalArgumentException("Kernel buffers must be non-null");
		validateArguments(inputA.length, inputB.length, outputElements,
			widthA, widthB, globalWorkSize, localWorkSize);
		int activeSnpColumns = Math.toIntExact(globalWorkSize[0]);
		int activeTraits = Math.toIntExact(globalWorkSize[1]);

		long inputABytes = Math.multiplyExact((long) inputA.length, Sizeof.DOUBLE);
		long inputBBytes = Math.multiplyExact((long) inputB.length, Sizeof.DOUBLE);
		long outputBytes = Math.multiplyExact((long) outputElements, Sizeof.DOUBLE);
		try {
			long phaseStart = profilingEnabled ? System.nanoTime() : 0;
			selectDevice();
			ensureInputABuffer(inputABytes);
			ensureInputBBuffer(inputBBytes);
			ensureOutputBuffer(outputBytes);
			double[] output = new double[outputElements];
			long setupNanos = elapsed(phaseStart);

			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(inputABuffer, Pointer.to(inputA), inputABytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemcpy(inputBBuffer, Pointer.to(inputB), inputBBytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemset(outputBuffer, 0, outputBytes);
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long uploadNanos = elapsed(phaseStart);

			Pointer alpha = Pointer.to(new double[] { normalization });
			Pointer beta = Pointer.to(new double[] { 0.0 });
			/*
			 * The Java arrays are row-major. Reinterpreting them as column-major lets
			 * cuBLAS compute C^T = B^T A^T without allocating transpose buffers.
			 */
			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCublas2.cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
				activeSnpColumns, activeTraits, widthA,
				alpha, inputBBuffer, widthB,
				inputABuffer, widthA,
				beta, outputBuffer, widthB);
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long computeNanos = elapsed(phaseStart);

			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(Pointer.to(output), outputBuffer, outputBytes, cudaMemcpyDeviceToHost);
			long downloadNanos = elapsed(phaseStart);
			lastExecutionMetrics = profilingEnabled
				? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
					inputABytes + inputBBytes, outputBytes)
				: GpuExecutionMetrics.EMPTY;
			return output;
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CUDA/cuBLAS eQTL execution failed on " + device.getName(), e);
		}
	}

	@Override
	public synchronized float[] executeFloatKernel(
		float[] inputA,
		float[] inputB,
		int outputElements,
		long localMemoryBytes,
		int widthA,
		int widthB,
		long[] globalWorkSize,
		long[] localWorkSize) {

		ensureOpen();
		if (!kernelReady)
			throw new GpuException("No CUDA eQTL operation has been prepared for " + device.getName());
		if (precision != GpuPrecision.FP32)
			throw new GpuException("The CUDA operation was prepared for " + precision.optionName());
		if (inputA == null || inputB == null)
			throw new IllegalArgumentException("Kernel buffers must be non-null");
		validateArguments(inputA.length, inputB.length, outputElements,
			widthA, widthB, globalWorkSize, localWorkSize);
		int activeSnpColumns = Math.toIntExact(globalWorkSize[0]);
		int activeTraits = Math.toIntExact(globalWorkSize[1]);

		long inputABytes = Math.multiplyExact((long) inputA.length, Sizeof.FLOAT);
		long inputBBytes = Math.multiplyExact((long) inputB.length, Sizeof.FLOAT);
		long outputBytes = Math.multiplyExact((long) outputElements, Sizeof.FLOAT);
		try {
			long phaseStart = profilingEnabled ? System.nanoTime() : 0;
			selectDevice();
			ensureInputABuffer(inputABytes);
			ensureInputBBuffer(inputBBytes);
			ensureOutputBuffer(outputBytes);
			float[] output = new float[outputElements];
			long setupNanos = elapsed(phaseStart);

			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(inputABuffer, Pointer.to(inputA), inputABytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemcpy(inputBBuffer, Pointer.to(inputB), inputBBytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemset(outputBuffer, 0, outputBytes);
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long uploadNanos = elapsed(phaseStart);

			Pointer alpha = Pointer.to(new float[] { (float) normalization });
			Pointer beta = Pointer.to(new float[] { 0.0f });
			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCublas2.cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
				activeSnpColumns, activeTraits, widthA,
				alpha, inputBBuffer, widthB,
				inputABuffer, widthA,
				beta, outputBuffer, widthB);
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long computeNanos = elapsed(phaseStart);

			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(Pointer.to(output), outputBuffer, outputBytes, cudaMemcpyDeviceToHost);
			long downloadNanos = elapsed(phaseStart);
			lastExecutionMetrics = profilingEnabled
				? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
					inputABytes + inputBBytes, outputBytes)
				: GpuExecutionMetrics.EMPTY;
			return output;
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CUDA/cuBLAS FP32 eQTL execution failed on " + device.getName(), e);
		}
	}

	@Override
	public synchronized double[] residualizeDoubleRows(double[] rowMajorValues, double[] rowMajorQ,
		int rowCount, int sampleCount, int covariateRank) {
		ensureOpen();
		validateResidualizationArguments(rowMajorValues == null ? 0 : rowMajorValues.length,
			rowMajorQ == null ? 0 : rowMajorQ.length, rowCount, sampleCount, covariateRank);
		long valuesBytes = Math.multiplyExact((long) rowMajorValues.length, Sizeof.DOUBLE);
		long qBytes = Math.multiplyExact((long) rowMajorQ.length, Sizeof.DOUBLE);
		long coefficientBytes = Math.multiplyExact((long) rowCount * covariateRank, Sizeof.DOUBLE);
		try {
			long phaseStart = profilingEnabled ? System.nanoTime() : 0;
			selectDevice();
			ensureResidualizationBuffers(valuesBytes, qBytes, coefficientBytes);
			double[] output = new double[rowMajorValues.length];
			long setupNanos = elapsed(phaseStart);

			boolean uploadProjection = !isCachedProjection(rowMajorQ, GpuPrecision.FP64,
				sampleCount, covariateRank);
			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(residualValuesBuffer, Pointer.to(rowMajorValues), valuesBytes, cudaMemcpyHostToDevice);
			if (uploadProjection) {
				JCuda.cudaMemcpy(residualQBuffer, Pointer.to(rowMajorQ), qBytes, cudaMemcpyHostToDevice);
				rememberProjection(rowMajorQ, GpuPrecision.FP64, sampleCount, covariateRank);
			}
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long uploadNanos = elapsed(phaseStart);

			Pointer one = Pointer.to(new double[] { 1.0 });
			Pointer zero = Pointer.to(new double[] { 0.0 });
			Pointer minusOne = Pointer.to(new double[] { -1.0 });
			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			// T^T = Q^T Y^T, where row-major host arrays are column-major transposes.
			JCublas2.cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
				covariateRank, rowCount, sampleCount,
				one, residualQBuffer, covariateRank,
				residualValuesBuffer, sampleCount,
				zero, residualCoefficientsBuffer, covariateRank);
			// Y_res^T = Y^T - Q T^T, updating the value buffer in place.
			JCublas2.cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N,
				sampleCount, rowCount, covariateRank,
				minusOne, residualQBuffer, covariateRank,
				residualCoefficientsBuffer, covariateRank,
				one, residualValuesBuffer, sampleCount);
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long computeNanos = elapsed(phaseStart);

			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(Pointer.to(output), residualValuesBuffer, valuesBytes, cudaMemcpyDeviceToHost);
			long downloadNanos = elapsed(phaseStart);
			lastExecutionMetrics = profilingEnabled
				? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
					valuesBytes + (uploadProjection ? qBytes : 0), valuesBytes)
				: GpuExecutionMetrics.EMPTY;
			return output;
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CUDA/cuBLAS FP64 residualization failed on " + device.getName(), e);
		}
	}

	@Override
	public synchronized float[] residualizeFloatRows(float[] rowMajorValues, float[] rowMajorQ,
		int rowCount, int sampleCount, int covariateRank) {
		ensureOpen();
		validateResidualizationArguments(rowMajorValues == null ? 0 : rowMajorValues.length,
			rowMajorQ == null ? 0 : rowMajorQ.length, rowCount, sampleCount, covariateRank);
		long valuesBytes = Math.multiplyExact((long) rowMajorValues.length, Sizeof.FLOAT);
		long qBytes = Math.multiplyExact((long) rowMajorQ.length, Sizeof.FLOAT);
		long coefficientBytes = Math.multiplyExact((long) rowCount * covariateRank, Sizeof.FLOAT);
		try {
			long phaseStart = profilingEnabled ? System.nanoTime() : 0;
			selectDevice();
			ensureResidualizationBuffers(valuesBytes, qBytes, coefficientBytes);
			float[] output = new float[rowMajorValues.length];
			long setupNanos = elapsed(phaseStart);

			boolean uploadProjection = !isCachedProjection(rowMajorQ, GpuPrecision.FP32,
				sampleCount, covariateRank);
			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(residualValuesBuffer, Pointer.to(rowMajorValues), valuesBytes, cudaMemcpyHostToDevice);
			if (uploadProjection) {
				JCuda.cudaMemcpy(residualQBuffer, Pointer.to(rowMajorQ), qBytes, cudaMemcpyHostToDevice);
				rememberProjection(rowMajorQ, GpuPrecision.FP32, sampleCount, covariateRank);
			}
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long uploadNanos = elapsed(phaseStart);

			Pointer one = Pointer.to(new float[] { 1.0f });
			Pointer zero = Pointer.to(new float[] { 0.0f });
			Pointer minusOne = Pointer.to(new float[] { -1.0f });
			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCublas2.cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
				covariateRank, rowCount, sampleCount,
				one, residualQBuffer, covariateRank,
				residualValuesBuffer, sampleCount,
				zero, residualCoefficientsBuffer, covariateRank);
			JCublas2.cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N,
				sampleCount, rowCount, covariateRank,
				minusOne, residualQBuffer, covariateRank,
				residualCoefficientsBuffer, covariateRank,
				one, residualValuesBuffer, sampleCount);
			if (profilingEnabled)
				JCuda.cudaDeviceSynchronize();
			long computeNanos = elapsed(phaseStart);

			phaseStart = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(Pointer.to(output), residualValuesBuffer, valuesBytes, cudaMemcpyDeviceToHost);
			long downloadNanos = elapsed(phaseStart);
			lastExecutionMetrics = profilingEnabled
				? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
					valuesBytes + (uploadProjection ? qBytes : 0), valuesBytes)
				: GpuExecutionMetrics.EMPTY;
			return output;
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CUDA/cuBLAS FP32 residualization failed on " + device.getName(), e);
		}
	}

	private long elapsed(long startedAtNanos) {
		return profilingEnabled ? System.nanoTime() - startedAtNanos : 0;
	}

	private static void validateArguments(int inputALength, int inputBLength, int outputElements,
			int widthA, int widthB, long[] globalWorkSize, long[] localWorkSize) {
		if (inputALength <= 0 || inputBLength <= 0 || outputElements <= 0) {
			throw new IllegalArgumentException("Kernel buffers must be non-empty");
		}
		if (widthA <= 0 || widthB <= 0) {
			throw new IllegalArgumentException("Matrix widths must be positive");
		}
		if (globalWorkSize == null || localWorkSize == null
				|| globalWorkSize.length != 2 || localWorkSize.length != 2) {
			throw new IllegalArgumentException("The eQTL operation requires two-dimensional work sizes");
		}
		if (globalWorkSize[0] <= 0 || globalWorkSize[1] <= 0
				|| globalWorkSize[0] > Integer.MAX_VALUE || globalWorkSize[1] > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Active work dimensions must be positive 32-bit values");
		}
		int activeSnpColumns = Math.toIntExact(globalWorkSize[0]);
		int activeTraits = Math.toIntExact(globalWorkSize[1]);
		long requiredInputA = Math.multiplyExact((long) activeTraits, widthA);
		long requiredInputB = Math.multiplyExact((long) widthA, widthB);
		long requiredOutput = Math.multiplyExact((long) activeTraits, widthB);
		if (activeSnpColumns > widthB || requiredInputA > inputALength
				|| requiredInputB > inputBLength || requiredOutput > outputElements)
			throw new IllegalArgumentException("Work dimensions exceed the supplied eQTL matrix buffers");
	}

	private static void validateResidualizationArguments(int valuesLength, int qLength,
		int rowCount, int sampleCount, int covariateRank) {
		if (rowCount <= 0 || sampleCount <= 0 || covariateRank <= 0 || covariateRank > sampleCount)
			throw new IllegalArgumentException("Invalid residualization matrix dimensions");
		if (valuesLength != Math.multiplyExact(rowCount, sampleCount)
			|| qLength != Math.multiplyExact(sampleCount, covariateRank))
			throw new IllegalArgumentException("Residualization buffers do not match their matrix dimensions");
	}

	private void ensureResidualizationBuffers(long valuesBytes, long qBytes, long coefficientBytes) {
		if (residualValuesBuffer == null || residualValuesCapacity < valuesBytes) {
			residualValuesBuffer = replaceBuffer(residualValuesBuffer, valuesBytes);
			residualValuesCapacity = valuesBytes;
		}
		if (residualQBuffer == null || residualQCapacity < qBytes) {
			residualQBuffer = replaceBuffer(residualQBuffer, qBytes);
			residualQCapacity = qBytes;
			residualProjectionSource = null;
		}
		if (residualCoefficientsBuffer == null || residualCoefficientsCapacity < coefficientBytes) {
			residualCoefficientsBuffer = replaceBuffer(residualCoefficientsBuffer, coefficientBytes);
			residualCoefficientsCapacity = coefficientBytes;
		}
	}

	private boolean isCachedProjection(Object source, GpuPrecision projectionPrecision,
		int sampleCount, int covariateRank) {
		return residualProjectionSource == source && residualProjectionPrecision == projectionPrecision
			&& residualProjectionSamples == sampleCount && residualProjectionRank == covariateRank;
	}

	private void rememberProjection(Object source, GpuPrecision projectionPrecision,
		int sampleCount, int covariateRank) {
		residualProjectionSource = source;
		residualProjectionPrecision = projectionPrecision;
		residualProjectionSamples = sampleCount;
		residualProjectionRank = covariateRank;
	}

	@Override
	public synchronized void releaseResidualizationResources() {
		ensureOpen();
		selectDevice();
		if (residualValuesBuffer != null) JCuda.cudaFree(residualValuesBuffer);
		if (residualQBuffer != null) JCuda.cudaFree(residualQBuffer);
		if (residualCoefficientsBuffer != null) JCuda.cudaFree(residualCoefficientsBuffer);
		residualValuesBuffer = residualQBuffer = residualCoefficientsBuffer = null;
		residualValuesCapacity = residualQCapacity = residualCoefficientsCapacity = 0;
		residualProjectionSource = null;
		residualProjectionPrecision = null;
		residualProjectionSamples = residualProjectionRank = 0;
	}

	private void ensureInputABuffer(long requiredBytes) {
		if (inputABuffer != null && inputACapacity >= requiredBytes) {
			return;
		}
		inputABuffer = replaceBuffer(inputABuffer, requiredBytes);
		inputACapacity = requiredBytes;
	}

	private void ensureInputBBuffer(long requiredBytes) {
		if (inputBBuffer != null && inputBCapacity >= requiredBytes) {
			return;
		}
		inputBBuffer = replaceBuffer(inputBBuffer, requiredBytes);
		inputBCapacity = requiredBytes;
	}

	private void ensureOutputBuffer(long requiredBytes) {
		if (outputBuffer != null && outputCapacity >= requiredBytes) {
			return;
		}
		outputBuffer = replaceBuffer(outputBuffer, requiredBytes);
		outputCapacity = requiredBytes;
	}

	private Pointer replaceBuffer(Pointer current, long requiredBytes) {
		if (current != null) {
			JCuda.cudaFree(current);
		}
		Pointer replacement = new Pointer();
		JCuda.cudaMalloc(replacement, requiredBytes);
		return replacement;
	}

	private void selectDevice() {
		JCuda.cudaSetDevice(device.ordinal());
	}

	private void ensureOpen() {
		if (closed) {
			throw new GpuException("GPU context is closed");
		}
	}

	@Override
	public synchronized void close() {
		if (closed) {
			return;
		}
		try {
			selectDevice();
			if (residualValuesBuffer != null) JCuda.cudaFree(residualValuesBuffer);
			if (residualQBuffer != null) JCuda.cudaFree(residualQBuffer);
			if (residualCoefficientsBuffer != null) JCuda.cudaFree(residualCoefficientsBuffer);
			if (inputABuffer != null) JCuda.cudaFree(inputABuffer);
			if (inputBBuffer != null) JCuda.cudaFree(inputBBuffer);
			if (outputBuffer != null) JCuda.cudaFree(outputBuffer);
			if (handle != null) JCublas2.cublasDestroy(handle);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("Could not release the CUDA resources for " + device.getName(), e);
		} finally {
			inputABuffer = inputBBuffer = outputBuffer = null;
			residualValuesBuffer = residualQBuffer = residualCoefficientsBuffer = null;
			handle = null;
			closed = true;
		}
	}
}
