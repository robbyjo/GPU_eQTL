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
	private double normalization;
	private GpuPrecision precision = GpuPrecision.FP64;
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
			selectDevice();
			ensureInputABuffer(inputABytes);
			ensureInputBBuffer(inputBBytes);
			ensureOutputBuffer(outputBytes);
			JCuda.cudaMemcpy(inputABuffer, Pointer.to(inputA), inputABytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemcpy(inputBBuffer, Pointer.to(inputB), inputBBytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemset(outputBuffer, 0, outputBytes);

			Pointer alpha = Pointer.to(new double[] { normalization });
			Pointer beta = Pointer.to(new double[] { 0.0 });
			/*
			 * The Java arrays are row-major. Reinterpreting them as column-major lets
			 * cuBLAS compute C^T = B^T A^T without allocating transpose buffers.
			 */
			JCublas2.cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
				activeSnpColumns, activeTraits, widthA,
				alpha, inputBBuffer, widthB,
				inputABuffer, widthA,
				beta, outputBuffer, widthB);

			double[] output = new double[outputElements];
			JCuda.cudaMemcpy(Pointer.to(output), outputBuffer, outputBytes, cudaMemcpyDeviceToHost);
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
			selectDevice();
			ensureInputABuffer(inputABytes);
			ensureInputBBuffer(inputBBytes);
			ensureOutputBuffer(outputBytes);
			JCuda.cudaMemcpy(inputABuffer, Pointer.to(inputA), inputABytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemcpy(inputBBuffer, Pointer.to(inputB), inputBBytes, cudaMemcpyHostToDevice);
			JCuda.cudaMemset(outputBuffer, 0, outputBytes);

			Pointer alpha = Pointer.to(new float[] { (float) normalization });
			Pointer beta = Pointer.to(new float[] { 0.0f });
			JCublas2.cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
				activeSnpColumns, activeTraits, widthA,
				alpha, inputBBuffer, widthB,
				inputABuffer, widthA,
				beta, outputBuffer, widthB);

			float[] output = new float[outputElements];
			JCuda.cudaMemcpy(Pointer.to(output), outputBuffer, outputBytes, cudaMemcpyDeviceToHost);
			return output;
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CUDA/cuBLAS FP32 eQTL execution failed on " + device.getName(), e);
		}
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
			if (inputABuffer != null) JCuda.cudaFree(inputABuffer);
			if (inputBBuffer != null) JCuda.cudaFree(inputBBuffer);
			if (outputBuffer != null) JCuda.cudaFree(outputBuffer);
			if (handle != null) JCublas2.cublasDestroy(handle);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("Could not release the CUDA resources for " + device.getName(), e);
		} finally {
			inputABuffer = inputBBuffer = outputBuffer = null;
			handle = null;
			closed = true;
		}
	}
}
