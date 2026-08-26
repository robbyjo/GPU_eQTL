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
import gov.nih.gpu.GpuPatternStatisticsPlan;
import gov.nih.gpu.GpuPatternStatisticsResult;
import gov.nih.gpu.GpuPatternStatisticsSupport;
import gov.nih.gpu.GpuPrecision;

import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUcontext;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import jcuda.jcublas.JCublas2;
import jcuda.jcublas.cublasHandle;
import jcuda.nvrtc.JNvrtc;
import jcuda.nvrtc.nvrtcProgram;
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
	private CUmodule patternModule;
	private CUfunction patternFinalizeFunction;
	private Pointer patternUpperBuffer;
	private Pointer patternSumsBuffer;
	private Pointer patternObservedBuffer;
	private Pointer patternCompactBuffer;
	private long patternUpperCapacity;
	private long patternSumsCapacity;
	private long patternObservedCapacity;
	private long patternCompactCapacity;
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
	public synchronized GpuPatternStatisticsResult computePatternStatisticsDouble(
		double[] aggregateInputs, int paddedSamples, int activeVariants, int variantCapacity,
		GpuPatternStatisticsPlan plan, boolean meanFill, int patternsPerBatch,
		int workGroupSize) {
		ensureOpen();
		if (!kernelReady || precision != GpuPrecision.FP64)
			throw new GpuException("The CUDA FP64 association operation must be compiled first");
		GpuPatternStatisticsSupport.validate(aggregateInputs, paddedSamples, activeVariants,
			variantCapacity, plan, patternsPerBatch, workGroupSize);
		if (plan.rank() > 64)
			throw new GpuException("CUDA compact pattern finalization supports at most 64 design columns");
		try {
			selectDevice();
			ensurePatternKernel();
			int rank = plan.rank();
			int tripledCapacity = Math.multiplyExact(variantCapacity, 3);
			int maximumRows = GpuPatternStatisticsSupport.roundUp(
				patternsPerBatch * plan.rowsPerPattern(), workGroupSize);
			double[] designMasks = new double[Math.multiplyExact(maximumRows, paddedSamples)];
			double[] upper = new double[Math.multiplyExact(patternsPerBatch, rank * rank)];
			double[] designSums = new double[Math.multiplyExact(patternsPerBatch, rank)];
			int[] observedCounts = new int[patternsPerBatch];
			double[] compact = new double[Math.multiplyExact(
				Math.multiplyExact(plan.patternCount(), activeVariants),
				GpuPatternStatisticsResult.VALUES_PER_CELL)];
			long aggregateBytes = Math.multiplyExact((long) aggregateInputs.length, Sizeof.DOUBLE);
			long maximumDesignBytes = Math.multiplyExact((long) designMasks.length, Sizeof.DOUBLE);
			long maximumProductBytes = Math.multiplyExact((long) maximumRows * tripledCapacity,
				Sizeof.DOUBLE);
			long maximumUpperBytes = Math.multiplyExact((long) upper.length, Sizeof.DOUBLE);
			long maximumSumsBytes = Math.multiplyExact((long) designSums.length, Sizeof.DOUBLE);
			long maximumObservedBytes = Math.multiplyExact((long) observedCounts.length, Sizeof.INT);
			long compactBytes = Math.multiplyExact((long) compact.length, Sizeof.DOUBLE);

			long setupStarted = profilingEnabled ? System.nanoTime() : 0;
			ensureInputABuffer(maximumDesignBytes);
			ensureInputBBuffer(aggregateBytes);
			ensureOutputBuffer(maximumProductBytes);
			patternUpperBuffer = ensurePatternBuffer(patternUpperBuffer,
				maximumUpperBytes, patternUpperCapacity);
			patternUpperCapacity = Math.max(patternUpperCapacity, maximumUpperBytes);
			patternSumsBuffer = ensurePatternBuffer(patternSumsBuffer,
				maximumSumsBytes, patternSumsCapacity);
			patternSumsCapacity = Math.max(patternSumsCapacity, maximumSumsBytes);
			patternObservedBuffer = ensurePatternBuffer(patternObservedBuffer,
				maximumObservedBytes, patternObservedCapacity);
			patternObservedCapacity = Math.max(patternObservedCapacity, maximumObservedBytes);
			patternCompactBuffer = ensurePatternBuffer(patternCompactBuffer,
				compactBytes, patternCompactCapacity);
			patternCompactCapacity = Math.max(patternCompactCapacity, compactBytes);
			long setupNanos = elapsed(setupStarted);

			long uploadNanos = 0;
			long computeNanos = 0;
			long uploadedBytes = aggregateBytes;
			long phaseStarted = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(inputBBuffer, Pointer.to(aggregateInputs), aggregateBytes,
				cudaMemcpyHostToDevice);
			if (profilingEnabled) JCuda.cudaDeviceSynchronize();
			uploadNanos += elapsed(phaseStarted);
			Pointer alpha = Pointer.to(new double[] {normalization});
			Pointer beta = Pointer.to(new double[] {0.0});
			for (int first = 0; first < plan.patternCount(); first += patternsPerBatch) {
				int count = Math.min(patternsPerBatch, plan.patternCount() - first);
				int rowCapacity = GpuPatternStatisticsSupport.roundUp(
					count * plan.rowsPerPattern(), workGroupSize);
				plan.fillBatch(first, count, paddedSamples, rowCapacity, designMasks,
					upper, designSums, observedCounts);
				long designBytes = Math.multiplyExact((long) rowCapacity * paddedSamples,
					Sizeof.DOUBLE);
				long upperBytes = Math.multiplyExact((long) count * rank * rank, Sizeof.DOUBLE);
				long sumsBytes = Math.multiplyExact((long) count * rank, Sizeof.DOUBLE);
				long observedBytes = Math.multiplyExact((long) count, Sizeof.INT);
				phaseStarted = profilingEnabled ? System.nanoTime() : 0;
				JCuda.cudaMemcpy(inputABuffer, Pointer.to(designMasks), designBytes,
					cudaMemcpyHostToDevice);
				JCuda.cudaMemcpy(patternUpperBuffer, Pointer.to(upper), upperBytes,
					cudaMemcpyHostToDevice);
				JCuda.cudaMemcpy(patternSumsBuffer, Pointer.to(designSums), sumsBytes,
					cudaMemcpyHostToDevice);
				JCuda.cudaMemcpy(patternObservedBuffer, Pointer.to(observedCounts), observedBytes,
					cudaMemcpyHostToDevice);
				if (profilingEnabled) JCuda.cudaDeviceSynchronize();
				uploadNanos += elapsed(phaseStarted);
				uploadedBytes += designBytes + upperBytes + sumsBytes + observedBytes;

				phaseStarted = profilingEnabled ? System.nanoTime() : 0;
				JCublas2.cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
					activeVariants * 3, rowCapacity, paddedSamples,
					alpha, inputBBuffer, tripledCapacity,
					inputABuffer, paddedSamples,
					beta, outputBuffer, tripledCapacity);
				launchPatternFinalizer(rank, tripledCapacity, activeVariants, count,
					first, meanFill);
				JCuda.cudaDeviceSynchronize();
				computeNanos += elapsed(phaseStarted);
			}
			phaseStarted = profilingEnabled ? System.nanoTime() : 0;
			JCuda.cudaMemcpy(Pointer.to(compact), patternCompactBuffer, compactBytes,
				cudaMemcpyDeviceToHost);
			long downloadNanos = elapsed(phaseStarted);
			lastExecutionMetrics = profilingEnabled
				? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos,
					downloadNanos, uploadedBytes, compactBytes)
				: GpuExecutionMetrics.EMPTY;
			return new GpuPatternStatisticsResult(plan.patternCount(), activeVariants, compact);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("CUDA compact pattern-statistics execution failed on "
				+ device.getName(), e);
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

	private void ensurePatternKernel() {
		if (patternFinalizeFunction != null) return;
		JCudaDriver.cuInit(0);
		CUcontext current = new CUcontext();
		JCudaDriver.cuCtxGetCurrent(current);
		nvrtcProgram program = new nvrtcProgram();
		JNvrtc.nvrtcCreateProgram(program, patternKernelSource(),
			"gpu_eqtl_pattern_statistics.cu", 0, null, null);
		try {
			String architecture = "--gpu-architecture=compute_" + device.computeMajor()
				+ device.computeMinor();
			try {
				JNvrtc.nvrtcCompileProgram(program, 1, new String[] {architecture});
			} catch (RuntimeException e) {
				long[] logSize = new long[1];
				JNvrtc.nvrtcGetProgramLogSize(program, logSize);
				String[] log = new String[1];
				JNvrtc.nvrtcGetProgramLog(program, log);
				throw new GpuException("NVRTC pattern-statistics compilation failed: "
					+ (log[0] == null ? "no compiler log" : log[0]), e);
			}
			String[] ptx = new String[1];
			JNvrtc.nvrtcGetPTX(program, ptx);
			patternModule = new CUmodule();
			JCudaDriver.cuModuleLoadData(patternModule, ptx[0]);
			patternFinalizeFunction = new CUfunction();
			JCudaDriver.cuModuleGetFunction(patternFinalizeFunction, patternModule,
				"finalize_pattern_statistics");
		} finally {
			JNvrtc.nvrtcDestroyProgram(program);
		}
	}

	private void launchPatternFinalizer(int rank, int productWidth, int variants,
		int patterns, int firstPattern, boolean meanFill) {
		int threads = 128;
		Pointer parameters = Pointer.to(
			Pointer.to(outputBuffer), Pointer.to(patternUpperBuffer),
			Pointer.to(patternSumsBuffer), Pointer.to(patternObservedBuffer),
			Pointer.to(patternCompactBuffer), Pointer.to(new int[] {rank}),
			Pointer.to(new int[] {productWidth}), Pointer.to(new int[] {variants}),
			Pointer.to(new int[] {firstPattern}), Pointer.to(new int[] {meanFill ? 1 : 0}));
		JCudaDriver.cuLaunchKernel(patternFinalizeFunction,
			(variants + threads - 1) / threads, patterns, 1,
			threads, 1, 1, 0, null, parameters, null);
	}

	private Pointer ensurePatternBuffer(Pointer current, long requiredBytes, long capacity) {
		return current != null && capacity >= requiredBytes
			? current : replaceBuffer(current, requiredBytes);
	}

	private static String patternKernelSource() {
		return "extern \"C\" __global__ void finalize_pattern_statistics("
			+ "const double* products, const double* upper, const double* design_sums, "
			+ "const int* observed, double* compact, int rank, int product_width, "
			+ "int variants, int first_pattern, int mean_fill) {\n"
			+ "  int variant = blockIdx.x * blockDim.x + threadIdx.x;\n"
			+ "  int pattern = blockIdx.y;\n"
			+ "  if (variant >= variants) return;\n"
			+ "  int rows_per_pattern = rank + 1;\n"
			+ "  int mask = (pattern * rows_per_pattern) * product_width + 3 * variant;\n"
			+ "  double called = products[mask];\n"
			+ "  double replacement = mean_fill ? (called > 0.0 ? products[mask + 1] / called : (0.0 / 0.0)) : 0.0;\n"
			+ "  double missing = (double)observed[pattern] - called;\n"
			+ "  double sum_squares = products[mask + 2] + missing * replacement * replacement;\n"
			+ "  double projection = 0.0; double q[64];\n"
			+ "  for (int row = 0; row < rank; row++) {\n"
			+ "    int product = (pattern * rows_per_pattern + 1 + row) * product_width + 3 * variant;\n"
			+ "    double value = products[product + 1] + replacement * (design_sums[pattern * rank + row] - products[product]);\n"
			+ "    for (int previous = 0; previous < row; previous++) value -= upper[(pattern * rank + previous) * rank + row] * q[previous];\n"
			+ "    value /= upper[(pattern * rank + row) * rank + row];\n"
			+ "    q[row] = value; projection += value * value;\n"
			+ "  }\n"
			+ "  int target = ((first_pattern + pattern) * variants + variant) * 3;\n"
			+ "  compact[target] = replacement; compact[target + 1] = sum_squares - projection; compact[target + 2] = sum_squares;\n"
			+ "}\n";
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
			if (patternModule != null) JCudaDriver.cuModuleUnload(patternModule);
			if (patternUpperBuffer != null) JCuda.cudaFree(patternUpperBuffer);
			if (patternSumsBuffer != null) JCuda.cudaFree(patternSumsBuffer);
			if (patternObservedBuffer != null) JCuda.cudaFree(patternObservedBuffer);
			if (patternCompactBuffer != null) JCuda.cudaFree(patternCompactBuffer);
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
			patternUpperBuffer = patternSumsBuffer = patternObservedBuffer = patternCompactBuffer = null;
			patternFinalizeFunction = null;
			patternModule = null;
			residualValuesBuffer = residualQBuffer = residualCoefficientsBuffer = null;
			handle = null;
			closed = true;
		}
	}
}
