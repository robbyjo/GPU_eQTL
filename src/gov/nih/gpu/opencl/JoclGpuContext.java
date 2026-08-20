/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.opencl;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuExecutionMetrics;
import gov.nih.gpu.GpuException;
import gov.nih.gpu.GpuPrecision;

import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_program;

import java.nio.charset.StandardCharsets;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_READ_WRITE;
import static org.jocl.CL.CL_MEM_WRITE_ONLY;
import static org.jocl.CL.CL_PROGRAM_BUILD_LOG;
import static org.jocl.CL.CL_TRUE;
import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.clBuildProgram;
import static org.jocl.CL.clCreateBuffer;
import static org.jocl.CL.clCreateCommandQueue;
import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clCreateKernel;
import static org.jocl.CL.clCreateProgramWithSource;
import static org.jocl.CL.clEnqueueNDRangeKernel;
import static org.jocl.CL.clEnqueueReadBuffer;
import static org.jocl.CL.clEnqueueWriteBuffer;
import static org.jocl.CL.clFinish;
import static org.jocl.CL.clGetProgramBuildInfo;
import static org.jocl.CL.clReleaseCommandQueue;
import static org.jocl.CL.clReleaseContext;
import static org.jocl.CL.clReleaseKernel;
import static org.jocl.CL.clReleaseMemObject;
import static org.jocl.CL.clReleaseProgram;
import static org.jocl.CL.clSetKernelArg;

/** JOCL implementation with reusable device buffers for repeated matrix blocks. */
final class JoclGpuContext implements GpuContext {
	private static final Pattern DATA_TYPE = Pattern.compile(
		"(?m)^\\s*#define\\s+DATATYPE\\s+(float|double)\\s*$");
	private final JoclGpuDevice device;
	private cl_context context;
	private cl_command_queue queue;
	private cl_program program;
	private cl_kernel kernel;
	private cl_mem inputABuffer;
	private cl_mem inputBBuffer;
	private cl_mem outputBuffer;
	private long inputACapacity;
	private long inputBCapacity;
	private long outputCapacity;
	private cl_program residualProgram;
	private cl_kernel residualCoefficientKernel;
	private cl_kernel residualApplyKernel;
	private cl_mem residualValuesBuffer;
	private cl_mem residualQBuffer;
	private cl_mem residualCoefficientsBuffer;
	private long residualValuesCapacity;
	private long residualQCapacity;
	private long residualCoefficientsCapacity;
	private Object residualProjectionSource;
	private GpuPrecision residualPrecision;
	private int residualProjectionSamples;
	private int residualProjectionRank;
	private GpuPrecision precision = GpuPrecision.FP64;
	private boolean profilingEnabled;
	private GpuExecutionMetrics lastExecutionMetrics = GpuExecutionMetrics.EMPTY;
	private boolean closed;

	JoclGpuContext(JoclGpuDevice device) {
		this.device = device;
		cl_context_properties properties = new cl_context_properties();
		properties.addProperty(CL_CONTEXT_PLATFORM, device.platformId());
		int[] status = new int[1];
		context = clCreateContext(properties, 1, new org.jocl.cl_device_id[] { device.deviceId() }, null, null, status);
		JoclGpuBackend.check(status[0], "clCreateContext");
		queue = clCreateCommandQueue(context, device.deviceId(), 0, status);
		if (status[0] != org.jocl.CL.CL_SUCCESS) {
			clReleaseContext(context);
			context = null;
			JoclGpuBackend.check(status[0], "clCreateCommandQueue");
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
		if (source == null)
			throw new IllegalArgumentException("Kernel source must not be null");
		Matcher dataType = DATA_TYPE.matcher(source);
		String expectedType = precision == GpuPrecision.FP32 ? "float" : "double";
		if (!dataType.find() || !expectedType.equals(dataType.group(1)))
			throw new GpuException("Kernel DATATYPE does not match " + precision.optionName());
		releaseProgramAndKernel();
		int[] status = new int[1];
		program = clCreateProgramWithSource(context, 1, new String[] { source }, null, status);
		JoclGpuBackend.check(status[0], "clCreateProgramWithSource");
		status[0] = clBuildProgram(program, 1, new org.jocl.cl_device_id[] { device.deviceId() }, null, null, null);
		if (status[0] != org.jocl.CL.CL_SUCCESS) {
			String buildLog = getBuildLog(program);
			int buildStatus = status[0];
			releaseProgramAndKernel();
			throw new GpuException("OpenCL kernel build failed on " + device.getName() + ": "
				+ org.jocl.CL.stringFor_errorCode(buildStatus) + " (" + buildStatus + ")\n" + buildLog);
		}
		kernel = clCreateKernel(program, kernelName, status);
		JoclGpuBackend.check(status[0], "clCreateKernel");
		this.precision = precision;
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
		if (kernel == null) {
			throw new GpuException("No GPU kernel has been compiled for " + device.getName());
		}
		if (precision != GpuPrecision.FP64)
			throw new GpuException("The OpenCL kernel was compiled for " + precision.optionName());
		if (inputA == null || inputB == null || outputElements <= 0) {
			throw new IllegalArgumentException("Kernel buffers must be non-null and non-empty");
		}
		if (globalWorkSize == null || localWorkSize == null || globalWorkSize.length != 2 || localWorkSize.length != 2) {
			throw new IllegalArgumentException("The eQTL kernel requires two-dimensional work sizes");
		}

		long inputABytes = Math.multiplyExact((long) inputA.length, Sizeof.cl_double);
		long inputBBytes = Math.multiplyExact((long) inputB.length, Sizeof.cl_double);
		long outputBytes = Math.multiplyExact((long) outputElements, Sizeof.cl_double);
		long phaseStart = profilingEnabled ? System.nanoTime() : 0;
		inputABuffer = ensureBuffer(inputABuffer, inputABytes, inputACapacity, CL_MEM_READ_ONLY);
		inputACapacity = Math.max(inputACapacity, inputABytes);
		inputBBuffer = ensureBuffer(inputBBuffer, inputBBytes, inputBCapacity, CL_MEM_READ_ONLY);
		inputBCapacity = Math.max(inputBCapacity, inputBBytes);
		outputBuffer = ensureBuffer(outputBuffer, outputBytes, outputCapacity, CL_MEM_WRITE_ONLY);
		outputCapacity = Math.max(outputCapacity, outputBytes);
		double[] output = new double[outputElements];
		long setupNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, inputABuffer, CL_TRUE, 0, inputABytes,
			Pointer.to(inputA), 0, null, null), "clEnqueueWriteBuffer(inputA)");
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, inputBBuffer, CL_TRUE, 0, inputBBytes,
			Pointer.to(inputB), 0, null, null), "clEnqueueWriteBuffer(inputB)");
		long uploadNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(outputBuffer)), "clSetKernelArg(0)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 1, Sizeof.cl_mem, Pointer.to(inputABuffer)), "clSetKernelArg(1)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(inputBBuffer)), "clSetKernelArg(2)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 3, localMemoryBytes, null), "clSetKernelArg(3)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 4, localMemoryBytes, null), "clSetKernelArg(4)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 5, Sizeof.cl_int, Pointer.to(new int[] { widthA })), "clSetKernelArg(5)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 6, Sizeof.cl_int, Pointer.to(new int[] { widthB })), "clSetKernelArg(6)");
		JoclGpuBackend.check(clEnqueueNDRangeKernel(queue, kernel, 2, null, globalWorkSize, localWorkSize,
			0, null, null), "clEnqueueNDRangeKernel");
		JoclGpuBackend.check(clFinish(queue), "clFinish");
		long computeNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueReadBuffer(queue, outputBuffer, CL_TRUE, 0, outputBytes,
			Pointer.to(output), 0, null, null), "clEnqueueReadBuffer");
		long downloadNanos = elapsed(phaseStart);
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
				inputABytes + inputBBytes, outputBytes)
			: GpuExecutionMetrics.EMPTY;
		return output;
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
		if (kernel == null)
			throw new GpuException("No GPU kernel has been compiled for " + device.getName());
		if (precision != GpuPrecision.FP32)
			throw new GpuException("The OpenCL kernel was compiled for " + precision.optionName());
		if (inputA == null || inputB == null || outputElements <= 0)
			throw new IllegalArgumentException("Kernel buffers must be non-null and non-empty");
		if (globalWorkSize == null || localWorkSize == null
				|| globalWorkSize.length != 2 || localWorkSize.length != 2)
			throw new IllegalArgumentException("The eQTL kernel requires two-dimensional work sizes");

		long inputABytes = Math.multiplyExact((long) inputA.length, Sizeof.cl_float);
		long inputBBytes = Math.multiplyExact((long) inputB.length, Sizeof.cl_float);
		long outputBytes = Math.multiplyExact((long) outputElements, Sizeof.cl_float);
		long phaseStart = profilingEnabled ? System.nanoTime() : 0;
		inputABuffer = ensureBuffer(inputABuffer, inputABytes, inputACapacity, CL_MEM_READ_ONLY);
		inputACapacity = Math.max(inputACapacity, inputABytes);
		inputBBuffer = ensureBuffer(inputBBuffer, inputBBytes, inputBCapacity, CL_MEM_READ_ONLY);
		inputBCapacity = Math.max(inputBCapacity, inputBBytes);
		outputBuffer = ensureBuffer(outputBuffer, outputBytes, outputCapacity, CL_MEM_WRITE_ONLY);
		outputCapacity = Math.max(outputCapacity, outputBytes);
		float[] output = new float[outputElements];
		long setupNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, inputABuffer, CL_TRUE, 0, inputABytes,
			Pointer.to(inputA), 0, null, null), "clEnqueueWriteBuffer(inputA)");
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, inputBBuffer, CL_TRUE, 0, inputBBytes,
			Pointer.to(inputB), 0, null, null), "clEnqueueWriteBuffer(inputB)");
		long uploadNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(outputBuffer)), "clSetKernelArg(0)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 1, Sizeof.cl_mem, Pointer.to(inputABuffer)), "clSetKernelArg(1)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(inputBBuffer)), "clSetKernelArg(2)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 3, localMemoryBytes, null), "clSetKernelArg(3)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 4, localMemoryBytes, null), "clSetKernelArg(4)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 5, Sizeof.cl_int, Pointer.to(new int[] { widthA })), "clSetKernelArg(5)");
		JoclGpuBackend.check(clSetKernelArg(kernel, 6, Sizeof.cl_int, Pointer.to(new int[] { widthB })), "clSetKernelArg(6)");
		JoclGpuBackend.check(clEnqueueNDRangeKernel(queue, kernel, 2, null, globalWorkSize, localWorkSize,
			0, null, null), "clEnqueueNDRangeKernel");
		JoclGpuBackend.check(clFinish(queue), "clFinish");
		long computeNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueReadBuffer(queue, outputBuffer, CL_TRUE, 0, outputBytes,
			Pointer.to(output), 0, null, null), "clEnqueueReadBuffer");
		long downloadNanos = elapsed(phaseStart);
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
				inputABytes + inputBBytes, outputBytes)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	@Override
	public synchronized double[] residualizeDoubleRows(double[] rowMajorValues, double[] rowMajorQ,
		int rowCount, int sampleCount, int covariateRank) {
		ensureOpen();
		validateResidualizationArguments(rowMajorValues == null ? 0 : rowMajorValues.length,
			rowMajorQ == null ? 0 : rowMajorQ.length, rowCount, sampleCount, covariateRank);
		ensureResidualProgram(GpuPrecision.FP64);
		long valuesBytes = Math.multiplyExact((long) rowMajorValues.length, Sizeof.cl_double);
		long qBytes = Math.multiplyExact((long) rowMajorQ.length, Sizeof.cl_double);
		long coefficientBytes = Math.multiplyExact((long) rowCount * covariateRank, Sizeof.cl_double);
		long phaseStart = profilingEnabled ? System.nanoTime() : 0;
		ensureResidualBuffers(valuesBytes, qBytes, coefficientBytes);
		double[] output = new double[rowMajorValues.length];
		long setupNanos = elapsed(phaseStart);

		boolean uploadProjection = !isCachedProjection(rowMajorQ, GpuPrecision.FP64,
			sampleCount, covariateRank);
		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, residualValuesBuffer, CL_TRUE, 0,
			valuesBytes, Pointer.to(rowMajorValues), 0, null, null), "clEnqueueWriteBuffer(residualValues)");
		if (uploadProjection) {
			JoclGpuBackend.check(clEnqueueWriteBuffer(queue, residualQBuffer, CL_TRUE, 0,
				qBytes, Pointer.to(rowMajorQ), 0, null, null), "clEnqueueWriteBuffer(residualQ)");
			rememberProjection(rowMajorQ, GpuPrecision.FP64, sampleCount, covariateRank);
		}
		long uploadNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		executeResidualKernels(rowCount, sampleCount, covariateRank);
		long computeNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueReadBuffer(queue, residualValuesBuffer, CL_TRUE, 0,
			valuesBytes, Pointer.to(output), 0, null, null), "clEnqueueReadBuffer(residualValues)");
		long downloadNanos = elapsed(phaseStart);
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
				valuesBytes + (uploadProjection ? qBytes : 0), valuesBytes)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	@Override
	public synchronized float[] residualizeFloatRows(float[] rowMajorValues, float[] rowMajorQ,
		int rowCount, int sampleCount, int covariateRank) {
		ensureOpen();
		validateResidualizationArguments(rowMajorValues == null ? 0 : rowMajorValues.length,
			rowMajorQ == null ? 0 : rowMajorQ.length, rowCount, sampleCount, covariateRank);
		ensureResidualProgram(GpuPrecision.FP32);
		long valuesBytes = Math.multiplyExact((long) rowMajorValues.length, Sizeof.cl_float);
		long qBytes = Math.multiplyExact((long) rowMajorQ.length, Sizeof.cl_float);
		long coefficientBytes = Math.multiplyExact((long) rowCount * covariateRank, Sizeof.cl_float);
		long phaseStart = profilingEnabled ? System.nanoTime() : 0;
		ensureResidualBuffers(valuesBytes, qBytes, coefficientBytes);
		float[] output = new float[rowMajorValues.length];
		long setupNanos = elapsed(phaseStart);

		boolean uploadProjection = !isCachedProjection(rowMajorQ, GpuPrecision.FP32,
			sampleCount, covariateRank);
		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, residualValuesBuffer, CL_TRUE, 0,
			valuesBytes, Pointer.to(rowMajorValues), 0, null, null), "clEnqueueWriteBuffer(residualValues)");
		if (uploadProjection) {
			JoclGpuBackend.check(clEnqueueWriteBuffer(queue, residualQBuffer, CL_TRUE, 0,
				qBytes, Pointer.to(rowMajorQ), 0, null, null), "clEnqueueWriteBuffer(residualQ)");
			rememberProjection(rowMajorQ, GpuPrecision.FP32, sampleCount, covariateRank);
		}
		long uploadNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		executeResidualKernels(rowCount, sampleCount, covariateRank);
		long computeNanos = elapsed(phaseStart);

		phaseStart = profilingEnabled ? System.nanoTime() : 0;
		JoclGpuBackend.check(clEnqueueReadBuffer(queue, residualValuesBuffer, CL_TRUE, 0,
			valuesBytes, Pointer.to(output), 0, null, null), "clEnqueueReadBuffer(residualValues)");
		long downloadNanos = elapsed(phaseStart);
		lastExecutionMetrics = profilingEnabled
			? new GpuExecutionMetrics(setupNanos, uploadNanos, computeNanos, downloadNanos,
				valuesBytes + (uploadProjection ? qBytes : 0), valuesBytes)
			: GpuExecutionMetrics.EMPTY;
		return output;
	}

	private long elapsed(long startedAtNanos) {
		return profilingEnabled ? System.nanoTime() - startedAtNanos : 0;
	}

	private cl_mem ensureBuffer(cl_mem current, long requiredBytes, long currentCapacity, long flags) {
		if (current != null && currentCapacity >= requiredBytes) {
			return current;
		}
		if (current != null) {
			clReleaseMemObject(current);
		}
		int[] status = new int[1];
		cl_mem replacement = clCreateBuffer(context, flags, requiredBytes, null, status);
		JoclGpuBackend.check(status[0], "clCreateBuffer");
		return replacement;
	}

	private String getBuildLog(cl_program targetProgram) {
		long[] size = new long[1];
		JoclGpuBackend.check(clGetProgramBuildInfo(targetProgram, device.deviceId(), CL_PROGRAM_BUILD_LOG,
			0, null, size), "clGetProgramBuildInfo");
		byte[] log = new byte[(int) size[0]];
		JoclGpuBackend.check(clGetProgramBuildInfo(targetProgram, device.deviceId(), CL_PROGRAM_BUILD_LOG,
			log.length, Pointer.to(log), null), "clGetProgramBuildInfo");
		int length = log.length;
		while (length > 0 && log[length - 1] == 0) {
			length--;
		}
		return new String(log, 0, length, StandardCharsets.UTF_8);
	}

	private static void validateResidualizationArguments(int valuesLength, int qLength,
		int rowCount, int sampleCount, int covariateRank) {
		if (rowCount <= 0 || sampleCount <= 0 || covariateRank <= 0 || covariateRank > sampleCount)
			throw new IllegalArgumentException("Invalid residualization matrix dimensions");
		if (valuesLength != Math.multiplyExact(rowCount, sampleCount)
			|| qLength != Math.multiplyExact(sampleCount, covariateRank))
			throw new IllegalArgumentException("Residualization buffers do not match their matrix dimensions");
	}

	private void ensureResidualProgram(GpuPrecision requestedPrecision) {
		if (residualProgram != null && residualPrecision == requestedPrecision)
			return;
		releaseResidualProgram();
		String type = requestedPrecision == GpuPrecision.FP32 ? "float" : "double";
		String extension = requestedPrecision == GpuPrecision.FP64
			? "#if defined(cl_khr_fp64)\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
				+ "#elif defined(cl_amd_fp64)\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n"
			: "";
		String source = extension + "#define DATATYPE " + type + "\n"
			+ "__kernel void projection_coefficients(__global const DATATYPE* values, "
			+ "__global const DATATYPE* q, __global DATATYPE* coefficients, "
			+ "const int samples, const int rank) {\n"
			+ "  int column = get_global_id(0), row = get_global_id(1);\n"
			+ "  if (column >= rank) return;\n"
			+ "  DATATYPE sum = (DATATYPE)0;\n"
			+ "  for (int sample = 0; sample < samples; sample++) "
			+ "sum += values[row * samples + sample] * q[sample * rank + column];\n"
			+ "  coefficients[row * rank + column] = sum;\n}\n"
			+ "__kernel void apply_projection(__global DATATYPE* values, "
			+ "__global const DATATYPE* q, __global const DATATYPE* coefficients, "
			+ "const int samples, const int rank) {\n"
			+ "  int sample = get_global_id(0), row = get_global_id(1);\n"
			+ "  if (sample >= samples) return;\n"
			+ "  DATATYPE fitted = (DATATYPE)0;\n"
			+ "  for (int column = 0; column < rank; column++) "
			+ "fitted += coefficients[row * rank + column] * q[sample * rank + column];\n"
			+ "  values[row * samples + sample] -= fitted;\n}\n";
		int[] status = new int[1];
		residualProgram = clCreateProgramWithSource(context, 1, new String[] { source }, null, status);
		JoclGpuBackend.check(status[0], "clCreateProgramWithSource(residualization)");
		status[0] = clBuildProgram(residualProgram, 1,
			new org.jocl.cl_device_id[] { device.deviceId() }, null, null, null);
		if (status[0] != org.jocl.CL.CL_SUCCESS) {
			String buildLog = getBuildLog(residualProgram);
			int buildStatus = status[0];
			releaseResidualProgram();
			throw new GpuException("OpenCL residualization kernel build failed on " + device.getName()
				+ ": " + org.jocl.CL.stringFor_errorCode(buildStatus) + " (" + buildStatus + ")\n" + buildLog);
		}
		residualCoefficientKernel = clCreateKernel(residualProgram, "projection_coefficients", status);
		JoclGpuBackend.check(status[0], "clCreateKernel(projection_coefficients)");
		residualApplyKernel = clCreateKernel(residualProgram, "apply_projection", status);
		JoclGpuBackend.check(status[0], "clCreateKernel(apply_projection)");
		residualPrecision = requestedPrecision;
		residualProjectionSource = null;
	}

	private void ensureResidualBuffers(long valuesBytes, long qBytes, long coefficientBytes) {
		residualValuesBuffer = ensureBuffer(residualValuesBuffer, valuesBytes,
			residualValuesCapacity, CL_MEM_READ_WRITE);
		residualValuesCapacity = Math.max(residualValuesCapacity, valuesBytes);
		if (residualQBuffer == null || residualQCapacity < qBytes) {
			residualQBuffer = ensureBuffer(residualQBuffer, qBytes, residualQCapacity, CL_MEM_READ_ONLY);
			residualQCapacity = qBytes;
			residualProjectionSource = null;
		}
		residualCoefficientsBuffer = ensureBuffer(residualCoefficientsBuffer, coefficientBytes,
			residualCoefficientsCapacity, CL_MEM_READ_WRITE);
		residualCoefficientsCapacity = Math.max(residualCoefficientsCapacity, coefficientBytes);
	}

	private void executeResidualKernels(int rowCount, int sampleCount, int covariateRank) {
		JoclGpuBackend.check(clSetKernelArg(residualCoefficientKernel, 0, Sizeof.cl_mem,
			Pointer.to(residualValuesBuffer)), "clSetKernelArg(projection values)");
		JoclGpuBackend.check(clSetKernelArg(residualCoefficientKernel, 1, Sizeof.cl_mem,
			Pointer.to(residualQBuffer)), "clSetKernelArg(projection Q)");
		JoclGpuBackend.check(clSetKernelArg(residualCoefficientKernel, 2, Sizeof.cl_mem,
			Pointer.to(residualCoefficientsBuffer)), "clSetKernelArg(projection coefficients)");
		JoclGpuBackend.check(clSetKernelArg(residualCoefficientKernel, 3, Sizeof.cl_int,
			Pointer.to(new int[] { sampleCount })), "clSetKernelArg(projection samples)");
		JoclGpuBackend.check(clSetKernelArg(residualCoefficientKernel, 4, Sizeof.cl_int,
			Pointer.to(new int[] { covariateRank })), "clSetKernelArg(projection rank)");
		JoclGpuBackend.check(clEnqueueNDRangeKernel(queue, residualCoefficientKernel, 2, null,
			new long[] { covariateRank, rowCount }, null, 0, null, null),
			"clEnqueueNDRangeKernel(projection_coefficients)");

		JoclGpuBackend.check(clSetKernelArg(residualApplyKernel, 0, Sizeof.cl_mem,
			Pointer.to(residualValuesBuffer)), "clSetKernelArg(apply values)");
		JoclGpuBackend.check(clSetKernelArg(residualApplyKernel, 1, Sizeof.cl_mem,
			Pointer.to(residualQBuffer)), "clSetKernelArg(apply Q)");
		JoclGpuBackend.check(clSetKernelArg(residualApplyKernel, 2, Sizeof.cl_mem,
			Pointer.to(residualCoefficientsBuffer)), "clSetKernelArg(apply coefficients)");
		JoclGpuBackend.check(clSetKernelArg(residualApplyKernel, 3, Sizeof.cl_int,
			Pointer.to(new int[] { sampleCount })), "clSetKernelArg(apply samples)");
		JoclGpuBackend.check(clSetKernelArg(residualApplyKernel, 4, Sizeof.cl_int,
			Pointer.to(new int[] { covariateRank })), "clSetKernelArg(apply rank)");
		JoclGpuBackend.check(clEnqueueNDRangeKernel(queue, residualApplyKernel, 2, null,
			new long[] { sampleCount, rowCount }, null, 0, null, null),
			"clEnqueueNDRangeKernel(apply_projection)");
		JoclGpuBackend.check(clFinish(queue), "clFinish(residualization)");
	}

	private boolean isCachedProjection(Object source, GpuPrecision projectionPrecision,
		int sampleCount, int covariateRank) {
		return residualProjectionSource == source && residualPrecision == projectionPrecision
			&& residualProjectionSamples == sampleCount && residualProjectionRank == covariateRank;
	}

	private void rememberProjection(Object source, GpuPrecision projectionPrecision,
		int sampleCount, int covariateRank) {
		residualProjectionSource = source;
		residualPrecision = projectionPrecision;
		residualProjectionSamples = sampleCount;
		residualProjectionRank = covariateRank;
	}

	@Override
	public synchronized void releaseResidualizationResources() {
		ensureOpen();
		releaseResidualBuffers();
		releaseResidualProgram();
	}

	private void releaseResidualBuffers() {
		if (residualValuesBuffer != null) clReleaseMemObject(residualValuesBuffer);
		if (residualQBuffer != null) clReleaseMemObject(residualQBuffer);
		if (residualCoefficientsBuffer != null) clReleaseMemObject(residualCoefficientsBuffer);
		residualValuesBuffer = residualQBuffer = residualCoefficientsBuffer = null;
		residualValuesCapacity = residualQCapacity = residualCoefficientsCapacity = 0;
		residualProjectionSource = null;
		residualProjectionSamples = residualProjectionRank = 0;
	}

	private void releaseResidualProgram() {
		if (residualCoefficientKernel != null) clReleaseKernel(residualCoefficientKernel);
		if (residualApplyKernel != null) clReleaseKernel(residualApplyKernel);
		if (residualProgram != null) clReleaseProgram(residualProgram);
		residualCoefficientKernel = residualApplyKernel = null;
		residualProgram = null;
		residualPrecision = null;
	}

	private void ensureOpen() {
		if (closed) {
			throw new GpuException("GPU context is closed");
		}
	}

	private void releaseProgramAndKernel() {
		if (kernel != null) {
			clReleaseKernel(kernel);
			kernel = null;
		}
		if (program != null) {
			clReleaseProgram(program);
			program = null;
		}
	}

	@Override
	public synchronized void close() {
		if (closed) {
			return;
		}
		closed = true;
		releaseResidualBuffers();
		releaseResidualProgram();
		releaseProgramAndKernel();
		if (inputABuffer != null) clReleaseMemObject(inputABuffer);
		if (inputBBuffer != null) clReleaseMemObject(inputBBuffer);
		if (outputBuffer != null) clReleaseMemObject(outputBuffer);
		inputABuffer = inputBBuffer = outputBuffer = null;
		if (queue != null) clReleaseCommandQueue(queue);
		if (context != null) clReleaseContext(context);
		queue = null;
		context = null;
	}
}
