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
import gov.nih.gpu.GpuException;

import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_program;

import java.nio.charset.StandardCharsets;

import static org.jocl.CL.CL_MEM_READ_ONLY;
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
	public synchronized void compileKernel(String source, String kernelName) {
		ensureOpen();
		releaseProgramAndKernel();
		int[] status = new int[1];
		program = clCreateProgramWithSource(context, 1, new String[] { source }, null, status);
		JoclGpuBackend.check(status[0], "clCreateProgramWithSource");
		status[0] = clBuildProgram(program, 1, new org.jocl.cl_device_id[] { device.deviceId() }, null, null, null);
		if (status[0] != org.jocl.CL.CL_SUCCESS) {
			String buildLog = getBuildLog();
			int buildStatus = status[0];
			releaseProgramAndKernel();
			throw new GpuException("OpenCL kernel build failed on " + device.getName() + ": "
				+ org.jocl.CL.stringFor_errorCode(buildStatus) + " (" + buildStatus + ")\n" + buildLog);
		}
		kernel = clCreateKernel(program, kernelName, status);
		JoclGpuBackend.check(status[0], "clCreateKernel");
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
		if (inputA == null || inputB == null || outputElements <= 0) {
			throw new IllegalArgumentException("Kernel buffers must be non-null and non-empty");
		}
		if (globalWorkSize == null || localWorkSize == null || globalWorkSize.length != 2 || localWorkSize.length != 2) {
			throw new IllegalArgumentException("The eQTL kernel requires two-dimensional work sizes");
		}

		long inputABytes = Math.multiplyExact((long) inputA.length, Sizeof.cl_double);
		long inputBBytes = Math.multiplyExact((long) inputB.length, Sizeof.cl_double);
		long outputBytes = Math.multiplyExact((long) outputElements, Sizeof.cl_double);
		inputABuffer = ensureBuffer(inputABuffer, inputABytes, inputACapacity, CL_MEM_READ_ONLY);
		inputACapacity = Math.max(inputACapacity, inputABytes);
		inputBBuffer = ensureBuffer(inputBBuffer, inputBBytes, inputBCapacity, CL_MEM_READ_ONLY);
		inputBCapacity = Math.max(inputBCapacity, inputBBytes);
		outputBuffer = ensureBuffer(outputBuffer, outputBytes, outputCapacity, CL_MEM_WRITE_ONLY);
		outputCapacity = Math.max(outputCapacity, outputBytes);

		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, inputABuffer, CL_TRUE, 0, inputABytes,
			Pointer.to(inputA), 0, null, null), "clEnqueueWriteBuffer(inputA)");
		JoclGpuBackend.check(clEnqueueWriteBuffer(queue, inputBBuffer, CL_TRUE, 0, inputBBytes,
			Pointer.to(inputB), 0, null, null), "clEnqueueWriteBuffer(inputB)");

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

		double[] output = new double[outputElements];
		JoclGpuBackend.check(clEnqueueReadBuffer(queue, outputBuffer, CL_TRUE, 0, outputBytes,
			Pointer.to(output), 0, null, null), "clEnqueueReadBuffer");
		return output;
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

	private String getBuildLog() {
		long[] size = new long[1];
		JoclGpuBackend.check(clGetProgramBuildInfo(program, device.deviceId(), CL_PROGRAM_BUILD_LOG,
			0, null, size), "clGetProgramBuildInfo");
		byte[] log = new byte[(int) size[0]];
		JoclGpuBackend.check(clGetProgramBuildInfo(program, device.deviceId(), CL_PROGRAM_BUILD_LOG,
			log.length, Pointer.to(log), null), "clGetProgramBuildInfo");
		int length = log.length;
		while (length > 0 && log[length - 1] == 0) {
			length--;
		}
		return new String(log, 0, length, StandardCharsets.UTF_8);
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
