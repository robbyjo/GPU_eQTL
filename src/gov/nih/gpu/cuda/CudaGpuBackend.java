/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.cuda;

import gov.nih.gpu.GpuBackend;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuException;

import jcuda.jcublas.JCublas2;
import jcuda.driver.JCudaDriver;
import jcuda.nvrtc.JNvrtc;
import jcuda.runtime.JCuda;
import jcuda.runtime.cudaDeviceProp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** NVIDIA CUDA Runtime and cuBLAS implementation of the GPU contracts. */
public final class CudaGpuBackend implements GpuBackend {
	private static final String BINDING_VERSION = "12.6.0";

	static {
		JCuda.setExceptionsEnabled(true);
		JCublas2.setExceptionsEnabled(true);
		JCudaDriver.setExceptionsEnabled(true);
		JNvrtc.setExceptionsEnabled(true);
	}

	@Override
	public String getName() {
		return "cuda";
	}

	@Override
	public String getRuntimeDescription() {
		return "JCuda " + BINDING_VERSION + " using the system CUDA Runtime and cuBLAS";
	}

	@Override
	public List<GpuDevice> discoverGpuDevices() {
		try {
			JCuda.initialize();
			JCublas2.initialize();
			int[] count = new int[1];
			JCuda.cudaGetDeviceCount(count);
			if (count[0] == 0) {
				return Collections.emptyList();
			}

			int[] driverVersion = new int[1];
			int[] runtimeVersion = new int[1];
			JCuda.cudaDriverGetVersion(driverVersion);
			JCuda.cudaRuntimeGetVersion(runtimeVersion);
			List<GpuDevice> devices = new ArrayList<GpuDevice>(count[0]);
			for (int ordinal = 0; ordinal < count[0]; ordinal++) {
				cudaDeviceProp properties = new cudaDeviceProp();
				JCuda.cudaGetDeviceProperties(properties, ordinal);
				devices.add(new CudaGpuDevice(ordinal, properties,
					formatCudaVersion(driverVersion[0]), formatCudaVersion(runtimeVersion[0])));
			}
			return Collections.unmodifiableList(devices);
		} catch (RuntimeException | LinkageError e) {
			throw new GpuException("The NVIDIA CUDA Runtime or cuBLAS could not be initialized. "
				+ "Install a compatible NVIDIA driver and CUDA runtime, or select the OpenCL backend.", e);
		}
	}

	private static String formatCudaVersion(int encoded) {
		if (encoded <= 0) {
			return "unknown";
		}
		int major = encoded / 1000;
		int minor = (encoded % 1000) / 10;
		int patch = encoded % 10;
		return patch == 0 ? major + "." + minor : major + "." + minor + "." + patch;
	}
}
