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

import jcuda.runtime.cudaDeviceProp;

import static jcuda.runtime.cudaComputeMode.cudaComputeModeProhibited;

/** Metadata and context factory for one NVIDIA CUDA device. */
final class CudaGpuDevice implements GpuDevice {
	private final int ordinal;
	private final String name;
	private final String driverVersion;
	private final String computeApiVersion;
	private final boolean available;
	private final boolean doublePrecision;
	private final boolean unifiedMemory;
	private final long totalGlobalMemory;
	private final long maxWorkGroupSize;
	private final long[] maxWorkItemSizes;

	CudaGpuDevice(int ordinal, cudaDeviceProp properties, String driverApiVersion, String runtimeVersion) {
		this.ordinal = ordinal;
		name = properties.getName();
		driverVersion = "CUDA driver API " + driverApiVersion;
		computeApiVersion = "CUDA Runtime " + runtimeVersion + " (compute capability "
			+ properties.major + "." + properties.minor + ")";
		available = properties.computeMode != cudaComputeModeProhibited;
		doublePrecision = properties.major > 1 || (properties.major == 1 && properties.minor >= 3);
		unifiedMemory = properties.integrated != 0;
		totalGlobalMemory = properties.totalGlobalMem;
		maxWorkGroupSize = properties.maxThreadsPerBlock;
		maxWorkItemSizes = new long[properties.maxThreadsDim.length];
		for (int i = 0; i < properties.maxThreadsDim.length; i++) {
			maxWorkItemSizes[i] = properties.maxThreadsDim[i];
		}
	}

	int ordinal() {
		return ordinal;
	}

	@Override
	public String getBackendName() { return "cuda"; }

	@Override
	public String getPlatformName() { return "NVIDIA CUDA"; }

	@Override
	public String getName() { return name; }

	@Override
	public String getVendor() { return "NVIDIA Corporation"; }

	@Override
	public String getDriverVersion() { return driverVersion; }

	@Override
	public String getComputeApiVersion() { return computeApiVersion; }

	@Override
	public boolean isAvailable() { return available; }

	@Override
	public boolean isCompilerAvailable() { return true; }

	@Override
	public boolean supportsDoublePrecision() { return doublePrecision; }

	@Override
	public boolean hasUnifiedMemory() { return unifiedMemory; }

	@Override
	public long getMaxAllocationBytes() { return totalGlobalMemory; }

	@Override
	public long getMaxWorkGroupSize() { return maxWorkGroupSize; }

	@Override
	public long[] getMaxWorkItemSizes() { return maxWorkItemSizes.clone(); }

	@Override
	public GpuContext openContext() {
		if (!available) {
			throw new GpuException("CUDA device " + name + " is in prohibited compute mode");
		}
		return new CudaGpuContext(this);
	}
}