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

import org.jocl.cl_device_id;
import org.jocl.cl_platform_id;

import static org.jocl.CL.CL_DEVICE_AVAILABLE;
import static org.jocl.CL.CL_DEVICE_COMPILER_AVAILABLE;
import static org.jocl.CL.CL_DEVICE_DOUBLE_FP_CONFIG;
import static org.jocl.CL.CL_DRIVER_VERSION;
import static org.jocl.CL.CL_DEVICE_EXTENSIONS;
import static org.jocl.CL.CL_DEVICE_HOST_UNIFIED_MEMORY;
import static org.jocl.CL.CL_DEVICE_GLOBAL_MEM_SIZE;
import static org.jocl.CL.CL_DEVICE_MAX_MEM_ALLOC_SIZE;
import static org.jocl.CL.CL_DEVICE_MAX_WORK_GROUP_SIZE;
import static org.jocl.CL.CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS;
import static org.jocl.CL.CL_DEVICE_MAX_WORK_ITEM_SIZES;
import static org.jocl.CL.CL_DEVICE_NAME;
import static org.jocl.CL.CL_DEVICE_VENDOR;
import static org.jocl.CL.CL_DEVICE_VERSION;
import static org.jocl.CL.CL_PLATFORM_NAME;

final class JoclGpuDevice implements GpuDevice {
	private final cl_platform_id platformId;
	private final cl_device_id deviceId;
	private final String platformName;
	private final String name;
	private final String vendor;
	private final String driverVersion;
	private final String computeApiVersion;
	private final String extensions;
	private final boolean available;
	private final boolean compilerAvailable;
	private final boolean doublePrecision;
	private final boolean unifiedMemory;
	private final long globalMemoryBytes;
	private final long maxAllocationBytes;
	private final long maxWorkGroupSize;
	private final long[] maxWorkItemSizes;

	JoclGpuDevice(cl_platform_id platformId, cl_device_id deviceId) {
		this.platformId = platformId;
		this.deviceId = deviceId;
		platformName = JoclGpuBackend.getPlatformString(platformId, CL_PLATFORM_NAME);
		name = JoclGpuBackend.getDeviceString(deviceId, CL_DEVICE_NAME);
		vendor = JoclGpuBackend.getDeviceString(deviceId, CL_DEVICE_VENDOR);
		driverVersion = JoclGpuBackend.getDeviceString(deviceId, CL_DRIVER_VERSION);
		computeApiVersion = JoclGpuBackend.getDeviceString(deviceId, CL_DEVICE_VERSION);
		extensions = JoclGpuBackend.getDeviceString(deviceId, CL_DEVICE_EXTENSIONS);
		available = JoclGpuBackend.getDeviceInt(deviceId, CL_DEVICE_AVAILABLE) != 0;
		compilerAvailable = JoclGpuBackend.getDeviceInt(deviceId, CL_DEVICE_COMPILER_AVAILABLE) != 0;
		long fpConfig;
		try {
			fpConfig = JoclGpuBackend.getDeviceLong(deviceId, CL_DEVICE_DOUBLE_FP_CONFIG);
		} catch (RuntimeException ignored) {
			fpConfig = 0;
		}
		doublePrecision = fpConfig != 0 || hasExtension("cl_khr_fp64") || hasExtension("cl_amd_fp64");
		boolean hostUnified;
		try {
			hostUnified = JoclGpuBackend.getDeviceInt(deviceId, CL_DEVICE_HOST_UNIFIED_MEMORY) != 0;
		} catch (RuntimeException ignored) {
			hostUnified = false;
		}
		unifiedMemory = hostUnified;
		globalMemoryBytes = JoclGpuBackend.getDeviceLong(deviceId, CL_DEVICE_GLOBAL_MEM_SIZE);
		maxAllocationBytes = JoclGpuBackend.getDeviceLong(deviceId, CL_DEVICE_MAX_MEM_ALLOC_SIZE);
		maxWorkGroupSize = JoclGpuBackend.getDeviceSizeT(deviceId, CL_DEVICE_MAX_WORK_GROUP_SIZE);
		int dimensions = JoclGpuBackend.getDeviceInt(deviceId, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);
		maxWorkItemSizes = JoclGpuBackend.getDeviceSizeTArray(deviceId, CL_DEVICE_MAX_WORK_ITEM_SIZES, dimensions);
	}

	cl_platform_id platformId() {
		return platformId;
	}

	cl_device_id deviceId() {
		return deviceId;
	}

	private boolean hasExtension(String extension) {
		for (String candidate : extensions.split("\\s+")) {
			if (extension.equals(candidate)) {
				return true;
			}
		}
		return false;
	}

	@Override
	public String getBackendName() { return "opencl"; }

	@Override
	public String getPlatformName() { return platformName; }

	@Override
	public String getName() { return name; }

	@Override
	public String getVendor() { return vendor; }

	@Override
	public String getDriverVersion() { return driverVersion; }

	@Override
	public String getComputeApiVersion() { return computeApiVersion; }

	@Override
	public boolean isAvailable() { return available; }

	@Override
	public boolean isCompilerAvailable() { return compilerAvailable; }

	@Override
	public boolean supportsDoublePrecision() { return doublePrecision; }

	@Override
	public boolean hasUnifiedMemory() { return unifiedMemory; }

	@Override
	public long getGlobalMemoryBytes() { return globalMemoryBytes; }

	@Override
	public long getMaxAllocationBytes() { return maxAllocationBytes; }

	@Override
	public long getMaxWorkGroupSize() { return maxWorkGroupSize; }

	@Override
	public long[] getMaxWorkItemSizes() { return maxWorkItemSizes.clone(); }

	@Override
	public GpuContext openContext() { return new JoclGpuContext(this); }
}
