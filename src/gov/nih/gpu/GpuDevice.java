/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

/** Vendor-neutral metadata for one compute device. */
public interface GpuDevice {
	String getBackendName();

	String getPlatformName();

	String getName();

	String getVendor();

	String getDriverVersion();

	String getComputeApiVersion();

	boolean isAvailable();

	boolean isCompilerAvailable();

	boolean supportsDoublePrecision();

	boolean hasUnifiedMemory();

	long getGlobalMemoryBytes();

	long getMaxAllocationBytes();

	long getMaxWorkGroupSize();

	long[] getMaxWorkItemSizes();

	GpuContext openContext();
}
