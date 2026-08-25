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

/** Metadata and context factory for the host CPU. */
final class CpuDevice implements GpuDevice {
	private final CpuBackend backend;
	private final String name;
	private final long heapBytes;
	private volatile String engineDescription = "native BLAS auto-detection; Java fallback available";

	CpuDevice(CpuBackend backend) {
		this.backend = backend;
		name = System.getProperty("os.arch", "unknown architecture") + " host CPU ("
			+ Runtime.getRuntime().availableProcessors() + " logical processors)";
		heapBytes = Runtime.getRuntime().maxMemory();
	}

	void recordEngine(CpuMatrixEngine engine) {
		engineDescription = engine.description();
	}

	@Override
	public String getBackendName() { return "cpu"; }

	@Override
	public String getPlatformName() { return "Java host process"; }

	@Override
	public String getName() { return name; }

	@Override
	public String getVendor() { return System.getProperty("java.vendor", "unknown"); }

	@Override
	public String getDriverVersion() {
		return System.getProperty("os.name", "unknown") + " " + System.getProperty("os.version", "");
	}

	@Override
	public String getComputeApiVersion() { return engineDescription; }

	@Override
	public boolean isAvailable() { return true; }

	@Override
	public boolean isCompilerAvailable() { return true; }

	@Override
	public boolean supportsDoublePrecision() { return true; }

	@Override
	public boolean hasUnifiedMemory() { return true; }

	@Override
	public long getGlobalMemoryBytes() { return heapBytes; }

	@Override
	public long getMaxAllocationBytes() { return Math.max(1, heapBytes / 2); }

	@Override
	public long getMaxWorkGroupSize() { return 1; }

	@Override
	public long[] getMaxWorkItemSizes() { return new long[] { 1, 1, 1 }; }

	@Override
	public GpuContext openContext() {
		CpuMatrixEngine engine = backend.openEngine();
		recordEngine(engine);
		return new CpuContext(this, engine);
	}
}
