/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/** Entry point for backend selection, device filtering, and diagnostics. */
public final class GpuRuntime {
	public static final int DEFAULT_ALIGNMENT = 64;
	private static final String BACKEND_PROPERTY = "eqtl.gpu.backend";

	private final GpuBackend backend;

	public GpuRuntime(GpuBackend backend) {
		if (backend == null) {
			throw new IllegalArgumentException("backend must not be null");
		}
		this.backend = backend;
	}

	public static GpuRuntime createDefault() {
		String requested = System.getProperty(BACKEND_PROPERTY, "auto").trim();
		if (requested.isEmpty() || "auto".equalsIgnoreCase(requested)) {
			return new GpuRuntime(new AutoGpuBackend());
		}
		if ("cuda".equalsIgnoreCase(requested)) {
			return new GpuRuntime(new CudaGpuBackend());
		}
		if ("opencl".equalsIgnoreCase(requested) || "jocl".equalsIgnoreCase(requested)) {
			return new GpuRuntime(new JoclGpuBackend());
		}
		throw new GpuException("Unsupported GPU backend '" + requested
			+ "'. Available backends: auto, cuda, opencl (JOCL)");
	}

	public GpuBackend getBackend() {
		return backend;
	}

	public List<GpuDevice> getGpuDevices(boolean onlyAvailable, boolean requireDoublePrecision) {
		List<GpuDevice> result = new ArrayList<GpuDevice>();
		for (GpuDevice device : backend.discoverGpuDevices()) {
			if (onlyAvailable && (!device.isAvailable() || !device.isCompilerAvailable())) {
				continue;
			}
			if (requireDoublePrecision && !device.supportsDoublePrecision()) {
				continue;
			}
			result.add(device);
		}
		return Collections.unmodifiableList(result);
	}

	public GpuContext[] openGpuContexts(boolean onlyAvailable, boolean requireDoublePrecision) {
		List<GpuDevice> devices = getGpuDevices(onlyAvailable, requireDoublePrecision);
		List<GpuContext> contexts = new ArrayList<GpuContext>(devices.size());
		try {
			for (GpuDevice device : devices) {
				contexts.add(device.openContext());
			}
			return contexts.toArray(new GpuContext[contexts.size()]);
		} catch (RuntimeException e) {
			for (GpuContext context : contexts) {
				context.close();
			}
			throw e;
		}
	}

	public String describeGpuDevices() {
		StringBuilder text = new StringBuilder();
		text.append("Backend: ").append(backend.getName()).append('\n');
		text.append("Runtime: ").append(backend.getRuntimeDescription()).append('\n');
		List<GpuDevice> devices = backend.discoverGpuDevices();
		if (devices.isEmpty()) {
			text.append("No GPU devices were reported.\n");
			return text.toString();
		}
		for (int i = 0; i < devices.size(); i++) {
			GpuDevice device = devices.get(i);
			text.append('\n').append("Device #").append(i + 1).append(": ").append(device.getName()).append('\n');
			text.append("  Platform: ").append(device.getPlatformName()).append('\n');
			text.append("  Vendor: ").append(device.getVendor()).append('\n');
			text.append("  Driver: ").append(device.getDriverVersion()).append('\n');
			text.append("  Compute API: ").append(device.getComputeApiVersion()).append('\n');
			text.append("  Available: ").append(device.isAvailable()).append('\n');
			text.append("  Compiler available: ").append(device.isCompilerAvailable()).append('\n');
			text.append("  Double precision: ").append(device.supportsDoublePrecision()).append('\n');
			text.append("  Unified host memory: ").append(device.hasUnifiedMemory()).append('\n');
			text.append("  Maximum allocation (bytes): ").append(device.getMaxAllocationBytes()).append('\n');
			text.append("  Maximum work-group size: ").append(device.getMaxWorkGroupSize()).append('\n');
		}
		return text.toString();
	}
}
