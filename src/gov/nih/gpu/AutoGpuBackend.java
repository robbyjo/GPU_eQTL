/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.cpu.CpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;

/**
 * Discovers every usable vendor backend while avoiding duplicate NVIDIA
 * devices. CUDA is preferred over OpenCL for an NVIDIA device because the
 * CUDA implementation uses vendor-tuned cuBLAS SGEMM/DGEMM.
 */
public final class AutoGpuBackend implements GpuBackend {
	private final List<GpuBackend> candidates;
	private List<GpuDevice> devices;
	private String runtimeDescription;

	public AutoGpuBackend() {
		this(Arrays.<GpuBackend>asList(new CudaGpuBackend(), new JoclGpuBackend(), new CpuBackend()));
	}

	AutoGpuBackend(List<GpuBackend> candidates) {
		if (candidates == null || candidates.isEmpty()) {
			throw new IllegalArgumentException("At least one GPU backend candidate is required");
		}
		this.candidates = List.copyOf(candidates);
	}

	@Override
	public String getName() {
		return "auto";
	}

	@Override
	public String getRuntimeDescription() {
		discoverGpuDevices();
		return runtimeDescription;
	}

	@Override
	public synchronized List<GpuDevice> discoverGpuDevices() {
		if (devices != null) {
			return devices;
		}

		List<GpuDevice> discovered = new ArrayList<GpuDevice>();
		List<String> details = new ArrayList<String>();
		Set<String> preferredCudaDevices = new HashSet<String>();
		Set<String> usableBackends = new HashSet<String>();
		for (GpuBackend candidate : candidates) {
			List<GpuDevice> candidateDevices;
			try {
				candidateDevices = candidate.discoverGpuDevices();
				details.add(candidate.getName() + ": " + candidate.getRuntimeDescription());
			} catch (RuntimeException | LinkageError e) {
				details.add(candidate.getName() + " unavailable (" + conciseMessage(e) + ")");
				continue;
			}

			for (GpuDevice device : candidateDevices) {
				boolean usable = isUsable(device);
				String key = deviceKey(device);
				if ("opencl".equalsIgnoreCase(candidate.getName()) && preferredCudaDevices.contains(key)) {
					continue;
				}
				discovered.add(device);
				if (usable) {
					usableBackends.add(candidate.getName().toLowerCase(Locale.ROOT));
					if ("cuda".equalsIgnoreCase(candidate.getName())) {
						preferredCudaDevices.add(key);
					}
				}
			}
		}

		String usable = usableBackends.isEmpty()
			? "no compute backend is currently usable"
			: "usable backend(s): " + String.join(", ", usableBackends);
		runtimeDescription = "Automatic GPU-first discovery with CPU fallback; " + usable
			+ ". " + String.join("; ", details);
		devices = Collections.unmodifiableList(discovered);
		return devices;
	}

	private static boolean isUsable(GpuDevice device) {
		return device.isAvailable() && device.isCompilerAvailable();
	}

	private static String deviceKey(GpuDevice device) {
		return normalize(device.getVendor()) + "|" + normalize(device.getName());
	}

	private static String normalize(String value) {
		return value == null ? "" : value.replaceAll("\\s+", " ").trim().toLowerCase(Locale.ROOT);
	}

	private static String conciseMessage(Throwable error) {
		String message = error.getMessage();
		return message == null || message.isBlank() ? error.getClass().getSimpleName() : message;
	}
}
