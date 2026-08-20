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

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

class AutoGpuBackendTest {
	@Test
	void cudaReplacesOnlyItsDuplicateOpenClDevice() {
		FakeDevice cudaNvidia = new FakeDevice("cuda", "NVIDIA Corporation", "RTX test", true);
		FakeDevice openclNvidia = new FakeDevice("opencl", "NVIDIA Corporation", "RTX test", true);
		FakeDevice openclIntel = new FakeDevice("opencl", "Intel", "Iris test", false);
		AutoGpuBackend backend = new AutoGpuBackend(List.of(
			new FakeBackend("cuda", List.of(cudaNvidia), null),
			new FakeBackend("opencl", List.of(openclNvidia, openclIntel), null)));

		List<GpuDevice> devices = backend.discoverGpuDevices();

		assertEquals(2, devices.size());
		assertSame(cudaNvidia, devices.get(0));
		assertSame(openclIntel, devices.get(1));
		assertEquals(List.of(cudaNvidia), new GpuRuntime(backend).getGpuDevices(true, true));
	}

	@Test
	void openClRemainsAvailableWhenCudaCannotInitialize() {
		FakeDevice opencl = new FakeDevice("opencl", "AMD", "Radeon test", true);
		AutoGpuBackend backend = new AutoGpuBackend(List.of(
			new FakeBackend("cuda", List.of(), new GpuException("CUDA missing")),
			new FakeBackend("opencl", List.of(opencl), null)));

		assertEquals(List.of(opencl), backend.discoverGpuDevices());
		assertTrue(backend.getRuntimeDescription().contains("cuda unavailable"));
	}

	@Test
	void fp32CudaStillReplacesItsNonFp64OpenClDuplicate() {
		FakeDevice cudaNvidia = new FakeDevice("cuda", "NVIDIA Corporation", "old RTX", false);
		FakeDevice openclNvidia = new FakeDevice("opencl", "NVIDIA Corporation", "old RTX", false);
		AutoGpuBackend backend = new AutoGpuBackend(List.of(
			new FakeBackend("cuda", List.of(cudaNvidia), null),
			new FakeBackend("opencl", List.of(openclNvidia), null)));

		assertEquals(List.of(cudaNvidia), backend.discoverGpuDevices());
		assertEquals(List.of(cudaNvidia), new GpuRuntime(backend).getGpuDevices(true, false));
	}

	@Test
	void systemPropertySelectsEachConcreteBackend() {
		String previous = System.getProperty("eqtl.gpu.backend");
		try {
			System.clearProperty("eqtl.gpu.backend");
			assertInstanceOf(AutoGpuBackend.class, GpuRuntime.createDefault().getBackend());
			System.setProperty("eqtl.gpu.backend", "cuda");
			assertInstanceOf(CudaGpuBackend.class, GpuRuntime.createDefault().getBackend());
			System.setProperty("eqtl.gpu.backend", "opencl");
			assertInstanceOf(JoclGpuBackend.class, GpuRuntime.createDefault().getBackend());
		} finally {
			if (previous == null) {
				System.clearProperty("eqtl.gpu.backend");
			} else {
				System.setProperty("eqtl.gpu.backend", previous);
			}
		}
	}

	private static final class FakeBackend implements GpuBackend {
		private final String name;
		private final List<GpuDevice> devices;
		private final RuntimeException failure;

		private FakeBackend(String name, List<GpuDevice> devices, RuntimeException failure) {
			this.name = name;
			this.devices = devices;
			this.failure = failure;
		}

		@Override
		public String getName() { return name; }

		@Override
		public String getRuntimeDescription() { return "test runtime"; }

		@Override
		public List<GpuDevice> discoverGpuDevices() {
			if (failure != null) throw failure;
			return devices;
		}
	}

	private static final class FakeDevice implements GpuDevice {
		private final String backend;
		private final String vendor;
		private final String name;
		private final boolean fp64;

		private FakeDevice(String backend, String vendor, String name, boolean fp64) {
			this.backend = backend;
			this.vendor = vendor;
			this.name = name;
			this.fp64 = fp64;
		}

		@Override
		public String getBackendName() { return backend; }

		@Override
		public String getPlatformName() { return "test"; }

		@Override
		public String getName() { return name; }

		@Override
		public String getVendor() { return vendor; }

		@Override
		public String getDriverVersion() { return "1"; }

		@Override
		public String getComputeApiVersion() { return "1"; }

		@Override
		public boolean isAvailable() { return true; }

		@Override
		public boolean isCompilerAvailable() { return true; }

		@Override
		public boolean supportsDoublePrecision() { return fp64; }

		@Override
		public boolean hasUnifiedMemory() { return false; }

		@Override
		public long getGlobalMemoryBytes() { return 1; }

		@Override
		public long getMaxAllocationBytes() { return 1; }

		@Override
		public long getMaxWorkGroupSize() { return 1; }

		@Override
		public long[] getMaxWorkItemSizes() { return new long[] { 1, 1 }; }

		@Override
		public GpuContext openContext() { throw new UnsupportedOperationException(); }
	}
}
