/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.gpu.cpu;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuDevice;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.util.Locale;

import static org.junit.jupiter.api.Assertions.assertTrue;

class CpuBackendTest {
	@Test
	void bundledOneMklLoadsWhenProfileIsEnabled() {
		Assumptions.assumeTrue(Boolean.getBoolean("eqtl.test.mkl"),
			"The current Maven profile does not bundle oneMKL");
		String previous = System.getProperty(CpuBackend.BLAS_PROPERTY);
		try {
			System.setProperty(CpuBackend.BLAS_PROPERTY, "mkl");
			GpuDevice device = new CpuBackend().discoverGpuDevices().get(0);
			try (GpuContext ignored = device.openContext()) {
				assertTrue(device.getComputeApiVersion().contains("Math Kernel Library"));
			}
		} finally {
			restoreProperty(CpuBackend.BLAS_PROPERTY, previous);
		}
	}

	@Test
	void bundledOpenBlasLoadsOnSupportedDesktopPlatform() {
		Assumptions.assumeTrue(bundledPlatform(), "No bundled OpenBLAS classifier for this platform");
		String previous = System.getProperty(CpuBackend.BLAS_PROPERTY);
		try {
			System.setProperty(CpuBackend.BLAS_PROPERTY, "openblas");
			GpuDevice device = new CpuBackend().discoverGpuDevices().get(0);
			try (GpuContext ignored = device.openContext()) {
				assertTrue(device.getComputeApiVersion().startsWith("OpenBLAS 0.3.34"));
			}
		} finally {
			restoreProperty(CpuBackend.BLAS_PROPERTY, previous);
		}
	}

	@Test
	void javaFallbackCanBeSelectedWithoutNativeLoading() {
		String previous = System.getProperty(CpuBackend.BLAS_PROPERTY);
		try {
			System.setProperty(CpuBackend.BLAS_PROPERTY, "java");
			GpuDevice device = new CpuBackend().discoverGpuDevices().get(0);
			try (GpuContext ignored = device.openContext()) {
				assertTrue(device.getComputeApiVersion().startsWith("portable Java"));
			}
		} finally {
			restoreProperty(CpuBackend.BLAS_PROPERTY, previous);
		}
	}

	private static boolean bundledPlatform() {
		String os = System.getProperty("os.name", "").toLowerCase(Locale.ROOT);
		String arch = System.getProperty("os.arch", "").toLowerCase(Locale.ROOT);
		boolean x64 = arch.equals("amd64") || arch.equals("x86_64");
		boolean arm64 = arch.equals("aarch64") || arch.equals("arm64");
		return (os.contains("windows") && x64)
			|| ((os.contains("linux") || os.contains("mac")) && (x64 || arm64));
	}

	private static void restoreProperty(String name, String value) {
		if (value == null) System.clearProperty(name);
		else System.setProperty(name, value);
	}
}
