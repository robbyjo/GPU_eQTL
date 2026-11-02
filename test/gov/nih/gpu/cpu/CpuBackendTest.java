/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.gpu.cpu;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuDevice;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.util.Locale;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.TimeUnit;

import static org.junit.jupiter.api.Assertions.assertEquals;
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
	void openBlasGemmAndThreadControlUseOneNativeLibrary() throws Exception {
		Assumptions.assumeTrue(bundledPlatform(), "No bundled OpenBLAS classifier for this platform");
		String java = Path.of(System.getProperty("java.home"), "bin",
			System.getProperty("os.name").startsWith("Windows") ? "java.exe" : "java").toString();
		ProcessBuilder builder = new ProcessBuilder(java, "-XX:ActiveProcessorCount=2",
			"-Xmx128m", "-cp", System.getProperty("java.class.path"),
			OpenBlasThreadProbe.class.getName()).redirectErrorStream(true);
		builder.environment().put("OPENBLAS_NUM_THREADS", "1");
		Process process = builder.start();
		try {
			assertTrue(process.waitFor(60, TimeUnit.SECONDS), "Isolated OpenBLAS probe timed out");
			String output = new String(process.getInputStream().readAllBytes(), StandardCharsets.UTF_8);
			assertEquals(0, process.exitValue(), output);
			assertTrue(output.contains("Single OpenBLAS library verified"), output);
		} finally {
			if (process.isAlive()) process.destroyForcibly();
		}
	}

	/** A fresh JVM keeps previously loaded test backends from masking duplicate libraries. */
	public static final class OpenBlasThreadProbe {
		public static void main(String[] args) {
			OpenBlasMatrixEngine engine = new OpenBlasMatrixEngine(2);
			var libraries = org.bytedeco.javacpp.Loader.getLoadedLibraries();
			if (!libraries.containsKey("jniopenblas_nolapack") || libraries.containsKey("jniopenblas"))
				throw new AssertionError("GEMM/thread control loaded different OpenBLAS libraries: "
					+ libraries.keySet().stream().filter(name -> name.contains("openblas")).toList());
			if (!engine.description().contains("requested 2 BLAS threads"))
				throw new AssertionError("Description must distinguish requested threads: " + engine.description());
			double[] left = {1, 2, 3, 4, 5, 6};
			double[] right = {2, 1, 0, 3, 4, -1};
			double[] product = new double[4];
			engine.multiplyDouble(left, right, product, 2, 2, 3, 3, 2, 2, 1, 0);
			double[] expected = {14, 4, 32, 13};
			for (int i = 0; i < product.length; i++)
				if (!Double.isFinite(product[i]) || Math.abs(product[i] - expected[i]) > 1e-12)
					throw new AssertionError("FP64 product " + i + ": " + product[i]);
			System.out.println("Single OpenBLAS library verified; " + engine.description());
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
