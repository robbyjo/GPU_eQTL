/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import gov.nih.eqtl.QeQTLAnalysis;
import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class GpuKernelIntegrationTest {
	@Test
	void automaticBackendMatchesCpuReferenceWhenFp64GpuIsPresent() {
		assertRealEqtlOperation(new AutoGpuBackend());
	}

	@Test
	void cudaBackendMatchesCpuReferenceWhenAvailable() {
		assertRealEqtlOperation(new CudaGpuBackend());
	}

	@Test
	void openClBackendMatchesCpuReferenceWhenAvailable() {
		assertRealEqtlOperation(new JoclGpuBackend());
	}

	@Test
	void automaticBackendFp32MatchesCpuReferenceWhenGpuIsPresent() {
		assertRealEqtlOperationFp32(new AutoGpuBackend());
	}

	@Test
	void cudaBackendFp32MatchesCpuReferenceWhenAvailable() {
		assertRealEqtlOperationFp32(new CudaGpuBackend());
	}

	@Test
	void openClBackendFp32MatchesCpuReferenceWhenAvailable() {
		assertRealEqtlOperationFp32(new JoclGpuBackend());
	}

	private static void assertRealEqtlOperation(GpuBackend backend) {
		List<GpuDevice> devices;
		try {
			devices = new GpuRuntime(backend).getGpuDevices(true, true);
		} catch (RuntimeException | LinkageError e) {
			Assumptions.abort("No usable " + backend.getName() + " runtime: " + e.getMessage());
			return;
		}
		Assumptions.assumeFalse(devices.isEmpty(), "No available FP64 device for " + backend.getName());

		int blockSize = 16;
		int rows = 64;
		int capacity = 32;
		int actualTraits = 3;
		int actualSnps = 5;
		double[] expression = new double[capacity * rows];
		double[] snps = new double[rows * capacity];
		for (int trait = 0; trait < actualTraits; trait++) {
			for (int row = 0; row < rows; row++) {
				expression[trait * rows + row] = (row - 31.5) * (trait + 1) + (row % (trait + 3));
			}
		}
		for (int row = 0; row < rows; row++) {
			for (int snp = 0; snp < actualSnps; snp++) {
				snps[row * capacity + snp] = (row % (snp + 4)) - (snp + 1.5);
			}
		}

		String line = System.lineSeparator();
		String source = "#define BLOCK_SIZE 16" + line
			+ "#define DATATYPE double" + line
			+ "#if defined(cl_khr_fp64)" + line
			+ "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" + line
			+ "#elif defined(cl_amd_fp64)" + line
			+ "#pragma OPENCL EXTENSION cl_amd_fp64 : enable" + line
			+ "#endif" + line
			+ "#define N_MIN_1 63" + line
			+ QeQTLAnalysis.eqtlReal;

		try (GpuContext context = devices.get(0).openContext()) {
			context.compileKernel(source, "eqtlReal");
			long localMemoryBytes = (long) (blockSize + 1) * (4L * blockSize) * Double.BYTES;
			double[] result = context.executeDoubleKernel(
				expression,
				snps,
				capacity * capacity,
				localMemoryBytes,
				rows,
				capacity,
				new long[] { blockSize, blockSize },
				new long[] { blockSize, blockSize });

			for (int trait = 0; trait < actualTraits; trait++) {
				for (int snp = 0; snp < actualSnps; snp++) {
					double expected = 0.0;
					for (int row = 0; row < rows; row++) {
						expected += expression[trait * rows + row] * snps[row * capacity + snp];
					}
					expected /= rows - 1.0;
					assertEquals(expected, result[trait * capacity + snp], 1e-11,
						backend.getName() + " trait " + trait + ", SNP " + snp);
				}
			}
			assertEquals(0.0, result[actualSnps], 0.0);
			assertEquals(0.0, result[actualTraits * capacity], 0.0);
		}
	}

	private static void assertRealEqtlOperationFp32(GpuBackend backend) {
		List<GpuDevice> devices;
		try {
			devices = new GpuRuntime(backend).getGpuDevices(true, false);
		} catch (RuntimeException | LinkageError e) {
			Assumptions.abort("No usable " + backend.getName() + " runtime: " + e.getMessage());
			return;
		}
		Assumptions.assumeFalse(devices.isEmpty(), "No available FP32 device for " + backend.getName());

		int blockSize = 16;
		int rows = 64;
		int capacity = 32;
		int actualTraits = 3;
		int actualSnps = 5;
		float[] expression = new float[capacity * rows];
		float[] snps = new float[rows * capacity];
		for (int trait = 0; trait < actualTraits; trait++)
			for (int row = 0; row < rows; row++)
				expression[trait * rows + row] = (float) ((row - 31.5) * (trait + 1) + (row % (trait + 3)));
		for (int row = 0; row < rows; row++)
			for (int snp = 0; snp < actualSnps; snp++)
				snps[row * capacity + snp] = (float) ((row % (snp + 4)) - (snp + 1.5));

		String line = System.lineSeparator();
		String source = "#define BLOCK_SIZE 16" + line
			+ "#define DATATYPE float" + line
			+ "#define N_MIN_1 63" + line
			+ QeQTLAnalysis.eqtlReal;

		try (GpuContext context = devices.get(0).openContext()) {
			context.compileKernel(source, "eqtlReal", GpuPrecision.FP32);
			long localMemoryBytes = (long) (blockSize + 1) * (4L * blockSize) * Float.BYTES;
			float[] result = context.executeFloatKernel(expression, snps, capacity * capacity,
				localMemoryBytes, rows, capacity,
				new long[] { blockSize, blockSize }, new long[] { blockSize, blockSize });

			for (int trait = 0; trait < actualTraits; trait++) {
				for (int snp = 0; snp < actualSnps; snp++) {
					double expected = 0.0;
					for (int row = 0; row < rows; row++)
						expected += (double) expression[trait * rows + row] * snps[row * capacity + snp];
					expected /= rows - 1.0;
					double tolerance = Math.max(2e-5, Math.abs(expected) * 5e-6);
					assertEquals(expected, result[trait * capacity + snp], tolerance,
						backend.getName() + " FP32 trait " + trait + ", SNP " + snp);
				}
			}
		}
	}
}
