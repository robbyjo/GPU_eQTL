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
import gov.nih.gpu.cpu.CpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;
import gov.nih.jama.QRDecomposition;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;

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
	void cpuBackendFp64MatchesScalarReference() {
		assertRealEqtlOperation(new CpuBackend());
	}

	@Test
	void portableJavaCpuEngineMatchesScalarReference() {
		String previous = System.getProperty(CpuBackend.BLAS_PROPERTY);
		try {
			System.setProperty(CpuBackend.BLAS_PROPERTY, "java");
			assertRealEqtlOperation(new CpuBackend());
		} finally {
			restoreProperty(CpuBackend.BLAS_PROPERTY, previous);
		}
	}

	@Test
	void oneMklCpuEngineFp64MatchesScalarReferenceWhenBundled() {
		withOneMkl(() -> assertRealEqtlOperation(new CpuBackend()));
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

	@Test
	void cpuBackendFp32MatchesScalarReference() {
		assertRealEqtlOperationFp32(new CpuBackend());
	}

	@Test
	void oneMklCpuEngineFp32MatchesScalarReferenceWhenBundled() {
		withOneMkl(() -> assertRealEqtlOperationFp32(new CpuBackend()));
	}

	@Test
	void cudaFp64ResidualizationMatchesCpuReferenceWhenAvailable() {
		assertResidualization(new CudaGpuBackend(), GpuPrecision.FP64);
	}

	@Test
	void openClFp64ResidualizationMatchesCpuReferenceWhenAvailable() {
		assertResidualization(new JoclGpuBackend(), GpuPrecision.FP64);
	}

	@Test
	void cpuFp64ResidualizationMatchesScalarReference() {
		assertResidualization(new CpuBackend(), GpuPrecision.FP64);
	}

	@Test
	void oneMklCpuEngineResidualizationMatchesScalarReferenceWhenBundled() {
		withOneMkl(() -> assertResidualization(new CpuBackend(), GpuPrecision.FP64));
	}

	@Test
	void cudaFp32ResidualizationMatchesCpuReferenceWhenAvailable() {
		assertResidualization(new CudaGpuBackend(), GpuPrecision.FP32);
	}

	@Test
	void openClFp32ResidualizationMatchesCpuReferenceWhenAvailable() {
		assertResidualization(new JoclGpuBackend(), GpuPrecision.FP32);
	}

	@Test
	void cpuFp32ResidualizationMatchesScalarReference() {
		assertResidualization(new CpuBackend(), GpuPrecision.FP32);
	}

	@Test
	void cudaCompactPatternStatisticsMatchExplicitQrWhenAvailable() {
		assertCompactPatternStatistics(new CudaGpuBackend());
	}

	@Test
	void openClCompactPatternStatisticsMatchExplicitQrWhenAvailable() {
		assertCompactPatternStatistics(new JoclGpuBackend());
	}

	@Test
	void cpuCompactPatternStatisticsMatchExplicitQr() {
		assertCompactPatternStatistics(new CpuBackend());
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

	private static void assertResidualization(GpuBackend backend, GpuPrecision precision) {
		List<GpuDevice> devices;
		try {
			devices = new GpuRuntime(backend).getGpuDevices(true, precision.requiresDoublePrecision());
		} catch (RuntimeException | LinkageError e) {
			Assumptions.abort("No usable " + backend.getName() + " runtime: " + e.getMessage());
			return;
		}
		Assumptions.assumeFalse(devices.isEmpty(), "No available " + precision.optionName()
			+ " device for " + backend.getName());

		int rows = 3;
		int samples = 8;
		int rank = 2;
		double inverseRootN = 1.0 / Math.sqrt(samples);
		double centeredNorm = Math.sqrt(42.0);
		double[] q = new double[samples * rank];
		for (int sample = 0; sample < samples; sample++) {
			q[sample * rank] = inverseRootN;
			q[sample * rank + 1] = (sample - 3.5) / centeredNorm;
		}
		double[] values = new double[] {
			2, 3, 5, 7, 11, 13, 17, 19,
			8, 5, 9, 2, 6, 5, 3, 5,
			1, 4, 1, 4, 2, 1, 3, 5
		};
		double[] expected = residualizeOnCpu(values, q, rows, samples, rank);

		try (GpuContext context = devices.get(0).openContext()) {
			context.setProfilingEnabled(true);
			boolean hostOnly = "cpu".equalsIgnoreCase(context.getDevice().getBackendName());
			if (precision == GpuPrecision.FP64) {
				double[] actual = context.residualizeDoubleRows(values, q, rows, samples, rank);
				assertArrayEquals(expected, actual, 2e-12, backend.getName() + " FP64 residualization");
				assertEquals(hostOnly ? 0 : (long) (values.length + q.length) * Double.BYTES,
					context.getLastExecutionMetrics().uploadedBytes());
				context.residualizeDoubleRows(values, q, rows, samples, rank);
				assertEquals(hostOnly ? 0 : (long) values.length * Double.BYTES,
					context.getLastExecutionMetrics().uploadedBytes(), "Q should be uploaded once per context");
			} else {
				float[] values32 = toFloat(values);
				float[] q32 = toFloat(q);
				float[] actual = context.residualizeFloatRows(values32, q32, rows, samples, rank);
				for (int i = 0; i < actual.length; i++)
					assertEquals(expected[i], actual[i], 2e-5,
						backend.getName() + " FP32 residualization element " + i);
				assertEquals(hostOnly ? 0 : (long) (values32.length + q32.length) * Float.BYTES,
					context.getLastExecutionMetrics().uploadedBytes());
				context.residualizeFloatRows(values32, q32, rows, samples, rank);
				assertEquals(hostOnly ? 0 : (long) values32.length * Float.BYTES,
					context.getLastExecutionMetrics().uploadedBytes(), "Q should be uploaded once per context");
			}
			assertArrayEquals(new double[] {
				2, 3, 5, 7, 11, 13, 17, 19,
				8, 5, 9, 2, 6, 5, 3, 5,
				1, 4, 1, 4, 2, 1, 3, 5 }, values, 0.0,
				"Residualization must not mutate the host input");
		}
	}

	private static double[] residualizeOnCpu(double[] values, double[] q,
		int rows, int samples, int rank) {
		double[] result = values.clone();
		for (int row = 0; row < rows; row++) {
			double[] coefficients = new double[rank];
			for (int column = 0; column < rank; column++)
				for (int sample = 0; sample < samples; sample++)
					coefficients[column] += values[row * samples + sample] * q[sample * rank + column];
			for (int sample = 0; sample < samples; sample++)
				for (int column = 0; column < rank; column++)
					result[row * samples + sample] -= coefficients[column] * q[sample * rank + column];
		}
		return result;
	}

	private static void assertCompactPatternStatistics(GpuBackend backend) {
		List<GpuDevice> devices;
		try {
			devices = new GpuRuntime(backend).getGpuDevices(true, true);
		} catch (RuntimeException | LinkageError e) {
			Assumptions.abort("No usable " + backend.getName() + " runtime: " + e.getMessage());
			return;
		}
		Assumptions.assumeFalse(devices.isEmpty(),
			"No available FP64 device for " + backend.getName());
		double[][] design = {
			{1, -3.5}, {1, -2.5}, {1, -1.5}, {1, -0.5},
			{1, 0.5}, {1, 1.5}, {1, 2.5}, {1, 3.5}
		};
		int[][] observed = {
			{0, 1, 2, 3, 4, 5, 6, 7},
			{0, 1, 2, 4, 5, 6, 7}
		};
		double[][] genotypes = {
			{0, 1, 2, Double.NaN, 1, 2, 0, 1.5},
			{2, 1, 0, 1, Double.NaN, 0.5, 1.5, 2},
			{0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2}
		};
		double[][] sums = new double[observed.length][design[0].length];
		double[][][] upper = new double[observed.length][][];
		for (int pattern = 0; pattern < observed.length; pattern++) {
			double[][] selected = new double[observed[pattern].length][design[0].length];
			for (int row = 0; row < observed[pattern].length; row++) {
				selected[row] = design[observed[pattern][row]].clone();
				for (int column = 0; column < design[0].length; column++)
					sums[pattern][column] += selected[row][column];
			}
			upper[pattern] = new QRDecomposition(selected).getR().getArray();
		}
		GpuPatternStatisticsPlan plan = new GpuPatternStatisticsPlan(design,
			new int[] {3, 9}, observed, sums, upper);
		int paddedSamples = 64;
		int variantCapacity = 16;
		double[] aggregate = new double[paddedSamples * variantCapacity * 3];
		for (int sample = 0; sample < design.length; sample++) {
			for (int variant = 0; variant < genotypes.length; variant++) {
				double value = genotypes[variant][sample];
				if (Double.isFinite(value)) {
					int base = sample * variantCapacity * 3 + variant * 3;
					aggregate[base] = 1;
					aggregate[base + 1] = value;
					aggregate[base + 2] = value * value;
				}
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
			+ "#define N_MIN_1 1" + line
			+ QeQTLAnalysis.eqtlReal;
		try (GpuContext context = devices.get(0).openContext()) {
			context.setProfilingEnabled(true);
			context.compileKernel(source, "eqtlReal", GpuPrecision.FP64);
			GpuPatternStatisticsResult actual = context.computePatternStatisticsDouble(
				aggregate, paddedSamples, genotypes.length, variantCapacity,
				plan, true, 1, 16);
			for (int pattern = 0; pattern < observed.length; pattern++) {
				double[][] selectedDesign = new double[observed[pattern].length][design[0].length];
				for (int row = 0; row < observed[pattern].length; row++)
					selectedDesign[row] = design[observed[pattern][row]].clone();
				double[][] q = new QRDecomposition(selectedDesign).getQ().getArray();
				for (int variant = 0; variant < genotypes.length; variant++) {
					double mean = 0;
					int called = 0;
					for (int sample : observed[pattern]) {
						double value = genotypes[variant][sample];
						if (Double.isFinite(value)) { mean += value; called++; }
					}
					mean /= called;
					double[] filled = new double[observed[pattern].length];
					double filledSquares = 0;
					for (int row = 0; row < observed[pattern].length; row++) {
						double value = genotypes[variant][observed[pattern][row]];
						filled[row] = Double.isFinite(value) ? value : mean;
						filledSquares += filled[row] * filled[row];
					}
					double[] residual = residualizeOnCpu(filled, flatten(q), 1,
						filled.length, q[0].length);
					double residualSquares = 0;
					for (double value : residual) residualSquares += value * value;
					assertEquals(mean, actual.value(pattern, variant,
						GpuPatternStatisticsResult.REPLACEMENT), 2e-12,
						backend.getName() + " pattern " + pattern + " variant " + variant
							+ " replacement " + java.util.Arrays.toString(actual.values()));
					assertEquals(filledSquares, actual.value(pattern, variant,
						GpuPatternStatisticsResult.FILLED_SUM_SQUARES), 2e-10);
					assertEquals(residualSquares, actual.value(pattern, variant,
						GpuPatternStatisticsResult.RESIDUAL_SUM_SQUARES), 2e-10);
					assertEquals(called, actual.value(pattern, variant,
						GpuPatternStatisticsResult.CALLED_COUNT), 2e-12);
					assertEquals(mean * called, actual.value(pattern, variant,
						GpuPatternStatisticsResult.DOSAGE_SUM), 2e-12);
				}
			}
			boolean hostOnly = "cpu".equalsIgnoreCase(context.getDevice().getBackendName());
			assertEquals(hostOnly ? 0 : (long) actual.values().length * Double.BYTES,
				context.getLastExecutionMetrics().downloadedBytes());
		}
	}

	private static double[] flatten(double[][] values) {
		double[] result = new double[values.length * values[0].length];
		for (int row = 0; row < values.length; row++)
			System.arraycopy(values[row], 0, result, row * values[row].length,
				values[row].length);
		return result;
	}

	private static float[] toFloat(double[] values) {
		float[] result = new float[values.length];
		for (int i = 0; i < values.length; i++)
			result[i] = (float) values[i];
		return result;
	}

	private static void restoreProperty(String name, String value) {
		if (value == null) System.clearProperty(name);
		else System.setProperty(name, value);
	}

	private static void withOneMkl(Runnable check) {
		Assumptions.assumeTrue(Boolean.getBoolean("eqtl.test.mkl"),
			"The current Maven profile does not bundle oneMKL");
		String previous = System.getProperty(CpuBackend.BLAS_PROPERTY);
		try {
			System.setProperty(CpuBackend.BLAS_PROPERTY, "mkl");
			check.run();
		} finally {
			restoreProperty(CpuBackend.BLAS_PROPERTY, previous);
		}
	}
}
