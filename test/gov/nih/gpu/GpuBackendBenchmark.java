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

import java.util.Arrays;
import java.util.List;

/**
 * Manual, transfer-inclusive backend comparison. This is deliberately not a
 * unit test: timings are indicative and should be repeated with representative
 * sample counts and tile sizes on the deployment machine.
 */
public final class GpuBackendBenchmark {
	private static final int WARMUP_RUNS = 3;
	private static final int MEASURED_RUNS = 9;
	private static volatile double sink;

	private GpuBackendBenchmark() { }

	public static void main(String[] arguments) {
		int[][] shapes;
		if (arguments.length == 2) {
			shapes = new int[][] { { Integer.parseInt(arguments[0]), Integer.parseInt(arguments[1]) } };
		} else if (arguments.length == 0) {
			shapes = new int[][] { { 128, 2048 }, { 512, 1024 }, { 2048, 512 } };
		} else {
			throw new IllegalArgumentException("Usage: GpuBackendBenchmark [sampleRows tileWidth]");
		}

		for (int[] shape : shapes) {
			run(new CudaGpuBackend(), shape[0], shape[1]);
			run(new JoclGpuBackend(), shape[0], shape[1]);
			run(new CpuBackend(), shape[0], shape[1]);
		}
		if (sink == Double.NEGATIVE_INFINITY) {
			System.out.println("Unreachable result guard");
		}
	}

	private static void run(GpuBackend backend, int rows, int width) {
		if (rows <= 1 || rows % 64 != 0 || width <= 0 || width % 16 != 0) {
			throw new IllegalArgumentException("sampleRows must be a multiple of 64 and tileWidth a multiple of 16");
		}
		List<GpuDevice> devices;
		try {
			devices = new GpuRuntime(backend).getGpuDevices(true, true);
		} catch (RuntimeException | LinkageError e) {
			System.out.println(backend.getName() + ": unavailable (" + e.getMessage() + ")");
			return;
		}
		if (devices.isEmpty()) {
			System.out.println(backend.getName() + ": no FP64 device");
			return;
		}

		int blockSize = 16;
		double[] inputA = new double[Math.multiplyExact(width, rows)];
		double[] inputB = new double[Math.multiplyExact(rows, width)];
		for (int i = 0; i < inputA.length; i++) inputA[i] = ((i * 17L) % 101 - 50) / 51.0;
		for (int i = 0; i < inputB.length; i++) inputB[i] = ((i * 31L) % 97 - 48) / 49.0;
		String line = System.lineSeparator();
		String source = "#define BLOCK_SIZE 16" + line + "#define DATATYPE double" + line
			+ "#if defined(cl_khr_fp64)" + line + "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" + line
			+ "#elif defined(cl_amd_fp64)" + line + "#pragma OPENCL EXTENSION cl_amd_fp64 : enable" + line
			+ "#endif" + line + "#define N_MIN_1 " + (rows - 1) + line + QeQTLAnalysis.eqtlReal;
		long localBytes = (long) (blockSize + 1) * (4L * blockSize) * Double.BYTES;
		long[] global = { width, width };
		long[] local = { blockSize, blockSize };
		long[] samples = new long[MEASURED_RUNS];
		try (GpuContext context = devices.get(0).openContext()) {
			context.compileKernel(source, "eqtlReal");
			for (int i = 0; i < WARMUP_RUNS; i++) {
				double[] result = execute(context, inputA, inputB, rows, width, localBytes, global, local);
				sink += result[(i * 7919) % result.length];
			}
			for (int i = 0; i < samples.length; i++) {
				long start = System.nanoTime();
				double[] result = execute(context, inputA, inputB, rows, width, localBytes, global, local);
				samples[i] = System.nanoTime() - start;
				sink += result[(i * 7919) % result.length];
			}
		}
		Arrays.sort(samples);
		System.out.printf("%s on %s, rows=%d, tile=%dx%d: median %.3f ms, best %.3f ms (including copies)%n",
			backend.getName(), devices.get(0).getName(), rows, width, width,
			samples[samples.length / 2] / 1e6, samples[0] / 1e6);
	}

	private static double[] execute(GpuContext context, double[] inputA, double[] inputB,
			int rows, int width, long localBytes, long[] global, long[] local) {
		return context.executeDoubleKernel(inputA, inputB, Math.multiplyExact(width, width), localBytes,
			rows, width, global, local);
	}
}
