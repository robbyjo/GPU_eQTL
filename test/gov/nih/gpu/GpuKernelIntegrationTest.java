/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import gov.nih.eqtl.QeQTLAnalysis;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

class GpuKernelIntegrationTest {
    @Test
    void realEqtlKernelMatchesCpuDotProductWhenFp64GpuIsPresent() {
        List<GpuDevice> devices;
        try {
            devices = GpuRuntime.createDefault().getGpuDevices(true, true);
        } catch (RuntimeException | UnsatisfiedLinkError e) {
            Assumptions.abort("No usable OpenCL runtime: " + e.getMessage());
            return;
        }
        Assumptions.assumeFalse(devices.isEmpty(), "No available FP64 OpenCL GPU");

        int blockSize = 16;
        int rows = 64;
        int columns = 16;
        double[] expression = new double[blockSize * rows];
        double[] snps = new double[rows * columns];
        double expected = 0;
        for (int row = 0; row < rows; row++) {
            double expressionValue = row - 31.5;
            double snpValue = (row % 7) - 3;
            expression[row] = expressionValue;
            snps[row * columns] = snpValue;
            expected += expressionValue * snpValue;
        }
        expected /= rows - 1.0;

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
                blockSize * columns,
                localMemoryBytes,
                rows,
                columns,
                new long[] { blockSize, blockSize },
                new long[] { blockSize, blockSize });

            assertEquals(expected, result[0], 1e-12);
            assertEquals(0.0, result[1], 0.0);
            assertEquals(0.0, result[columns], 0.0);
        }
    }
}
