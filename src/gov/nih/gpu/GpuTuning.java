/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import java.util.List;

/** Deterministic device-memory and worker-count recommendations. */
public final class GpuTuning {
    private static final long MIB = 1024L * 1024L;
    private static final long MINIMUM_VRAM_HEADROOM = 256L * MIB;
    private static final long TARGET_OUTPUT_ALLOCATION = 1024L * MIB;
    private static final int KERNEL_TILE = 16;
    private static final int LARGE_BLOCK_GRANULARITY = 512;

    public record BlockRecommendation(int blockSize, long estimatedDeviceBytes,
        String limitingDevice) { }

    private GpuTuning() { }

    public static BlockRecommendation recommendBlockSize(List<? extends GpuDevice> devices,
        int sampleCount, int snpCount, int traitCount, GpuPrecision precision) {
        if (devices == null || devices.isEmpty())
            throw new IllegalArgumentException("At least one GPU is required for automatic block sizing");
        if (sampleCount <= 1 || snpCount <= 0 || traitCount <= 0)
            throw new IllegalArgumentException("Matrix dimensions must be positive");
        if (precision == null)
            throw new IllegalArgumentException("precision must not be null");

        long paddedSamples = roundUp(sampleCount, GpuRuntime.DEFAULT_ALIGNMENT);
        long matrixLimit = Math.min(roundUp(snpCount, KERNEL_TILE), roundUp(traitCount, KERNEL_TILE));
        long candidate = Math.min(matrixLimit, (long) Math.sqrt(Integer.MAX_VALUE));
        String limitingDevice = devices.get(0).getName();

        for (GpuDevice device : devices) {
            long totalBytes = device.getGlobalMemoryBytes();
            long maxAllocation = device.getMaxAllocationBytes();
            if (totalBytes <= 0)
                totalBytes = maxAllocation;
            if (maxAllocation <= 0)
                maxAllocation = totalBytes;
            if (totalBytes <= 0 || maxAllocation <= 0)
                throw new IllegalArgumentException("GPU does not report usable memory limits: " + device.getName());

            long headroom = Math.max(MINIMUM_VRAM_HEADROOM, totalBytes / 8);
            long usableBytes = totalBytes > headroom ? totalBytes - headroom : totalBytes / 2;
            long outputBudget = Math.min(maxAllocation, TARGET_OUTPUT_ALLOCATION);
            long outputLimit = floorSqrt(outputBudget / precision.bytes());
            long inputLimit = maxAllocation / Math.multiplyExact(paddedSamples, precision.bytes());
            long totalElements = usableBytes / precision.bytes();
            long totalLimit = (long) Math.floor(Math.sqrt((double) paddedSamples * paddedSamples
                + totalElements) - paddedSamples);
            long deviceLimit = Math.min(outputLimit, Math.min(inputLimit, totalLimit));
            if (deviceLimit < candidate) {
                candidate = deviceLimit;
                limitingDevice = device.getName();
            }
        }

        int granularity = candidate >= LARGE_BLOCK_GRANULARITY
            ? LARGE_BLOCK_GRANULARITY : KERNEL_TILE;
        int blockSize = (int) (candidate - candidate % granularity);
        if (blockSize < KERNEL_TILE)
            throw new IllegalArgumentException("GPU memory limits cannot hold the minimum 16-row tile");
        long estimated = estimateDeviceBytes(blockSize, sampleCount, precision);
        return new BlockRecommendation(blockSize, estimated, limitingDevice);
    }

    public static int recommendThreadCount(int cpuCores, int deviceCount,
        long requiredIterations, boolean streamed) {
        if (cpuCores <= 0 || deviceCount <= 0 || requiredIterations <= 0)
            throw new IllegalArgumentException("CPU, GPU, and iteration counts must be positive");
        int cpuBudget = Math.max(1, cpuCores - 1);
        int workersPerDevice = streamed ? 1 : 2;
        long deviceBudget = Math.multiplyExact((long) deviceCount, workersPerDevice);
        return (int) Math.max(1, Math.min(requiredIterations, Math.min(cpuBudget, deviceBudget)));
    }

    public static long estimateDeviceBytes(int blockSize, int sampleCount, GpuPrecision precision) {
        long paddedSamples = roundUp(sampleCount, GpuRuntime.DEFAULT_ALIGNMENT);
        long elements = Math.addExact(Math.multiplyExact(2L * blockSize, paddedSamples),
            Math.multiplyExact((long) blockSize, blockSize));
        return Math.multiplyExact(elements, precision.bytes());
    }

    private static long floorSqrt(long value) {
        if (value <= 0)
            return 0;
        long root = (long) Math.sqrt(value);
        while (root + 1 <= value / (root + 1))
            root++;
        while (root > value / root)
            root--;
        return root;
    }

    private static long roundUp(long value, int multiple) {
        long remainder = value % multiple;
        return remainder == 0 ? value : Math.addExact(value, multiple - remainder);
    }
}
