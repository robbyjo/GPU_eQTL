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
    private static final int PIPELINE_WORKERS_PER_DEVICE = 4;
    private static final int FULL_MEMORY_WORKERS_PER_DEVICE = 2;
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

        long desiredJobs = Math.min((long) snpCount,
            Math.multiplyExact((long) devices.size(), PIPELINE_WORKERS_PER_DEVICE));
        long concurrencyLimit = (snpCount + desiredJobs - 1) / desiredJobs;
        int concurrencyGranularity = concurrencyLimit >= LARGE_BLOCK_GRANULARITY
            ? LARGE_BLOCK_GRANULARITY : KERNEL_TILE;
        concurrencyLimit -= concurrencyLimit % concurrencyGranularity;
        concurrencyLimit = Math.max(KERNEL_TILE, concurrencyLimit);
        if (concurrencyLimit < candidate) {
            candidate = concurrencyLimit;
            limitingDevice = "workload concurrency";
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
        return recommendThreadCount(cpuCores, deviceCount, requiredIterations, streamed,
            Long.MAX_VALUE, 1);
    }

    public static int recommendThreadCount(int cpuCores, int deviceCount,
        long requiredIterations, boolean streamed, long availableHeapBytes,
        long estimatedWorkerBytes) {
        if (cpuCores <= 0 || deviceCount <= 0 || requiredIterations <= 0)
            throw new IllegalArgumentException("CPU, GPU, and iteration counts must be positive");
        if (availableHeapBytes <= 0 || estimatedWorkerBytes <= 0)
            throw new IllegalArgumentException("Heap and per-worker memory estimates must be positive");
        int cpuBudget = Math.max(1, cpuCores - 1);
        int workersPerDevice = streamed ? PIPELINE_WORKERS_PER_DEVICE : FULL_MEMORY_WORKERS_PER_DEVICE;
        long deviceBudget = Math.multiplyExact((long) deviceCount, workersPerDevice);
        long heapHeadroom = Math.max(MINIMUM_VRAM_HEADROOM, availableHeapBytes / 4);
        long usableHeap = availableHeapBytes > heapHeadroom
            ? availableHeapBytes - heapHeadroom : Math.max(1, availableHeapBytes / 2);
        long heapBudget = Math.max(1, usableHeap / estimatedWorkerBytes);
        return (int) Math.max(1, Math.min(requiredIterations,
            Math.min(cpuBudget, Math.min(deviceBudget, heapBudget))));
    }

    public static long estimateStreamedWorkerBytes(int sampleCount, int genotypeRows,
        int expressionRows, GpuPrecision precision) {
        if (sampleCount <= 1 || genotypeRows <= 0 || expressionRows <= 0 || precision == null)
            throw new IllegalArgumentException("Streamed worker dimensions and precision must be positive");
        long paddedSamples = roundUp(sampleCount, GpuRuntime.DEFAULT_ALIGNMENT);
        long preparedElements = Math.multiplyExact((long) genotypeRows + expressionRows, sampleCount);
        long packedElements = Math.multiplyExact((long) genotypeRows + expressionRows, paddedSamples);
        long resultElements = Math.multiplyExact((long) genotypeRows, expressionRows);
        long preparedBytes = Math.multiplyExact(preparedElements, Double.BYTES);
        long packedBytes = Math.multiplyExact(packedElements, precision.bytes());
        long resultBytes = Math.multiplyExact(resultElements, precision.bytes());
        return Math.addExact(1L * MIB, Math.addExact(preparedBytes,
            Math.addExact(packedBytes, resultBytes)));
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
