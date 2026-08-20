/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.gpu;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.List;

import org.junit.jupiter.api.Test;

class GpuTuningTest {
    @Test
    void blockSizeUsesPrecisionAndTheSmallestSelectedDevice() {
        FakeDevice large = new FakeDevice("large", 8L << 30, 8L << 30);
        FakeDevice smallAllocation = new FakeDevice("small", 4L << 30, 512L << 20);

        assertEquals(11264, GpuTuning.recommendBlockSize(List.of(large),
            2005, 250000, 136840, GpuPrecision.FP64).blockSize());
        assertEquals(16384, GpuTuning.recommendBlockSize(List.of(large),
            2005, 250000, 136840, GpuPrecision.FP32).blockSize());
        assertEquals(8192, GpuTuning.recommendBlockSize(List.of(large, smallAllocation),
            2005, 250000, 136840, GpuPrecision.FP64).blockSize());
    }

    @Test
    void threadRecommendationUsesAllGpusWithoutUsingEveryCpuCore() {
        assertEquals(2, GpuTuning.recommendThreadCount(32, 2, 100, true));
        assertEquals(4, GpuTuning.recommendThreadCount(32, 2, 100, false));
        assertEquals(1, GpuTuning.recommendThreadCount(32, 4, 1, false));
        assertEquals(3, GpuTuning.recommendThreadCount(4, 4, 100, false));
    }

    private static final class FakeDevice implements GpuDevice {
        private final String name;
        private final long globalMemory;
        private final long maxAllocation;

        private FakeDevice(String name, long globalMemory, long maxAllocation) {
            this.name = name;
            this.globalMemory = globalMemory;
            this.maxAllocation = maxAllocation;
        }

        @Override public String getBackendName() { return "fake"; }
        @Override public String getPlatformName() { return "fake"; }
        @Override public String getName() { return name; }
        @Override public String getVendor() { return "fake"; }
        @Override public String getDriverVersion() { return "1"; }
        @Override public String getComputeApiVersion() { return "1"; }
        @Override public boolean isAvailable() { return true; }
        @Override public boolean isCompilerAvailable() { return true; }
        @Override public boolean supportsDoublePrecision() { return true; }
        @Override public boolean hasUnifiedMemory() { return false; }
        @Override public long getGlobalMemoryBytes() { return globalMemory; }
        @Override public long getMaxAllocationBytes() { return maxAllocation; }
        @Override public long getMaxWorkGroupSize() { return 1024; }
        @Override public long[] getMaxWorkItemSizes() { return new long[] { 1024, 1024 }; }
        @Override public GpuContext openContext() { throw new UnsupportedOperationException(); }
    }
}
