/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;

class GpuRuntimeTest {
    @Test
    void discoveryFiltersUnavailableCompilerlessAndNonFp64Devices() {
        FakeDevice usable = new FakeDevice(true, true, true);
        GpuRuntime runtime = new GpuRuntime(new FakeBackend(Arrays.asList(
            usable,
            new FakeDevice(false, true, true),
            new FakeDevice(true, false, true),
            new FakeDevice(true, true, false))));

        List<GpuDevice> devices = runtime.getGpuDevices(true, true);

        assertEquals(1, devices.size());
        assertSame(usable, devices.get(0));
    }

    private static final class FakeBackend implements GpuBackend {
        private final List<GpuDevice> devices;

        private FakeBackend(List<GpuDevice> devices) {
            this.devices = devices;
        }

        @Override
        public String getName() { return "fake"; }

        @Override
        public String getRuntimeDescription() { return "test runtime"; }

        @Override
        public List<GpuDevice> discoverGpuDevices() { return devices; }
    }

    private static final class FakeDevice implements GpuDevice {
        private final boolean available;
        private final boolean compilerAvailable;
        private final boolean doublePrecision;

        private FakeDevice(boolean available, boolean compilerAvailable, boolean doublePrecision) {
            this.available = available;
            this.compilerAvailable = compilerAvailable;
            this.doublePrecision = doublePrecision;
        }

        @Override
        public String getBackendName() { return "fake"; }

        @Override
        public String getPlatformName() { return "test"; }

        @Override
        public String getName() { return "test device"; }

        @Override
        public String getVendor() { return "test vendor"; }

        @Override
        public String getDriverVersion() { return "1"; }

        @Override
        public String getComputeApiVersion() { return "1"; }

        @Override
        public boolean isAvailable() { return available; }

        @Override
        public boolean isCompilerAvailable() { return compilerAvailable; }

        @Override
        public boolean supportsDoublePrecision() { return doublePrecision; }

        @Override
        public boolean hasUnifiedMemory() { return false; }

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
