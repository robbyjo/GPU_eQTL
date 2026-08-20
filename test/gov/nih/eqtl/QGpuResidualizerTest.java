/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.jupiter.api.Test;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuPrecision;

class QGpuResidualizerTest {
    @Test
    void concurrentBlocksUseExclusiveBorrowedContextsAndReleaseProjectionBuffers() throws Exception {
        CyclicBarrier barrier = new CyclicBarrier(2);
        FakeContext first = new FakeContext("first", barrier);
        FakeContext second = new FakeContext("second", barrier);
        double[][] q = {
            { 0.5 }, { 0.5 }, { 0.5 }, { 0.5 }
        };
        QGpuResidualizer residualizer = new QGpuResidualizer(
            new GpuContext[] { first, second }, q, GpuPrecision.FP64, new QeQTLProfiler(false));
        ExecutorService workers = Executors.newFixedThreadPool(2);
        try {
            Future<double[][]> a = workers.submit(() -> residualizer.residualize(
                new double[][] { { 1, 2, 3, 4 } }, q, "A"));
            Future<double[][]> b = workers.submit(() -> residualizer.residualize(
                new double[][] { { 4, 6, 8, 10 } }, q, "B"));
            assertArrayEquals(new double[] { -1.5, -0.5, 0.5, 1.5 }, a.get()[0], 0.0);
            assertArrayEquals(new double[] { -3, -1, 1, 3 }, b.get()[0], 0.0);
        } finally {
            workers.shutdownNow();
            residualizer.close();
        }
        assertEquals(1, first.calls.get());
        assertEquals(1, second.calls.get());
        assertEquals(1, first.releases.get());
        assertEquals(1, second.releases.get());
        assertEquals(0, first.closes.get(), "Borrowed contexts remain open for the eQTL scheduler");
        assertEquals(0, second.closes.get(), "Borrowed contexts remain open for the eQTL scheduler");
    }

    private static final class FakeContext implements GpuContext {
        private final FakeDevice device;
        private final CyclicBarrier barrier;
        private final AtomicInteger calls = new AtomicInteger();
        private final AtomicInteger releases = new AtomicInteger();
        private final AtomicInteger closes = new AtomicInteger();

        private FakeContext(String name, CyclicBarrier barrier) {
            device = new FakeDevice(name, this);
            this.barrier = barrier;
        }

        @Override public GpuDevice getDevice() { return device; }
        @Override public void compileKernel(String source, String kernelName) { }
        @Override public double[] executeDoubleKernel(double[] a, double[] b, int elements,
            long localBytes, int widthA, int widthB, long[] global, long[] local) {
            throw new UnsupportedOperationException();
        }
        @Override public double[] residualizeDoubleRows(double[] values, double[] q,
            int rows, int samples, int rank) {
            calls.incrementAndGet();
            try {
                barrier.await();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
            double[] result = values.clone();
            for (int row = 0; row < rows; row++) {
                double coefficient = 0;
                for (int sample = 0; sample < samples; sample++)
                    coefficient += values[row * samples + sample] * q[sample * rank];
                for (int sample = 0; sample < samples; sample++)
                    result[row * samples + sample] -= coefficient * q[sample * rank];
            }
            return result;
        }
        @Override public void releaseResidualizationResources() { releases.incrementAndGet(); }
        @Override public void close() { closes.incrementAndGet(); }
    }

    private record FakeDevice(String name, GpuContext context) implements GpuDevice {
        @Override public String getBackendName() { return "fake"; }
        @Override public String getPlatformName() { return "test"; }
        @Override public String getName() { return name; }
        @Override public String getVendor() { return "test"; }
        @Override public String getDriverVersion() { return "test"; }
        @Override public String getComputeApiVersion() { return "test"; }
        @Override public boolean isAvailable() { return true; }
        @Override public boolean isCompilerAvailable() { return true; }
        @Override public boolean supportsDoublePrecision() { return true; }
        @Override public boolean hasUnifiedMemory() { return false; }
        @Override public long getGlobalMemoryBytes() { return 1L << 30; }
        @Override public long getMaxAllocationBytes() { return 1L << 29; }
        @Override public long getMaxWorkGroupSize() { return 256; }
        @Override public long[] getMaxWorkItemSizes() { return new long[] { 256, 256, 64 }; }
        @Override public GpuContext openContext() { return context; }
    }
}
