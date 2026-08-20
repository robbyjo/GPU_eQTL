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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class GpuContextPoolTest {
    @Test
    void collectionConstructorTracksEveryContext() {
        FakeContext first = new FakeContext();
        FakeContext second = new FakeContext();
        GpuContextPool pool = new GpuContextPool(Arrays.asList(first, second));

        assertEquals(2, pool.getNumContexts());
        assertEquals(2, pool.getAllContexts().size());
        assertSame(first, pool.reserveContext());
        assertSame(second, pool.reserveContext());
    }

    @Test
    void releasedContextCanBeReservedAgain() {
        FakeContext context = new FakeContext();
        GpuContextPool pool = new GpuContextPool(new GpuContext[] { context });

        assertSame(context, pool.reserveContext());
        pool.releaseContext(context);
        pool.releaseContext(context);
        assertSame(context, pool.reserveContext());
    }

    @Test
    void closeClosesContextsAndRejectsReservations() {
        FakeContext context = new FakeContext();
        GpuContextPool pool = new GpuContextPool(new GpuContext[] { context });

        assertFalse(context.closed);
        pool.close();
        assertTrue(context.closed);
        assertThrows(GpuException.class, pool::reserveContext);
    }

    private static final class FakeContext implements GpuContext {
        private boolean closed;

        @Override
        public GpuDevice getDevice() { return null; }

        @Override
        public void compileKernel(String source, String kernelName) { }

        @Override
        public double[] executeDoubleKernel(double[] inputA, double[] inputB, int outputElements,
                long localMemoryBytes, int widthA, int widthB, long[] globalWorkSize, long[] localWorkSize) {
            return new double[outputElements];
        }

        @Override
        public void close() { closed = true; }
    }
}
