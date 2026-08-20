/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

/** Host-observed timings for the most recent exclusive GPU-context execution. */
public record GpuExecutionMetrics(long bufferSetupNanoseconds, long uploadNanoseconds,
    long computeNanoseconds, long downloadNanoseconds, long uploadedBytes, long downloadedBytes) {

    public static final GpuExecutionMetrics EMPTY = new GpuExecutionMetrics(0, 0, 0, 0, 0, 0);

    public GpuExecutionMetrics {
        if (bufferSetupNanoseconds < 0 || uploadNanoseconds < 0 || computeNanoseconds < 0
            || downloadNanoseconds < 0 || uploadedBytes < 0 || downloadedBytes < 0)
            throw new IllegalArgumentException("GPU execution metrics must not be negative");
    }
}
