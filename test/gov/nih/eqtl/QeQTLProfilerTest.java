/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QeQTLProfilerTest {
    @Test
    void disabledProfilerDoesNotRecord() {
        QeQTLProfiler profiler = new QeQTLProfiler(false);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_COMPUTE, 10, 2, 3);
        assertEquals(0, profiler.measurement(QeQTLProfiler.Phase.GPU_COMPUTE).calls());
    }

    @Test
    void measurementsAreAggregatedAndWrittenAsCsv(@TempDir Path directory) throws Exception {
        QeQTLProfiler profiler = new QeQTLProfiler(true);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_UPLOAD, 2_000_000, 1, 1024);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_UPLOAD, 4_000_000, 2, 2048);

        QeQTLProfiler.Measurement measurement = profiler.measurement(QeQTLProfiler.Phase.GPU_UPLOAD);
        assertEquals(2, measurement.calls());
        assertEquals(6_000_000, measurement.nanoseconds());
        assertEquals(3, measurement.units());
        assertEquals(3072, measurement.bytes());

        Path output = directory.resolve("profile.csv");
        profiler.writeCsv(output);
        String csv = Files.readString(output);
        assertTrue(csv.startsWith("phase,calls,total_nanoseconds"));
        assertTrue(csv.contains("gpu_upload,2,6000000,0.006000000,3.000000,3,3072"));
    }
}
