/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.EnumMap;
import java.util.Locale;
import java.util.concurrent.atomic.LongAdder;

/** Low-overhead, opt-in phase measurements for an eQTL analysis run. */
public final class QeQTLProfiler {
    public enum Phase {
        METADATA_AND_ALIGNMENT("metadata_alignment"),
        VARIANT_QC("variant_qc"),
        CACHE_SIGNATURES("cache_signatures"),
        GENOTYPE_CACHE_OPEN_OR_BUILD("genotype_cache_open_or_build"),
        EXPRESSION_CACHE_OPEN_OR_BUILD("expression_cache_open_or_build"),
        TRAIT_CACHE_MEMORY_LOAD("trait_cache_memory_load"),
        GENOTYPE_CACHE_READ("genotype_cache_read"),
        EXPRESSION_CACHE_READ("expression_cache_read"),
        GENOTYPE_PACK("genotype_pack"),
        EXPRESSION_PACK("expression_pack"),
		GPU_RESIDUALIZATION_UPLOAD("gpu_residualization_upload"),
		GPU_RESIDUALIZATION_BUFFER_SETUP("gpu_residualization_buffer_setup"),
		GPU_RESIDUALIZATION_COMPUTE("gpu_residualization_compute"),
		GPU_RESIDUALIZATION_DOWNLOAD("gpu_residualization_download"),
        GPU_CONTEXT_WAIT("gpu_context_wait"),
        GPU_BUFFER_SETUP("gpu_buffer_setup"),
        GPU_UPLOAD("gpu_upload"),
        GPU_COMPUTE("gpu_compute"),
        GPU_DOWNLOAD("gpu_download"),
        CPU_RESULTS_AND_WRITE("cpu_results_and_write"),
        KERNEL_COMPILE("kernel_compile"),
        OUTPUT_ASSEMBLY("output_assembly"),
        ANALYSIS_WALL("analysis_wall");

        private final String label;

        Phase(String label) {
            this.label = label;
        }

        public String label() {
            return label;
        }
    }

    public record Measurement(long calls, long nanoseconds, long units, long bytes) { }

    private static final class Metric {
        private final LongAdder calls = new LongAdder();
        private final LongAdder nanoseconds = new LongAdder();
        private final LongAdder units = new LongAdder();
        private final LongAdder bytes = new LongAdder();

        private Measurement snapshot() {
            return new Measurement(calls.sum(), nanoseconds.sum(), units.sum(), bytes.sum());
        }
    }

    private final boolean enabled;
    private final EnumMap<Phase, Metric> metrics = new EnumMap<>(Phase.class);

    public QeQTLProfiler(boolean enabled) {
        this.enabled = enabled;
        for (Phase phase : Phase.values())
            metrics.put(phase, new Metric());
    }

    public boolean isEnabled() {
        return enabled;
    }

    public long start() {
        return enabled ? System.nanoTime() : 0L;
    }

    public void record(Phase phase, long startedAtNanos) {
        record(phase, startedAtNanos, 0, 0);
    }

    public void record(Phase phase, long startedAtNanos, long units, long bytes) {
        if (!enabled)
            return;
        recordElapsed(phase, System.nanoTime() - startedAtNanos, units, bytes);
    }

    public void recordElapsed(Phase phase, long nanoseconds, long units, long bytes) {
        if (!enabled)
            return;
        if (nanoseconds < 0 || units < 0 || bytes < 0)
            throw new IllegalArgumentException("Profile measurements must not be negative");
        Metric metric = metrics.get(phase);
        metric.calls.increment();
        metric.nanoseconds.add(nanoseconds);
        metric.units.add(units);
        metric.bytes.add(bytes);
    }

    public Measurement measurement(Phase phase) {
        return metrics.get(phase).snapshot();
    }

    public void printSummary(PrintStream output) {
        if (!enabled)
            return;
        output.println("Profiling summary (phase totals may overlap when workers run concurrently):");
        output.printf(Locale.ROOT, "%-34s %8s %13s %13s %14s %12s%n",
            "Phase", "Calls", "Total (s)", "Average (ms)", "Units", "MiB");
        for (Phase phase : Phase.values()) {
            Measurement value = measurement(phase);
            if (value.calls() == 0)
                continue;
            double seconds = value.nanoseconds() / 1_000_000_000.0;
            double averageMillis = value.nanoseconds() / 1_000_000.0 / value.calls();
            double mebibytes = value.bytes() / (1024.0 * 1024.0);
            output.printf(Locale.ROOT, "%-34s %8d %13.6f %13.3f %14d %12.2f%n",
                phase.label(), value.calls(), seconds, averageMillis, value.units(), mebibytes);
        }
    }

    public void writeCsv(Path path) throws IOException {
        if (!enabled || path == null)
            return;
        Path normalized = path.toAbsolutePath().normalize();
        Path parent = normalized.getParent();
        if (parent != null)
            Files.createDirectories(parent);
        try (BufferedWriter writer = Files.newBufferedWriter(normalized, StandardCharsets.UTF_8)) {
            writer.write("phase,calls,total_nanoseconds,total_seconds,average_milliseconds,units,bytes");
            writer.newLine();
            for (Phase phase : Phase.values()) {
                Measurement value = measurement(phase);
                writer.write(String.format(Locale.ROOT, "%s,%d,%d,%.9f,%.6f,%d,%d",
                    phase.label(), value.calls(), value.nanoseconds(),
                    value.nanoseconds() / 1_000_000_000.0,
                    value.calls() == 0 ? 0.0 : value.nanoseconds() / 1_000_000.0 / value.calls(),
                    value.units(), value.bytes()));
                writer.newLine();
            }
        }
    }
}
