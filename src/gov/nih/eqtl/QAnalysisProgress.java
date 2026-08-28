/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.util.Locale;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.function.Consumer;
import java.util.function.LongSupplier;

/** Thread-safe, rate-limited progress reporting for association comparisons. */
final class QAnalysisProgress implements AutoCloseable {
    static final long DEFAULT_REPORT_INTERVAL_NANOS = TimeUnit.SECONDS.toNanos(15);

    private final String label;
    private long totalComparisons;
    private final long initialComparisons;
    private final boolean dynamicTotal;
    private final LongSupplier clock;
    private final Consumer<String> output;
    private final long reportIntervalNanos;
    private final long started;
    private long lastReport;
    private long completedComparisons;
    private boolean finalReportWritten;
    private ScheduledExecutorService reporter;

    QAnalysisProgress(String label, long totalComparisons, long initialComparisons) {
        this(label, totalComparisons, initialComparisons, System::nanoTime,
            System.out::println, DEFAULT_REPORT_INTERVAL_NANOS, false);
        startReporter();
    }

    static QAnalysisProgress dynamic(String label, long completedComparisons) {
        QAnalysisProgress progress = new QAnalysisProgress(label, completedComparisons,
            completedComparisons, System::nanoTime, System.out::println,
            DEFAULT_REPORT_INTERVAL_NANOS, true);
        progress.startReporter();
        return progress;
    }

    private void startReporter() {
        reporter = Executors.newSingleThreadScheduledExecutor(action -> {
            Thread thread = new Thread(action, "gpu-eqtl-association-progress");
            thread.setDaemon(true);
            return thread;
        });
        reporter.scheduleAtFixedRate(this::heartbeat, reportIntervalNanos,
            reportIntervalNanos, TimeUnit.NANOSECONDS);
    }

    QAnalysisProgress(String label, long totalComparisons, long initialComparisons,
        LongSupplier clock, Consumer<String> output, long reportIntervalNanos) {
        this(label, totalComparisons, initialComparisons, clock, output,
            reportIntervalNanos, false);
    }

    QAnalysisProgress(String label, long totalComparisons, long initialComparisons,
        LongSupplier clock, Consumer<String> output, long reportIntervalNanos,
        boolean dynamicTotal) {
        if (label == null || label.isBlank())
            throw new IllegalArgumentException("Progress label is required");
        if (totalComparisons < 0 || (!dynamicTotal && totalComparisons == 0))
            throw new IllegalArgumentException("Total association comparisons must be positive");
        if (initialComparisons < 0 || initialComparisons > totalComparisons)
            throw new IllegalArgumentException("Initial association progress is out of range");
        if (clock == null || output == null || reportIntervalNanos < 0)
            throw new IllegalArgumentException("Progress clock, output, and interval must be valid");
        this.label = label;
        this.totalComparisons = totalComparisons;
        this.initialComparisons = initialComparisons;
        this.completedComparisons = initialComparisons;
        this.dynamicTotal = dynamicTotal;
        this.clock = clock;
        this.output = output;
        this.reportIntervalNanos = reportIntervalNanos;
        started = clock.getAsLong();
        lastReport = started;
        output.accept(dynamicTotal
            ? dynamicStartMessage(label, initialComparisons)
            : startMessage(label, totalComparisons, initialComparisons));
    }

    synchronized void registerComparisons(long count) {
        if (!dynamicTotal)
            throw new IllegalStateException("Cannot register work on fixed-total progress");
        if (count < 0)
            throw new IllegalArgumentException("Registered association comparisons cannot be negative");
        totalComparisons = Math.addExact(totalComparisons, count);
    }

    synchronized void addComparisons(long count) {
        if (count <= 0)
            throw new IllegalArgumentException("Completed association comparisons must be positive");
        if (completedComparisons > totalComparisons - count)
            throw new IllegalStateException("Association progress exceeded its declared total");
        completedComparisons += count;
        long now = clock.getAsLong();
        boolean complete = !dynamicTotal && completedComparisons == totalComparisons;
        if (complete || now - lastReport >= reportIntervalNanos)
            report(now, complete);
    }

    synchronized void heartbeat() {
        if (!finalReportWritten) {
            long now = clock.getAsLong();
            if (now - lastReport >= reportIntervalNanos)
                report(now, false);
        }
    }

    synchronized void complete() {
        if (completedComparisons != totalComparisons)
            throw new IllegalStateException("Association completed after " + completedComparisons
                + " of " + totalComparisons + " declared comparisons");
        if (!finalReportWritten)
            report(clock.getAsLong(), true);
    }

    @Override
    public synchronized void close() {
        stopReporter();
    }

    private void report(long now, boolean complete) {
        output.accept(dynamicTotal && !complete
            ? dynamicProgressMessage(label, completedComparisons, totalComparisons,
                completedComparisons - initialComparisons, Math.max(0, now - started))
            : progressMessage(label, completedComparisons, totalComparisons,
                completedComparisons - initialComparisons, Math.max(0, now - started), complete));
        lastReport = now;
        finalReportWritten |= complete;
        if (complete)
            stopReporter();
    }

    private void stopReporter() {
        if (reporter != null) {
            reporter.shutdownNow();
            reporter = null;
        }
    }

    static String startMessage(String label, long totalComparisons, long initialComparisons) {
        String message = label + " started: " + String.format(Locale.ROOT, "%,d", totalComparisons)
            + " variant-trait comparison(s)";
        if (initialComparisons > 0)
            message += "; resumed " + String.format(Locale.ROOT, "%,d", initialComparisons)
                + " completed comparison(s)";
        return message + "; progress approximately every 15 seconds.";
    }

    static String dynamicStartMessage(String label, long initialComparisons) {
        String message = label + " started; exact active comparison total will be finalized "
            + "from pattern-specific variant QC";
        if (initialComparisons > 0)
            message += "; resumed " + String.format(Locale.ROOT, "%,d", initialComparisons)
                + " completed active comparison(s)";
        return message + "; progress approximately every 15 seconds.";
    }

    static String dynamicProgressMessage(String label, long completedComparisons,
        long discoveredComparisons, long completedThisRun, long elapsedNanos) {
        double elapsedSeconds = elapsedNanos / 1_000_000_000.0;
        double rate = elapsedSeconds > 0 ? completedThisRun / elapsedSeconds : 0;
        StringBuilder message = new StringBuilder(label).append(" progress: ")
            .append(String.format(Locale.ROOT, "%,d completed / %,d discovered active comparisons",
                completedComparisons, discoveredComparisons))
            .append("; total still being determined; elapsed=").append(duration(elapsedSeconds));
        if (rate > 0)
            message.append("; rate=").append(String.format(Locale.ROOT, "%,.0f", rate)).append("/s");
        return message.toString();
    }

    static String progressMessage(String label, long completedComparisons, long totalComparisons,
        long completedThisRun, long elapsedNanos, boolean complete) {
        double percent = totalComparisons == 0 && complete ? 100.0
            : 100.0 * completedComparisons / totalComparisons;
        double elapsedSeconds = elapsedNanos / 1_000_000_000.0;
        double rate = elapsedSeconds > 0 ? completedThisRun / elapsedSeconds : 0;
        StringBuilder message = new StringBuilder(label)
            .append(complete ? " complete: " : " progress: ")
            .append(String.format(Locale.ROOT, "%,d/%,d (%.1f%%)",
                completedComparisons, totalComparisons, percent))
            .append(" comparisons; elapsed=").append(duration(elapsedSeconds));
        if (rate > 0) {
            message.append("; rate=").append(String.format(Locale.ROOT, "%,.0f", rate)).append("/s");
            if (!complete)
                message.append("; ETA=").append(duration((totalComparisons - completedComparisons) / rate));
        }
        return message.toString();
    }

    private static String duration(double seconds) {
        if (!Double.isFinite(seconds) || seconds < 0)
            return "unknown";
        long rounded = Math.round(seconds);
        long hours = rounded / 3600;
        long minutes = (rounded % 3600) / 60;
        long remainingSeconds = rounded % 60;
        if (hours > 0)
            return String.format(Locale.ROOT, "%dh%02dm%02ds", hours, minutes, remainingSeconds);
        if (minutes > 0)
            return String.format(Locale.ROOT, "%dm%02ds", minutes, remainingSeconds);
        return remainingSeconds + "s";
    }
}
