/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import org.junit.jupiter.api.Test;

class QAnalysisProgressTest {
    @Test
    void rateLimitsIntermediateUpdatesAndAlwaysReportsCompletion() {
        AtomicLong clock = new AtomicLong();
        List<String> messages = new ArrayList<>();
        QAnalysisProgress progress = new QAnalysisProgress("Association", 100, 0,
            clock::get, messages::add, TimeUnit.SECONDS.toNanos(10));

        clock.set(TimeUnit.SECONDS.toNanos(1));
        progress.addComparisons(25);
        assertEquals(1, messages.size());

        clock.set(TimeUnit.SECONDS.toNanos(10));
        progress.addComparisons(25);
        assertEquals(2, messages.size());
        assertTrue(messages.get(1).contains("50/100 (50.0%)"));
        assertTrue(messages.get(1).contains("rate=5/s"));
        assertTrue(messages.get(1).contains("ETA=10s"));

        clock.set(TimeUnit.SECONDS.toNanos(20));
        progress.addComparisons(50);
        progress.complete();
        assertEquals(3, messages.size());
        assertTrue(messages.get(2).contains("Association complete: 100/100 (100.0%)"));
    }

    @Test
    void resumedWorkIsShownButExcludedFromCurrentRunRate() {
        String start = QAnalysisProgress.startMessage("Association", 1_000, 400);
        assertTrue(start.contains("resumed 400"));
        String update = QAnalysisProgress.progressMessage("Association", 500, 1_000,
            100, TimeUnit.SECONDS.toNanos(10), false);
        assertTrue(update.contains("rate=10/s"));
        assertTrue(update.contains("ETA=50s"));
    }

    @Test
    void heartbeatReportsEvenWhenNoTileHasFinished() {
        AtomicLong clock = new AtomicLong();
        List<String> messages = new ArrayList<>();
        QAnalysisProgress progress = new QAnalysisProgress("Association", 100, 0,
            clock::get, messages::add, TimeUnit.SECONDS.toNanos(15));
        clock.set(TimeUnit.SECONDS.toNanos(15));
        progress.heartbeat();
        progress.close();
        assertEquals(2, messages.size());
        assertTrue(messages.get(1).contains("0/100 (0.0%)"));
        assertTrue(messages.get(1).contains("elapsed=15s"));
    }

    @Test
    void dynamicProgressCountsOnlyDiscoveredActiveComparisons() {
        AtomicLong clock = new AtomicLong();
        List<String> messages = new ArrayList<>();
        QAnalysisProgress progress = new QAnalysisProgress("Pattern association", 20, 20,
            clock::get, messages::add, TimeUnit.SECONDS.toNanos(10), true);
        assertTrue(messages.get(0).contains("exact active comparison total"));
        assertTrue(messages.get(0).contains("resumed 20"));

        progress.registerComparisons(12);
        clock.set(TimeUnit.SECONDS.toNanos(10));
        progress.addComparisons(5);
        assertTrue(messages.get(1).contains("25 completed / 32 discovered active comparisons"));
        assertTrue(messages.get(1).contains("total still being determined"));
        progress.addComparisons(7);
        progress.complete();
        assertTrue(messages.get(2).contains("32/32 (100.0%)"));
    }
}
