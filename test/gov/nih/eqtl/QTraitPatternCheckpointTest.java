/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QTraitPatternCheckpointTest {
    @Test
    void commitsResumesAndAssemblesPatternGroupsInOrder(@TempDir Path directory)
        throws Exception {
        Path root = directory.resolve("checkpoint");
        QTraitPatternCheckpoint first = QTraitPatternCheckpoint.open(
            root, "signature", 2, false, true);
        Files.writeString(first.groupOutput(1), "header\ngroup1\n");
        first.commitResultFromOutput(1);
        first.commitQcLine(1, "qc1");
        first.commitEmptyResult(0);
        first.commitQcLine(0, "qc0");

        QTraitPatternCheckpoint resumed = QTraitPatternCheckpoint.open(
            root, "signature", 2, true, true);
        assertEquals(2, resumed.completedResultCount());
        assertTrue(resumed.isQcComplete(0));
        assertTrue(resumed.isQcComplete(1));
        Path results = directory.resolve("results.csv");
        Path qc = directory.resolve("qc.tsv");
        resumed.assembleResults(results, "header");
        resumed.assembleQc(qc, "qc-header");
        assertEquals("header" + System.lineSeparator() + "group1" + System.lineSeparator(),
            Files.readString(results));
        assertEquals("qc-header" + System.lineSeparator() + "qc0" + System.lineSeparator()
            + "qc1" + System.lineSeparator(), Files.readString(qc));
        assertThrows(IOException.class, () -> QTraitPatternCheckpoint.open(
            root, "different", 2, true, true));
    }

    @Test
    void successfulNonRetainedCheckpointCleansOnlyItsKnownDirectory(@TempDir Path directory)
        throws Exception {
        Path root = directory.resolve("checkpoint");
        QTraitPatternCheckpoint checkpoint = QTraitPatternCheckpoint.open(
            root, "signature", 1, false, false);
        checkpoint.commitEmptyResult(0);
        checkpoint.assembleResults(directory.resolve("result.csv"), "header");
        checkpoint.finishSuccess();
        assertFalse(Files.exists(root));
        assertTrue(Files.isRegularFile(directory.resolve("result.csv")));
    }
}
