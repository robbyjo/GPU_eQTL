/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
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

class QAnalysisCheckpointTest {
    @Test
    void completedPartsAreAssembledInBlockOrderAndCleaned(@TempDir Path directory) throws Exception {
        Path checkpointDirectory = directory.resolve("checkpoint");
        QAnalysisCheckpoint checkpoint = QAnalysisCheckpoint.open(checkpointDirectory,
            "signature", 3, false, false);
        checkpoint.writeBlock(1, writer -> writer.write("block1\n"));
        checkpoint.writeBlock(0, writer -> writer.write("block0\n"));
        checkpoint.writeBlock(2, writer -> { });
        assertEquals(3, checkpoint.completedCount());

        Path output = directory.resolve("result.csv");
        checkpoint.assemble(output, "header");
        assertEquals("header" + System.lineSeparator() + "block0\nblock1\n", Files.readString(output));
        assertFalse(Files.exists(checkpointDirectory));
    }

    @Test
    void matchingCheckpointResumesAndMismatchedCheckpointStops(@TempDir Path directory) throws Exception {
        Path checkpointDirectory = directory.resolve("checkpoint");
        QAnalysisCheckpoint first = QAnalysisCheckpoint.open(checkpointDirectory,
            "signature", 2, false, true);
        first.writeBlock(0, writer -> writer.write("first\n"));

        QAnalysisCheckpoint resumed = QAnalysisCheckpoint.open(checkpointDirectory,
            "signature", 2, true, true);
        assertTrue(resumed.isComplete(0));
        assertFalse(resumed.isComplete(1));
        assertEquals(1, resumed.completedCount());
        assertThrows(IOException.class, () -> QAnalysisCheckpoint.open(checkpointDirectory,
            "different", 2, true, true));
        assertThrows(IOException.class, () -> QAnalysisCheckpoint.open(checkpointDirectory,
            "signature", 2, false, true));
    }
}
