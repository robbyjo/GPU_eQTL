/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledResultCheckpointTest {
    @Test void partsResumeWithDigestAndRejectMutation(@TempDir Path root)throws Exception {
        Path dir=root.resolve("checkpoint");
        try(var parts=QPooledResultCheckpoint.open(dir,"signature",2,false)) {
            parts.writeBlock(0,w->w.write("a,1\n"));assertTrue(parts.isComplete(0));assertFalse(parts.isComplete(1));
        }
        try(var parts=QPooledResultCheckpoint.open(dir,"signature",2,true)) {
            assertTrue(parts.isComplete(0));parts.writeBlock(1,w->w.write("b,2\n"));
            parts.assemble(root.resolve("out.csv"),"ID,Value");assertEquals(3,Files.readAllLines(root.resolve("out.csv")).size());
            Files.writeString(dir.resolve("block-00000000.part"),"tampered\n");assertThrows(java.io.IOException.class,()->parts.isComplete(0));
        }
    }
    @Test void incompletePartsAreRetainedAndPairedCleanupIsSafe(@TempDir Path root)throws Exception {
        Path dir=root.resolve("checkpoint");
        try(var parts=QPooledResultCheckpoint.open(dir,"signature",1,false)) {
            Files.writeString(dir.resolve("block-00000000.part"),"uncommitted\n");assertFalse(parts.isComplete(0));
            parts.writeBlock(0,w->w.write("valid\n"));
            try(var retained=Files.list(root.resolve("checkpoint.uncommitted"))){assertEquals(1,retained.count());}
            parts.assemble(root.resolve("out.csv"),"header");parts.removeCompleted();assertFalse(Files.exists(dir));
        }
    }
    @Test void missingCommittedPartAndTruncatedDigestIndexFail(@TempDir Path root)throws Exception {
        Path dir=root.resolve("checkpoint");
        try(var parts=QPooledResultCheckpoint.open(dir,"signature",1,false)) {
            parts.writeBlock(0,w->w.write("valid\n"));Files.delete(dir.resolve("block-00000000.part"));
            assertThrows(java.io.IOException.class,()->parts.isComplete(0));
        }
        Files.write(dir.resolve("trait-digests.bin"),new byte[]{1,2});assertThrows(java.io.IOException.class,()->QPooledResultCheckpoint.open(dir,"signature",1,true));
    }
}
