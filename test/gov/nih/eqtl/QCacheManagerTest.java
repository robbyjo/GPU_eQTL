/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.FileTime;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.time.Instant;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QCacheManagerTest {
    @Test
    void inventoriesAndRecoverablyPrunesOnlyUnlockedKnownArtifacts(@TempDir Path root)
        throws Exception {
        Path stale = root.resolve("matrix.qcache");
        Path locked = root.resolve("raw.qraw");
        Path lock = root.resolve("raw.qraw.lock");
        Path unrelated = root.resolve("notes.txt");
        Files.writeString(stale, "prepared-cache");
        Files.writeString(locked, "raw-cache");
        Files.writeString(lock, "active");
        Files.writeString(unrelated, "keep me");
        FileTime old = FileTime.from(Instant.now().minusSeconds(10L * 24 * 60 * 60));
        Files.setLastModifiedTime(stale, old);
        Files.setLastModifiedTime(locked, old);
        Files.setLastModifiedTime(unrelated, old);

		Path dryReport = root.resolve("dry.tsv");
		Path applyReport = root.resolve("apply.tsv");
		try (FileChannel channel = FileChannel.open(lock, StandardOpenOption.WRITE);
			FileLock active = channel.lock()) {
			QCacheManager.Summary dry = QCacheManager.run(new QCacheManager.Options(
				root, dryReport, true, false, 5));
			assertEquals(2, dry.candidates());
			assertEquals(1, dry.locked());
			assertEquals(0, dry.moved());
			assertTrue(Files.exists(stale));
			assertTrue(Files.readString(dryReport).contains("would-move-to-trash"));

			QCacheManager.Summary applied = QCacheManager.run(new QCacheManager.Options(
				root, applyReport, true, true, 5));
			assertEquals(1, applied.moved());
		}
		assertFalse(Files.exists(stale));
		assertTrue(Files.exists(locked));
		assertTrue(Files.exists(unrelated));
        try (var trash = Files.walk(root.resolve(".trash"))) {
            assertTrue(trash.anyMatch(path -> path.getFileName().toString().equals("matrix.qcache")));
        }
        assertTrue(Files.readString(applyReport).contains("moved-to-trash"));
    }
}
