/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.FileTime;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QMissingnessScanCacheTest {
    @Test
    void exactScanIsReusedAndCorruptionIsRebuilt(@TempDir Path directory) throws Exception {
        Path matrix = directory.resolve("traits.csv");
        Files.writeString(matrix, "id,S1,S2,S3\n"
            + "t1,1,NA,3\n"
            + "t2,2,,4\n"
            + "t3,5,6,7\n"
            + "t4,.,8,9\n");
        QDelimitedMatrixSource source = new QDelimitedMatrixSource(matrix, new char[] {','}, "#");
        Path cacheDirectory = directory.resolve("cache");
        int[] order = {2, 0, 1};

        QMissingnessScan first = QMissingnessScan.scanOrLoad("trait", source, order,
            cacheDirectory, false);
        assertEquals(3, first.totalMissingValues());
        assertEquals(3, first.patterns().size());
        assertArrayEquals(new String[] {"S3", "S1", "S2"}, first.sampleIds());
        Path cache = Files.list(cacheDirectory).findFirst().orElseThrow();
        FileTime marker = FileTime.fromMillis(1_000_000);
        Files.setLastModifiedTime(cache, marker);

        QMissingnessScan reused = QMissingnessScan.scanOrLoad("trait", source, order,
            cacheDirectory, false);
        assertEquals(marker, Files.getLastModifiedTime(cache));
        assertEquivalent(first, reused);

        byte[] corrupted = Files.readAllBytes(cache);
        corrupted[corrupted.length / 2] ^= 0x40;
        Files.write(cache, corrupted);
        QMissingnessScan rebuilt = QMissingnessScan.scanOrLoad("trait", source, order,
            cacheDirectory, false);
        assertTrue(Files.getLastModifiedTime(cache).toMillis() > marker.toMillis());
        assertEquivalent(first, rebuilt);
    }

    @Test
    void sourceMetadataAndSelectedOrderChangeTheSignature(@TempDir Path directory)
        throws Exception {
        Path matrix = directory.resolve("traits.csv");
        Files.writeString(matrix, "id,S1,S2\nt1,1,NA\n");
        QDelimitedMatrixSource first = new QDelimitedMatrixSource(matrix, new char[] {','}, "#");
        String original = QMissingnessScan.cacheSignature(first, new int[] {0, 1});
        String reordered = QMissingnessScan.cacheSignature(first, new int[] {1, 0});
        assertTrue(!original.equals(reordered));

        Files.writeString(matrix, "id,S1,S2\nt1,1,NA\nt2,2,3\n");
        QDelimitedMatrixSource changed = new QDelimitedMatrixSource(matrix, new char[] {','}, "#");
        String updated = QMissingnessScan.cacheSignature(changed, new int[] {0, 1});
        assertTrue(!original.equals(updated));
    }

    private static void assertEquivalent(QMissingnessScan expected,
        QMissingnessScan actual) {
        assertEquals(expected.rowCount(), actual.rowCount());
        assertEquals(expected.totalMissingValues(), actual.totalMissingValues());
        assertArrayEquals(expected.sampleIds(), actual.sampleIds());
        assertArrayEquals(expected.sampleMissingValues(), actual.sampleMissingValues());
        assertEquals(expected.missingRows(), actual.missingRows());
        assertEquals(expected.patterns().size(), actual.patterns().size());
        for (int i = 0; i < expected.patterns().size(); i++) {
            assertEquals(expected.patterns().get(i).id(), actual.patterns().get(i).id());
            assertEquals(expected.patterns().get(i).missingSamples(),
                actual.patterns().get(i).missingSamples());
            assertArrayEquals(expected.patterns().get(i).rowIndices(),
                actual.patterns().get(i).rowIndices());
        }
        for (long row = 0; row < expected.rowCount(); row++)
            assertEquals(expected.rowMean(row), actual.rowMean(row));
    }
}
