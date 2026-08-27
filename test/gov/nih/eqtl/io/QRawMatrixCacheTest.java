/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QRawMatrixCacheTest {
    @TempDir Path temporaryDirectory;

    @Test
    void preservesAlignedOrderSupportsSubsetsAndDetectsRowCorruption() throws Exception {
        Path csv = temporaryDirectory.resolve("raw.csv");
        Files.writeString(csv, ",S1,S2,S3\nv1,0,1,2\nv2,2,1,0\n");
        QDelimitedMatrixSource source = new QDelimitedMatrixSource(csv, cCommonDelimiter, "#");
        int[] aligned = {2, 0};
        String signature = QRawMatrixCache.signature(source, aligned);
        QRawMatrixCache cache = QRawMatrixCache.openOrBuild(temporaryDirectory.resolve("cache"),
            signature, source, aligned, 1, false);
        QRawMatrixCache reopened = QRawMatrixCache.openIfPresent(
            temporaryDirectory.resolve("cache"), signature, source);
        assertEquals(cache.path(), reopened.path());
        assertNull(QRawMatrixCache.openIfPresent(temporaryDirectory.resolve("absent"),
            signature, source));
        assertArrayEquals(new String[] {"S3", "S1"}, cache.metadata().sampleIds());
        try (QMatrixRowSource.BlockReader reader = cache.open(new int[] {1})) {
            QMatrixRowSource.Block block = reader.readBlock(10);
            assertArrayEquals(new String[] {"v1", "v2"}, block.rowIds());
            assertArrayEquals(new double[] {0}, block.values()[0], 0);
            assertArrayEquals(new double[] {2}, block.values()[1], 0);
            assertNull(reader.readBlock(1));
        }
        QMatrixRowSource.Block selected = cache.readSelected(Set.of("v2"),
            new int[] {1, 0});
        assertArrayEquals(new String[] {"v2"}, selected.rowIds());
        assertArrayEquals(new double[] {2, 0}, selected.values()[0], 0);
        assertEquals(1, cache.indexedReadStatistics().selectionCalls());
        assertEquals(1, cache.indexedReadStatistics().selectedRows());
        assertEquals(2, cache.indexedReadStatistics().indexedRows());
        assertEquals(16, cache.indexedReadStatistics().numericBytesRead());
        Path indexPath = cache.path().resolveSibling(
            cache.path().getFileName().toString() + ".idx");
        assertTrue(Files.isRegularFile(indexPath));
        QRawMatrixCache indexedReopen = QRawMatrixCache.open(cache.path(), signature);
        assertArrayEquals(new String[] {"v1"},
            indexedReopen.readSelected(Set.of("v1"), null).rowIds());
        assertTrue(indexedReopen.indexedReadStatistics().persistentIndexReused());

        try (RandomAccessFile file = new RandomAccessFile(indexPath.toFile(), "rw")) {
            file.seek(file.length() - 1);
            file.writeByte(file.readByte() ^ 1);
        }
        QRawMatrixCache rebuiltIndex = QRawMatrixCache.open(cache.path(), signature);
        assertArrayEquals(new String[] {"v2"},
            rebuiltIndex.readSelected(Set.of("v2"), null).rowIds());
        assertTrue(!rebuiltIndex.indexedReadStatistics().persistentIndexReused());

        try (RandomAccessFile file = new RandomAccessFile(cache.path().toFile(), "rw")) {
            long position = file.length() - 6;
            file.seek(position);
            file.writeByte(file.readByte() ^ 1);
        }
        QRawMatrixCache corrupt = QRawMatrixCache.open(cache.path(), signature);
        try (QMatrixRowSource.BlockReader reader = corrupt.open(null)) {
            assertThrows(IOException.class, () -> reader.readBlock(10));
        }
    }
}
