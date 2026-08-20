/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QDelimitedMatrixSourceTest {
    @Test
    void readsMetadataAndReordersEachBlock() throws Exception {
        Path file = Path.of("test/resources/eqtl-reference/genotype.csv");
        QDelimitedMatrixSource source = new QDelimitedMatrixSource(file, cCommonDelimiter, "#");
        assertEquals(3, source.metadata().rowCount());
        assertEquals(8, source.metadata().columnCount());
        assertArrayEquals(new String[] { "G3", "G1", "G8", "G2", "G7", "G4", "G6", "G5" },
            source.metadata().sampleIds());

        int[] canonicalOrder = { 1, 3, 0, 5, 7, 6, 4, 2 };
        try (QDelimitedMatrixSource.BlockReader reader = source.open(canonicalOrder)) {
            QDelimitedMatrixSource.Block first = reader.readBlock(2);
            assertEquals(0, first.rowOffset());
            assertArrayEquals(new String[] { "rs1", "rs2" }, first.rowIds());
            assertArrayEquals(new double[] { 0.1, 0.4, 1.7, 1.4, 1.9, 0.8, 0.2, 1.1 },
                first.values()[0]);
            assertEquals(1, reader.readBlock(2).rowCount());
            assertEquals(null, reader.readBlock(2));
        }
    }

    @Test
    void rejectsDuplicateSampleAndRowIdentifiers(@TempDir Path temporaryDirectory) throws Exception {
        Path duplicateSamples = temporaryDirectory.resolve("duplicate-samples.csv");
        Files.writeString(duplicateSamples, ",S1,S1\nr1,0,1\n");
        IOException sampleError = assertThrows(IOException.class,
            () -> new QDelimitedMatrixSource(duplicateSamples, cCommonDelimiter, "#"));
        assertTrue(sampleError.getMessage().contains("Duplicate sample identifier"));

        Path duplicateRows = temporaryDirectory.resolve("duplicate-rows.csv");
        Files.writeString(duplicateRows, ",S1,S2\nr1,0,1\nr1,1,0\n");
        IOException rowError = assertThrows(IOException.class,
            () -> new QDelimitedMatrixSource(duplicateRows, cCommonDelimiter, "#"));
        assertTrue(rowError.getMessage().contains("Duplicate row identifier"));
    }
}
