/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QDataType;
import gov.nih.eqtl.QMissingValuePolicy;

class QMissingnessPolicyTest {
    @Test
    void scansPatternsAndAppliesStaticPolicies(@TempDir Path directory) throws Exception {
        QDelimitedMatrixSource source = source(directory);
        QMissingnessScan scan = QMissingnessScan.scan("predictor", source, null);

        assertEquals(1, scan.totalMissingValues());
        assertEquals(2, scan.patterns().size());
        assertEquals(1.0, scan.rowMean(1), 0);
        assertArrayEquals(new long[] {0, 1, 0, 0}, scan.sampleMissingValues());

        QMatrixRowSource mean = new QPolicyMatrixSource(source, scan, QMissingValuePolicy.MEAN);
        try (QMatrixRowSource.BlockReader reader = mean.open(null)) {
            assertArrayEquals(new double[] {0, 1, 1, 2}, reader.readBlock(3).values()[1], 0);
        }

        QMatrixRowSource zero = new QPolicyMatrixSource(source, scan, QMissingValuePolicy.ZERO);
        try (QMatrixRowSource.BlockReader reader = zero.open(null)) {
            assertArrayEquals(new double[] {0, 0, 1, 2}, reader.readBlock(3).values()[1], 0);
        }

        QMatrixRowSource excluded = new QPolicyMatrixSource(source, scan, QMissingValuePolicy.EXCLUDE_ROW);
        assertEquals(2, excluded.metadata().rowCount());
        try (QMatrixRowSource.BlockReader reader = excluded.open(null)) {
            assertArrayEquals(new String[] {"1:100:A:G", "1:120:G:A"}, reader.readBlock(3).rowIds());
        }

        QMatrixRowSource error = new QPolicyMatrixSource(source, scan, QMissingValuePolicy.ERROR);
        try (QMatrixRowSource.BlockReader reader = error.open(null)) {
            assertThrows(IllegalArgumentException.class, () -> reader.readBlock(3));
        }
    }

    @Test
    void localPatternUsesNearestFlankingGenotypeDonor(@TempDir Path directory) throws Exception {
        QDelimitedMatrixSource source = source(directory);
        QMissingnessScan scan = QMissingnessScan.scan("predictor", source, null);
        QMatrixRowSource imputed = new QLocalPatternImputedSource(source, scan, 1);

        try (QMatrixRowSource.BlockReader reader = imputed.open(null)) {
            QMatrixRowSource.Block block = reader.readBlock(3);
            assertArrayEquals(new double[] {0, 0, 1, 2}, block.values()[1], 0,
                "S2 matches S1 at both flanks, so the local donor dosage is 0 rather than row mean 1");
        }
    }

    @Test
    void writesAuditableMissingnessReport(@TempDir Path directory) throws Exception {
        QDelimitedMatrixSource predictor = source(directory);
        Path traitsFile = directory.resolve("traits.csv");
        Files.writeString(traitsFile, ",S1,S2,S3,S4\nt1,1,2,NA,4\n");
        QDelimitedMatrixSource traits = new QDelimitedMatrixSource(traitsFile, cCommonDelimiter, "#");
        Path report = directory.resolve("result.missingness.tsv");

        QMissingnessReport.write(report,
            QMissingnessScan.scan("predictor", predictor, null), QDataType.GENOTYPE,
            QMissingValuePolicy.LOCAL_PATTERN,
            QMissingnessScan.scan("trait", traits, null), QDataType.PROTEOMICS,
            QMissingValuePolicy.PATTERN, null, "complete-samples",
            new String[] {"S1", "S2", "S3", "S4"});

        String text = Files.readString(report);
        assertTrue(text.contains("SUMMARY\tpredictor"));
        assertTrue(text.contains("1:110:C:T"));
        assertTrue(text.contains("SUMMARY\ttrait"));
        assertTrue(text.contains("proteomics"));
        assertTrue(text.contains("local-pattern"));
    }

    private QDelimitedMatrixSource source(Path directory) throws Exception {
        Path file = directory.resolve("predictor.csv");
        Files.writeString(file, ",S1,S2,S3,S4\n"
            + "1:100:A:G,0,0,2,2\n"
            + "1:110:C:T,0,NA,1,2\n"
            + "1:120:G:A,0,0,1,2\n");
        return new QDelimitedMatrixSource(file, cCommonDelimiter, "#");
    }
}
