/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.FileTime;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QeQTLPreprocessor;
import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.jama.QRDecomposition;

class QBinaryMatrixCacheTest {
    @Test
    void indexedPreparedRowsRoundTripExactlyAndCacheIsReused(@TempDir Path directory) throws Exception {
        QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/genotype.csv"), cCommonDelimiter, "#");
        QDelimitedMatrixSource expression = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/expression.csv"), cCommonDelimiter, "#");
        QCovariateTable covariates = QCovariateTable.load(
            Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");
        QSampleAlignment alignment = covariates.align(genotype.metadata().sampleIds(),
            expression.metadata().sampleIds(), "genotype_id", "expression_id");
        double[][] covariateQ = new QRDecomposition(covariates.buildModelMatrix(
            new String[] { "Age", "Batch" }, null).values()).getQ().getArray();
        String signature = QBinaryMatrixCache.signature("Genotype", genotype,
            alignment.genotypeColumnOrder(), covariateQ);

        Path cachePath;
        PreparedBlock cached;
        try (QBinaryMatrixCache cache = QBinaryMatrixCache.openOrBuild(directory, "Genotype",
            signature, genotype, alignment.genotypeColumnOrder(), covariateQ, 2, false)) {
            cachePath = cache.path();
            assertEquals(3, cache.rowCount());
            assertEquals(8, cache.sampleCount());
            cached = cache.readBlock(1, 2);
        }
        assertArrayEquals(new String[] { "rs2", "rs3" }, cached.rowIds());

        PreparedBlock direct;
        try (QDelimitedMatrixSource.BlockReader reader = genotype.open(alignment.genotypeColumnOrder())) {
            QDelimitedMatrixSource.Block all = reader.readBlock(3);
            QDelimitedMatrixSource.Block tail = new QDelimitedMatrixSource.Block(1,
                new String[] { all.rowIds()[1], all.rowIds()[2] },
                new double[][] { all.values()[1], all.values()[2] });
            direct = QeQTLPreprocessor.prepare(tail, covariateQ, "Genotype");
        }
        for (int row = 0; row < 2; row++) {
            assertArrayEquals(direct.values()[row], cached.values()[row], 0.0);
            assertEquals(direct.standardDeviations()[row], cached.standardDeviations()[row], 0.0);
        }

        FileTime modified = Files.getLastModifiedTime(cachePath);
        try (QBinaryMatrixCache reused = QBinaryMatrixCache.openOrBuild(directory, "Genotype",
            signature, genotype, alignment.genotypeColumnOrder(), covariateQ, 2, false)) {
            assertEquals(cachePath, reused.path());
        }
        assertEquals(modified, Files.getLastModifiedTime(cachePath));

        double[][] changedQ = new double[covariateQ.length][];
        for (int row = 0; row < covariateQ.length; row++)
            changedQ[row] = covariateQ[row].clone();
        changedQ[0][0] = Math.nextUp(changedQ[0][0]);
        assertNotEquals(signature, QBinaryMatrixCache.signature("Genotype", genotype,
            alignment.genotypeColumnOrder(), changedQ));
    }
}
