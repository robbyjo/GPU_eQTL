/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.nio.file.Path;
import java.nio.file.Files;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.jama.QRDecomposition;

class QCovariateTableTest {
    @Test
    void bridgesDifferentIdsReordersAndAutomaticallyEncodesTextFactors() throws Exception {
        QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/genotype.csv"), cCommonDelimiter, "#");
        QDelimitedMatrixSource expression = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/expression.csv"), cCommonDelimiter, "#");
        QCovariateTable covariates = QCovariateTable.load(
            Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");

        QSampleAlignment alignment = covariates.align(genotype.metadata().sampleIds(),
            expression.metadata().sampleIds(), null, null);
        assertEquals("genotype_id", alignment.genotypeIdColumn());
        assertEquals("expression_id", alignment.expressionIdColumn());
        assertArrayEquals(new int[] { 1, 3, 0, 5, 7, 6, 4, 2 }, alignment.genotypeColumnOrder());
        assertArrayEquals(new int[] { 2, 0, 6, 4, 3, 7, 5, 1 }, alignment.expressionColumnOrder());
        assertEquals(8, alignment.genotypeReorderedCount());
        assertEquals(8, alignment.expressionReorderedCount());

        QCovariateTable.ModelMatrix model = covariates.buildModelMatrix(
            new String[] { "Age", "Batch" }, null);
        assertArrayEquals(new String[] { "(Intercept)", "Age", "Batch[B]" }, model.columnNames());
        assertArrayEquals(new String[] { "Batch" }, model.automaticFactors());
        assertArrayEquals(new double[] { 1, 20, 0 }, model.values()[0]);
        assertArrayEquals(new double[] { 1, 25, 1 }, model.values()[2]);
        QRDecomposition qr = new QRDecomposition(model.values());
        assertEquals(3, qr.getRank());
        assertEquals(3, qr.getNumColumns());
    }

    @Test
    void failsWhenAnExplicitBridgeColumnDoesNotMatch() throws Exception {
        QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/genotype.csv"), cCommonDelimiter, "#");
        QDelimitedMatrixSource expression = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/expression.csv"), cCommonDelimiter, "#");
        QCovariateTable covariates = QCovariateTable.load(
            Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");
        assertThrows(IllegalArgumentException.class, () -> covariates.align(
            genotype.metadata().sampleIds(), expression.metadata().sampleIds(), "subject_id", "expression_id"));
    }

	@Test
	void inspectsSelectedCovariatesAndBuildsOnlyCompleteSamples(@TempDir Path directory) throws Exception {
		Path file = directory.resolve("covariates.csv");
		Files.writeString(file, "id,age,batch,unused\nS1,40,A,NA\nS2,NA,B,x\nS3,60,A,y\n");
		QCovariateTable table = QCovariateTable.load(file, cCommonDelimiter, "#");
		QCovariateTable.Missingness missing = table.inspectMissingness(new String[] {"age", "batch"});
		assertEquals(1, missing.totalMissing());
		assertArrayEquals(new boolean[] {true, false, true}, missing.completeRows());
		assertEquals(2, missing.completeSampleCount());
		QCovariateTable.ModelMatrix model = table.buildModelMatrix(
			new String[] {"age"}, null, missing.completeRows());
		assertEquals(2, model.values().length);
		assertArrayEquals(new double[] {1, 40}, model.values()[0]);
		assertArrayEquals(new double[] {1, 60}, model.values()[1]);
	}
}
