/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.FileTime;
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QeQTLPreprocessor;
import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.jama.QRDecomposition;
import gov.nih.gpu.GpuPrecision;

class QBinaryMatrixCacheTest {
	@Test
	void concurrentPreparationPreservesCacheRowOrder(@TempDir Path directory) throws Exception {
		QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
			Path.of("test/resources/eqtl-reference/genotype.csv"), cCommonDelimiter, "#");
		QDelimitedMatrixSource expression = new QDelimitedMatrixSource(
			Path.of("test/resources/eqtl-reference/expression.csv"), cCommonDelimiter, "#");
		QCovariateTable covariates = QCovariateTable.load(
			Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");
		QSampleAlignment alignment = covariates.align(genotype.metadata().sampleIds(),
			expression.metadata().sampleIds(), "genotype_id", "expression_id");
		double[][] q = new QRDecomposition(covariates.buildModelMatrix(
			new String[] { "Age", "Batch" }, null).values()).getQ().getArray();
		CyclicBarrier barrier = new CyclicBarrier(2);
		AtomicInteger active = new AtomicInteger();
		AtomicInteger maximumActive = new AtomicInteger();
		QeQTLPreprocessor.Residualizer concurrent = new QeQTLPreprocessor.Residualizer() {
			@Override
			public double[][] residualize(double[][] values, double[][] projection, String matrixName) {
				int current = active.incrementAndGet();
				maximumActive.accumulateAndGet(current, Math::max);
				try {
					barrier.await();
					double[][] result = new double[values.length][values[0].length];
					for (int row = 0; row < values.length; row++) {
						System.arraycopy(values[row], 0, result[row], 0, values[row].length);
						for (int column = 0; column < projection[0].length; column++) {
							double coefficient = 0;
							for (int sample = 0; sample < projection.length; sample++)
								coefficient += values[row][sample] * projection[sample][column];
							for (int sample = 0; sample < projection.length; sample++)
								result[row][sample] -= coefficient * projection[sample][column];
						}
					}
					return result;
				} catch (Exception e) {
					throw new RuntimeException(e);
				} finally {
					active.decrementAndGet();
				}
			}

			@Override public int concurrency() { return 2; }
			@Override public String cacheSignatureTag() { return "test-concurrent"; }
		};
		String signature = QBinaryMatrixCache.signature("Genotype", genotype,
			alignment.genotypeColumnOrder(), q, concurrent.cacheSignatureTag());
		try (QBinaryMatrixCache cache = QBinaryMatrixCache.openOrBuild(directory, "Genotype",
			signature, genotype, alignment.genotypeColumnOrder(), q, 2, false, concurrent)) {
			assertArrayEquals(new String[] { "rs1", "rs2", "rs3" }, cache.readBlock(0, 3).rowIds());
		}
		assertEquals(2, maximumActive.get());
	}

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
            assertNotEquals(QBinaryMatrixCache.analysisSignature(cache, cache,
                    16, 16, 0, 4, 0.0, false, false, GpuPrecision.FP64),
                QBinaryMatrixCache.analysisSignature(cache, cache,
                    16, 16, 0, 4, 0.0, false, false, GpuPrecision.FP32));
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
		assertNotEquals(signature, QBinaryMatrixCache.signature("Genotype", genotype,
			alignment.genotypeColumnOrder(), covariateQ, "gpu-projection-v1-fp64-cuda"));

        try (RandomAccessFile mutable = new RandomAccessFile(cachePath.toFile(), "rw")) {
            mutable.seek(20); // Stable v1 header position of the index offset.
            long indexOffset = mutable.readLong();
            mutable.seek(indexOffset + Integer.BYTES + Long.BYTES);
            long firstRowOffset = mutable.readLong();
            mutable.seek(firstRowOffset + Integer.BYTES);
            int firstIdentifierByte = mutable.readUnsignedByte();
            mutable.seek(firstRowOffset + Integer.BYTES);
            mutable.writeByte(firstIdentifierByte ^ 1);
        }
        try (QBinaryMatrixCache corrupt = QBinaryMatrixCache.open(cachePath, "Genotype", signature)) {
            assertThrows(IOException.class, () -> corrupt.readBlock(0, 1));
        }
    }
}
