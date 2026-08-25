/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QMissingValuePolicy;

class QPatternVariantSourceTest {
    @TempDir Path temporaryDirectory;

    @Test
    void cachesExactSubsetStatisticsAndUsesPatternMeanForMissingDosage() throws Exception {
        QMatrixRowSource raw = source();
        int[] columns = {0, 1, 3};
        QPatternVariantSource first = QPatternVariantSource.openOrBuild(
            temporaryDirectory.resolve("cache"), "pattern-1", raw, columns,
            QMissingValuePolicy.MEAN, false, 0, 0, false);

        assertEquals(4, first.summary().inputVariants());
        assertEquals(2, first.summary().includedVariants());
        assertEquals(2, first.summary().monomorphicVariants());
        assertEquals(2, first.summary().missingGenotypes());
        assertTrue(!first.summary().reused());
        try (QMatrixRowSource.BlockReader reader = first.open(columns)) {
            QMatrixRowSource.Block block = reader.readBlock(10);
            assertArrayEquals(new String[] {"v1", "v2"}, block.rowIds());
            assertArrayEquals(new double[] {0, 1, 0.5}, block.values()[0], 0);
            assertArrayEquals(new double[] {0, 0, 2}, block.values()[1], 0);
            assertNull(reader.readBlock(1));
        }

        QPatternVariantSource reused = QPatternVariantSource.openOrBuild(
            temporaryDirectory.resolve("cache"), "pattern-1", raw, columns,
            QMissingValuePolicy.MEAN, false, 0, 0, false);
        assertTrue(reused.summary().reused());
        assertEquals(first.summary().includedVariants(), reused.summary().includedVariants());
    }

    @Test
    void appliesOptionalPatternSpecificMacFilterAfterSampleRemoval() throws Exception {
        int[] columns = {0, 1, 3};
        QPatternVariantSource filtered = QPatternVariantSource.openOrBuild(
            temporaryDirectory.resolve("filtered"), "pattern-2", source(), columns,
            QMissingValuePolicy.MEAN, true, 0, 2, false);
        assertEquals(1, filtered.summary().includedVariants());
        assertEquals(1, filtered.summary().belowMinimumMac());
        try (QMatrixRowSource.BlockReader reader = filtered.open(columns)) {
            QMatrixRowSource.Block block = reader.readBlock(10);
            assertArrayEquals(new String[] {"v2"}, block.rowIds());
            assertArrayEquals(new double[] {0, 0, 2}, block.values()[0], 0);
        }
    }

    @Test
    void zeroFillUsesTheFilledPatternVarianceRatherThanObservedOnlyVariance() throws Exception {
        int[] columns = {0, 1, 3};
        QPatternVariantSource zero = QPatternVariantSource.openOrBuild(
            temporaryDirectory.resolve("zero"), "pattern-zero", source(), columns,
            QMissingValuePolicy.ZERO, false, 0, 0, false);
        assertEquals(3, zero.summary().includedVariants());
        try (QMatrixRowSource.BlockReader reader = zero.open(columns)) {
            QMatrixRowSource.Block block = reader.readBlock(10);
            assertArrayEquals(new String[] {"v1", "v2", "v4"}, block.rowIds());
            assertArrayEquals(new double[] {1, 1, 0}, block.values()[2], 0);
        }
    }

    @Test
    void alignedVcfCandidatesAreReadjustedExactlyForTraitPattern() throws Exception {
        QVariantMatrixSource.Options options = new QVariantMatrixSource.Options(
            QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
            QVariantMatrixSource.MissingPolicy.PRESERVE,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 2, null, 2, null,
            null, null, null, QGenomicRegions.Coordinates.ONE_BASED,
            QVariantMatrixSource.FrequencyScope.PATTERN);
        QVariantMatrixSource variants = new QVariantMatrixSource(
            Path.of("test/resources/variant-reference/genotype.vcf"), options);
        assertEquals(5, variants.metadata().rowCount());
        String rawSignature = QRawMatrixCache.signature(variants, null);
        try (QRawMatrixCache raw = QRawMatrixCache.openOrBuild(
            temporaryDirectory.resolve("raw"), rawSignature, variants, null, 2, false)) {
            // Source order is S2,S1,S3,S4; remove S1 for the second trait pattern.
            int[] observed = {0, 2, 3};
            QPatternVariantSource pattern = QPatternVariantSource.openOrBuild(
                temporaryDirectory.resolve("pattern"), "missing-s1", raw, observed,
                QMissingValuePolicy.MEAN, true, 0, 2, false);
            assertEquals(5, pattern.summary().inputVariants());
            assertEquals(3, pattern.summary().includedVariants());
            assertEquals(2, pattern.summary().monomorphicVariants());
            try (QMatrixRowSource.BlockReader reader = pattern.open(observed)) {
                assertArrayEquals(new String[] {"1:130:T:C", "1:140:A:C", "1:150:G:T"},
                    reader.readBlock(10).rowIds());
            }
        }
    }

    private QMatrixRowSource source() throws IOException {
        Path path = temporaryDirectory.resolve("raw-source.bin");
        if (!Files.exists(path))
            Files.write(path, new byte[] {1});
        return new QMatrixRowSource() {
            private final String[] samples = {"S1", "S2", "S3", "S4"};
            private final String[] ids = {"v1", "v2", "v3", "v4"};
            private final double[][] rows = {
                {0, 1, 2, Double.NaN},
                {0, 0, 2, 2},
                {1, 1, 2, 1},
                {1, 1, 2, Double.NaN}
            };

            @Override
            public Metadata metadata() {
                return new Metadata(path, rows.length, samples.length, samples, "test-raw-v1");
            }

            @Override
            public BlockReader open(int[] columnOrder) {
                int[] selected = columnOrder == null ? new int[] {0, 1, 2, 3} : columnOrder.clone();
                return new BlockReader() {
                    private boolean read;

                    @Override
                    public Block readBlock(int maximumRows) {
                        if (read)
                            return null;
                        read = true;
                        double[][] values = new double[rows.length][selected.length];
                        for (int row = 0; row < rows.length; row++)
                            for (int column = 0; column < selected.length; column++)
                                values[row][column] = rows[row][selected[column]];
                        return new Block(0, ids.clone(), values);
                    }

                    @Override public void close() { }
                };
            }
        };
    }
}
