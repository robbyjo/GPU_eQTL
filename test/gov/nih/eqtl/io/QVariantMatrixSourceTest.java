/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

class QVariantMatrixSourceTest {
    @TempDir Path temporaryDirectory;

    @Test
    void readsGzippedVcfFiltersVariantsAndPreservesHeaderSampleOrder() throws Exception {
        Path gzipped = gzipFixture();
        Path qc = temporaryDirectory.resolve("variant-qc.tsv");
        QVariantMatrixSource source = new QVariantMatrixSource(gzipped,
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
                QVariantMatrixSource.MissingPolicy.MEAN, 0, 2, qc));

        assertEquals(QVariantMatrixSource.GenotypeField.DS, source.selectedField());
        assertArrayEquals(new String[] {"S2", "S1", "S3", "S4"}, source.metadata().sampleIds());
        assertEquals(4, source.metadata().rowCount());
        assertEquals(7, source.summary().inputRecords());
        assertEquals(1, source.summary().monomorphicVariants());
        assertEquals(1, source.summary().singletons());
        assertEquals(1, source.summary().doubletons());
        assertEquals(1, source.summary().multiallelicRecords());

        try (QMatrixRowSource.BlockReader reader = source.open(new int[] {1, 0, 2, 3})) {
            QMatrixRowSource.Block first = reader.readBlock(2);
            assertArrayEquals(new String[] {"1:120:G:A", "1:130:T:C"}, first.rowIds());
            assertArrayEquals(new double[] {2, 0, 0, 0}, first.values()[0], 0);
            assertArrayEquals(new double[] {1, 0, 1, 2}, first.values()[1], 0);

            QMatrixRowSource.Block second = reader.readBlock(2);
            assertArrayEquals(new String[] {"1:140:A:C", "1:150:G:T"}, second.rowIds());
            assertArrayEquals(new double[] {1, 0, 1, 2}, second.values()[0], 0);
            assertArrayEquals(new double[] {0.8, 0.2, 1.4, 2}, second.values()[1], 0);
            assertNull(reader.readBlock(1));
        }

        List<String> qcLines = Files.readAllLines(qc);
        assertEquals(8, qcLines.size());
        assertTrue(qcLines.get(1).contains("rsMono\tA\tG"));
        assertTrue(qcLines.get(2).endsWith("singleton\tfalse\tmac_below_minimum"));
        assertTrue(qcLines.get(4).contains("LowQual\tDS"));
        assertTrue(qcLines.get(7).endsWith("multiallelic\tfalse\tmultiallelic"));
    }

    @Test
    void missingDosageIsFatalByDefault() throws Exception {
        Path gzipped = gzipFixture();
        IOException error = assertThrows(IOException.class, () -> new QVariantMatrixSource(gzipped,
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.DS,
                QVariantMatrixSource.MissingPolicy.ERROR, 0, 0, null)));
        assertTrue(error.getMessage().contains("--genotype-missing"));
    }

    @Test
    void alignedSubsetControlsVariantFrequencyAndMonomorphicQc() throws Exception {
        Path qc = temporaryDirectory.resolve("subset-qc.tsv");
        QVariantMatrixSource source = QVariantMatrixSource.openForAlignment(gzipFixture(),
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.DS,
                QVariantMatrixSource.MissingPolicy.PRESERVE, 0, 0, qc));
        assertEquals(-1, source.metadata().rowCount());
        assertThrows(IllegalStateException.class, source::summary);

        source.selectAnalysisSamples(new int[] {0, 2, 3}); // Exclude S1, the only ALT carrier at 110/120.
        assertEquals(3, source.analysisSampleCount());
        assertEquals(3, source.metadata().rowCount());
        assertEquals(3, source.summary().monomorphicVariants());
        try (QMatrixRowSource.BlockReader reader = source.open(new int[] {0, 2, 3})) {
            assertArrayEquals(new String[] {"1:130:T:C", "1:140:A:C", "1:150:G:T"},
                reader.readBlock(10).rowIds());
        }
        String qcText = Files.readString(qc);
        assertTrue(qcText.contains("1:120:G:A\t1\t120\trsDouble"));
        assertTrue(qcText.contains("\tmonomorphic\tfalse\tmonomorphic"));
        assertTrue(qcText.contains("1:140:A:C\t1\t140\trsMissing"));
    }

	@Test
	void preservedMissingDosageRemainsVisibleToTheCommonQcScanner() throws Exception {
		QVariantMatrixSource source = new QVariantMatrixSource(gzipFixture(),
			options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.DS,
				QVariantMatrixSource.MissingPolicy.PRESERVE, 0, 0, null));
		QMissingnessScan scan = QMissingnessScan.scan("predictor", source, null);
		assertTrue(scan.totalMissingValues() > 0);
		assertTrue(scan.missingRows().stream().anyMatch(row -> row.rowId().equals("1:140:A:C")));
	}

    @Test
    void readsBcf22VersionWithTheRestrictedDiploidReader() throws Exception {
        Path bcf22 = makeBcf22Fixture();
        QVariantMatrixSource source = new QVariantMatrixSource(bcf22,
            options(QVariantMatrixSource.Format.BCF, QVariantMatrixSource.GenotypeField.AUTO,
                QVariantMatrixSource.MissingPolicy.MEAN, 0.2, 0, null));

        assertEquals(QVariantMatrixSource.GenotypeField.DS, source.selectedField());
        assertArrayEquals(new String[] {"S2", "S1", "S3", "S4"}, source.metadata().sampleIds());
        assertEquals(4, source.metadata().rowCount());
        try (QMatrixRowSource.BlockReader reader = source.open(null)) {
            QMatrixRowSource.Block block = reader.readBlock(10);
            assertArrayEquals(new String[] {"1:120:G:A", "1:130:T:C", "1:140:A:C", "1:150:G:T"},
                block.rowIds());
            assertArrayEquals(new double[] {0, 2, 0, 0}, block.values()[0], 0);
            assertArrayEquals(new double[] {0.2, 0.8, 1.4, 2}, block.values()[3], 1e-6);
        }
    }

    @Test
    void exactHweCalculationIsSymmetricAndBounded() {
        double first = QVariantMatrixSource.hardyWeinbergExact(48, 40, 12);
        double reversed = QVariantMatrixSource.hardyWeinbergExact(12, 40, 48);
        assertEquals(first, reversed, 0);
        assertTrue(first >= 0 && first <= 1);
        assertEquals(1, QVariantMatrixSource.hardyWeinbergExact(4, 0, 0), 0);
    }

    @Test
    void formatsProgressWithCountsRatePercentageAndEta() {
        String knownTotal = QVariantMatrixSource.progressMessage("Variant input", 250_000,
            1_000_000, 200_000, 600_000_000_000L, false);
        assertTrue(knownTotal.contains("250,000 / 1,000,000 records (25.0%)"));
        assertTrue(knownTotal.contains("retained=200,000"));
        assertTrue(knownTotal.contains("elapsed=10.0 min"));
        assertTrue(knownTotal.contains("rate=416.7 records/s"));
        assertTrue(knownTotal.contains("ETA=30.0 min"));

        String unknownTotal = QVariantMatrixSource.progressMessage("Variant QC", 250_000,
            -1, 200_000, 60_000_000_000L, false);
        assertTrue(unknownTotal.contains("250,000 records scanned"));
        assertTrue(!unknownTotal.contains("%"));
        assertTrue(!unknownTotal.contains("ETA="));
    }

    @Test
    void readsExternalBcftoolsFixtureWhenConfigured() throws Exception {
        String configured = System.getProperty("eqtl.test.bcf");
        assumeTrue(configured != null && !configured.isBlank(),
            "Set -Deqtl.test.bcf to an independently generated BCF 2.2 file");
        QVariantMatrixSource source = new QVariantMatrixSource(Path.of(configured),
            options(QVariantMatrixSource.Format.BCF, QVariantMatrixSource.GenotypeField.GT,
                QVariantMatrixSource.MissingPolicy.MEAN, 0, 0, null));
        assertTrue(source.metadata().rowCount() > 0);
        assertTrue(source.metadata().columnCount() > 0);
        try (QMatrixRowSource.BlockReader reader = source.open(null)) {
            assertTrue(reader.readBlock(2).rowCount() > 0);
        }
    }

    private QVariantMatrixSource.Options options(QVariantMatrixSource.Format format,
        QVariantMatrixSource.GenotypeField field, QVariantMatrixSource.MissingPolicy missing,
        double minMaf, double minMac, Path qc) {
        return new QVariantMatrixSource.Options(format, field, missing,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, minMaf, minMac, qc);
    }

    private Path gzipFixture() throws IOException {
        Path output = temporaryDirectory.resolve("genotype.vcf.gz");
        try (InputStream input = fixture(); OutputStream raw = Files.newOutputStream(output);
             GZIPOutputStream gzip = new GZIPOutputStream(raw)) {
            input.transferTo(gzip);
        }
        return output;
    }

    private Path makeBcf22Fixture() throws Exception {
        Path bcf21 = temporaryDirectory.resolve("source.bcf");
        try (VCFIterator input = new VCFIteratorBuilder().open(fixturePath());
             VariantContextWriter writer = new VariantContextWriterBuilder()
                 .setOutputFile(bcf21.toFile())
                 .unsetOption(Options.INDEX_ON_THE_FLY)
                 .setOption(Options.FORCE_BCF)
                 .build()) {
            writer.writeHeader(input.getHeader());
            while (input.hasNext()) {
                VariantContext variant = input.next();
                writer.add(variant);
            }
        }

        byte[] uncompressed = Files.readAllBytes(bcf21);
        assertEquals('B', uncompressed[0]);
        assertEquals(2, uncompressed[3]);
        assertEquals(1, uncompressed[4]);
        uncompressed[4] = 2;
        Path bcf22 = temporaryDirectory.resolve("genotype.bcf");
        try (GZIPOutputStream output = new GZIPOutputStream(Files.newOutputStream(bcf22))) {
            output.write(uncompressed);
        }
        return bcf22;
    }

    private InputStream fixture() {
        try {
            return Files.newInputStream(fixturePath());
        } catch (Exception e) {
            throw new IllegalStateException(e);
        }
    }

    private Path fixturePath() {
        return Path.of("test", "resources", "variant-reference", "genotype.vcf")
            .toAbsolutePath().normalize();
    }
}
