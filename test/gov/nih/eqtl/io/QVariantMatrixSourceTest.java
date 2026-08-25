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

import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
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
        assertTrue(qcLines.get(2).contains(
            "singleton\tfalse\tmac_below_minimum\t.\taligned\tmac_below_minimum"));
        assertTrue(qcLines.get(4).contains("LowQual\tDS"));
        assertTrue(qcLines.get(7).contains("multiallelic\tfalse\tmultiallelic\t.\taligned\t."));
    }

    @Test
    void indexedRegionsMergeQueriesPreserveOverlappingSetsAndReportEmptySets() throws Exception {
        IndexedFixture fixture = makeIndexedVcfFixture();
        Path regions = temporaryDirectory.resolve("regions.tsv");
        Files.writeString(regions,
            "setA\t1\t115\t135\nsetB\t1\t125\t145\nempty\t1\t900\t950\n");
        Path qc = temporaryDirectory.resolve("indexed-regions.tsv");
        QVariantMatrixSource.Options options = new QVariantMatrixSource.Options(
            QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
            QVariantMatrixSource.MissingPolicy.MEAN,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 0, qc, 2, null,
            fixture.index(), null, regions, QGenomicRegions.Coordinates.ONE_BASED,
            QVariantMatrixSource.FrequencyScope.ALIGNED);

        QVariantMatrixSource source = new QVariantMatrixSource(fixture.vcf(), options);
        assertEquals(fixture.index(), source.variantIndex());
        assertEquals(2, source.regions().queryIntervals().size());
        assertEquals(3, source.summary().inputRecords());
        assertEquals(3, source.metadata().rowCount());
        assertEquals(3, source.summary().regionSets());
        assertEquals(List.of("empty"), source.summary().emptyRegionSets());
        try (QMatrixRowSource.BlockReader reader = source.open(null)) {
            assertArrayEquals(new String[] {"1:120:G:A", "1:130:T:C", "1:140:A:C"},
                reader.readBlock(10).rowIds());
        }
        List<String> qcLines = Files.readAllLines(qc);
        assertTrue(qcLines.get(0).endsWith(
            "region_sets\tfrequency_scope\taligned_frequency_reason"));
        assertTrue(qcLines.stream().anyMatch(line -> line.startsWith("1:130:T:C\t")
            && line.contains("\tsetA;setB\taligned\t.")));
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
    void queriesIndexedBcfWithTheBcf22CompatibleCodec() throws Exception {
        IndexedFixture fixture = makeIndexedBcfFixture();
        Path qc = temporaryDirectory.resolve("indexed-bcf.tsv");
        QVariantMatrixSource.Options options = new QVariantMatrixSource.Options(
            QVariantMatrixSource.Format.BCF, QVariantMatrixSource.GenotypeField.AUTO,
            QVariantMatrixSource.MissingPolicy.MEAN,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 0, qc, 2, null,
            fixture.index(), "set=1:125-145", null, QGenomicRegions.Coordinates.ONE_BASED,
            QVariantMatrixSource.FrequencyScope.ALIGNED);
        QVariantMatrixSource source = new QVariantMatrixSource(fixture.vcf(), options);
        assertEquals(2, source.summary().inputRecords());
        try (QMatrixRowSource.BlockReader reader = source.open(null)) {
            assertArrayEquals(new String[] {"1:130:T:C", "1:140:A:C"},
                reader.readBlock(10).rowIds());
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
    void sequentialAndParallelQcProduceIdenticalOrderedResults() throws Exception {
        Path sequentialQc = temporaryDirectory.resolve("sequential-qc.tsv");
        Path parallelQc = temporaryDirectory.resolve("parallel-qc.tsv");
        Path genotype = gzipFixture();
        QVariantMatrixSource sequential = new QVariantMatrixSource(genotype,
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
                QVariantMatrixSource.MissingPolicy.MEAN, 0, 2, sequentialQc, 1));
        QVariantMatrixSource parallel = new QVariantMatrixSource(genotype,
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
                QVariantMatrixSource.MissingPolicy.MEAN, 0, 2, parallelQc, 4));

        assertEquals(1, sequential.qcThreadCount());
        assertEquals(4, parallel.qcThreadCount());
        assertEquals(sequential.summary(), parallel.summary());
        assertArrayEquals(Files.readAllBytes(sequentialQc), Files.readAllBytes(parallelQc));
        try (QMatrixRowSource.BlockReader sequentialReader = sequential.open(new int[] {1, 0, 2, 3});
             QMatrixRowSource.BlockReader parallelReader = parallel.open(new int[] {1, 0, 2, 3})) {
            QMatrixRowSource.Block sequentialBlock = sequentialReader.readBlock(10);
            QMatrixRowSource.Block parallelBlock = parallelReader.readBlock(10);
            assertArrayEquals(sequentialBlock.rowIds(), parallelBlock.rowIds());
            for (int row = 0; row < sequentialBlock.rowCount(); row++)
                assertArrayEquals(sequentialBlock.values()[row], parallelBlock.values()[row], 0);
        }
    }

    @Test
    void interruptedQcResumesFromDurableOrderedPartsAndCompletedScanIsReusable() throws Exception {
        Path genotype = gzipFixture();
        Path resumedQc = temporaryDirectory.resolve("resumed-qc.tsv");
        Path checkpoint = temporaryDirectory.resolve("resumed-checkpoint");
        String batchProperty = "eqtl.variant.qc.checkpoint.records";
        String failProperty = "eqtl.test.variant.qc.fail.after";
        System.setProperty(batchProperty, "2");
        System.setProperty(failProperty, "4");
        try {
            QVariantMatrixSource.Options resumableOptions = new QVariantMatrixSource.Options(
                QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
                QVariantMatrixSource.MissingPolicy.MEAN,
                QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 2, resumedQc, 3, checkpoint);
            IOException interrupted = assertThrows(IOException.class,
                () -> new QVariantMatrixSource(genotype, resumableOptions));
            assertTrue(interrupted.getMessage().contains("Injected variant-QC interruption"));
            assertTrue(!Files.exists(resumedQc));

            System.clearProperty(failProperty);
            QVariantMatrixSource resumed = new QVariantMatrixSource(genotype, resumableOptions);
            assertEquals(4, resumed.resumedQcRecords());
            assertEquals(7, resumed.summary().inputRecords());
            assertEquals(4, resumed.summary().includedVariants());

            Path freshQc = temporaryDirectory.resolve("fresh-qc.tsv");
            QVariantMatrixSource fresh = new QVariantMatrixSource(genotype,
                new QVariantMatrixSource.Options(QVariantMatrixSource.Format.VCF,
                    QVariantMatrixSource.GenotypeField.AUTO, QVariantMatrixSource.MissingPolicy.MEAN,
                    QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 2, freshQc, 3,
                    temporaryDirectory.resolve("fresh-checkpoint")));
            assertEquals(fresh.summary(), resumed.summary());
            assertArrayEquals(Files.readAllBytes(freshQc), Files.readAllBytes(resumedQc));

            Files.delete(resumedQc);
            QVariantMatrixSource reused = new QVariantMatrixSource(genotype, resumableOptions);
            assertEquals(7, reused.resumedQcRecords());
            assertTrue(Files.isRegularFile(resumedQc));
            assertArrayEquals(Files.readAllBytes(freshQc), Files.readAllBytes(resumedQc));
        } finally {
            System.clearProperty(batchProperty);
            System.clearProperty(failProperty);
        }
    }

    @Test
    void interruptedIndexedQcSeeksDirectlyToDurableBoundary() throws Exception {
        IndexedFixture fixture = makeIndexedVcfFixture();
        Path qc = temporaryDirectory.resolve("indexed-resume.tsv");
        Path checkpoint = temporaryDirectory.resolve("indexed-resume-checkpoint");
        String batchProperty = "eqtl.variant.qc.checkpoint.records";
        String failProperty = "eqtl.test.variant.qc.fail.after";
        QVariantMatrixSource.Options options = new QVariantMatrixSource.Options(
            QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
            QVariantMatrixSource.MissingPolicy.MEAN,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 0, qc, 2, checkpoint,
            fixture.index(), "all=1:1-1000", null, QGenomicRegions.Coordinates.ONE_BASED,
            QVariantMatrixSource.FrequencyScope.ALIGNED);
        System.setProperty(batchProperty, "2");
        System.setProperty(failProperty, "4");
        try {
            assertThrows(IOException.class, () -> new QVariantMatrixSource(fixture.vcf(), options));
            System.clearProperty(failProperty);
            QVariantMatrixSource resumed = new QVariantMatrixSource(fixture.vcf(), options);
            assertEquals(4, resumed.resumedQcRecords());
            assertTrue(resumed.indexedResumeUsed());
            assertEquals(7, resumed.summary().inputRecords());
            assertTrue(Files.isRegularFile(qc));
        } finally {
            System.clearProperty(batchProperty);
            System.clearProperty(failProperty);
        }
    }

    @Test
    void rejectsIndexMutationBetweenAlignmentAndQc() throws Exception {
        IndexedFixture fixture = makeIndexedVcfFixture();
        QVariantMatrixSource.Options options = new QVariantMatrixSource.Options(
            QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.AUTO,
            QVariantMatrixSource.MissingPolicy.MEAN,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 0, null, 1, null,
            fixture.index(), "all=1:1-1000", null, QGenomicRegions.Coordinates.ONE_BASED,
            QVariantMatrixSource.FrequencyScope.ALIGNED);
        QVariantMatrixSource source = QVariantMatrixSource.openForAlignment(fixture.vcf(), options);
        Files.write(fixture.index(), new byte[] {0}, java.nio.file.StandardOpenOption.APPEND);
        IOException error = assertThrows(IOException.class,
            () -> source.selectAnalysisSamples(new int[] {0, 1, 2, 3}));
        assertTrue(error.getMessage().contains("index changed"));
    }

    @Test
    void excludedHeaderSamplesAreNotEvaluatedForDosageQc() throws Exception {
        Path invalidExcludedDosage = temporaryDirectory.resolve("invalid-excluded-dosage.vcf");
        String fixture = Files.readString(fixturePath());
        String validRecord = "1\t110\trsSingle\tC\tT\t.\tPASS\t.\tGT:DS\t0/0:0\t0/1:1\t0/0:0\t0/0:0";
        String invalidRecord = validRecord.replace("0/1:1", "0/1:9");
        assertTrue(fixture.contains(validRecord));
        Files.writeString(invalidExcludedDosage, fixture.replace(validRecord, invalidRecord));

        QVariantMatrixSource subset = QVariantMatrixSource.openForAlignment(invalidExcludedDosage,
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.DS,
                QVariantMatrixSource.MissingPolicy.PRESERVE, 0, 0, null, 3));
        subset.selectAnalysisSamples(new int[] {0, 2, 3});
        assertEquals(7, subset.summary().inputRecords());

        QVariantMatrixSource allSamples = QVariantMatrixSource.openForAlignment(invalidExcludedDosage,
            options(QVariantMatrixSource.Format.VCF, QVariantMatrixSource.GenotypeField.DS,
                QVariantMatrixSource.MissingPolicy.PRESERVE, 0, 0, null, 3));
        IOException error = assertThrows(IOException.class,
            () -> allSamples.selectAnalysisSamples(new int[] {0, 1, 2, 3}));
        assertTrue(error.getMessage().contains("Invalid DS value '9'"));
    }

    @Test
    void automaticQcThreadRecommendationRetainsVariantParallelism() {
        assertEquals(1, QVariantMatrixSource.recommendQcThreads(1));
        assertEquals(4, QVariantMatrixSource.recommendQcThreads(4));
        int automatic = QVariantMatrixSource.recommendQcThreads(0);
        assertTrue(automatic >= 1 && automatic <= 16);
        if (Runtime.getRuntime().availableProcessors() > 1)
            assertTrue(automatic >= 2);
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
        return options(format, field, missing, minMaf, minMac, qc, 2);
    }

    private QVariantMatrixSource.Options options(QVariantMatrixSource.Format format,
        QVariantMatrixSource.GenotypeField field, QVariantMatrixSource.MissingPolicy missing,
        double minMaf, double minMac, Path qc, int qcThreads) {
        return new QVariantMatrixSource.Options(format, field, missing,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, minMaf, minMac, qc, qcThreads);
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

    private record IndexedFixture(Path vcf, Path index) { }

    private IndexedFixture makeIndexedVcfFixture() throws Exception {
        Path output = temporaryDirectory.resolve("indexed.vcf.gz");
        try (VCFIterator input = new VCFIteratorBuilder().open(fixturePath());
             VariantContextWriter writer = new VariantContextWriterBuilder()
                 .setOutputPath(output)
                 .setReferenceDictionary(input.getHeader().getSequenceDictionary())
                 .setIndexCreator(new TabixIndexCreator(input.getHeader().getSequenceDictionary(),
                     TabixFormat.VCF))
                 .setOption(Options.INDEX_ON_THE_FLY)
                 .build()) {
            writer.writeHeader(input.getHeader());
            while (input.hasNext())
                writer.add(input.next());
        }
        Path index = Path.of(output.toString() + ".tbi");
        assertTrue(Files.isRegularFile(index));
        return new IndexedFixture(output, index);
    }

    private IndexedFixture makeIndexedBcfFixture() throws Exception {
        Path output = temporaryDirectory.resolve("indexed.bcf");
        try (VCFIterator input = new VCFIteratorBuilder().open(fixturePath());
             VariantContextWriter writer = new VariantContextWriterBuilder()
                 .setOutputPath(output)
                 .setReferenceDictionary(input.getHeader().getSequenceDictionary())
                 .setOption(Options.FORCE_BCF)
                 .setOption(Options.INDEX_ON_THE_FLY)
                 .build()) {
            writer.writeHeader(input.getHeader());
            while (input.hasNext())
                writer.add(input.next());
        }
        Path index = Path.of(output.toString() + ".idx");
        assertTrue(Files.isRegularFile(index));
        byte[] bcf = Files.readAllBytes(output);
        assertEquals(2, bcf[3]);
        assertEquals(1, bcf[4]);
        bcf[4] = 2;
        Files.write(output, bcf);
        return new IndexedFixture(output, index);
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
