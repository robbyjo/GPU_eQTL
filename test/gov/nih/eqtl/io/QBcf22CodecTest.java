/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.io;

import static org.junit.jupiter.api.Assertions.*;

import java.io.ByteArrayInputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Base64;
import java.util.List;
import java.util.function.UnaryOperator;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFIteratorBuilder;

class QBcf22CodecTest {
    @TempDir Path temporary;

    @Test
    void nativeSparseTagDictionaryPreservesAllDecodedFieldsAndSampleOrder() throws Exception {
        assertVariantsMatchVcf("sparse-dictionary");
    }

    @Test
    void nativeSparseContigDictionaryUsesEncodedIndicesRatherThanHeaderOrder() throws Exception {
        assertVariantsMatchVcf("sparse-contigs");
    }

    @Test
    void alignedGtAndDsQcMatchVcfIncludingNonPassAndMissingCalls() throws Exception {
        for (var field : List.of(QVariantMatrixSource.GenotypeField.GT, QVariantMatrixSource.GenotypeField.DS)) {
            Path bcf = materialize("sparse-dictionary");
            Path vcfQc = temporary.resolve(field + "-vcf.tsv"), bcfQc = temporary.resolve(field + "-bcf.tsv");
            var vcf = QVariantMatrixSource.openForAlignment(fixture("sparse-dictionary.vcf"),
                options(QVariantMatrixSource.Format.VCF, field, vcfQc));
            var binary = QVariantMatrixSource.openForAlignment(bcf,
                options(QVariantMatrixSource.Format.BCF, field, bcfQc));
            int[] selected = {2, 0, 1};
            vcf.selectAnalysisSamples(selected); binary.selectAnalysisSamples(selected);
            assertEquals(3, binary.metadata().rowCount());
            assertTrue(binary.metadata().cacheSignatureTag().endsWith(";" + QVariantMatrixSource.BCF_READER_SIGNATURE));
            assertFalse(vcf.metadata().cacheSignatureTag().contains(QVariantMatrixSource.BCF_READER_SIGNATURE));
            assertEquals(vcf.summary(), binary.summary());
            assertEquals(Files.readString(vcfQc), Files.readString(bcfQc));
            for (int pass = 0; pass < 2; pass++) {
                try (var a = vcf.open(selected); var b = binary.open(selected)) {
                    var expected = a.readBlock(10); var actual = b.readBlock(10);
                    assertArrayEquals(expected.rowIds(), actual.rowIds());
                    for (int row = 0; row < 3; row++)
                        assertArrayEquals(expected.values()[row], actual.values()[row], 0);
                    assertNull(a.readBlock(1)); assertNull(b.readBlock(1));
                }
            }
        }
    }

    @Test
    void indexedTraversalKeepsOriginalHeaderOffsetsAndSparseContigIdentities() throws Exception {
        Path bcf = temporary.resolve("indexed.bcf");
        Files.write(bcf, unpack("sparse-contigs"));
        var index = IndexFactory.createDynamicIndex(bcf.toFile(), new QBcf22Codec());
        Path indexPath = temporary.resolve("indexed.bcf.idx"); index.write(indexPath);
        var options = new QVariantMatrixSource.Options(QVariantMatrixSource.Format.BCF,
            QVariantMatrixSource.GenotypeField.GT, QVariantMatrixSource.MissingPolicy.PRESERVE,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 0, null, 1, null,
            indexPath, "1:150-250", null, QGenomicRegions.Coordinates.ONE_BASED,
            QVariantMatrixSource.FrequencyScope.ALIGNED);
        var source = new QVariantMatrixSource(bcf, options);
        assertEquals(1, source.metadata().rowCount());
        try (var rows = source.open(null)) {
            assertArrayEquals(new String[]{"1:200:C:T"}, rows.readBlock(2).rowIds());
        }
    }

    @Test
    void correctedBcfReaderDoesNotReuseAnOlderRawCacheIdentity() throws Exception {
        var source = new QVariantMatrixSource(materialize("sparse-dictionary"),
            options(QVariantMatrixSource.Format.BCF, QVariantMatrixSource.GenotypeField.GT, null));
        QMatrixRowSource legacyIdentity = new QMatrixRowSource() {
            @Override public Metadata metadata() {
                var m = source.metadata();
                return new Metadata(m.path(), m.rowCount(), m.columnCount(), m.sampleIds(),
                    m.cacheSignatureTag().replace(";" + QVariantMatrixSource.BCF_READER_SIGNATURE, ""));
            }
            @Override public BlockReader open(int[] columns) throws java.io.IOException { return source.open(columns); }
        };
        Path cache = temporary.resolve("raw-cache");
        String oldSignature = QRawMatrixCache.signature(legacyIdentity, null);
        String correctedSignature = QRawMatrixCache.signature(source, null);
        assertNotEquals(oldSignature, correctedSignature);
        try (var ignored = QRawMatrixCache.openOrBuild(cache, oldSignature, legacyIdentity, null, 2, false)) { }
        assertNull(QRawMatrixCache.openIfPresent(cache, correctedSignature, source));
        try (var rebuilt = QRawMatrixCache.openOrBuild(cache, correctedSignature, source, null, 2, false);
                var rows = rebuilt.open(null)) {
            assertArrayEquals(new String[]{"1:100:A:G", "1:200:C:T", "2:300:G:A"}, rows.readBlock(4).rowIds());
        }
    }

    @Test
    void implicitIndicesFollowTheLargestPrecedingIndexAndSharedIdsKeepTheirSlot() throws Exception {
        // Preserve the native record IDs while making selected header IDs implicit.
        // An unused explicit predecessor advances the next implicit slot, regardless
        // of how many names have been seen. A repeated namespace reuses its first ID.
        byte[] bytes = changeHeader(unpack("sparse-dictionary"), h -> h
            .replace("##INFO=<ID=AF", "##INFO=<ID=Unused,Number=0,Type=Flag,Description=\"Unused\",IDX=52>\n##INFO=<ID=AF")
            .replace(",IDX=53>", ">")
            .replace("second namespace\",IDX=8>", "second namespace\">")
            .replace(",IDX=1>", ">"));
        var rows = decode(bytes);
        assertEquals(3, rows.size());
        assertEquals("2", rows.get(2).getContig());
        assertEquals(0.5, rows.get(0).getAttributeAsDouble("AF", -1), 1e-7);
        assertTrue(rows.get(1).hasAttribute("LowQual"));
    }

    @Test
    void rejectsConflictingDuplicateAndInvalidExplicitIndices() throws Exception {
        for (UnaryOperator<String> change : List.<UnaryOperator<String>>of(
                h -> h.replace("IDX=31>", "IDX=27>"),
                h -> h.replace("second namespace\",IDX=8>", "second namespace\",IDX=9>"),
                h -> h.replace("IDX=53>", "IDX=-1>"),
                h -> h.replace(",IDX=8>", ">").replace("second namespace\">", "second namespace\",IDX=8>"),
                h -> h.replace("IDX=53>", "IDX=2147483646>"),
                h -> h.replace("IDX=53>", "IDX=2147483648>"),
                h -> h.replace("IDX=53>", "IDX=x>"))) {
            assertThrows(TribbleException.class, () -> decode(changeHeader(unpack("sparse-dictionary"), change)));
        }
    }

    @Test
    void rejectsUndeclaredInfoFormatFilterAndContigSlotsInsteadOfMislabeling() throws Exception {
        for (String token : List.of("IDX=53>", "IDX=27>", "IDX=8>", "IDX=1>")) {
            byte[] bad = changeHeader(unpack("sparse-dictionary"), h -> h.replace(token, "IDX=98>"));
            var error = assertThrows(TribbleException.class, () -> decode(bad));
            assertTrue(error.getMessage().contains("Undeclared BCF"), error.toString());
        }
    }

    @Test
    void verySparseHeaderDoesNotAllocateAnArrayUpToTheLargestIndex() throws Exception {
        byte[] bytes = changeHeader(unpack("sparse-dictionary"), h -> h.replace("IDX=53>", "IDX=1000000000>"));
        // The record still refers to 53, so it must fail as undeclared, not OOM.
        var error = assertThrows(TribbleException.class, () -> decode(bytes));
        assertTrue(error.getMessage().contains("dictionary index 53"));
    }

    private void assertVariantsMatchVcf(String name) throws Exception {
        List<VariantContext> binary = decode(unpack(name));
        try (var reader = new VCFIteratorBuilder().open(fixture("sparse-dictionary.vcf"))) {
            assertEquals(3, binary.size());
            for (var actual : binary) {
                var expected = reader.next();
                assertEquals(expected.getContig(), actual.getContig());
                assertEquals(expected.getStart(), actual.getStart());
                assertEquals(expected.getID(), actual.getID());
                assertEquals(expected.getAlleles(), actual.getAlleles());
                assertEquals(expected.getFilters(), actual.getFilters());
                assertEquals(expected.getAttributeAsDouble("AF", -1), actual.getAttributeAsDouble("AF", -1), 1e-7);
                assertEquals(expected.hasAttribute("LowQual"), actual.hasAttribute("LowQual"));
                for (String sample : List.of("S2", "S1", "S3", "S4")) {
                    var a = expected.getGenotype(sample); var b = actual.getGenotype(sample);
                    assertEquals(a.getAlleles(), b.getAlleles());
                    assertEquals(a.isPhased(), b.isPhased());
                    Object expectedDs = a.getExtendedAttribute("DS"), actualDs = b.getExtendedAttribute("DS");
                    if (expectedDs == null || expectedDs.toString().equals("."))
                        assertTrue(actualDs == null || actualDs.toString().equals("."));
                    else
                        assertEquals(Double.parseDouble(expectedDs.toString()), Double.parseDouble(actualDs.toString()), 0);
                }
            }
            assertFalse(reader.hasNext());
        }
    }

    private List<VariantContext> decode(byte[] bytes) throws Exception {
        var result = new ArrayList<VariantContext>();
        var codec = new QBcf22Codec();
        try (var input = new PositionalBufferedStream(new ByteArrayInputStream(bytes))) {
            var header = codec.readHeader(input);
            assertEquals(input.getPosition(), header.getHeaderEnd());
            while (!codec.isDone(input)) {
                var row = codec.decode(input);
                row.getGenotypes().size(); // Materialize lazy FORMAT data while its reader owns the codec.
                result.add(row);
            }
        }
        return result;
    }

    private byte[] unpack(String name) throws Exception {
        byte[] compressed = Base64.getMimeDecoder().decode(Files.readString(fixture(name + ".bcf.base64")));
        try (var input = new GZIPInputStream(new ByteArrayInputStream(compressed))) { return input.readAllBytes(); }
    }
    private Path materialize(String name) throws Exception {
        Path path = temporary.resolve(name + ".bcf"); Files.write(path, unpack(name)); return path;
    }
    private byte[] changeHeader(byte[] bytes, UnaryOperator<String> change) {
        int length = ByteBuffer.wrap(bytes, 5, 4).order(ByteOrder.LITTLE_ENDIAN).getInt();
        String old = new String(bytes, 9, length, StandardCharsets.UTF_8);
        String updated = change.apply(old); assertNotEquals(old, updated, "Fixture mutation must take effect");
        byte[] header = updated.getBytes(StandardCharsets.UTF_8);
        var out = ByteBuffer.allocate(bytes.length - length + header.length).order(ByteOrder.LITTLE_ENDIAN);
        out.put(bytes, 0, 5).putInt(header.length).put(header).put(bytes, 9 + length, bytes.length - 9 - length);
        return out.array();
    }
    private Path fixture(String name) { return Path.of("test/resources/variant-reference", name); }
    private QVariantMatrixSource.Options options(QVariantMatrixSource.Format format,
            QVariantMatrixSource.GenotypeField field, Path qc) {
        return new QVariantMatrixSource.Options(format, field, QVariantMatrixSource.MissingPolicy.PRESERVE,
            QVariantMatrixSource.MultiallelicPolicy.EXCLUDE, 0, 0, qc, 4);
    }
}
