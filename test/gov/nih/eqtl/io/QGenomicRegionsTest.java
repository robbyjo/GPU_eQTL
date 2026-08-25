/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl.io;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

class QGenomicRegionsTest {
    @TempDir Path temporaryDirectory;

    @Test
    void parsesBedHeaderResolvesAliasesAndOrdersQueriesBySourceContig() throws Exception {
        Path file = temporaryDirectory.resolve("regions.tsv");
        Files.writeString(file, "set_id\tchrom\tstart\tend\nsetA\tchr1\t99\t110\n");
        QGenomicRegions regions = QGenomicRegions.load("setM=chrM:5-6", file,
            QGenomicRegions.Coordinates.BED, List.of("MT", "1"));
        assertEquals(List.of("setM", "setA"), regions.setIds());
        assertEquals(List.of(
            new QGenomicRegions.QueryInterval("MT", 5, 6),
            new QGenomicRegions.QueryInterval("1", 100, 110)), regions.queryIntervals());
        assertEquals("setA", regions.memberships(variant("1", 105)));
        assertEquals(List.of("setM"), regions.emptySetIds(Set.of("setA")));
    }

    @Test
    void rejectsAmbiguousCaseAndChrAlias() {
        assertThrows(IOException.class, () -> QGenomicRegions.load("x=CHR1:1-2", null,
            QGenomicRegions.Coordinates.ONE_BASED, List.of("1", "chr1")));
    }

    private static VariantContext variant(String contig, int position) {
        return new VariantContextBuilder("test", contig, position, position,
            List.of(Allele.create("A", true), Allele.create("G"))).make();
    }
}
