/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.gpu.GpuPrecision;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

class QSlidingWindowAnalysisTest {
    @AfterEach
    void restoreState() {
        QeQTLAnalysis.config = null;
        QeQTLAnalysis.profiler = null;
        QeQTLAnalysis.mContexts = null;
        QeQTLAnalysis.simplifyResult = false;
        QeQTLAnalysis.rsqOnly = false;
        QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        QeQTLAnalysis.residualizationMode = QResidualizationMode.AUTO;
    }

    @Test
    void csvVcfAndBcfSlidingWindowsMatchAndCsvResumeIsByteIdentical(
        @TempDir Path directory) throws Exception {
        Path csv = directory.resolve("variants.csv");
        Files.writeString(csv, ",S2,S1,S3,S4\n"
            + "1:110:C:T,0,1,0,0\n"
            + "1:120:G:A,0,2,0,0\n"
            + "1:130:T:C,0,1,1,2\n"
            + "1:140:A:C,0,,1,2\n"
            + "1:150:G:T,0.25,0.75,1.5,2.0\n");
        Path expression = Path.of("test/resources/variant-reference/expression.csv")
            .toAbsolutePath().normalize();
        Path vcfFixture = Path.of("test/resources/variant-reference/genotype.vcf")
            .toAbsolutePath().normalize();
        Path vcf = directory.resolve("genotype.vcf");
        Files.writeString(vcf, Files.readString(vcfFixture).replace(
            "\t0/0:0.2\t0/1:0.8\t0/1:1.4\t1/1:2.0",
            "\t0/0:0.25\t0/1:0.75\t0/1:1.5\t1/1:2.0"));
        Path bcf = makeBcf22Fixture(directory, vcf);

        for (String method : List.of("burden", "skat", "skat-o")) {
            Path csvRoot = directory.resolve(method + "-csv");
            Path csvOutput = directory.resolve(method + "-csv.csv");
            run(csvRoot, csv, "csv", expression, csvOutput, method, false, true);
            byte[] initial = Files.readAllBytes(csvOutput);
            Path initialAudit = Path.of(csvOutput + ".sets.tsv");

            run(csvRoot, csv, "csv", expression, csvOutput, method, true, true);
            assertArrayEquals(initial, Files.readAllBytes(csvOutput), method + " resume");

            Path vcfOutput = directory.resolve(method + "-vcf.csv");
            run(directory.resolve(method + "-vcf"), vcf, "vcf", expression,
                vcfOutput, method, false, false);
            assertArrayEquals(initial, Files.readAllBytes(vcfOutput), method + " CSV/VCF");
            assertArrayEquals(Files.readAllBytes(initialAudit),
                Files.readAllBytes(Path.of(vcfOutput + ".sets.tsv")), method + " audit");

            Path bcfOutput = directory.resolve(method + "-bcf.csv");
            run(directory.resolve(method + "-bcf"), bcf, "bcf", expression,
                bcfOutput, method, false, false);
            assertArrayEquals(initial, Files.readAllBytes(bcfOutput), method + " CSV/BCF");
            assertArrayEquals(Files.readAllBytes(initialAudit),
                Files.readAllBytes(Path.of(bcfOutput + ".sets.tsv")), method + " BCF audit");
            assertTrue(Files.readAllLines(csvOutput).size() > 1);
            assertEquals(8, Files.readAllLines(initialAudit).size());
        }
    }

    private static Path makeBcf22Fixture(Path directory, Path vcf) throws Exception {
        Path bcf21 = directory.resolve("source.bcf");
        try (VCFIterator input = new VCFIteratorBuilder().open(vcf);
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
        Path bcf22 = directory.resolve("genotype.bcf");
        try (GZIPOutputStream output = new GZIPOutputStream(Files.newOutputStream(bcf22))) {
            output.write(uncompressed);
        }
        return bcf22;
    }

    private static void run(Path root, Path genotype, String format, Path expression,
        Path output, String method, boolean resume, boolean keepCheckpoints) throws Exception {
        Files.createDirectories(root);
        Map<String, String> values = new LinkedHashMap<>();
        values.put("ini.path", root.toString() + java.io.File.separator);
        values.put("genotype_file", genotype.toString());
        values.put("expression_file", expression.toString());
        values.put("output_file", output.toString());
        values.put("genotype_format", format);
        values.put("analysis", method);
        values.put("window_size", "30");
        values.put("window_stride", "10");
        values.put("set_block_size", "2");
        values.put("set_degenerate", "skip");
        values.put("predictor_missing", "mean");
        values.put("trait_missing", "error");
        values.put("min_mac", "0");
        values.put("min_maf", "0");
        values.put("variant_qc_threads", "1");
        values.put("genotype_block_rows", "2");
        values.put("expression_block_rows", "1");
        values.put("cache_dir", root.resolve("cache").toString());
        values.put("checkpoint_dir", root.resolve("checkpoint").toString());
        values.put("precision", "fp64");
        values.put("skat_o_simulations", "99");
        values.put("skat_o_seed", "73129");
        if (resume) values.put("resume", "true");
        if (keepCheckpoints) values.put("keep_checkpoints", "true");
        QeQTLAnalysis.config = new QeQTLAnalysisConfig(values);
        QeQTLAnalysis.profiler = new QeQTLProfiler(false);
        QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        QeQTLAnalysis.residualizationMode = QResidualizationMode.CPU;
        QeQTLAnalysis.simplifyResult = false;
        QeQTLAnalysis.rsqOnly = false;
        QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(), genotype.toString(),
            expression.toString(), null, null, null, "none", Double.NaN, 0, true,
            2, 1, format);
    }
}
