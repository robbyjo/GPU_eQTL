/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;

class QRelatedSetPipelineTest {
    record Fixture(Path g,Path y,Path cov,Path grm,Path sets) { }
    private Fixture fixture(Path root) throws Exception {
        var lines = Files.readAllLines(Path.of("test/resources/related-set-reference/data.csv"));
        StringBuilder g = new StringBuilder("ID"), y = new StringBuilder("ID"), cov = new StringBuilder("ID,x\n");
        for (int i = 1; i < lines.size(); i++) g.append(',').append(lines.get(i).split(",")[0]);
        for (int i = lines.size()-1; i >= 1; i--) y.append(',').append(lines.get(i).split(",")[0]);
        for (int j = 0; j < 3; j++) {
            g.append("\n1:").append(100+j).append(":A:G");
            for (int i = 1; i < lines.size(); i++) g.append(',').append(lines.get(i).split(",")[3+j].replace("NA",""));
        }
        for (String trait : List.of("complete","pattern","same-mask-new-trait")) {
            y.append('\n').append(trait);
            for (int i = lines.size()-1; i >= 1; i--) {
                y.append(',');
                if (!trait.equals("complete") && Set.of(3,7,11,20,29).contains(i)) continue;
                double value = Double.parseDouble(lines.get(i).split(",")[2]);
                y.append(trait.equals("same-mask-new-trait") ? value*value : value);
            }
        }
        for (int i = 1; i < lines.size(); i++) {
            String[] row = lines.get(i).split(","); cov.append(row[0]).append(',').append(row[1]).append('\n');
        }
        Path gp = root.resolve("g.csv"), yp = root.resolve("y.csv"), cp = root.resolve("cov.csv"), rp = root.resolve("grm.tsv"), sp = root.resolve("sets.tsv");
        Files.writeString(gp,g.append('\n')); Files.writeString(yp,y.append('\n')); Files.writeString(cp,cov);
        Files.copy(Path.of("test/resources/related-set-reference/grm.tsv"),rp);
        Files.writeString(sp,"SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\tWEIGHT\n"
            +"all\t1:100:A:G\tA\tG\tG\t1\nall\t1:101:A:G\tA\tG\tG\t.7\nall\t1:102:A:G\tA\tG\tG\t1.3\n"
            +"single\t1:100:A:G\tA\tG\tG\t1\n");
        return new Fixture(gp,yp,cp,rp,sp);
    }
    private QeQTLAnalysisConfig config(Fixture f,Path output,String method,boolean related) {
        Map<String,String> values = new LinkedHashMap<>();
        values.put("genotype_file",f.g().toString()); values.put("expression_file",f.y().toString()); values.put("output_file",output.toString());
        values.put("covariate_file",f.cov().toString()); values.put("covariate_fixed","x");
        values.put("genotype_id_column","ID"); values.put("expression_id_column","ID");
        values.put("variant_sets",f.sets().toString()); values.put("analysis",method); values.put("set_block_size","1");
        values.put("trait_missing","pattern"); values.put("predictor_missing","mean"); values.put("keep_checkpoints","true");
        values.put("skat_o_method","integration"); values.put("precision","fp64");
        if (related) values.put("relationship_matrix",f.grm().toString());
        return new QeQTLAnalysisConfig(values);
    }
    private void run(Fixture f,QeQTLAnalysisConfig config) throws Exception {
        QeQTLAnalysis.config = config; QeQTLAnalysis.profiler = new QeQTLProfiler(false); QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        String old = System.getProperty(CpuBackend.BLAS_PROPERTY); System.setProperty(CpuBackend.BLAS_PROPERTY,"java");
        try (var context = new CpuBackend().discoverGpuDevices().get(0).openContext()) {
            QeQTLAnalysis.mContexts = new GpuContext[]{context};
            QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(),f.g().toString(),f.y().toString(),f.cov().toString(),
                new String[]{"x"},null,"none",Double.NaN,0,true,16,1,"csv");
        } finally {
            if (old == null) System.clearProperty(CpuBackend.BLAS_PROPERTY); else System.setProperty(CpuBackend.BLAS_PROPERTY,old);
            QeQTLAnalysis.config = null; QeQTLAnalysis.mContexts = null;
        }
    }
    @Test void allMethodsRelatedAndFixedDeleteExactlyAndResumeWithAudits(@TempDir Path root) throws Exception {
        Fixture f = fixture(root);
        var references = Files.readAllLines(Path.of("test/resources/related-set-reference/scores.csv"));
        for (boolean related : new boolean[]{false,true}) for (String method : List.of("burden","skat","skat-o")) {
            Path output = root.resolve(method+"-"+related+".csv"); var config = config(f,output,method,related);
            run(f,config); byte[] before = Files.readAllBytes(output), audit = Files.readAllBytes(Path.of(output+".sets.tsv"));
            var rows = Files.readAllLines(output); assertEquals(7,rows.size());
            for (int trait = 0; trait < 3; trait++) {
                String[] row = rows.get(1+trait*2).split(","), ref = references.get(1+trait*2+(related ? 1 : 0)).split(",");
                assertEquals("all",row[0]); assertEquals(ref[0],row[1]); assertEquals(ref[2],row[4]);
                if (method.equals("skat")) assertEquals(Double.parseDouble(ref[14]),Double.parseDouble(row[6]),3e-6*Math.max(1,Double.parseDouble(ref[14])));
                if (method.equals("burden") && related) {
                    assertEquals("1",row[5]); assertEquals(Double.parseDouble(ref[13]),Double.parseDouble(row[6]),3e-6);
                }
            }
            // Recovery from an interrupted tile: preserve all other parts, ignore a .partial.
            Path checkpoints = Path.of(output+".checkpoint");
            Files.delete(checkpoints.resolve("block-00000003.part")); Files.writeString(checkpoints.resolve("block-00000003.partial"),"incomplete");
            config.set("resume","true"); run(f,config);
            assertArrayEquals(before,Files.readAllBytes(output)); assertArrayEquals(audit,Files.readAllBytes(Path.of(output+".sets.tsv")));
            config.set("skat_o_relative_tolerance","0.01");
            assertThrows(java.io.IOException.class,()->run(f,config)); assertArrayEquals(before,Files.readAllBytes(output));
        }
    }
    @Test void masksAndGrmContentBindResumeAndValidationDoesNotFit(@TempDir Path root) throws Exception {
        Fixture f = fixture(root); Path output = root.resolve("out.csv"); var config = config(f,output,"burden",true);
        run(f,config); byte[] before = Files.readAllBytes(output); config.set("resume","true");
        var stamp = Files.getLastModifiedTime(f.y()); String original = Files.readString(f.y());
        Files.writeString(f.y(),original.replace("pattern,","patterN,")); Files.setLastModifiedTime(f.y(),stamp);
        assertThrows(java.io.IOException.class,()->run(f,config));
        Files.writeString(f.y(),original); Files.setLastModifiedTime(f.y(),stamp);
        var grmStamp = Files.getLastModifiedTime(f.grm()); String grm = Files.readString(f.grm());
        Files.writeString(f.grm(),grm.replace("0.5","0.4")); Files.setLastModifiedTime(f.grm(),grmStamp);
        assertThrows(java.io.IOException.class,()->run(f,config)); assertArrayEquals(before,Files.readAllBytes(output));
        config.set("validate_only","true"); config.set("reml_max_iterations","1"); run(f,config);
        assertArrayEquals(before,Files.readAllBytes(output));
    }
    @Test void unestimableMasksSkipOnlyByPolicyAndFitFailuresAreFatal(@TempDir Path root) throws Exception {
        Fixture f = fixture(root); Path output = root.resolve("unestimable.csv"); var config = config(f,output,"skat",true);
        String original = Files.readString(f.y());
        Files.writeString(f.y(),original.lines().findFirst().orElseThrow()+"\nmissing,0"+",".repeat(47)+"\n");
        assertThrows(IllegalArgumentException.class,()->run(f,config));
        assertFalse(Files.exists(Path.of(output+".checkpoint")));
        config.set("unestimable_trait_patterns","skip"); run(f,config);
        assertEquals(1,Files.readAllLines(output).size()); assertTrue(Files.readString(Path.of(output+".sets.tsv")).contains("skipped-unestimable-mask"));
        Files.writeString(f.y(),original);
        var failed = config(f,root.resolve("failed-fit.csv"),"skat",true);
        failed.set("reml_max_iterations","1"); failed.set("set_degenerate","skip"); failed.set("unestimable_trait_patterns","skip");
        var error = assertThrows(IllegalArgumentException.class,()->run(f,failed));
        assertTrue(error.getMessage().contains("converge")); assertFalse(Files.exists(root.resolve("failed-fit.csv")));
    }
    @Test void patternFrequencyIsExplicitAndUsesOnlyObservedSamples(@TempDir Path root) throws Exception {
        Fixture f = fixture(root);
        for (String scope : List.of("aligned","pattern")) {
            Path output = root.resolve(scope+".csv"); var config = config(f,output,"burden",false);
            config.set("frequency_scope",scope); config.set("set_max_maf","0.05"); config.set("set_degenerate","skip");
            run(f,config);
            var lines = Files.readAllLines(output);
            assertEquals(scope.equals("aligned") ? 1 : 5,lines.size());
            if (scope.equals("pattern")) {
                assertTrue(lines.get(1).startsWith("all,pattern,burden,1,43,"));
                assertTrue(lines.get(3).startsWith("all,same-mask-new-trait,burden,1,43,"));
            }
        }
    }
}
