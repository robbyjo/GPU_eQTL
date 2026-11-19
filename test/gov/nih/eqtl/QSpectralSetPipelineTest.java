/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QSpectralSetPipelineTest {
    private Map<String,String> fixture(Path root)throws Exception {
        StringBuilder g=new StringBuilder("ID"),y=new StringBuilder("ID"),cov=new StringBuilder("ID,x\n"),ped=new StringBuilder("ID\tSIRE\tDAM\n");
        int n=48;for(int i=0;i<n;i++){g.append(",s").append(i);y.append(",s").append(i);cov.append("s").append(i).append(',').append((i-24.)/17).append('\n');ped.append("s").append(i).append("\tf").append(i/4).append("\tm").append(i/4).append('\n');}
        for(int j=0;j<n/4;j++)ped.append("f").append(j).append("\t0\t0\n").append("m").append(j).append("\t0\t0\n");
        for(int variant=0;variant<3;variant++) {
            g.append("\n1:").append(100+variant).append(":A:G");
            for(int i=0;i<n;i++)g.append(',').append((i+variant*3)%7==0?1:0);
        }
        for(int trait=0;trait<3;trait++) {
            y.append("\nt").append(trait);
            for(int i=0;i<n;i++) {
                y.append(',');if(trait>0&&Set.of(3,7,11,20,29).contains(i))continue;
                double value=.3+.02*i+.1*Math.sin(i)+.25*Math.sin(i/4.);
                y.append(trait==2?value*value:value);
            }
        }
        Files.writeString(root.resolve("g.csv"),g.append('\n'));Files.writeString(root.resolve("y.csv"),y.append('\n'));
        Files.writeString(root.resolve("cov.csv"),cov);Files.writeString(root.resolve("ped.tsv"),ped);
        Files.writeString(root.resolve("sets.tsv"),"SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\tWEIGHT\nall\t1:100:A:G\tA\tG\tG\t1\nall\t1:101:A:G\tA\tG\tG\t.7\nall\t1:102:A:G\tA\tG\tG\t1.3\nsingle\t1:100:A:G\tA\tG\tG\t1\n");
        Map<String,String> config=new LinkedHashMap<>();
        config.put("genotype_file",root.resolve("g.csv").toString());config.put("expression_file",root.resolve("y.csv").toString());
        config.put("covariate_file",root.resolve("cov.csv").toString());config.put("covariate_fixed","x");
        config.put("pedigree_file",root.resolve("ped.tsv").toString());config.put("mixed_model","score");
        config.put("genotype_id_column","ID");config.put("expression_id_column","ID");
        config.put("variant_sets",root.resolve("sets.tsv").toString());config.put("analysis","skat");config.put("set_block_size","1");
        config.put("trait_missing","pattern");config.put("predictor_missing","error");config.put("keep_checkpoints","true");config.put("reml_max_iterations","1000");
        return config;
    }
    private void run(QeQTLAnalysisConfig config,String engine)throws Exception {
        String old=System.getProperty(CpuBackend.BLAS_PROPERTY);System.setProperty(CpuBackend.BLAS_PROPERTY,engine);
        QeQTLAnalysis.config=config;QeQTLAnalysis.profiler=new QeQTLProfiler(false);QeQTLAnalysis.gpuPrecision=GpuPrecision.FP64;
        try(var context=new CpuBackend().discoverGpuDevices().get(0).openContext()) {
            QeQTLAnalysis.mContexts=new GpuContext[]{context};
            QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(),config.get("genotype_file"),config.get("expression_file"),config.get("covariate_file"),
                new String[]{"x"},null,"none",Double.NaN,0,true,16,1,"csv");
        }finally {
            if(old==null)System.clearProperty(CpuBackend.BLAS_PROPERTY);else System.setProperty(CpuBackend.BLAS_PROPERTY,old);
            QeQTLAnalysis.config=null;QeQTLAnalysis.mContexts=null;
        }
    }
    @Test void pipelineMatchesLegacyNativeAndPortableAndResumesOnlyCompatibleFits(@TempDir Path root)throws Exception {
        var base=fixture(root);Path legacy=root.resolve("legacy.csv");base.put("output_file",legacy.toString());
        var config=new QeQTLAnalysisConfig(base);run(config,"java");List<String> reference=Files.readAllLines(legacy);
        assertEquals(7,reference.size());
        for(String engine:List.of("java","openblas")) {
            Path out=root.resolve("spectral-"+engine+".csv");base.put("output_file",out.toString());base.put("set_null_fitter","spectral");
            config=new QeQTLAnalysisConfig(base);run(config,engine);List<String> actual=Files.readAllLines(out);
            assertEquals(reference.size(),actual.size());
            for(int i=1;i<actual.size();i++) {
                String[] a=actual.get(i).split(","),b=reference.get(i).split(",");
                assertArrayEquals(Arrays.copyOf(b,6),Arrays.copyOf(a,6));
                assertEquals(Double.parseDouble(b[6]),Double.parseDouble(a[6]),Math.abs(Double.parseDouble(b[6]))*2e-5+1e-9);
                assertEquals(Double.parseDouble(b[11]),Double.parseDouble(a[11]),1e-4);
                assertTrue(actual.get(i).contains("fit-method=spectral-reml"));
            }
            byte[] before=Files.readAllBytes(out),audit=Files.readAllBytes(Path.of(out+".sets.tsv"));
            config.set("resume","true");run(config,engine);assertArrayEquals(before,Files.readAllBytes(out));assertArrayEquals(audit,Files.readAllBytes(Path.of(out+".sets.tsv")));
            config.set("set_null_fitter","legacy");final var incompatible=config;
            assertThrows(java.io.IOException.class,()->run(incompatible,engine));
            assertArrayEquals(before,Files.readAllBytes(out));
        }
    }
    @Test void unsupportedSpectralRequestsFailBeforeAnalysis() {
        var config=new QeQTLAnalysisConfig(Map.of("analysis","skat","mixed_model","none","set_null_fitter","spectral"));
        assertThrows(IllegalArgumentException.class,config::requireSupportedAnalysisModel);
    }
}
