/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;
import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;

class QPooledScorePipelineTest {
    @Test void cliExactMasksCisAuditAndResumeMatchDenseReference(@TempDir Path root)throws Exception {
        var d=QPooledPedigreeModelTest.data();StringBuilder g=new StringBuilder("ID"),y=new StringBuilder("ID");
        for(String id:d.ids())g.append(',').append(id);
        for(int i=95;i>=0;i--)y.append(',').append(d.ids()[i]);
        g.append("\n1:100:A:G");for(double value:d.g())g.append(',').append(value);
        g.append("\n1:110:A:C");for(int i=0;i<96;i++)g.append(",2");
        g.append("\n2:100:G:T");for(double value:d.g())g.append(',').append(value);
        for(String name:List.of("full","missing")) {
            y.append('\n').append(name);
            for(int i=95;i>=0;i--){y.append(',');if(name.equals("missing")&&Set.of(2,7,18,41,52,80,93).contains(i))continue;y.append(d.y()[i]);}
        }
        Files.writeString(root.resolve("g.csv"),g.append('\n'));Files.writeString(root.resolve("y.csv"),y.append('\n'));
        Files.writeString(root.resolve("positions.tsv"),"ID\tCHROM\tPOS\nfull\tchr1\t105\nmissing\t1\t105\n");
        List<String> args=new ArrayList<>(List.of("--genotype",root.resolve("g.csv").toString(),"--expression",root.resolve("y.csv").toString(),
            "--output",root.resolve("out.csv").toString(),"--covariates",QPooledPedigreeModelTest.ROOT.resolve("data.csv").toAbsolutePath().toString(),
            "--genotype-id-column","ID","--expression-id-column","ID","--fixed-covariates","WHI,CARDIA,x,FHS_batch,WHI_x2,CARDIA_year2",
            "--pedigree",QPooledPedigreeModelTest.ROOT.resolve("pedigree.tsv").toAbsolutePath().toString(),"--mixed-model","score","--mixed-null-fitter","pedigree-block",
            "--residual-group-column","Group","--pedigree-group","FHS","--trait-missing","pattern","--predictor-missing","error","--reml-max-iterations","1000",
            "--trait-positions",root.resolve("positions.tsv").toString(),"--mixed-cis-window","5","--keep-checkpoints","--genotype-block-rows","1"));
        var config=QeQTLCommandLine.parse(args.toArray(String[]::new)).config();
        config.set("validate_only","true");run(config);
        assertFalse(Files.exists(root.resolve("out.csv")));assertFalse(Files.exists(root.resolve("out.csv.pooled-checkpoint")));
        Files.writeString(root.resolve("positions.tsv"),"ID\tCHROM\tPOS\nfull\t1\t105\n");
        assertThrows(java.io.IOException.class,()->run(config));
        Files.writeString(root.resolve("positions.tsv"),"ID\tCHROM\tPOS\nfull\tchr1\t105\nmissing\t1\t105\n");
        config.set("validate_only","false");run(config);
        Path output=root.resolve("out.csv"),audit=Path.of(output+".pooled.tsv");byte[] original=Files.readAllBytes(output),originalAudit=Files.readAllBytes(audit);
        var lines=Files.readAllLines(output);assertEquals(5,lines.size());
        assertTrue(lines.get(1).contains("\"full\",96,"));assertTrue(lines.get(3).contains("\"missing\",89,"));
        var expected=Files.readAllLines(QPooledPedigreeModelTest.ROOT.resolve("expected.csv"));
        for(int i=0;i<2;i++)assertEquals(Double.parseDouble(expected.get(i+1).split(",")[9]),Double.parseDouble(lines.get(1+i*2).split(",")[7]),2e-5);
        assertTrue(lines.get(2).contains("monomorphic-or-collinear"));assertEquals(3,Files.readAllLines(audit).size());
        assertTrue(Files.readAllLines(audit).get(2).contains("CARDIA=22;FHS=44;WHI=23"));
        config.set("resume","true");run(config);assertArrayEquals(original,Files.readAllBytes(output));assertArrayEquals(originalAudit,Files.readAllBytes(audit));
        Files.writeString(root.resolve("positions.tsv"),"ID\tCHROM\tPOS\nfull\t1\t106\nmissing\t1\t105\n");
        assertThrows(java.io.IOException.class,()->run(config));assertArrayEquals(original,Files.readAllBytes(output));
    }
    private void run(QeQTLAnalysisConfig config)throws Exception {
        QeQTLAnalysis.config=config;QeQTLAnalysis.profiler=new QeQTLProfiler(false);QeQTLAnalysis.gpuPrecision=GpuPrecision.FP64;
        try(var context=new CpuBackend().discoverGpuDevices().get(0).openContext()) {
            QeQTLAnalysis.mContexts=new GpuContext[]{context};
            QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(),config.getGenotypeFilename(),config.getExpressionFilename(),config.getCovariateFilename(),config.getFixedCovariates(),null,"none",Double.NaN,0,true,16,1,"csv");
        }finally{QeQTLAnalysis.config=null;QeQTLAnalysis.mContexts=null;}
    }
    @Test void unsafeOrIgnoredCombinationsAreRejected()throws Exception {
        for(String[] args:List.of(new String[]{"--mixed-null-fitter","pedigree-block"},new String[]{"--residual-group-column","Group"},new String[]{"--mixed-null-fitter","unknown"}))
            assertThrows(IllegalArgumentException.class,()->QeQTLCommandLine.parse(args));
    }
    @Test void cisWindowsIncludeBothEndpointsAndRejectMissingExactIds(@TempDir Path root)throws Exception {
        Path positions=root.resolve("positions.tsv");
        Files.writeString(positions,"ID\tCHROM\tPOS\nt.1\tchr1\t9223372036854775806\n");
        var index=new QPooledScoreRunner.CisIndex(positions,new String[]{"1:9223372036854775807:A:C","1:9223372036854775805:A:T","2:10:A:G"},1);
        assertArrayEquals(new int[]{0,1},index.rows("t.1"));
        assertThrows(java.io.IOException.class,()->index.rows("t"));
    }
    @Test void cudaProductsMatchScalar()throws Exception{check(new CudaGpuBackend());}
    @Test void openClProductsMatchScalar()throws Exception{check(new JoclGpuBackend());}
    @Test void nativeCpuProductsMatchScalar()throws Exception{check(new CpuBackend());}
    private void check(GpuBackend backend)throws Exception {
        List<GpuDevice> devices;try{devices=new GpuRuntime(backend).getGpuDevices(true,true);}catch(RuntimeException|LinkageError error){Assumptions.abort("Runtime unavailable: "+error);return;}
        Assumptions.assumeFalse(devices.isEmpty(),"No eligible device");var d=QPooledPedigreeModelTest.data();
        var geometry=QPooledPedigreeModelTest.model(d).geometry(java.util.stream.IntStream.range(0,96).toArray());var fit=geometry.fit(d.y());
        try(var pool=new GpuContextPool(new GpuContext[]{devices.get(0).openContext()})) {
            QPooledScoreRunner.compile(pool);
            for(int count:new int[]{1,17,33}) {
                double[][] g=new double[count][];for(int a=0;a<count;a++){double[] row=d.g().clone();for(int i=0;i<row.length;i++)row[i]+=.01*Math.sin((i+1)*(a+1));g[a]=geometry.rotate(row);}
                double[][] actual=QPooledScoreRunner.products(g,fit,pool);
                for(int a=0;a<count;a++){assertArrayEquals(fit.products(g[a]),actual[a],2e-10);assertEquals(fit.score(g[a],fit.products(g[a])).log10P(),fit.score(g[a],actual[a]).log10P(),1e-9);}
            }
            System.out.println("POOLED_SCORE_PRODUCTS_PASS "+backend.getName()+" "+devices.get(0).getName());
        }
    }
}
