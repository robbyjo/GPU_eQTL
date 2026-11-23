/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledBatchPipelineTest {
    @Test void exactMasksThresholdsAuditResumeAndIdentityMatchOriginal(@TempDir Path root)throws Exception {
        var d=QPooledPedigreeModelTest.data();StringBuilder g=new StringBuilder("ID"),y=new StringBuilder("ID");
        for(String id:d.ids())g.append(',').append(id);
        for(int i=95;i>=0;i--)y.append(',').append(d.ids()[i]);
        for(int a=0;a<17;a++){g.append('\n').append("v"+a);for(int i=0;i<96;i++)g.append(',').append(a==0?2:a==1?d.x()[i][1]:d.g()[i]+.02*Math.sin((i+1)*(a+1)));}
        for(int a=0;a<9;a++){y.append('\n').append("t"+a);for(int i=95;i>=0;i--){y.append(',');if(Set.of(3,4,7).contains(a)&&Set.of(2,7,18,41,52,80,93).contains(i))continue;y.append(d.y()[i]*(1+a/10.)+.02*(a%3)*Math.sin(i*.3));}}
        Files.writeString(root.resolve("g.csv"),g.append('\n'));Files.writeString(root.resolve("y.csv"),y.append('\n'));
        for(String threshold:List.of("1","0.1")) {
            var reference=config(root,"reference"+threshold,threshold,1);run(reference);
            var candidate=config(root,"batch"+threshold,threshold,3);run(candidate);
            Path output=Path.of(candidate.getOutputFilename()),audit=Path.of(output+".pooled.tsv");
            compare(Path.of(reference.getOutputFilename()),output);
            assertEquals(Files.readString(Path.of(reference.getOutputFilename()+".pooled.tsv")),Files.readString(audit));
            byte[] resultBytes=Files.readAllBytes(output),auditBytes=Files.readAllBytes(audit);
            candidate.set("resume","true");run(candidate);
            assertArrayEquals(resultBytes,Files.readAllBytes(output));assertArrayEquals(auditBytes,Files.readAllBytes(audit));
            Path digest=Path.of(output+".pooled-checkpoint/audit/trait-digests.bin");
            byte[] index=Files.readAllBytes(digest);int nameBytes=((index[4]&255)<<8)|(index[5]&255);
            Arrays.fill(index,10+nameBytes,10+nameBytes+32,(byte)0);Files.write(digest,index);
            run(candidate);
            assertArrayEquals(resultBytes,Files.readAllBytes(output));assertArrayEquals(auditBytes,Files.readAllBytes(audit));
            assertTrue(Files.isDirectory(Path.of(output+".pooled-checkpoint.uncommitted/audit")));
            candidate.set("mixed_trait_batch_size","2");assertThrows(java.io.IOException.class,()->run(candidate));
            assertArrayEquals(resultBytes,Files.readAllBytes(output));
            candidate.set("mixed_trait_batch_size","1");assertThrows(java.io.IOException.class,()->run(candidate));
        }
        // Corrupted rotated rows must fail without regeneration or final output.
        Path cache;
        try(var files=Files.list(root.resolve("cache/centered-rotated-v1"))){cache=files.filter(p->p.toString().endsWith(".qraw")).findFirst().orElseThrow();}
        byte[] bytes=Files.readAllBytes(cache);bytes[bytes.length-8]^=1;Files.write(cache,bytes);
        var corrupted=config(root,"corrupted","1",3);
        assertThrows(java.io.IOException.class,()->run(corrupted));
        assertFalse(Files.exists(Path.of(corrupted.getOutputFilename())));
    }
    @Test void invalidBatchOptionsFailBeforeAssociation() {
        for(String value:List.of("0","-1","257","NaN"))assertThrows(IllegalArgumentException.class,()->QeQTLCommandLine.parse(new String[]{"--mixed-trait-batch-size",value}));
        assertThrows(IllegalArgumentException.class,()->QeQTLCommandLine.parse(new String[]{"--mixed-trait-batch-size","2"}));
    }
    private static QeQTLAnalysisConfig config(Path root,String name,String p,int batch)throws Exception {
        return QeQTLCommandLine.parse(new String[]{"--genotype",root.resolve("g.csv").toString(),"--expression",root.resolve("y.csv").toString(),"--output",root.resolve(name+".csv").toString(),
            "--covariates",QPooledPedigreeModelTest.ROOT.resolve("data.csv").toAbsolutePath().toString(),"--genotype-id-column","ID","--expression-id-column","ID","--fixed-covariates","WHI,CARDIA,x,FHS_batch,WHI_x2,CARDIA_year2",
            "--pedigree",QPooledPedigreeModelTest.ROOT.resolve("pedigree.tsv").toAbsolutePath().toString(),"--mixed-model","score","--mixed-null-fitter","pedigree-block","--residual-group-column","Group","--pedigree-group","FHS",
            "--trait-missing","pattern","--predictor-missing","error","--reml-max-iterations","1000","--genotype-block-rows","7","--keep-checkpoints","--cache-dir",root.resolve("cache").toString(),"--threshold","pval",p,"--mixed-trait-batch-size",Integer.toString(batch)}).config();
    }
    private static void compare(Path first,Path second)throws Exception {
        var a=Files.readAllLines(first);var b=Files.readAllLines(second);assertEquals(a.size(),b.size());assertEquals(a.get(0),b.get(0));
        for(int i=1;i<a.size();i++){
            String[] x=a.get(i).split(","),y=b.get(i).split(",");assertEquals(x.length,y.length);
            for(int j:new int[]{0,1,2,3,8,9,10,14})assertEquals(x[j],y[j]);
            for(int j:new int[]{4,5,6,7,11,13}){double expected=Double.parseDouble(x[j]),actual=Double.parseDouble(y[j]);if(Double.isNaN(expected))assertTrue(Double.isNaN(actual));else assertEquals(expected,actual,1e-8*Math.max(1,Math.abs(expected)));}
        }
    }
    private static void run(QeQTLAnalysisConfig config)throws Exception {
        QeQTLAnalysis.config=config;QeQTLAnalysis.profiler=new QeQTLProfiler(false);QeQTLAnalysis.gpuPrecision=GpuPrecision.FP64;
        try(var context=new CpuBackend().discoverGpuDevices().get(0).openContext()){
            QeQTLAnalysis.mContexts=new GpuContext[]{context};
            QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(),config.getGenotypeFilename(),config.getExpressionFilename(),config.getCovariateFilename(),config.getFixedCovariates(),null,"none",Double.NaN,0,true,16,1,"csv");
        }finally{QeQTLAnalysis.config=null;QeQTLAnalysis.mContexts=null;}
    }
}
