/* Copyright 2026 Roby Joehanes; GNU GPL version 3. Public synthetic reference only. */
package gov.nih.eqtl;
import java.nio.file.*;
import java.util.*;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;
import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;

class QPooledBatchProductsTest {
    @org.junit.jupiter.api.Test void heapGuardReclaimsOnceAndStillRejectsInsufficientCapacity() {
        var reclaimed=new java.util.concurrent.atomic.AtomicInteger();
        QPooledBatchProducts.requireHeap(100,()->400,reclaimed::incrementAndGet);
        org.junit.jupiter.api.Assertions.assertEquals(0,reclaimed.get());
        QPooledBatchProducts.requireHeap(100,()->reclaimed.get()==0?150:400,reclaimed::incrementAndGet);
        org.junit.jupiter.api.Assertions.assertEquals(1,reclaimed.get());
        org.junit.jupiter.api.Assertions.assertThrows(IllegalArgumentException.class,()->QPooledBatchProducts.requireHeap(100,()->150,reclaimed::incrementAndGet));
        org.junit.jupiter.api.Assertions.assertEquals(2,reclaimed.get());
    }
    @org.junit.jupiter.api.io.TempDir Path temporary;
    @org.junit.jupiter.api.Test void portableJavaMatchesScalar() throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas","java");
        try{check("cpu");}finally{if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    @org.junit.jupiter.api.Test void cpuMatchesScalar() throws Exception {check("cpu");}
    @org.junit.jupiter.api.Test void cudaMatchesScalar() throws Exception {check("cuda");}
    @org.junit.jupiter.api.Test void openClMatchesScalar() throws Exception {check("opencl");}
    private void check(String backend) throws Exception {
        Path root=QPooledPedigreeModelTest.ROOT;
        var lines=Files.readAllLines(root.resolve("data.csv"));int n=lines.size()-1;
        double[][] x=new double[n][7];double[] y=new double[n],g=new double[n];String[] groups=new String[n],ids=new String[n];
        for(int i=0;i<n;i++){String[] v=lines.get(i+1).split(",");ids[i]=v[0];groups[i]=v[1];y[i]=Double.parseDouble(v[2]);g[i]=Double.parseDouble(v[3]);x[i][0]=1;for(int j=0;j<6;j++)x[i][j+1]=Double.parseDouble(v[j+4]);}
        GpuBackend be=switch(backend){case "cuda"->new CudaGpuBackend();case "opencl"->new JoclGpuBackend();default->new CpuBackend();};
        List<GpuDevice> devices;try{devices=new GpuRuntime(be).getGpuDevices(true,true);}catch(RuntimeException|LinkageError error){org.junit.jupiter.api.Assumptions.abort("Runtime unavailable: "+error);return;}org.junit.jupiter.api.Assumptions.assumeFalse(devices.isEmpty(),"Required backend unavailable");
        long checked=0;double maximum=0;
        try(var pool=new GpuContextPool(new GpuContext[]{devices.get(0).openContext()})) {
            QPooledScoreRunner.compile(pool);
            for(boolean missing:new boolean[]{false,true}) {
                Set<Integer> absent=missing?Set.of(2,7,18,41,52,80,93):Set.of();int[] mask=java.util.stream.IntStream.range(0,n).filter(i->!absent.contains(i)).toArray();
                var model=new QPooledPedigreeModel(x,groups,"FHS",ids,QMixedModel.readPedigree(root.resolve("pedigree.tsv")),1000);var geometry=model.geometry(mask);
                var fits=new QPooledPedigreeModel.Fit[3];
                for(int t=0;t<fits.length;t++){double[] yy=new double[mask.length];for(int i=0;i<yy.length;i++)yy[i]=y[mask[i]]*(t+1)+.2*t*Math.sin(mask[i]*.37);fits[t]=geometry.fit(yy);}
                for(int rows:new int[]{1,17,33}) {
                    double[][] gg=new double[rows][];
                    for(int a=0;a<rows;a++){double[] v=new double[mask.length];for(int i=0;i<v.length;i++)v[i]=a==0?2:a==1?x[mask[i]][1]:a==2?x[mask[i]][1]+1e-7*g[mask[i]]:g[mask[i]]+.01*Math.sin((a+1)*(mask[i]+1));double mean=Arrays.stream(v).average().orElseThrow();for(int i=0;i<v.length;i++)v[i]-=mean;gg[a]=geometry.rotate(v);}
                    var plan=new QPooledBatchProducts(fits,pool);var products=plan.calculate(gg);
                    for(int a=0;a<rows;a++)for(int t=0;t<fits.length;t++) {
                        var expected=fits[t].score(gg[a],fits[t].products(gg[a]));var actual=plan.score(gg[a],products,a,t,Double.NEGATIVE_INFINITY);
                        if(actual==null||!expected.status().equals(actual.status()))throw new AssertionError("Status mismatch");
                        double[] ee={expected.beta(),expected.se(),expected.statistic(),expected.log10P()},aa={actual.beta(),actual.se(),actual.statistic(),actual.log10P()};
                        for(int j=0;j<ee.length;j++){if(Double.isNaN(ee[j])&&Double.isNaN(aa[j]))continue;double delta=Math.abs(ee[j]-aa[j])/Math.max(1,Math.abs(ee[j]));if(!(delta<1e-8))throw new AssertionError("Field "+j+" difference "+delta);maximum=Math.max(maximum,delta);}
                        for(double cutoff:new double[]{0,1,5,100}) {
                            var filtered=plan.score(gg[a],products,a,t,cutoff);
                            if(filtered==null&&expected.statistic()>=cutoff*(1+1e-10)&&expected.status().equals("ok"))throw new AssertionError("Lost retained association");
                        }
                        checked++;
                    }
                }
            }
            checkBoundary(pool);
        }
        System.out.println("BATCH_REFERENCE_PASS backend="+backend+" pairs="+checked+" maximumRelativeError="+maximum);
    }

    private void checkBoundary(GpuContextPool pool)throws Exception {
        int m=24,n=3*m;double[][] x=new double[n][3];double[] y=new double[n];String[] ids=new String[n],groups=new String[n];
        StringBuilder pedigree=new StringBuilder("ID\tSIRE\tDAM\n");
        for(int i=0;i<n;i++){
            int group=i/m,j=i%m;ids[i]="synthetic-"+i;groups[i]=new String[]{"FHS","WHI","CARDIA"}[group];
            x[i][0]=1;x[i][1]=group==1?1:0;x[i][2]=group==2?1:0;
            y[i]=group==1?0:Math.sin(j*.71)*(group==0?1:.3)+Math.cos(j*.23)*.1;
            if(group==0)pedigree.append(ids[i]).append("\t0\t0\n");
        }
        Path file=temporary.resolve("boundary.tsv");Files.writeString(file,pedigree);
        var geometry=new QPooledPedigreeModel(x,groups,"FHS",ids,QMixedModel.readPedigree(file),1000).geometry(java.util.stream.IntStream.range(0,n).toArray());
        var fit=geometry.fit(y);
        // Also exercise the maximum configured trait batch without fitting 256
        // duplicate nulls. All slots retain their own packed score columns.
        var fits=new QPooledPedigreeModel.Fit[256];Arrays.fill(fits,fit);var plan=new QPooledBatchProducts(fits,pool);
        for(double shift:new double[]{0,1,2}) {
            double[] g=new double[n];for(int i=0;i<n;i++)g[i]=i/m==1?shift:(i%3);
            double mean=Arrays.stream(g).average().orElseThrow();for(int i=0;i<n;i++)g[i]-=mean;
            double u=0,information=0;
            for(int group=0;group<3;group++){
                double gy=0,gg=0;for(int j=0;j<m;j++){gy+=y[group*m+j]/m;gg+=g[group*m+j]/m;}
                int component=geometry.names.indexOf(group==0?"total:FHS":group==1?"residual:WHI":"residual:CARDIA");
                double variance=fit.variance[component];
                for(int j=0;j<m;j++){double centered=g[group*m+j]-gg;u+=centered*(y[group*m+j]-gy)/variance;information+=centered*centered/variance;}
            }
            double[] rotated=geometry.rotate(g);var products=plan.calculate(new double[][]{rotated});
            for(int t=0;t<fits.length;t++){
                var actual=plan.score(rotated,products,0,t,Double.NEGATIVE_INFINITY);
                org.junit.jupiter.api.Assertions.assertEquals("ok",actual.status());
                org.junit.jupiter.api.Assertions.assertEquals(u/information,actual.beta(),1e-9);
                org.junit.jupiter.api.Assertions.assertEquals(1/Math.sqrt(information),actual.se(),1e-9);
                org.junit.jupiter.api.Assertions.assertEquals(u*u/information,actual.statistic(),1e-8);
                org.junit.jupiter.api.Assertions.assertEquals(jdistlib.ChiSquare.cumulative(u*u/information,1,false,true)/Math.log(10),actual.log10P(),1e-8);
            }
        }
    }
}

