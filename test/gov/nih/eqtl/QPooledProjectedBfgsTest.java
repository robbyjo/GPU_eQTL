/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledProjectedBfgsTest {
    @TempDir Path root;
    @Test void portableActiveBoundDoesNotContaminateFreeSteps() throws Exception { check("java"); }
    @Test void nativeActiveBoundDoesNotContaminateFreeSteps() throws Exception { check("openblas"); }

    private void check(String engine) throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            int m=64,n=3*m;Random random=new Random(2);
            double[][] x=new double[n][12];double[] y=new double[n];String[] groups=new String[n],ids=new String[n];
            StringBuilder pedigree=new StringBuilder("ID\tSIRE\tDAM\n");
            for(int f=0;f<m/4;f++)pedigree.append("father"+f+"\t0\t0\nmother"+f+"\t0\t0\n");
            for(int i=0;i<n;i++) {
                int group=i/m;groups[i]=new String[]{"FHS","WHI","CARDIA"}[group];ids[i]="synthetic-"+i;
                x[i][0]=1;x[i][1]=group==1?1:0;x[i][2]=group==2?1:0;
                for(int j=3;j<x[i].length;j++)x[i][j]=random.nextGaussian();
                y[i]=random.nextGaussian()*(group==0?.025:group==1?0:1);
                if(group==0)pedigree.append(ids[i]+"\tfather"+(i/4)+"\tmother"+(i/4)+"\n");
            }
            Path file=root.resolve("pedigree.tsv");Files.writeString(file,pedigree);
            var geometry=new QPooledPedigreeModel(x,groups,"FHS",ids,QMixedModel.readPedigree(file),100)
                .geometry(java.util.stream.IntStream.range(0,n).toArray());
            var fit=geometry.fit(y);
            assertTrue(fit.maximumScore<=1e-6);assertTrue(fit.iterations<=100);
            assertEquals(12,geometry.p);assertFalse(geometry.combinedIdentity);
            int fhs=geometry.names.indexOf("residual:FHS"),whi=geometry.names.indexOf("residual:WHI"),cardia=geometry.names.indexOf("residual:CARDIA");
            assertTrue(fit.variance[whi]<1e-8);
            // At the noiseless cohort's lower bound the shared slopes are pinned
            // to zero. Equal-size sibships give independent within/between-family
            // contrasts, an analytic REML limit unrelated to the BFGS machinery.
            double mean=0,within=0,between=0;for(int i=0;i<m;i++)mean+=y[i]/m;
            for(int f=0;f<m/4;f++) {
                double family=0;for(int j=0;j<4;j++)family+=y[4*f+j]/4;
                between+=4*(family-mean)*(family-mean);
                for(int j=0;j<4;j++)within+=Math.pow(y[4*f+j]-family,2);
            }
            double vw=within/48,vb=between/15,additive=Math.max(0,(vb-vw)/2);
            double residual=additive>0?vw-.5*additive:(within+between)/63;
            assertTrue(residual>0,"Fixture requires a positive residual limit");
            assertEquals(additive,fit.variance[0],2e-8);assertEquals(residual,fit.variance[fhs],2e-8);
            mean=0;for(int j=0;j<m;j++)mean+=y[2*m+j]/m;
            double sse=0;for(int j=0;j<m;j++)sse+=Math.pow(y[2*m+j]-mean,2);
            assertEquals(sse/63,fit.variance[cardia],2e-7);
            for(double scale:new double[]{1e-6,1e6}) {
                var scaled=geometry.fit(Arrays.stream(y).map(v->v*scale).toArray());
                assertTrue(scaled.maximumScore<=1e-6);
                for(int a=0;a<fit.variance.length;a++)assertEquals(fit.variance[a],scaled.variance[a]/scale/scale,2e-8);
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
}
