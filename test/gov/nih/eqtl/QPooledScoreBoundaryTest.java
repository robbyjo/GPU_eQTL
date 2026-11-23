/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledScoreBoundaryTest {
    @TempDir Path root;
    @Test void portableBoundaryScoresMatchIndependentWithinCohortGls() throws Exception { check("java"); }
    @Test void nativeBoundaryScoresMatchIndependentWithinCohortGls() throws Exception { check("openblas"); }
    private void check(String engine) throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            int m=24,n=3*m;double[][] x=new double[n][3];double[] y=new double[n];String[] ids=new String[n],groups=new String[n];
            StringBuilder pedigree=new StringBuilder("ID\tSIRE\tDAM\n");
            for(int i=0;i<n;i++) {
                int group=i/m,j=i%m;ids[i]="synthetic-"+i;groups[i]=new String[]{"FHS","WHI","CARDIA"}[group];
                x[i][0]=1;x[i][1]=group==1?1:0;x[i][2]=group==2?1:0;
                y[i]=group==1?0:Math.sin(j*.71)*(group==0?1:.3)+Math.cos(j*.23)*.1;
                if(group==0)pedigree.append(ids[i]).append("\t0\t0\n");
            }
            Path file=root.resolve("ped.tsv");Files.writeString(file,pedigree);
            var geometry=new QPooledPedigreeModel(x,groups,"FHS",ids,QMixedModel.readPedigree(file),1000).geometry(java.util.stream.IntStream.range(0,n).toArray());
            var fit=geometry.fit(y);assertTrue(fit.maximumScore<=1e-6);
            for(double shift:new double[]{0,1,2}) {
                double[] g=new double[n];for(int i=0;i<n;i++)g[i]=i/m==1?shift:(i%3);
                double mean=Arrays.stream(g).average().orElseThrow();for(int i=0;i<n;i++)g[i]-=mean;
                double u=0,information=0;
                for(int group=0;group<3;group++) {
                    double gy=0,gg=0;for(int j=0;j<m;j++){gy+=y[group*m+j]/m;gg+=g[group*m+j]/m;}
                    int component=geometry.names.indexOf(group==0?"total:FHS":group==1?"residual:WHI":"residual:CARDIA");
                    double variance=fit.variance[component];
                    for(int j=0;j<m;j++){double centered=g[group*m+j]-gg;u+=centered*(y[group*m+j]-gy)/variance;information+=centered*centered/variance;}
                }
                double[] rotated=geometry.rotate(g);var actual=fit.score(rotated,fit.products(rotated));
                assertEquals("ok",actual.status());
                assertEquals(u/information,actual.beta(),1e-9);
                assertEquals(1/Math.sqrt(information),actual.se(),1e-9);
                assertEquals(u*u/information,actual.statistic(),1e-8);
                assertEquals(jdistlib.ChiSquare.cumulative(u*u/information,1,false,true)/Math.log(10),actual.log10P(),1e-8);
            }
            double[] collinear=new double[n];for(int i=0;i<n;i++)collinear[i]=i/m;
            double[] rotated=geometry.rotate(collinear);
            assertEquals("monomorphic-or-collinear",fit.score(rotated,fit.products(rotated)).status());
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
}
