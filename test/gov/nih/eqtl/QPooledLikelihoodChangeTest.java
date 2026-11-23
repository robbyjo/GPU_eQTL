/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledLikelihoodChangeTest {
    @TempDir Path root;
    @Test void portableBoundaryLineSearchMatchesAnalyticLimit() throws Exception { boundary("java"); }
    @Test void nativeBoundaryLineSearchMatchesAnalyticLimit() throws Exception { boundary("openblas"); }
    @Test void portableTinyLikelihoodChangesMatchAnalyticReml() throws Exception { changes("java"); }
    @Test void nativeTinyLikelihoodChangesMatchAnalyticReml() throws Exception { changes("openblas"); }
    @Test void portableUnequalComponentChangesMatchIndependentDenseR() throws Exception { mixedChanges("java"); }
    @Test void nativeUnequalComponentChangesMatchIndependentDenseR() throws Exception { mixedChanges("openblas"); }
    private void mixedChanges(String engine) throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            var d=QPooledPedigreeModelTest.data();var model=QPooledPedigreeModelTest.model(d);
            var lines=Files.readAllLines(QPooledPedigreeModelTest.ROOT.resolve("relative-changes.csv"));
            for(String line:lines.subList(1,lines.size())) {
                String[] v=line.split(",");Set<Integer> missing=v[0].equals("missing")?Set.of(2,7,18,41,52,80,93):Set.of();
                int[] mask=java.util.stream.IntStream.range(0,d.y().length).filter(i->!missing.contains(i)).toArray();
                var geometry=model.geometry(mask);double[] y=Arrays.stream(mask).mapToDouble(i->d.y()[i]).toArray();
                double[] from=Arrays.stream(new double[]{.2,.03,.15,.34}).map(Math::log).toArray(),to=from.clone();
                for(int j=0;j<to.length;j++)to[j]+=Double.parseDouble(v[j+1]);
                double expected=Double.parseDouble(v[5]);
                assertEquals(expected,geometry.objectiveChange(from,to,y),1e-9,"forward "+line);
                assertEquals(-expected,geometry.objectiveChange(to,from,y),1e-9,"reverse "+line);
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    record Data(double[][] x,double[] y,double[] noise,String[] groups,String[] ids,Path pedigree) { }
    private Data data(int extra) throws Exception {
        int m=24,n=3*m;double[][] x=new double[n][3+extra];double[] y=new double[n],noise=new double[n];String[] groups=new String[n],ids=new String[n];
        StringBuilder ped=new StringBuilder("ID\tSIRE\tDAM\n");Random random=new Random(0);
        for(int i=0;i<n;i++) {
            int g=i/m;groups[i]=new String[]{"FHS","WHI","CARDIA"}[g];ids[i]="synthetic-"+i;
            x[i][0]=1;x[i][1]=g==1?1:0;x[i][2]=g==2?1:0;
            for(int j=3;j<3+extra;j++)x[i][j]=random.nextGaussian();
            noise[i]=random.nextGaussian()*(g==0?1:g==1?(extra>0?0:.7):.3);
            y[i]=noise[i]+(extra>0?.5*x[i][3]:0);
            if(g==0)ped.append(ids[i]).append("\t0\t0\n");
        }
        Path file=root.resolve("ped.tsv");Files.writeString(file,ped);return new Data(x,y,noise,groups,ids,file);
    }
    private QPooledPedigreeModel.Geometry geometry(Data d) throws Exception {
        return new QPooledPedigreeModel(d.x,d.groups,"FHS",d.ids,QMixedModel.readPedigree(d.pedigree),1000).geometry(java.util.stream.IntStream.range(0,d.y.length).toArray());
    }
    private void boundary(String engine) throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            Data d=data(5);var geometry=geometry(d);var fit=geometry.fit(d.y);
            assertTrue(fit.maximumScore<=1e-6);assertTrue(fit.iterations<=1000);assertEquals(8,geometry.p);
            // The noiseless WHI cohort identifies the five shared slopes. At its
            // lower variance bound, the other variances approach within-cohort SSE/23.
            double u=0,information=0;double[] g=new double[72];
            for(int group:new int[]{0,2}) {
                double mean=0;for(int j=0;j<24;j++)mean+=d.noise[group*24+j]/24;
                double sse=0;for(int j=0;j<24;j++)sse+=Math.pow(d.noise[group*24+j]-mean,2);
                double variance=sse/23;int component=geometry.names.indexOf(group==0?"total:FHS":"residual:CARDIA");
                assertEquals(variance,fit.variance[component],2e-7);
                for(int j=0;j<24;j++){double value=j%3-1;g[group*24+j]=value;u+=value*(d.noise[group*24+j]-mean)/variance;information+=value*value/variance;}
            }
            assertTrue(fit.variance[geometry.names.indexOf("residual:WHI")]<1e-8);
            double[] rotated=geometry.rotate(g);var score=fit.score(rotated,fit.products(rotated));
            assertEquals("ok",score.status());assertEquals(u/information,score.beta(),2e-7);
            assertEquals(1/Math.sqrt(information),score.se(),2e-7);assertEquals(u*u/information,score.statistic(),2e-6);
            assertEquals(jdistlib.ChiSquare.cumulative(u*u/information,1,false,true)/Math.log(10),score.log10P(),2e-6);
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    private void changes(String engine) throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            Data d=data(0);var geometry=geometry(d);double[] theta=new double[geometry.names.size()];Arrays.fill(theta,Math.log(.1));
            for(int group=0;group<3;group++) {
                double mean=0,sse=0;for(int j=0;j<24;j++)mean+=d.y[group*24+j]/24;
                for(int j=0;j<24;j++)sse+=Math.pow(d.y[group*24+j]-mean,2);
                int component=geometry.names.indexOf(new String[]{"total:FHS","residual:WHI","residual:CARDIA"}[group]);theta[component]=Math.log(sse/23);
            }
            for(double step:new double[]{-1e-8,1e-8,-.1,.2}) {
                double[] next=theta.clone();for(int a=1;a<next.length;a++)next[a]+=step;
                // At the analytic REML optimum: delta=-df/2*(1-exp(-h)-h).
                // Use its Taylor series for tiny h so the independent oracle does
                // not itself subtract nearly equal full log likelihoods.
                double expected=69*.5*(Math.abs(step)<1e-5?step*step*(.5-step/6+step*step/24):step+Math.expm1(-step));
                double actual=geometry.objectiveChange(theta,next,d.y);
                assertTrue(actual>0);assertEquals(expected,actual,Math.max(2e-21,expected*2e-6));
                assertEquals(-actual,geometry.objectiveChange(next,theta,d.y),Math.max(2e-21,expected*2e-6));
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
}
