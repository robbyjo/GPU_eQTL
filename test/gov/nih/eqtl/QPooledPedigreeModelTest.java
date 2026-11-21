/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;

class QPooledPedigreeModelTest {
    static final Path ROOT=Path.of("test/resources/pooled-reference");
    record Data(double[][] x,double[] y,double[] g,String[] group,String[] ids) { }
    static Data data()throws Exception {
        var lines=Files.readAllLines(ROOT.resolve("data.csv"));int n=lines.size()-1;
        double[][] x=new double[n][7];double[] y=new double[n],g=new double[n];String[] group=new String[n],ids=new String[n];
        for(int i=0;i<n;i++){String[] v=lines.get(i+1).split(",");ids[i]=v[0];group[i]=v[1];y[i]=Double.parseDouble(v[2]);g[i]=Double.parseDouble(v[3]);x[i][0]=1;for(int a=0;a<6;a++)x[i][a+1]=Double.parseDouble(v[a+4]);}
        return new Data(x,y,g,group,ids);
    }
    static QPooledPedigreeModel model(Data d)throws Exception{return new QPooledPedigreeModel(d.x,d.group,"FHS",d.ids,QMixedModel.readPedigree(ROOT.resolve("pedigree.tsv")),1000);}
    @Test void independentDenseRReferencesMatchWithExactObservedMasks()throws Exception {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas","java");
        try {
            Data d=data();var model=model(d);
            for(String line:Files.readAllLines(ROOT.resolve("expected.csv")).subList(1,3)) {
                String[] v=line.split(",");Set<Integer> missing=v[0].equals("missing")?Set.of(2,7,18,41,52,80,93):Set.of();
                int[] mask=java.util.stream.IntStream.range(0,d.y.length).filter(i->!missing.contains(i)).toArray();
                double[] y=Arrays.stream(mask).mapToDouble(i->d.y[i]).toArray(),g=Arrays.stream(mask).mapToDouble(i->d.g[i]).toArray();
                var geometry=model.geometry(mask);var fit=geometry.fit(y);double[] rotated=geometry.rotate(g);var score=fit.score(rotated,fit.products(rotated));
                assertEquals(Integer.parseInt(v[1]),geometry.n);assertEquals(7,geometry.p);
                double[] actual={fit.variance[0],fit.variance[1],fit.variance[2],fit.variance[3],score.beta(),score.se(),score.statistic(),score.log10P(),fit.logLik};
                for(int a=0;a<actual.length;a++)assertEquals(Double.parseDouble(v[a+2]),actual[a],2e-5,"case="+v[0]+" field="+a);
                assertTrue(fit.maximumScore<=1e-6);
                assertSame(geometry,model.geometry(mask));
            }
            assertEquals(2,model.geometryBuilds());
        }finally{if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    @Test void exactIdentityReductionRankAndResponseScalingAreExplicit()throws Exception {
        Data d=data();var model=model(d);int[] mask=java.util.stream.IntStream.range(0,96).filter(i->i>=48||i%4==0).toArray();
        var geometry=model.geometry(mask);assertTrue(geometry.combinedIdentity);
        double[] y=Arrays.stream(mask).mapToDouble(i->d.y[i]).toArray();var fit=geometry.fit(y);
        assertEquals(0,fit.variance[0]);
        double[] g=geometry.rotate(Arrays.stream(mask).mapToDouble(i->d.g[i]).toArray());var score=fit.score(g,fit.products(g));
        for(double scale:new double[]{1e-100,1e-6,1e6,1e100}) {
            double[] yy=Arrays.stream(y).map(v->v*scale).toArray();var scaled=geometry.fit(yy);var actual=scaled.score(g,scaled.products(g));
            assertEquals(score.log10P(),actual.log10P(),2e-6);assertEquals(score.beta(),actual.beta()/scale,2e-6);
        }
        var subset=model.geometry(java.util.stream.IntStream.range(48,96).toArray());assertTrue(subset.p<7);
        assertThrows(QPooledPedigreeModel.Unestimable.class,()->model.geometry(new int[]{1,2}));
    }
    @Test void projectedCovarianceConfoundingAndSingletonGroupsFailExplicitly()throws Exception {
        Data d=data();
        // One sibship plus cohort intercept: its additive kernel becomes half
        // the cohort residual kernel after fixed-effect projection.
        int[] mask={0,1,2,3,48,49,50,51};
        double[][] x=new double[96][2];for(int i=0;i<96;i++){x[i][0]=1;x[i][1]=i>=48?1:0;}
        var m=model(new Data(x,d.y,d.g,d.group,d.ids));
        assertThrows(QPooledPedigreeModel.Unestimable.class,()->m.geometry(mask));
        int[] singleton=java.util.stream.IntStream.range(0,49).toArray();
        assertThrows(QPooledPedigreeModel.Unestimable.class,()->m.geometry(singleton));
    }
}
