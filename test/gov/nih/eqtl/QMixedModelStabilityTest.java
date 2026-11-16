/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.reml.VarianceComponent;
import org.jlinalg.reml.RemlOptions;
import org.junit.jupiter.api.Test;

class QMixedModelStabilityTest {
    @Test void scoreRefinementMatchesIndependentCombinedPedigreeGroupReference() throws Exception {
        var data=QMixedModelTest.data();int n=data.y().length;
        var pedigree=QMixedModel.readPedigree(QMixedModelTest.ROOT.resolve("pedigree.tsv"));
        double[] relationship=new double[n*n];
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)relationship[i*n+j]=pedigree.relationship(data.ids()[i],data.ids()[j]);
        var table=gov.nih.eqtl.io.QCovariateTable.load(QMixedModelTest.ROOT.resolve("data.csv"),new char[]{','},"#");
        var components=List.of(new VarianceComponent("pedigree",n,relationship),
            VarianceComponent.randomIntercept("group",Arrays.asList(table.columnValues("group",null))),
            VarianceComponent.identity("residual",n));
        var expected=QMixedModelTest.model("combined","residual").fitNull(data.y()).fit();
        var actual=QRemlScoreRefinement.fit(data.y(),data.x(),components,RemlOptions.builder()
            .maximumIterations(200).initialVariances(.4,.4,.4).build());
        assertArrayEquals(expected.varianceComponents(),actual.varianceComponents(),3e-6);
        assertArrayEquals(expected.fixedEffects(),actual.fixedEffects(),3e-6);
        assertArrayEquals(expected.fixedEffectCovariance(),actual.fixedEffectCovariance(),3e-6);
        assertEquals(expected.restrictedLogLikelihood(),actual.restrictedLogLikelihood(),3e-6);
    }
    @Test void scoreRefinementHonorsActiveVarianceBounds() {
        int n=32,size=4;double[][] x=new double[n][1];double[] y=new double[n];
        var labels=new ArrayList<Integer>();
        for(int i=0;i<n;i++){x[i][0]=1;labels.add(i/size);y[i]=new double[]{-1,-.5,.5,1}[i%size];}
        var components=List.of(VarianceComponent.randomIntercept("family",labels),VarianceComponent.identity("residual",n));
        var fit=QRemlScoreRefinement.fit(y,x,components,RemlOptions.builder().maximumIterations(100)
            .varianceBounds(.01,10).initialVariances(.3,1).build());
        double g=fit.varianceComponents()[0],e=fit.varianceComponents()[1],v=e+size*g;
        assertEquals(.01,g,1e-15);
        // Each group mean is zero: the genetic score points outside the allowed
        // region, while the residual score must satisfy the original tolerance.
        double sg=-.5*g*size*7/v;
        double se=.5*e*(20/(e*e)-24/e-7/v);
        assertTrue(sg<0);assertTrue(Math.abs(se)<=1e-6);
    }
    @Test void plateauRefinementRequiresTheOriginalScoreToleranceAndBudget() {
        int n=32,size=4,groups=n/size;
        double[][] x=new double[n][1];for(var row:x)row[0]=1/Math.sqrt(n);
        var labels=new ArrayList<Integer>();for(int i=0;i<n;i++)labels.add(i/size);
        var components=List.of(VarianceComponent.randomIntercept("family",labels),VarianceComponent.identity("residual",n));
        var random=new Random(4);double[] family=new double[groups],y=new double[n],means=new double[groups];
        for(int i=0;i<groups;i++)family[i]=.4*random.nextGaussian();
        double mean=0;
        for(int i=0;i<n;i++){y[i]=family[i/size]+random.nextGaussian();mean+=y[i]/n;means[i/size]+=y[i]/size;}
        double within=0,between=0;
        for(int i=0;i<n;i++)within+=Math.pow(y[i]-means[i/size],2);
        for(double m:means)between+=size*Math.pow(m-mean,2);
        double residual=within/(n-groups),genetic=(between/(groups-1)-residual)/size;
        var fit=QRemlScoreRefinement.fit(y,x,components,RemlOptions.builder().maximumIterations(12)
            .varianceBounds(1e-10,1e10).initialVariances(genetic*(1+4e-7),residual*(1-4e-7)).build());
        double g=fit.varianceComponents()[0],e=fit.varianceComponents()[1],v=e+size*g;
        double betweenScore=between/(v*v)-(groups-1)/v;
        double sg=.5*g*size*betweenScore;
        double se=.5*e*(within/(e*e)-(n-groups)/e+betweenScore);
        assertTrue(Math.max(Math.abs(sg),Math.abs(se))<=1e-6,"Independent projected REML score check");
        assertEquals(genetic,g,genetic*2e-6);assertEquals(residual,e,residual*2e-6);
        assertTrue(fit.converged());assertTrue(fit.iterations()<=12);
        assertThrows(IllegalArgumentException.class,()->QRemlScoreRefinement.fit(y,x,components,
            RemlOptions.builder().maximumIterations(1).initialVariances(3,4).build()));
    }
    @Test void retryPreservesMultipleCovariateUnitsAndResponseOffsets() throws Exception {
        var data=QMixedModelTest.data();
        var baseline=QMixedModelTest.model("combined","residual").fitNull(data.y());
        var config=new QeQTLAnalysisConfig(new LinkedHashMap<>());
        config.set("pedigree_file",QMixedModelTest.ROOT.resolve("pedigree.tsv").toString());
        config.set("covariate_random","group");
        double[][] x=Arrays.stream(data.x()).map(double[]::clone).toArray(double[][]::new);
        double[] y=data.y().clone();
        for(int i=0;i<y.length;i++){x[i][1]*=1e6;y[i]=1+1e-10*y[i];}
        var table=gov.nih.eqtl.io.QCovariateTable.load(QMixedModelTest.ROOT.resolve("data.csv"),new char[]{','},"#");
        var shifted=QMixedModel.create(config,x,table,null,data.ids()).fitNull(y);
        assertTrue(shifted.fit().converged());
        for(int i=0;i<baseline.fit().varianceComponents().length;i++)
            assertEquals(baseline.fit().varianceComponents()[i],shifted.fit().varianceComponents()[i]/1e-20,2e-5);
        assertEquals(baseline.fit().fixedEffects()[0],(shifted.fit().fixedEffects()[0]-1)/1e-10,3e-6);
        assertEquals(baseline.fit().fixedEffects()[1],shifted.fit().fixedEffects()[1]/1e-16,3e-6);
        double[] expected=baseline.fit().fixedEffectCovariance(),actual=shifted.fit().fixedEffectCovariance();
        for(int i=0;i<2;i++)for(int j=0;j<2;j++)
            assertEquals(expected[i*2+j],actual[i*2+j]/(1e-20*Math.pow(1e-6,i+j)),2e-5);
        double[] a=baseline.center(data.g()),b=shifted.center(data.g());
        assertEquals(baseline.score(a,baseline.project(a)).statistic(),shifted.score(b,shifted.project(b)).statistic(),2e-5);
        assertTrue(shifted.fit().iterations()<=200);
    }
    @Test void offsetSmallVarianceNullMatchesBalancedAnovaAndDirectGls() {
        int n=32,size=4,groups=n/size;
        double[][] x=new double[n][1];for(var row:x)row[0]=1;
        var labels=new ArrayList<Integer>();for(int i=0;i<n;i++)labels.add(i/size);
        var random=new Random(4);double[] family=new double[groups],y=new double[n];
        for(int i=0;i<groups;i++)family[i]=.4*random.nextGaussian();
        for(int i=0;i<n;i++)y[i]=1+1e-10*(family[i/size]+random.nextGaussian());
        var model=new QMixedModel(x,List.of(VarianceComponent.randomIntercept("family",labels),
            VarianceComponent.identity("residual",n)),200,"residual");
        var state=model.fitNull(y);
        assertTrue(state.fit().converged());
        // Closed-form REML for a balanced one-way random-intercept model.
        // Center first so the reference does not subtract two numbers near one.
        double[] centered=new double[n],means=new double[groups];double mean=0;
        for(int i=0;i<n;i++){centered[i]=y[i]-y[0];mean+=centered[i]/n;means[i/size]+=centered[i]/size;}
        double within=0,between=0;
        for(int i=0;i<n;i++)within+=Math.pow(centered[i]-means[i/size],2);
        for(double value:means)between+=size*Math.pow(value-mean,2);
        double residual=within/(n-groups),genetic=(between/(groups-1)-residual)/size;
        assertTrue(genetic>0);
        assertEquals(genetic,state.fit().varianceComponents()[0],genetic*3e-5);
        assertEquals(residual,state.fit().varianceComponents()[1],residual*3e-5);
        assertEquals(y[0]+mean,state.fit().fixedEffects()[0],4*Math.ulp(1.0));
        double likelihood=-.5*((n-1)*Math.log(2*Math.PI)+(groups-1)*Math.log(residual+size*genetic)
            +(n-groups)*Math.log(residual)+Math.log(n)+within/residual+between/(residual+size*genetic));
        assertEquals(likelihood,state.fit().restrictedLogLikelihood(),2e-6);
        double[] g=new double[n],sums=new double[groups];double gm=0;
        for(int i=0;i<n;i++){g[i]=Math.sin(.7*i);gm+=g[i]/n;sums[i/size]+=g[i];}
        double u=0,information=0;
        for(int i=0;i<n;i++) {
            double py=(centered[i]-mean-genetic/(residual+size*genetic)*size*(means[i/size]-mean))/residual;
            double pg=(g[i]-gm-genetic/(residual+size*genetic)*(sums[i/size]-size*gm))/residual;
            u+=g[i]*py;information+=g[i]*pg;
        }
        double[] cg=state.center(g);
        assertEquals(u*u/information,state.score(cg,state.project(cg)).statistic(),3e-5);
    }
}
