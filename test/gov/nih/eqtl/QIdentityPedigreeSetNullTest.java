/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class QIdentityPedigreeSetNullTest {
    private static QMixedModel identityModel(double[][] x,int iterations) {
        return new QMixedModel(x,List.of(VarianceComponent.identity("additive-genetic",x.length),
            VarianceComponent.identity("residual",x.length)),iterations,"residual");
    }
    @Test void combinedVarianceAndScoresMatchClosedFormRemlAndDirectGls() {
        int n=32;double[][] x=new double[n][2];double[] y=new double[n],g=new double[n];
        double xm=0,ym=0,gm=0;
        for(int i=0;i<n;i++) {
            x[i][0]=1;x[i][1]=i-13.0;
            y[i]=.4+.03*x[i][1]+.2*Math.sin(.7*i);g[i]=Math.cos(.3*i)+.1*i;
            xm+=x[i][1]/n;ym+=y[i]/n;gm+=g[i]/n;
        }
        double sxx=0,sxy=0,sxg=0;
        for(int i=0;i<n;i++){double dx=x[i][1]-xm;sxx+=dx*dx;sxy+=dx*(y[i]-ym);sxg+=dx*(g[i]-gm);}
        double slope=sxy/sxx,intercept=ym-slope*xm,rss=0,cross=0,gg=0;
        for(int i=0;i<n;i++) {
            double residual=y[i]-intercept-slope*x[i][1],cg=g[i]-gm-sxg/sxx*(x[i][1]-xm);
            rss+=residual*residual;cross+=cg*residual;gg+=cg*cg;
        }
        double variance=rss/(n-2),likelihood=-.5*((n-2)*(Math.log(2*Math.PI)+Math.log(variance)+1)+Math.log(n*sxx));
        var state=identityModel(x,200).fitSetNull(y);
        assertTrue(state.combinedIdentityVariance());assertTrue(state.fit().converged());
        assertEquals(List.of("total-genetic-plus-residual"),state.fit().componentNames());
        assertEquals(variance,state.fit().varianceComponents()[0],variance*2e-7);
        assertArrayEquals(new double[]{intercept,slope},state.fit().fixedEffects(),2e-8);
        assertEquals(likelihood,state.fit().restrictedLogLikelihood(),2e-7);
        double[] centered=state.center(g);
        var score=state.score(centered,state.project(centered));
        assertEquals(cross/gg,score.beta(),2e-8);
        assertEquals(cross*cross/(variance*gg),score.statistic(),2e-7);
        double[] p=state.projection();
        for(int i=0;i<n;i++)for(int j=0;j<n;j++) {
            double expected=((i==j?1.:0.)-1./n-(x[i][1]-xm)*(x[j][1]-xm)/sxx)/variance;
            assertEquals(expected,p[i*n+j],Math.abs(expected)*3e-7+1e-10);
        }
    }
    @Test void identityDecisionUsesTheObservedPrincipalSubmatrix() {
        int n=36;double[][] x=new double[n][2];double[] covariance=new double[n*n],y=new double[n];
        for(int i=0;i<n;i++){x[i][0]=1;x[i][1]=i;covariance[i*n+i]=1;y[i]=.3+.2*Math.sin(i);}
        covariance[32*n+33]=covariance[33*n+32]=.5;
        var full=new QMixedModel(x,List.of(new VarianceComponent("additive-genetic",n,covariance),
            VarianceComponent.identity("residual",n)),200,"residual");
        int[] mask=java.util.stream.IntStream.range(0,32).toArray();
        var subset=full.subset(mask).fitSetNull(Arrays.copyOf(y,32));
        assertTrue(subset.combinedIdentityVariance());
        assertEquals(32,subset.fit().observations());
    }
    @Test void relatedAndAdditionalCovarianceComponentsKeepTheirExistingFits() throws Exception {
        var data=QMixedModelTest.data();
        for(String kind:List.of("pedigree","combined")) {
            var model=QMixedModelTest.model(kind,"residual");
            var expected=model.fitNull(data.y());var actual=model.fitSetNull(data.y());
            assertFalse(actual.combinedIdentityVariance());
            assertArrayEquals(expected.fit().varianceComponents(),actual.fit().varianceComponents(),0);
            assertArrayEquals(expected.projection(),actual.projection(),0);
        }
    }
    @Test void noApproximateMergeOrSilentNonconvergence() {
        int n=32;double[][] x=new double[n][2];double[] y=new double[n],near=new double[n*n];
        for(int i=0;i<n;i++){x[i][0]=1;x[i][1]=i;y[i]=3*i+Math.sin(.7*i);near[i*n+i]=1;}
        near[1]=near[n]=1e-12;
        var model=new QMixedModel(x,List.of(new VarianceComponent("additive-genetic",n,near),
            VarianceComponent.identity("residual",n)),200,"residual");
        assertThrows(IllegalArgumentException.class,()->model.fitSetNull(y));
        assertThrows(IllegalArgumentException.class,()->identityModel(x,200).fitNull(y));
        assertThrows(IllegalArgumentException.class,()->identityModel(x,1).fitSetNull(y));
        assertThrows(IllegalArgumentException.class,()->identityModel(x,200).fitSetNull(new double[n]));
        var wrongSize=new QMixedModel(x,List.of(VarianceComponent.identity("additive-genetic",n-1),
            VarianceComponent.identity("residual",n-1)),200,"residual");
        assertThrows(IllegalArgumentException.class,()->wrongSize.fitSetNull(y));
    }
}
