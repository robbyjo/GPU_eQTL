/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import org.junit.jupiter.api.Test;

class QSkatOmnibusIntegrationTest {
    @Test void correlatedRankTwoMatchesIndependentConditionalNormalQuadratureIntoSmallTails() throws Exception {
        var lines = Files.readAllLines(QRelatedSetTest.ROOT.resolve("skat-o-rank-two.csv"));
        double worst = 0;
        for (String line : lines.subList(1,lines.size())) {
            String[] v = line.split(","); double cross = Double.parseDouble(v[0]), a = Double.parseDouble(v[1]);
            var distribution = new QSetScoreDistribution(new double[][]{{1,cross},{cross,.7}});
            var r = distribution.omnibus(new double[]{a*.8,-a*.4},new double[]{0,.25,.5,.75,1},false,10000,17,1e-6,1048576);
            double expected = Double.parseDouble(v[3]); worst = Math.max(worst,Math.abs(r.p()/expected-1));
            assertEquals(Double.parseDouble(v[2]),r.minimumP(),Double.parseDouble(v[2])*1e-9);
            assertEquals(expected,r.p(),expected*2e-6,line);
            assertTrue(r.p() >= r.minimumP() && r.p() <= 5*r.minimumP());
        }
        System.out.println("SKATO_CONDITIONAL_R_REFERENCE worst relative error="+worst);
    }
    @Test void isotropicGridAdjustmentMatchesIndependentChiSquareIntegration() throws Exception {
        var lines = Files.readAllLines(QRelatedSetTest.ROOT.resolve("skat-o-isotropic.csv"));
        for (String line : lines.subList(1,lines.size())) {
            String[] v = line.split(","); int rank = Integer.parseInt(v[0]);
            double[][] covariance = new double[rank][rank]; for (int i = 0; i < rank; i++) covariance[i][i] = 1;
            double[] u = new double[rank]; u[0] = Double.parseDouble(v[2]); u[1] = -u[0];
            var r = new QSetScoreDistribution(covariance).omnibus(u,new double[]{0,1},false,10000,17,1e-6,1048576);
            double expected = Double.parseDouble(v[3]); assertEquals(expected,r.p(),expected*2e-6,line);
            assertTrue(r.p() > r.minimumP()*1.8);
        }
    }
    @Test void underresolvedIntegrationFailsAndScaleDoesNotChangeTail() {
        double[][] covariance = {{1,.25},{.25,.7}}; double[] u = {7.2,-3.6};
        var d = new QSetScoreDistribution(covariance);
        assertThrows(ArithmeticException.class,()->d.omnibus(u,new double[]{0,.5,1},false,99,17,1e-6,256));
        var expected = d.omnibus(u,new double[]{0,.5,1},false,99,17,1e-6,1048576);
        for (double scale : new double[]{1e-80,1e80}) {
            double[][] c = {{scale,.25*scale},{.25*scale,.7*scale}};
            var actual = new QSetScoreDistribution(c).omnibus(new double[]{u[0]*Math.sqrt(scale),u[1]*Math.sqrt(scale)},
                new double[]{0,.5,1},false,99,17,1e-6,1048576);
            assertEquals(expected.p(),actual.p(),expected.p()*1e-9);
        }
    }
    @Test void randomizedRadialRuleIsDeterministicAndReportsConvergence() {
        var d = new QSetScoreDistribution(new double[][]{{1,0,0},{0,.8,0},{0,0,.6}});
        var a = d.omnibus(new double[]{1,-2,1},new double[]{0,.5,1},false,99,19,.02,1048576);
        var c = d.omnibus(new double[]{1,-2,1},new double[]{0,.5,1},false,99,19,.02,1048576);
        assertEquals(a,c); assertTrue(a.method().contains("randomized-quadrature"));
    }
    @Test void nonisotropicRankThreeSmallTailsMatchCompQuadFormAndCubature() throws Exception {
        var lines = Files.readAllLines(QRelatedSetTest.ROOT.resolve("skat-o-rank-three.csv"));
        double worst = 0;
        for (String line : lines.subList(1,lines.size())) {
            String[] v = line.split(","); double score = Double.parseDouble(v[1]), expected = Double.parseDouble(v[2]);
            var d = new QSetScoreDistribution(new double[][]{{.8,0,0},{0,1,0},{0,0,1.2}});
            var a = d.omnibus(new double[]{score,-score,0},new double[]{0,.25,.5,.75,1},false,99,19,.005,1048576);
            assertEquals(Double.parseDouble(v[0]),a.minimumP(),Double.parseDouble(v[0])*3e-5);
            assertEquals(expected,a.p(),expected*.005,line);
            worst = Math.max(worst,Math.abs(a.p()/expected-1));
        }
        System.out.println("SKATO_COMPQUADFORM_CUBATURE worst relative error="+worst);
    }
}
