/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

class QQuadraticFormDistributionTest {
    @Test void weakEigenvaluesDoNotForceAnUnresolvableExtremeTailIntegral() throws Exception {
        // Independent positive convolution of chi-square(1) and epsilon times
        // chi-square(rank-1), evaluated at 70 and 110 decimal digits. The old
        // rank-only truncation bound fails its relative-error check here.
        var lines=java.nio.file.Files.readAllLines(java.nio.file.Path.of("test/resources/tail-remainder-reference.csv"));
        assertEquals(5,lines.size());
        for(String row:lines.subList(1,lines.size())) {
            String[] fields=row.split(",");
            double q=Double.parseDouble(fields[0]),epsilon=Double.parseDouble(fields[1]);
            int rank=Integer.parseInt(fields[2]);double expected=Double.parseDouble(fields[3]);
            double[] lambda=new double[rank];java.util.Arrays.fill(lambda,epsilon);lambda[rank-1]=1;
            var direct=QTiltedQuadraticForm.probability(q,lambda,false);
            assertTrue(direct.estimatedRelativeError()<=2e-11,row);
            assertEquals(expected,direct.probability(),expected*2e-11,row);
            assertEquals(expected,QQuadraticFormDistribution.survival(q,lambda).pValue(),expected*2e-11,row);
            assertEquals(q,QQuadraticFormDistribution.critical(lambda,expected),q*3e-10,row);
            for(double scale:new double[]{1e-100,1e100}) {
                double[] scaled=java.util.Arrays.stream(lambda).map(x->x*scale).toArray();
                assertEquals(expected,QQuadraticFormDistribution.survival(q*scale,scaled).pValue(),expected*2e-11,row);
            }
        }
    }

    @Test void spectrumBoundRetainsTheLawForTwoDominantEigenvalues() {
        double epsilon=1e-6;int rank=160;
        double[] lambda=new double[rank];java.util.Arrays.fill(lambda,epsilon);
        lambda[rank-1]=lambda[rank-2]=1;
        for(double q:new double[]{1,10,100,500,1100}) {
            // Condition on epsilon*chi-square(rank-2); the remaining pair is
            // exponential, giving an independent closed-form convolution.
            double expected=Math.exp(-q/2-(rank-2)/2.0*Math.log1p(-epsilon))
                *jdistlib.ChiSquare.cumulative((1-epsilon)*q/epsilon,rank-2,true,false)
                +jdistlib.ChiSquare.cumulative(q/epsilon,rank-2,false,false);
            var actual=QTiltedQuadraticForm.probability(q,lambda,false);
            assertEquals(expected,actual.probability(),expected*2e-11,"q="+q);
            assertTrue(actual.estimatedRelativeError()<=2e-11);
        }
    }

    @Test void broadSpectraMatchIndependentDaviesAndImhofReferences() throws Exception {
        var lines=java.nio.file.Files.readAllLines(java.nio.file.Path.of("test/resources/difficult-tail-reference.csv"));
        assertTrue(lines.size()>10);
        for(String row:lines.subList(1,lines.size())) {
            String[] fields=row.split(",");
            double q=Double.parseDouble(fields[0]),expected=Double.parseDouble(fields[2]);
            double[] lambda=java.util.Arrays.stream(fields[1].split(";")).mapToDouble(Double::parseDouble).toArray();
            try {
                assertEquals(expected,QSkatReference.mixturePValue(q,lambda).pValue(),Math.max(2e-11,expected*3e-8),row);
            } catch (RuntimeException exception) {
                throw new AssertionError("Independent reference: "+row,exception);
            }
        }
    }
    @Test void tailsAreMonotoneAndQuantilesInvertTheSameDistribution() {
        for (double[] spectrum : new double[][]{{.5,2},{.5,2,7},{.1,.5,1,2,4},{.01,1,10}}) {
            double previous = 1;
            for (double q : new double[]{0,.001,.1,1,5,10,20,50,100}) {
                double p = QSkatReference.mixturePValue(q,spectrum).pValue();
                assertTrue(p <= previous+1e-14,"Nonmonotone p-value at "+q); previous = p;
                if (q > 0 && p < .999 && p > 1e-10) {
                    double critical = QQuadraticFormDistribution.critical(spectrum,p);
                    assertEquals(q,critical,2e-9*Math.max(1,q));
                }
            }
        }
    }
    @Test void broadSpectraUseTheSameLawWithoutMomentApproximation() {
        var tail=QSkatReference.mixturePValue(10,new double[]{1e-12,1,10});
        assertTrue(tail.method().startsWith("tilted-characteristic-inversion"));
        double limiting=QSkatReference.mixturePValue(10,new double[]{1,10}).pValue();
        assertEquals(limiting,tail.pValue(),2e-11);
        // Independent convolution: chi-square(2) is exponential. Conditioning
        // on the small chi-square(1) gives this closed form, even in tiny tails.
        for(double epsilon:new double[]{1e-12,1e-6,.01})
            for(double q:new double[]{.001,.1,1,10,100,1000}) {
                double expected=Math.exp(-q/2)/Math.sqrt(1-epsilon)
                    *jdistlib.ChiSquare.cumulative((1-epsilon)*q/epsilon,1,true,false)
                    +jdistlib.ChiSquare.cumulative(q/epsilon,1,false,false);
                assertEquals(expected,QSkatReference.mixturePValue(q,new double[]{epsilon,1,1}).pValue(),
                    expected*2e-9,"epsilon="+epsilon+" q="+q);
            }
    }
    @Test void difficultSpectrumQuantilesAndScalingRemainConsistent() {
        double[] lambda={1e-12,.01,1};
        double previous=1;
        for(double q:new double[]{.001,.1,1,10,100,1000}) {
            double p=QSkatReference.mixturePValue(q,lambda).pValue();
            assertTrue(p<=previous);previous=p;
            for(double scale:new double[]{1e-100,1e100}) {
                double[] scaled=java.util.Arrays.stream(lambda).map(x->x*scale).toArray();
                assertEquals(p,QSkatReference.mixturePValue(q*scale,scaled).pValue(),p*2e-9);
            }
            if(p>1e-200)assertEquals(q,QQuadraticFormDistribution.critical(lambda,p),2e-8*Math.max(1,q));
        }
        for(double p:new double[]{1e-12,Math.nextDown(1.0)}) {
            double critical=QQuadraticFormDistribution.critical(lambda,p);
            assertEquals(p,QSkatReference.mixturePValue(critical,lambda).pValue(),Math.max(1e-15,p*2e-8));
        }
    }
    @Test void tinyLowerTailIsIntegratedDirectly() {
        // Inside the ellipsoid sum(lambda*z^2)<=q the Gaussian density is
        // between exp(-q/(2*min(lambda))) and one times its value at zero.
        // Its volume therefore bounds F(q) by this leading term with relative
        // difference <=q/(2*min(lambda)), below 5e-13 in these cases.
        double[] lambda={1e-12,1,1};
        for(double q:new double[]{1e-24,1e-28}) {
            double leading=Math.sqrt(2/Math.PI)/3*Math.pow(q,1.5)/Math.sqrt(1e-12);
            var actual=QTiltedQuadraticForm.probability(q,lambda,true);
            assertEquals(leading,actual.probability(),leading*2e-10);
            assertTrue(actual.estimatedRelativeError()<=2e-11);
        }
    }
    @Test void smallRankTwoTailUsesRelativeRatherThanAbsoluteQuadrature() {
        // For lambda=(1,2), the rank-two angular integral is independently
        // evaluated with mpmath at 80 decimal places in the reference script.
        double small = QSkatReference.mixturePValue(400,new double[]{1,2}).pValue();
        assertEquals(2.960962523383312411e-45,small,2.960962523383312411e-45*2e-10);
        assertEquals(2.540961256183794017e-219,QSkatReference.mixturePValue(2000,new double[]{1,2}).pValue(),2.540961256183794017e-219*2e-10);
        assertTrue(small > 0 && small < 1e-40);
        assertEquals(small,QSkatReference.mixturePValue(4e-98,new double[]{1e-100,2e-100}).pValue(),small*2e-11);
    }
}
