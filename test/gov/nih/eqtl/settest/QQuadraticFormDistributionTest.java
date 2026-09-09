/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

class QQuadraticFormDistributionTest {
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
    @Test void unresolvedSpectraFailWithoutMomentApproximation() {
        assertThrows(IllegalArgumentException.class,()->QSkatReference.mixturePValue(10,new double[]{1e-12,1,10}));
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
