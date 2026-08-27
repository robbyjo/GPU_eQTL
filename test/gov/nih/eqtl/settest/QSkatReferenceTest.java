/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

class QSkatReferenceTest {
    private static final double[][] VARIANTS = {
        {1, -1, 0, 0},
        {0, 0, 1, -1}
    };
    private static final double[] TRAIT = {1, -1, 1, -1};

    @Test
    void deterministicSkatMatchesEqualEigenvalueChiSquareReference() {
        QSkatReference.Result result = QSkatReference.calculate(VARIANTS,
            new double[] {1, 1}, TRAIT, 4);
        assertEquals(8, result.statistic(), 1e-12);
        assertArrayEquals(new double[] {2, 2}, result.eigenvalues(), 1e-12);
        assertEquals(Math.exp(-2), result.pValue(), 2e-7);
        assertEquals("exact-equal-eigenvalue-chi-square", result.pValueMethod());
    }

    @Test
    void singleKernelColumnUsesExactScaledChiSquare() {
        QSkatReference.Result result = QSkatReference.calculate(
            new double[][] {VARIANTS[0]}, new double[] {2}, TRAIT, 4);
        assertEquals("exact-scaled-chi-square", result.pValueMethod());
        assertTrue(result.pValue() > 0 && result.pValue() <= 1);
    }

    @Test
    void skatOIsSeededCorrelatedAndReportsAdjustedPInsteadOfMinimum() {
        double[] grid = {0, 0.5, 1};
        QSkatReference.OmnibusResult first = QSkatReference.calculateOmnibus(VARIANTS,
            new double[] {1, 1}, TRAIT, 4, grid, 999, 12345);
        QSkatReference.OmnibusResult second = QSkatReference.calculateOmnibus(VARIANTS,
            new double[] {1, 1}, TRAIT, 4, grid, 999, 12345);
        assertEquals(first.minimumComponentP(), second.minimumComponentP());
        assertEquals(first.adjustedP(), second.adjustedP());
        assertEquals(first.simulations(), second.simulations());
        assertEquals(first.seed(), second.seed());
        for (int component = 0; component < first.components().size(); component++) {
            assertEquals(first.components().get(component).rho(),
                second.components().get(component).rho());
            assertEquals(first.components().get(component).result().statistic(),
                second.components().get(component).result().statistic());
            assertEquals(first.components().get(component).result().pValue(),
                second.components().get(component).result().pValue());
            assertArrayEquals(first.components().get(component).result().eigenvalues(),
                second.components().get(component).result().eigenvalues());
        }
        assertEquals(3, first.components().size());
        assertEquals(QSkatReference.calculate(VARIANTS, new double[] {1, 1}, TRAIT, 4).pValue(),
            first.components().get(0).result().pValue(), 1e-12);
        assertTrue(first.adjustedP() >= 1.0 / 1000);
        assertNotEquals(first.minimumComponentP(), first.adjustedP());
        assertEquals("correlated-parametric-null", first.adjustmentMethod());
    }
}
