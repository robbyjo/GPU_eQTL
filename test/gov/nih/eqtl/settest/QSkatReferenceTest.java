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
    void nearlyCollinearVariantsAndExtremeWeightsRemainFinite() {
        double[][] nearlyCollinear = {
            {1, -1, 0.5, -0.5, 0.25, -0.25},
            {1.0001, -0.9999, 0.5001, -0.4999, 0.2501, -0.2499}
        };
        double[] trait = {0.5, -1, 1.5, -0.25, 0.75, -1.5};
        QSkatReference.Result collinear = QSkatReference.calculate(nearlyCollinear,
            new double[] {1, 1}, trait, 5);
        assertTrue(Double.isFinite(collinear.statistic()));
        assertTrue(Double.isFinite(collinear.pValue()));
        assertTrue(collinear.eigenvalues().length >= 1);

        QSkatReference.Result extreme = QSkatReference.calculate(VARIANTS,
            new double[] {1e-50, 1e50}, TRAIT, 4);
        assertTrue(Double.isFinite(extreme.statistic()));
        assertTrue(Double.isFinite(extreme.pValue()));
        assertEquals("exact-scaled-chi-square", extreme.pValueMethod());
    }

    @Test
    void namedSatterthwaiteFallbackIsFiniteWhenImhofBudgetIsExhausted() {
        QSkatReference.MixturePValue fallback = QSkatReference.mixturePValue(
            3.25, new double[] {0.5, 2.0, 7.0}, 0);
        assertEquals("satterthwaite-fallback", fallback.method());
        assertTrue(fallback.pValue() > 0 && fallback.pValue() <= 1);
        assertTrue(Double.isFinite(fallback.degreesOfFreedom()));
        assertTrue(Double.isFinite(fallback.scale()));
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

    @Test
    void skatOMonteCarloPValuesHonorDeclaredResolution() {
        double[] grid = {0, 0.5, 1};
        QSkatReference.OmnibusResult coarse = QSkatReference.calculateOmnibus(VARIANTS,
            new double[] {1, 1}, TRAIT, 4, grid, 9, 8172);
        QSkatReference.OmnibusResult fine = QSkatReference.calculateOmnibus(VARIANTS,
            new double[] {1, 1}, TRAIT, 4, grid, 999, 8172);
        assertEquals(Math.rint(coarse.adjustedP() * 10), coarse.adjustedP() * 10, 1e-12);
        assertEquals(Math.rint(fine.adjustedP() * 1000), fine.adjustedP() * 1000, 1e-12);
        assertTrue(coarse.adjustedP() >= 0.1);
        assertTrue(fine.adjustedP() >= 0.001);
        assertEquals(9, coarse.simulations());
        assertEquals(999, fine.simulations());
    }
}
