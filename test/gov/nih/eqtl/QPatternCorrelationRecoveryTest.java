/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class QPatternCorrelationRecoveryTest {
    private static final double REPORTED = -1.000000061191563;
    private static final int[] OBSERVED = {0, 1, 3, 4, 5, 6, 7, 8};

    @Test
    void reportedOvershootIsRecomputedRatherThanClampedToPerfectCorrelation() {
        double[][] design = design();
        double[] genotype = new double[10];
        double[] trait = new double[10];
        double[] u = {1, -2, 1, 0, 0, 0, 0, 0};
        double[] v = {0, 0, 0, 0, 0, 1, -2, 1};
        double sd = Math.sqrt(6.0 / 7);
        double noise = 0.0003;
        double traitSd = sd * Math.sqrt(1 + noise * noise);
        for (int i = 0; i < OBSERVED.length; i++) {
            int sample = OBSERVED[i];
            genotype[sample] = 10000 + 2 * design[sample][1] + u[i];
            trait[sample] = (u[i] + noise * v[i]) / traitSd;
        }
        // Neither excluded samples nor the sign of an overshoot may alter the mask.
        genotype[2] = Double.NaN;
        trait[9] = 1e100;
        for (double sign : new double[] {-1, 1}) {
            double[] signedTrait = trait.clone();
            for (int sample : OBSERVED) signedTrait[sample] *= sign;
            var result = QPatternCorrelationRecovery.verify(sign * -REPORTED,
                sd * (1 - 6e-8), genotype, signedTrait, 0, OBSERVED, design);
            double expected = sign / Math.sqrt(1 + noise * noise);
            assertEquals(expected, result.correlation(), 2e-15);
            assertEquals(sd, result.genotypeSd(), 2e-13);
            assertEquals(1, result.traitScale(), 2e-15);
            var statistics = QeQTLStatistics.calculate(result.correlation(),
                traitSd * result.traitScale(), result.genotypeSd(), 5, 0);
            assertEquals(1 / (1 + noise * noise), statistics.rSquared(), 3e-15);
            assertEquals(sign, statistics.effect(), 3e-13);
            assertTrue(Double.isFinite(statistics.tStatistic()));
            assertTrue(Double.isFinite(statistics.log10P()));
        }
    }

    @Test
    void selectedMaskAndMeanOrZeroReplacementMatchIndependentScalarProjection() {
        double[][] design = design();
        double[] genotype = {0, 1, 5000, Double.NaN, 1, 2, 0, 1.5, 2, -9000};
        for (double replacement : new double[] {7.5 / 7, 0}) {
            double[] residual = independentResidual(genotype, replacement, design);
            double sumSquares = 0;
            for (double value : residual) sumSquares += value * value;
            double sd = Math.sqrt(sumSquares / 7);
            double[] trait = new double[10];
            for (int i = 0; i < OBSERVED.length; i++)
                trait[OBSERVED[i]] = -residual[i] / sd;
            var result = QPatternCorrelationRecovery.verify(REPORTED, sd,
                genotype, trait, replacement, OBSERVED, design);
            assertEquals(-1, result.correlation(), 2e-15);
            assertEquals(sd, result.genotypeSd(), 2e-15);
            assertEquals(1, result.traitScale(), 2e-15);
        }
    }

    @Test
    void recoveryWindowRetainsOrdinaryValidationAndRejectsDisagreement() {
        for (double value : new double[] {0, 1, -1, 1 + 1e-8, -1 - 1e-8,
            Math.nextUp(1 + 1e-6), -1.01, Double.NaN, Double.POSITIVE_INFINITY})
            assertFalse(QPatternCorrelationRecovery.canVerify(value));
        assertTrue(QPatternCorrelationRecovery.canVerify(REPORTED));
        assertTrue(QPatternCorrelationRecovery.canVerify(1 + 1e-6));
        double[] genotype = {0, 1, 0, 2, 1, 2, 0, 1.5, 2, 0};
        double[] unrelatedTrait = {1, 0, 0, -1, 0, 0, 1, -2, 1, 0};
        assertThrows(IllegalArgumentException.class, () -> QPatternCorrelationRecovery.verify(
            REPORTED, 1, genotype, unrelatedTrait, 0, OBSERVED, design()));
        assertThrows(IllegalArgumentException.class, () -> QPatternCorrelationRecovery.verify(
            REPORTED, 0, genotype, unrelatedTrait, 0, OBSERVED, design()));
    }

    private static double[][] design() {
        double[][] result = new double[10][2];
        for (int sample = 0; sample < result.length; sample++) result[sample][0] = 1;
        for (int i = 0; i < OBSERVED.length; i++) result[OBSERVED[i]][1] = i - 3;
        return result;
    }

    // Closed-form intercept/slope projection, independent of the production QR.
    private static double[] independentResidual(double[] genotype, double replacement,
        double[][] design) {
        double[] result = new double[OBSERVED.length];
        double mean = 0;
        for (int i = 0; i < result.length; i++) {
            double value = genotype[OBSERVED[i]];
            result[i] = Double.isNaN(value) ? replacement : value;
            mean += result[i] / result.length;
        }
        double xx = 0, xy = 0;
        for (int i = 0; i < result.length; i++) {
            double x = design[OBSERVED[i]][1] - 0.5;
            xx += x * x;
            xy += x * (result[i] - mean);
        }
        for (int i = 0; i < result.length; i++)
            result[i] -= mean + (design[OBSERVED[i]][1] - 0.5) * xy / xx;
        return result;
    }
}
