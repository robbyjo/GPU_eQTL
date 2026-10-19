/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

import gov.nih.gpu.GpuPrecision;

class QeQTLStatisticsTest {
    @Test
    void reportedGenotypeOuterRoundoffIsClampedBeforeInference() {
        for (double sign : new double[] {1, -1}) {
            double correlation = QeQTLStatistics.validateCorrelation(
                sign * 1.0000000001232048, GpuPrecision.FP64);
            assertEquals(sign, correlation);
            var result = QeQTLStatistics.calculate(correlation, 2, 4, 20, 0);
            assertEquals(1, result.rSquared());
            assertEquals(sign * 0.5, result.effect());
            assertEquals(sign * Double.POSITIVE_INFINITY, result.tStatistic());
            assertEquals(Double.NEGATIVE_INFINITY, result.log10P());
        }
    }

    @Test
    void fp64RoundoffAllowanceHasAFiniteBoundary() {
        double boundary = 1.0 + 1e-8;
        for (double sign : new double[] {1, -1}) {
            assertEquals(sign, QeQTLStatistics.validateCorrelation(
                sign * boundary, GpuPrecision.FP64));
            assertThrows(IllegalArgumentException.class, () -> QeQTLStatistics.validateCorrelation(
                sign * Math.nextUp(boundary), GpuPrecision.FP64));
            assertThrows(IllegalArgumentException.class, () -> QeQTLStatistics.validateCorrelation(
                sign * 1.000001, GpuPrecision.FP64));
        }
    }

    @Test
    void inRangeCorrelationsAreUnchangedAndNonFiniteValuesAreRejected() {
        for (GpuPrecision precision : GpuPrecision.values()) {
            for (double value : new double[] {-1, -0.9999999999, -0.25, 0, 0.25, 0.9999999999, 1})
                assertEquals(value, QeQTLStatistics.validateCorrelation(value, precision));
            for (double value : new double[] {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY})
                assertThrows(IllegalArgumentException.class,
                    () -> QeQTLStatistics.validateCorrelation(value, precision));
        }
    }

    @Test
    void fp32RetainsItsExistingRoundoffAllowance() {
        for (double sign : new double[] {1, -1}) {
            assertEquals(sign, QeQTLStatistics.validateCorrelation(
                sign * (1 + 1e-4), GpuPrecision.FP32));
            assertThrows(IllegalArgumentException.class, () -> QeQTLStatistics.validateCorrelation(
                sign * Math.nextUp(1 + 1e-4), GpuPrecision.FP32));
        }
    }
}
