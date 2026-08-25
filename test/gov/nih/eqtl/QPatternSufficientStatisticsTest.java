/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static gov.nih.utils.QStatsUtils.calcStdDevAndStandardize;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

import gov.nih.jama.QRDecomposition;

class QPatternSufficientStatisticsTest {
    @Test
    void paddedTraitAndRStatisticsMatchExplicitQrResidualization() {
        double[][] design = {
            {1, 21, 0}, {1, 34, 1}, {1, 28, 0}, {1, 43, 1},
            {1, 37, 0}, {1, 52, 1}, {1, 46, 0}, {1, 31, 1}
        };
        int[] observed = {0, 1, 2, 4, 5, 6, 7};
        double[] genotype = {0, 1, 2, 0.5, 1, 2, 0, 1.5};
        double[] trait = {2.1, 1.8, 3.2, Double.NaN, 2.7, 4.0, 2.2, 3.5};
        double[][] selectedDesign = selectRows(design, observed);
        QRDecomposition qr = new QRDecomposition(selectedDesign);
        double[][] q = qr.getQ().getArray();

        double[] selectedGenotype = select(genotype, observed);
        double[] selectedTrait = select(trait, observed);
        double[] residualGenotype = residualize(selectedGenotype, q);
        double[] residualTrait = residualize(selectedTrait, q);
        double explicitGenotypeSd = calcStdDevAndStandardize(residualGenotype);
        double traitSd = calcStdDevAndStandardize(residualTrait);
        double explicitCorrelation = QeQTLPreprocessor.correlation(residualGenotype, residualTrait);

        double sumSquares = 0;
        double[] xTransposeGenotype = new double[selectedDesign[0].length];
        for (int sample = 0; sample < observed.length; sample++) {
            double value = selectedGenotype[sample];
            sumSquares += value * value;
            for (int covariate = 0; covariate < xTransposeGenotype.length; covariate++)
                xTransposeGenotype[covariate] += selectedDesign[sample][covariate] * value;
        }
        double residualSumSquares = QPatternSufficientStatistics.residualSumSquares(
            sumSquares, xTransposeGenotype, qr.getR().getArray());
        double sufficientGenotypeSd = QPatternSufficientStatistics.standardDeviation(
            residualSumSquares, observed.length);
        assertEquals(explicitGenotypeSd, sufficientGenotypeSd, 1e-12);

        double[] paddedStandardizedTrait = new double[trait.length];
        for (int sample = 0; sample < observed.length; sample++)
            paddedStandardizedTrait[observed[sample]] = residualTrait[sample];
        double rawDot = 0;
        for (int sample = 0; sample < genotype.length; sample++)
            rawDot += genotype[sample] * paddedStandardizedTrait[sample];
        double sufficientCorrelation = QPatternSufficientStatistics.correlation(
            rawDot, sufficientGenotypeSd, observed.length);
        assertEquals(explicitCorrelation, sufficientCorrelation, 1e-12);

        QeQTLStatistics.Result explicit = QeQTLStatistics.calculate(explicitCorrelation,
            traitSd, explicitGenotypeSd, observed.length - selectedDesign[0].length - 1, 0);
        QeQTLStatistics.Result sufficient = QeQTLStatistics.calculate(sufficientCorrelation,
            traitSd, sufficientGenotypeSd, observed.length - selectedDesign[0].length - 1, 0);
        assertEquals(explicit.rSquared(), sufficient.rSquared(), 1e-12);
        assertEquals(explicit.effect(), sufficient.effect(), 1e-12);
        assertEquals(explicit.tStatistic(), sufficient.tStatistic(), 1e-12);
        assertEquals(explicit.log10P(), sufficient.log10P(), 1e-12);
    }

    @Test
    void invalidOrZeroVarianceIsRejected() {
        assertThrows(IllegalArgumentException.class,
            () -> QPatternSufficientStatistics.standardDeviation(0, 10));
        assertThrows(IllegalArgumentException.class,
            () -> QPatternSufficientStatistics.residualSumSquares(1, new double[] {2},
                new double[][] {{1}}));
    }

    private static double[] residualize(double[] values, double[][] q) {
        double[] coefficients = new double[q[0].length];
        for (int sample = 0; sample < values.length; sample++)
            for (int column = 0; column < coefficients.length; column++)
                coefficients[column] += values[sample] * q[sample][column];
        double[] residuals = values.clone();
        for (int sample = 0; sample < values.length; sample++)
            for (int column = 0; column < coefficients.length; column++)
                residuals[sample] -= q[sample][column] * coefficients[column];
        return residuals;
    }

    private static double[][] selectRows(double[][] values, int[] rows) {
        double[][] result = new double[rows.length][];
        for (int i = 0; i < rows.length; i++)
            result[i] = values[rows[i]].clone();
        return result;
    }

    private static double[] select(double[] values, int[] rows) {
        double[] result = new double[rows.length];
        for (int i = 0; i < rows.length; i++)
            result[i] = values[rows[i]];
        return result;
    }
}
