/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

/** Exact complete-predictor statistics for a trait-specific sample/covariate subset. */
final class QPatternSufficientStatistics {
    private static final double NEGATIVE_ROUNDOFF_TOLERANCE = 1e-10;

    private QPatternSufficientStatistics() { }

    /**
     * Return {@code g' (I - QQ') g} from {@code g'g}, {@code X'g}, and the
     * upper-triangular {@code R} in the thin decomposition {@code X = QR}.
     */
    static double residualSumSquares(double sumSquares, double[] xTransposeGenotype,
        double[][] upperR) {
		return residualSumSquares(sumSquares, xTransposeGenotype, upperR,
			new double[xTransposeGenotype == null ? 0 : xTransposeGenotype.length]);
	}

	static double residualSumSquares(double sumSquares, double[] xTransposeGenotype,
		double[][] upperR, double[] qTransposeWorkspace) {
        if (!Double.isFinite(sumSquares) || sumSquares < 0
            || xTransposeGenotype == null || upperR == null
            || xTransposeGenotype.length == 0 || upperR.length != xTransposeGenotype.length
			|| qTransposeWorkspace == null
			|| qTransposeWorkspace.length < xTransposeGenotype.length)
            throw new IllegalArgumentException("Invalid pattern sufficient statistics");
        int rank = xTransposeGenotype.length;
        double projectionSumSquares = 0;
        for (int row = 0; row < rank; row++) {
            if (upperR[row] == null || upperR[row].length != rank
                || !Double.isFinite(upperR[row][row]) || upperR[row][row] == 0)
                throw new IllegalArgumentException("Invalid or singular pattern R matrix");
            double value = xTransposeGenotype[row];
            if (!Double.isFinite(value))
                throw new IllegalArgumentException("Non-finite X-transpose-genotype statistic");
            for (int previous = 0; previous < row; previous++)
                value -= upperR[previous][row] * qTransposeWorkspace[previous];
            value /= upperR[row][row];
			qTransposeWorkspace[row] = value;
            projectionSumSquares += value * value;
        }
		return validateResidualSumSquares(sumSquares, sumSquares - projectionSumSquares);
	}

	static double validateResidualSumSquares(double sumSquares, double residual) {
		if (!Double.isFinite(sumSquares) || sumSquares < 0 || !Double.isFinite(residual))
			throw new IllegalArgumentException("Invalid genotype residual sufficient statistics");
		double tolerance = Math.max(1.0, sumSquares) * NEGATIVE_ROUNDOFF_TOLERANCE;
		if (residual < -tolerance)
			throw new IllegalArgumentException(
				"Invalid negative genotype residual sum of squares " + residual);
		return Math.max(0, residual);
	}

    static double standardDeviation(double residualSumSquares, int sampleCount) {
        if (!Double.isFinite(residualSumSquares) || residualSumSquares <= 0 || sampleCount <= 1)
            throw new IllegalArgumentException("Pattern genotype variance is zero or invalid");
        return Math.sqrt(residualSumSquares / (sampleCount - 1.0));
    }

    /** Solve {@code (X'X) beta = X'y} through the QR upper triangle {@code R}. */
    static void leastSquaresCoefficients(double[] xTransposeResponse, double[][] upperR,
        double[] qTransposeWorkspace, double[] coefficients) {
        if (xTransposeResponse == null || upperR == null || qTransposeWorkspace == null
            || coefficients == null || xTransposeResponse.length == 0
            || upperR.length != xTransposeResponse.length
            || qTransposeWorkspace.length < xTransposeResponse.length
            || coefficients.length < xTransposeResponse.length)
            throw new IllegalArgumentException("Invalid pattern least-squares dimensions");
        int rank = xTransposeResponse.length;
        for (int row = 0; row < rank; row++) {
            if (upperR[row] == null || upperR[row].length != rank
                || !Double.isFinite(upperR[row][row]) || upperR[row][row] == 0)
                throw new IllegalArgumentException("Invalid or singular pattern R matrix");
            double value = xTransposeResponse[row];
            if (!Double.isFinite(value))
                throw new IllegalArgumentException("Non-finite X-transpose-response statistic");
            for (int previous = 0; previous < row; previous++)
                value -= upperR[previous][row] * qTransposeWorkspace[previous];
            qTransposeWorkspace[row] = value / upperR[row][row];
        }
        for (int row = rank - 1; row >= 0; row--) {
            double value = qTransposeWorkspace[row];
            for (int column = row + 1; column < rank; column++)
                value -= upperR[row][column] * coefficients[column];
            coefficients[row] = value / upperR[row][row];
        }
    }

    /** A standardized padded trait makes its raw genotype dot product an exact residual cross-product. */
    static double correlation(double rawGenotypeDotStandardizedTrait,
        double genotypeStandardDeviation, int sampleCount) {
        if (!Double.isFinite(rawGenotypeDotStandardizedTrait)
            || !Double.isFinite(genotypeStandardDeviation) || genotypeStandardDeviation <= 0
            || sampleCount <= 1)
            throw new IllegalArgumentException("Invalid pattern correlation statistics");
        return rawGenotypeDotStandardizedTrait
            / (genotypeStandardDeviation * (sampleCount - 1.0));
    }
}
