/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
package gov.nih.eqtl;

import static gov.nih.utils.QDataUtils.kUndefinedValue;
import static gov.nih.utils.QStatsUtils.calcStdDevAndStandardize;

import gov.nih.eqtl.io.QDelimitedMatrixSource.Block;
import gov.nih.utils.matrix.EMultiplicationMode;
import gov.nih.utils.matrix.QMatrixUtils;

/** Shared block preprocessing for full-memory and streamed analyses. */
public final class QeQTLPreprocessor {
    public record PreparedBlock(long rowOffset, String[] rowIds, double[][] values, double[] standardDeviations) { }

    private QeQTLPreprocessor() { }

    public static PreparedBlock prepare(Block block, double[][] covariateQ, String matrixName) {
        double[][] values = block.values();
        int sampleCount = values[0].length;
        for (int row = 0; row < values.length; row++) {
            if (values[row].length != sampleCount)
                throw new IllegalArgumentException(matrixName + " row '" + block.rowIds()[row] + "' has the wrong sample count");
            for (double value : values[row])
                if (value == kUndefinedValue || !Double.isFinite(value))
                    throw new IllegalArgumentException(matrixName + " row '" + block.rowIds()[row] + "' contains a missing or non-finite value");
        }

        if (covariateQ != null) {
            double[][] residuals = QMatrixUtils.parallelMatrixMultiplication(values, covariateQ, null, 1,
                values.length, sampleCount, EMultiplicationMode.XMinusXYYt);
            for (int row = 0; row < values.length; row++)
                System.arraycopy(residuals[row], 0, values[row], 0, sampleCount);
        }

        double[] standardDeviations = new double[values.length];
        for (int row = 0; row < values.length; row++) {
            standardDeviations[row] = calcStdDevAndStandardize(values[row]);
            if (!(standardDeviations[row] > 0) || !Double.isFinite(standardDeviations[row]))
                throw new IllegalArgumentException(matrixName + " row '" + block.rowIds()[row]
                    + "' has zero or invalid variance after covariate adjustment");
        }
        return new PreparedBlock(block.rowOffset(), block.rowIds(), values, standardDeviations);
    }

    public static double correlation(double[] standardizedX, double[] standardizedY) {
        if (standardizedX.length != standardizedY.length)
            throw new IllegalArgumentException("Sample counts differ");
        double sum = 0;
        for (int i = 0; i < standardizedX.length; i++)
            sum += standardizedX[i] * standardizedY[i];
        return sum / (standardizedX.length - 1);
    }
}
