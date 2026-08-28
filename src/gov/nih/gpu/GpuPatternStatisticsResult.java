/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

/** Compact pattern-by-variant statistics returned by a compute context. */
public record GpuPatternStatisticsResult(int patternCount, int variantCount, double[] values) {
    public static final int VALUES_PER_CELL = 5;
    public static final int REPLACEMENT = 0;
    public static final int RESIDUAL_SUM_SQUARES = 1;
    public static final int FILLED_SUM_SQUARES = 2;
    public static final int CALLED_COUNT = 3;
    public static final int DOSAGE_SUM = 4;

    public GpuPatternStatisticsResult {
        if (patternCount <= 0 || variantCount <= 0 || values == null
            || values.length != Math.multiplyExact(
                Math.multiplyExact(patternCount, variantCount), VALUES_PER_CELL))
            throw new IllegalArgumentException("Invalid compact pattern-statistics result");
    }

    public double value(int pattern, int variant, int field) {
        if (pattern < 0 || pattern >= patternCount || variant < 0 || variant >= variantCount
            || field < 0 || field >= VALUES_PER_CELL)
            throw new IndexOutOfBoundsException();
        return values[(pattern * variantCount + variant) * VALUES_PER_CELL + field];
    }
}
