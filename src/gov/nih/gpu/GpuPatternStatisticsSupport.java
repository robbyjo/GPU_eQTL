/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

/** Portable reference/fallback implementation of compact exact-pattern finalization. */
public final class GpuPatternStatisticsSupport {
    private GpuPatternStatisticsSupport() { }

    public static GpuPatternStatisticsResult calculate(GpuContext context,
        double[] aggregateInputs, int paddedSamples, int activeVariants, int variantCapacity,
        GpuPatternStatisticsPlan plan, boolean meanFill, int patternsPerBatch,
        int workGroupSize) {
        validate(aggregateInputs, paddedSamples, activeVariants, variantCapacity,
            plan, patternsPerBatch, workGroupSize);
        int tripledCapacity = Math.multiplyExact(variantCapacity, 3);
        double[] compact = new double[Math.multiplyExact(
            Math.multiplyExact(plan.patternCount(), activeVariants),
            GpuPatternStatisticsResult.VALUES_PER_CELL)];
        int rank = plan.rank();
        int maximumRows = roundUp(patternsPerBatch * plan.rowsPerPattern(), workGroupSize);
        double[] designMasks = new double[Math.multiplyExact(maximumRows, paddedSamples)];
        double[] upper = new double[Math.multiplyExact(patternsPerBatch, rank * rank)];
        double[] designSums = new double[Math.multiplyExact(patternsPerBatch, rank)];
        int[] observedCounts = new int[patternsPerBatch];
        double[] qTranspose = new double[rank];
        for (int first = 0; first < plan.patternCount(); first += patternsPerBatch) {
            int count = Math.min(patternsPerBatch, plan.patternCount() - first);
            int rowCapacity = roundUp(count * plan.rowsPerPattern(), workGroupSize);
            plan.fillBatch(first, count, paddedSamples, rowCapacity, designMasks,
                upper, designSums, observedCounts);
            long localBytes = (long) (workGroupSize + 1) * (4L * workGroupSize)
                * Double.BYTES;
            double[] products = context.executeDoubleKernel(designMasks, aggregateInputs,
                Math.multiplyExact(rowCapacity, tripledCapacity), localBytes, paddedSamples,
                tripledCapacity, new long[] {roundUp(activeVariants * 3, workGroupSize),
                    rowCapacity}, new long[] {workGroupSize, workGroupSize});
            finalizeBatch(products, compact, first, count, activeVariants, tripledCapacity,
                rank, upper, designSums, observedCounts, meanFill, qTranspose);
        }
        return new GpuPatternStatisticsResult(plan.patternCount(), activeVariants, compact);
    }

    public static void validate(double[] aggregateInputs, int paddedSamples,
        int activeVariants, int variantCapacity, GpuPatternStatisticsPlan plan,
        int patternsPerBatch, int workGroupSize) {
        if (aggregateInputs == null || plan == null || paddedSamples < plan.sampleCount()
            || activeVariants <= 0 || activeVariants > variantCapacity || variantCapacity <= 0
            || aggregateInputs.length != Math.multiplyExact(paddedSamples,
                Math.multiplyExact(variantCapacity, 3))
            || patternsPerBatch <= 0 || workGroupSize <= 0)
            throw new IllegalArgumentException("Invalid compact pattern-statistics operation");
    }

    public static int roundUp(int number, int multiple) {
        int remainder = number % multiple;
        return remainder == 0 ? number : number + multiple - remainder;
    }

    private static void finalizeBatch(double[] products, double[] compact, int firstPattern,
        int patternCount, int variants, int tripledCapacity, int rank, double[] upper,
        double[] designSums, int[] observedCounts, boolean meanFill, double[] qTranspose) {
        int rowsPerPattern = rank + 1;
        for (int pattern = 0; pattern < patternCount; pattern++) {
            int rowBase = pattern * rowsPerPattern;
            for (int variant = 0; variant < variants; variant++) {
                int maskBase = rowBase * tripledCapacity + 3 * variant;
                double called = products[maskBase];
                double sum = products[maskBase + 1];
                double replacement = meanFill && called > 0 ? sum / called : 0.0;
                if (meanFill && called <= 0) replacement = Double.NaN;
                double missing = observedCounts[pattern] - called;
                double sumSquares = products[maskBase + 2]
                    + missing * replacement * replacement;
                double projection = 0;
                for (int row = 0; row < rank; row++) {
                    int productBase = (rowBase + 1 + row) * tripledCapacity + 3 * variant;
                    double xCalled = products[productBase];
                    double xGenotype = products[productBase + 1];
                    double value = xGenotype + replacement
                        * (designSums[pattern * rank + row] - xCalled);
                    for (int previous = 0; previous < row; previous++)
                        value -= upper[(pattern * rank + previous) * rank + row]
                            * qTranspose[previous];
                    value /= upper[(pattern * rank + row) * rank + row];
                    qTranspose[row] = value;
                    projection += value * value;
                }
                int target = ((firstPattern + pattern) * variants + variant)
                    * GpuPatternStatisticsResult.VALUES_PER_CELL;
                compact[target + GpuPatternStatisticsResult.REPLACEMENT] = replacement;
                compact[target + GpuPatternStatisticsResult.RESIDUAL_SUM_SQUARES]
                    = sumSquares - projection;
                compact[target + GpuPatternStatisticsResult.FILLED_SUM_SQUARES] = sumSquares;
                compact[target + GpuPatternStatisticsResult.CALLED_COUNT] = called;
                compact[target + GpuPatternStatisticsResult.DOSAGE_SUM] = sum;
            }
        }
    }
}
