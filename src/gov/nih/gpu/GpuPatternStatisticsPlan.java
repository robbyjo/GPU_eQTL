/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import java.util.Arrays;

/** Immutable host description of the exact sample masks and covariate models used by a
 * genotype-outer pattern-statistics operation. */
public final class GpuPatternStatisticsPlan {
    private final double[][] design;
    private final int[] patternIds;
    private final int[][] observedSamples;
    private final double[][] designSums;
    private final double[][][] upperTriangles;
    private final int rank;

    public GpuPatternStatisticsPlan(double[][] design, int[] patternIds,
        int[][] observedSamples, double[][] designSums, double[][][] upperTriangles) {
        if (design == null || design.length == 0 || design[0] == null
            || design[0].length == 0 || patternIds == null || observedSamples == null
            || designSums == null || upperTriangles == null
            || patternIds.length == 0 || observedSamples.length != patternIds.length
            || designSums.length != patternIds.length
            || upperTriangles.length != patternIds.length)
            throw new IllegalArgumentException("Invalid pattern-statistics plan dimensions");
        rank = design[0].length;
        this.design = design;
        this.patternIds = patternIds.clone();
        this.observedSamples = observedSamples.clone();
        this.designSums = designSums.clone();
        this.upperTriangles = upperTriangles.clone();
        for (double[] row : design)
            if (row == null || row.length != rank)
                throw new IllegalArgumentException("Pattern-statistics design is ragged");
        for (int pattern = 0; pattern < patternIds.length; pattern++) {
            int[] observed = observedSamples[pattern];
            double[] sums = designSums[pattern];
            double[][] upper = upperTriangles[pattern];
            if (observed == null || observed.length == 0 || sums == null || sums.length != rank
                || upper == null || upper.length != rank)
                throw new IllegalArgumentException("Invalid pattern-statistics model " + pattern);
            int previous = -1;
            for (int sample : observed) {
                if (sample <= previous || sample < 0 || sample >= design.length)
                    throw new IllegalArgumentException("Pattern sample indices must be sorted and unique");
                previous = sample;
            }
            for (double[] row : upper)
                if (row == null || row.length != rank)
                    throw new IllegalArgumentException("Pattern R matrix is ragged");
        }
    }

    public int sampleCount() { return design.length; }
    public int rank() { return rank; }
    public int patternCount() { return patternIds.length; }
    public int patternId(int planIndex) { return patternIds[planIndex]; }
    public int rowsPerPattern() { return rank + 1; }

    /** Fill reusable batch arrays. Model arrays use local batch order. */
    public int fillBatch(int firstPattern, int maximumPatterns, int paddedSamples,
        int rowCapacity, double[] designMasks, double[] flattenedUpperTriangles,
        double[] flattenedDesignSums, int[] observedCounts) {
        if (firstPattern < 0 || firstPattern >= patternCount() || maximumPatterns <= 0
            || paddedSamples < sampleCount())
            throw new IllegalArgumentException("Invalid pattern batch bounds");
        int count = Math.min(maximumPatterns, patternCount() - firstPattern);
        int activeRows = Math.multiplyExact(count, rowsPerPattern());
        if (rowCapacity < activeRows
            || designMasks.length < Math.multiplyExact(rowCapacity, paddedSamples)
            || flattenedUpperTriangles.length < Math.multiplyExact(count, rank * rank)
            || flattenedDesignSums.length < Math.multiplyExact(count, rank)
            || observedCounts.length < count)
            throw new IllegalArgumentException("Pattern batch workspace is too small");
        Arrays.fill(designMasks, 0.0);
        for (int local = 0; local < count; local++) {
            int planIndex = firstPattern + local;
            int rowBase = local * rowsPerPattern();
            for (int sample : observedSamples[planIndex]) {
                designMasks[rowBase * paddedSamples + sample] = 1.0;
                for (int column = 0; column < rank; column++)
                    designMasks[(rowBase + 1 + column) * paddedSamples + sample]
                        = design[sample][column];
            }
            observedCounts[local] = observedSamples[planIndex].length;
            System.arraycopy(designSums[planIndex], 0, flattenedDesignSums,
                local * rank, rank);
            double[][] upper = upperTriangles[planIndex];
            for (int row = 0; row < rank; row++)
                System.arraycopy(upper[row], 0, flattenedUpperTriangles,
                    (local * rank + row) * rank, rank);
        }
        return count;
    }
}
