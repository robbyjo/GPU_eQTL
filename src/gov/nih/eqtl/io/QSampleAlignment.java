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
package gov.nih.eqtl.io;

/** Column permutations that put both matrices in one canonical sample order. */
public record QSampleAlignment(
    int[] genotypeColumnOrder,
    int[] expressionColumnOrder,
    String genotypeIdColumn,
    String expressionIdColumn,
    int genotypeReorderedCount,
    int expressionReorderedCount,
    QSampleAlignmentPolicy policy,
    String genotypeIdStripPrefix,
    String expressionIdStripPrefix,
    int genotypeIdsStripped,
    int expressionIdsStripped,
    int genotypeExtraSampleCount,
    int expressionExtraSampleCount) {

    public QSampleAlignment {
        genotypeColumnOrder = genotypeColumnOrder.clone();
        expressionColumnOrder = expressionColumnOrder.clone();
        if (genotypeColumnOrder.length != expressionColumnOrder.length)
            throw new IllegalArgumentException("Aligned predictor and trait sample counts differ");
        if (policy == null)
            throw new IllegalArgumentException("Sample alignment policy is required");
        genotypeIdStripPrefix = genotypeIdStripPrefix == null ? "" : genotypeIdStripPrefix;
        expressionIdStripPrefix = expressionIdStripPrefix == null ? "" : expressionIdStripPrefix;
        if (genotypeIdsStripped < 0 || expressionIdsStripped < 0
            || genotypeExtraSampleCount < 0 || expressionExtraSampleCount < 0)
            throw new IllegalArgumentException("Sample alignment audit counts must not be negative");
    }

    @Override
    public int[] genotypeColumnOrder() {
        return genotypeColumnOrder.clone();
    }

    @Override
    public int[] expressionColumnOrder() {
        return expressionColumnOrder.clone();
    }

    public int sampleCount() {
        return genotypeColumnOrder.length;
    }

    public QSampleAlignment select(boolean[] retainedSamples) {
        if (retainedSamples == null || retainedSamples.length != sampleCount())
            throw new IllegalArgumentException("Sample selection has the wrong length");
        int retained = 0;
        for (boolean keep : retainedSamples)
            if (keep) retained++;
        if (retained == 0)
            throw new IllegalArgumentException("Sample selection removed every aligned sample");
        int[] genotype = new int[retained];
        int[] expression = new int[retained];
        int output = 0;
        for (int sample = 0; sample < retainedSamples.length; sample++) {
            if (!retainedSamples[sample])
                continue;
            genotype[output] = genotypeColumnOrder[sample];
            expression[output] = expressionColumnOrder[sample];
            output++;
        }
        return new QSampleAlignment(genotype, expression, genotypeIdColumn, expressionIdColumn,
            reorderedCount(genotype), reorderedCount(expression), policy,
            genotypeIdStripPrefix, expressionIdStripPrefix,
            genotypeIdsStripped, expressionIdsStripped,
            genotypeExtraSampleCount, expressionExtraSampleCount);
    }

    private static int reorderedCount(int[] order) {
        int count = 0;
        for (int i = 0; i < order.length; i++)
            if (order[i] != i) count++;
        return count;
    }
}
