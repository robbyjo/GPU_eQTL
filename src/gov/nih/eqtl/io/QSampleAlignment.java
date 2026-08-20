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
    int expressionReorderedCount) {

    public QSampleAlignment {
        genotypeColumnOrder = genotypeColumnOrder.clone();
        expressionColumnOrder = expressionColumnOrder.clone();
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
}
