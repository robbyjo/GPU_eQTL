/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import java.util.LinkedHashSet;
import java.util.Set;

import gov.nih.eqtl.QeQTLPreprocessor;
import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QMatrixRowSource.Block;
import gov.nih.jama.QRDecomposition;

/** Reusable FP64 continuous-trait null model for deterministic CPU set-test references. */
public final class QSetTestNullModel {
    private final int sampleCount;
    private final int covariateRank;
    private final int residualDegreesOfFreedom;
    private final double[][] covariateQ;
    private final PreparedBlock traits;

    private QSetTestNullModel(int sampleCount, int covariateRank,
        int residualDegreesOfFreedom, double[][] covariateQ, PreparedBlock traits) {
        this.sampleCount = sampleCount;
        this.covariateRank = covariateRank;
        this.residualDegreesOfFreedom = residualDegreesOfFreedom;
        this.covariateQ = covariateQ;
        this.traits = traits;
    }

    public static QSetTestNullModel create(String[] traitIds, double[][] traitValues,
        double[][] covariateDesign) {
        if (traitIds == null || traitValues == null || covariateDesign == null
            || traitIds.length == 0 || traitIds.length != traitValues.length
            || covariateDesign.length == 0 || covariateDesign[0].length == 0)
            throw new IllegalArgumentException("Traits and a nonempty covariate design are required");
        int samples = covariateDesign.length;
        int columns = covariateDesign[0].length;
        double[][] design = new double[samples][columns];
        for (int sample = 0; sample < samples; sample++) {
            if (covariateDesign[sample] == null || covariateDesign[sample].length != columns)
                throw new IllegalArgumentException("Covariate design is ragged");
            for (int column = 0; column < columns; column++) {
                double value = covariateDesign[sample][column];
                if (!Double.isFinite(value))
                    throw new IllegalArgumentException("Covariate design contains a non-finite value");
                design[sample][column] = value;
            }
        }
        QRDecomposition qr = new QRDecomposition(design);
        if (!qr.isFullRank())
            throw new IllegalArgumentException("Set-test covariate design is rank deficient: rank "
                + qr.getRank() + " for " + qr.getNumColumns() + " columns");
        int dfe = samples - qr.getRank() - 1;
        if (dfe <= 0)
            throw new IllegalArgumentException("Set-test null model has non-positive residual DF: "
                + dfe);
        double[][] copiedTraits = new double[traitValues.length][];
        String[] copiedIds = traitIds.clone();
        Set<String> seenTraitIds = new LinkedHashSet<>();
        for (int trait = 0; trait < traitValues.length; trait++) {
            if (copiedIds[trait] == null || copiedIds[trait].isBlank())
                throw new IllegalArgumentException("Trait IDs must be nonblank");
            copiedIds[trait] = copiedIds[trait].trim();
            if (!seenTraitIds.add(copiedIds[trait]))
                throw new IllegalArgumentException("Duplicate set-test trait ID '"
                    + copiedIds[trait] + "'");
            if (traitValues[trait] == null || traitValues[trait].length != samples)
                throw new IllegalArgumentException("Trait '" + copiedIds[trait]
                    + "' has the wrong sample count");
            copiedTraits[trait] = traitValues[trait].clone();
        }
        double[][] q = qr.getQ().getArray();
        PreparedBlock prepared = QeQTLPreprocessor.prepare(
            new Block(0, copiedIds, copiedTraits), q, "Set-test trait");
        return new QSetTestNullModel(samples, qr.getRank(), dfe, q, prepared);
    }

    public int sampleCount() { return sampleCount; }
    public int covariateRank() { return covariateRank; }
    public int residualDegreesOfFreedom() { return residualDegreesOfFreedom; }
    public String[] traitIds() { return traits.rowIds().clone(); }

    double[][] covariateQ() { return covariateQ; }
    PreparedBlock traits() { return traits; }
}
