/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import gov.nih.eqtl.QMissingValuePolicy;

/** Explicit aligned-cohort masks and failure behavior shared by set tests. */
public record QSetTestPolicy(double minMaf, double maxMaf, double minMac, double maxMac,
    QMissingValuePolicy predictorMissing, FailurePolicy absentVariant,
    FailurePolicy degenerateSet) {

    public enum FailurePolicy {
        ERROR, SKIP;

        public static FailurePolicy parse(String value, FailurePolicy defaultValue) {
            if (value == null || value.isBlank()) return defaultValue;
            return switch (value.trim().toLowerCase(java.util.Locale.ROOT)) {
                case "error" -> ERROR;
                case "skip" -> SKIP;
                default -> throw new IllegalArgumentException(
                    "Set-test failure policy must be error or skip");
            };
        }
    }

    public QSetTestPolicy {
        if (!Double.isFinite(minMaf) || !Double.isFinite(maxMaf)
            || minMaf < 0 || maxMaf > 0.5 || minMaf > maxMaf)
            throw new IllegalArgumentException("Set-test MAF bounds must satisfy 0 <= min <= max <= 0.5");
        if (!Double.isFinite(minMac) || minMac < 0 || Double.isNaN(maxMac)
            || maxMac < minMac)
            throw new IllegalArgumentException("Set-test MAC bounds must satisfy 0 <= min <= max");
        if (predictorMissing != QMissingValuePolicy.ERROR
            && predictorMissing != QMissingValuePolicy.MEAN
            && predictorMissing != QMissingValuePolicy.ZERO)
            throw new IllegalArgumentException(
                "Set-test predictor missingness must be error, mean, or zero");
        if (absentVariant == null || degenerateSet == null)
            throw new IllegalArgumentException("Set-test failure policies are required");
    }

    /** No implicit rare-variant filter; callers must explicitly select a maximum MAF/MAC. */
    public static QSetTestPolicy unfiltered(QMissingValuePolicy predictorMissing) {
        return new QSetTestPolicy(0, 0.5, 0, Double.POSITIVE_INFINITY,
            predictorMissing, FailurePolicy.ERROR, FailurePolicy.ERROR);
    }

    public boolean includes(double maf, double mac) {
        return Double.isFinite(maf) && Double.isFinite(mac)
            && maf >= minMaf && maf <= maxMaf && mac >= minMac && mac <= maxMac;
    }
}
