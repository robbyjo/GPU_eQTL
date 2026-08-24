/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.util.Locale;

/** Explicit rule for relating matrix header samples to covariate rows. */
public enum QSampleAlignmentPolicy {
    STRICT("strict"),
    COVARIATE_SUBSET("covariate-subset");

    private final String optionName;

    QSampleAlignmentPolicy(String optionName) {
        this.optionName = optionName;
    }

    public String optionName() {
        return optionName;
    }

    public static QSampleAlignmentPolicy parse(String value) {
        if (value == null || value.isBlank())
            return STRICT;
        String normalized = value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
        for (QSampleAlignmentPolicy policy : values())
            if (policy.optionName.equals(normalized))
                return policy;
        throw new IllegalArgumentException("sample_alignment must be strict or covariate-subset");
    }
}
