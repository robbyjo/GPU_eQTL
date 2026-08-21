/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.util.Locale;

/** Explicit handling for missing numeric matrix values. */
public enum QMissingValuePolicy {
    ERROR("error"),
    MEAN("mean"),
    ZERO("zero"),
    LOCAL_PATTERN("local-pattern"),
    EXCLUDE_ROW("exclude-row"),
    PATTERN("pattern");

    private final String optionName;

    QMissingValuePolicy(String optionName) {
        this.optionName = optionName;
    }

    public String optionName() {
        return optionName;
    }

    public static QMissingValuePolicy parse(String value, QMissingValuePolicy defaultValue) {
        if (value == null || value.isBlank())
            return defaultValue;
        String normalized = value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
        if (normalized.equals("exclude") || normalized.equals("exclude-variant")
            || normalized.equals("exclude-trait"))
            normalized = "exclude-row";
        if (normalized.equals("dynamic") || normalized.equals("complete-case"))
            normalized = "pattern";
        for (QMissingValuePolicy policy : values())
            if (policy.optionName.equals(normalized))
                return policy;
        if (normalized.equals("flanking") || normalized.equals("neighbor-pattern")
            || normalized.equals("local-haplotype"))
            return LOCAL_PATTERN;
        throw new IllegalArgumentException("Missing-value policy must be error, mean, zero, local-pattern, exclude-row, or pattern: " + value);
    }
}
