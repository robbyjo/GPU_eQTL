/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import java.util.Locale;

enum QUnestimableTraitPolicy {
    ERROR, SKIP;

    static QUnestimableTraitPolicy parse(String value) {
        if (value == null || value.isBlank()) return ERROR;
        return switch (value.trim().toLowerCase(Locale.ROOT).replace('_', '-')) {
            case "error", "fail" -> ERROR;
            case "skip", "exclude", "audit-skip" -> SKIP;
            default -> throw new IllegalArgumentException(
                "unestimable_trait_patterns must be error or skip");
        };
    }

    String optionName() { return name().toLowerCase(Locale.ROOT); }
}
