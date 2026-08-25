/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import java.util.Locale;

enum QTraitPatternScheduler {
    AUTO, PATTERN, GENOTYPE;

    static QTraitPatternScheduler parse(String value) {
        if (value == null || value.isBlank()) return AUTO;
        return switch (value.trim().toLowerCase(Locale.ROOT).replace('_', '-')) {
            case "auto" -> AUTO;
            case "pattern", "pattern-outer" -> PATTERN;
            case "genotype", "genotype-outer", "variant", "variant-outer" -> GENOTYPE;
            default -> throw new IllegalArgumentException(
                "trait_pattern_scheduler must be auto, pattern, or genotype");
        };
    }

    String optionName() { return name().toLowerCase(Locale.ROOT); }
}
