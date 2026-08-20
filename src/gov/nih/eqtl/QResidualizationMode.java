/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.util.Locale;

/** Where fixed-effect covariate projection is applied. */
public enum QResidualizationMode {
    AUTO("auto"), GPU("gpu"), CPU("cpu");

    private final String optionName;

    QResidualizationMode(String optionName) {
        this.optionName = optionName;
    }

    public String optionName() {
        return optionName;
    }

    public static QResidualizationMode parse(String value) {
        if (value == null || value.isBlank())
            return AUTO;
        return switch (value.trim().toLowerCase(Locale.ROOT)) {
            case "auto" -> AUTO;
            case "gpu" -> GPU;
            case "cpu" -> CPU;
            default -> throw new IllegalArgumentException(
                "residualization must be auto, gpu, or cpu, not '" + value + "'");
        };
    }
}
