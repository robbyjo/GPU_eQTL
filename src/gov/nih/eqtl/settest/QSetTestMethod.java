/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import java.util.Locale;

/** User-visible analysis selection for ordinary single-variant and set tests. */
public enum QSetTestMethod {
    EQTL("eqtl"), BURDEN("burden"), SKAT("skat"), SKAT_O("skat-o");

    private final String optionName;

    QSetTestMethod(String optionName) { this.optionName = optionName; }

    public String optionName() { return optionName; }

    public boolean isSetTest() { return this != EQTL; }

    public static QSetTestMethod parse(String value) {
        if (value == null || value.isBlank()) return EQTL;
        String normalized = value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
        for (QSetTestMethod method : values())
            if (method.optionName.equals(normalized)) return method;
        throw new IllegalArgumentException("analysis must be eqtl, burden, skat, or skat-o");
    }
}
