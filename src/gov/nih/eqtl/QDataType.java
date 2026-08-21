/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.util.Locale;

/** Declared biological role of one numeric analysis matrix. */
public enum QDataType {
    GENOTYPE("genotype"),
    EXPRESSION("expression"),
    METHYLATION("methylation"),
    PROTEOMICS("proteomics"),
    CONTINUOUS("continuous");

    private final String optionName;

    QDataType(String optionName) {
        this.optionName = optionName;
    }

    public String optionName() {
        return optionName;
    }

    public boolean isGenotype() {
        return this == GENOTYPE;
    }

    public static QDataType parse(String value, QDataType defaultValue) {
        if (value == null || value.isBlank())
            return defaultValue;
        String normalized = value.trim().toLowerCase(Locale.ROOT);
        if (normalized.equals("protein") || normalized.equals("proteomic"))
            normalized = "proteomics";
        if (normalized.equals("dna-methylation") || normalized.equals("dnm"))
            normalized = "methylation";
        for (QDataType type : values())
            if (type.optionName.equals(normalized))
                return type;
        throw new IllegalArgumentException("Data type must be genotype, expression, methylation, proteomics, or continuous: " + value);
    }
}
