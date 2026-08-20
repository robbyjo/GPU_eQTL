/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import java.util.Locale;

/** Numeric precision used for the GPU matrix product. */
public enum GpuPrecision {
    FP64("fp64", Double.BYTES, true),
    FP32("fp32", Float.BYTES, false);

    private final String optionName;
    private final int bytes;
    private final boolean requiresDoublePrecision;

    GpuPrecision(String optionName, int bytes, boolean requiresDoublePrecision) {
        this.optionName = optionName;
        this.bytes = bytes;
        this.requiresDoublePrecision = requiresDoublePrecision;
    }

    public String optionName() {
        return optionName;
    }

    public int bytes() {
        return bytes;
    }

    public boolean requiresDoublePrecision() {
        return requiresDoublePrecision;
    }

    public static GpuPrecision parse(String value) {
        if (value == null || value.isBlank())
            return FP64;
        return switch (value.trim().toLowerCase(Locale.ROOT)) {
            case "fp64", "double", "64" -> FP64;
            case "fp32", "float", "single", "32" -> FP32;
            default -> throw new IllegalArgumentException(
                "precision must be fp64 or fp32, not '" + value + "'");
        };
    }
}
