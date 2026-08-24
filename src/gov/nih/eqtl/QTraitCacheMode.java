/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.util.Locale;

/** Residency policy for prepared trait rows in bounded-RAM analysis. */
public enum QTraitCacheMode {
	AUTO("auto"),
	MEMORY("memory"),
	DISK("disk");

	private final String optionName;

	QTraitCacheMode(String optionName) {
		this.optionName = optionName;
	}

	public String optionName() {
		return optionName;
	}

	public static QTraitCacheMode parse(String value) {
		if (value == null || value.isBlank())
			return AUTO;
		String normalized = value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
		for (QTraitCacheMode mode : values())
			if (mode.optionName.equals(normalized))
				return mode;
		throw new IllegalArgumentException("trait_cache must be auto, memory, or disk");
	}
}
