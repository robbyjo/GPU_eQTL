/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

import java.util.List;

/** Pluggable discovery and execution backend. */
public interface GpuBackend {
	String getName();

	String getRuntimeDescription();

	List<GpuDevice> discoverGpuDevices();
}
