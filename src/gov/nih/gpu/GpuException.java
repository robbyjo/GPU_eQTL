/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu;

/** Signals a GPU discovery, compilation, or execution failure. */
public class GpuException extends RuntimeException {
	private static final long serialVersionUID = 1L;

	public GpuException(String message) {
		super(message);
	}

	public GpuException(String message, Throwable cause) {
		super(message, cause);
	}
}
