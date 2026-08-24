/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.IOException;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;

/** Read-only FP64 rows after alignment, residualization, and standardization. */
public interface QPreparedMatrix {
    String signature();
    int sampleCount();
    long rowCount();
    PreparedBlock readBlock(long rowOffset, int maximumRows) throws IOException;
}
