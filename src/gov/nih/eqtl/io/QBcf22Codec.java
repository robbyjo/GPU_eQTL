/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package gov.nih.eqtl.io;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;

/**
 * HTSJDK's BCF decoder implements the BCF 2 record layout but rejects the
 * current 2.2 version byte. This narrowly extends its version policy. The
 * matrix reader separately restricts accepted records to biallelic diploid
 * GT/DS values, avoiding the variable-ploidy vector-end case that differs
 * between the old and current specifications.
 */
final class QBcf22Codec extends BCF2Codec {
    @Override
    protected void validateVersionCompatibility(BCFVersion supportedVersion, BCFVersion actualVersion) {
        if (actualVersion.getMajorVersion() != 2
            || (actualVersion.getMinorVersion() != 1 && actualVersion.getMinorVersion() != 2)) {
            throw new TribbleException("GPU eQTL supports BCF 2.1 and 2.2; found " + actualVersion);
        }
    }
}
