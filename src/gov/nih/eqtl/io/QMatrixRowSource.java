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

import java.io.IOException;
import java.nio.file.Path;

/** Metadata-first, re-readable matrix rows for bounded-RAM analysis. */
public interface QMatrixRowSource {
    record Metadata(Path path, long rowCount, int columnCount, String[] sampleIds,
        String cacheSignatureTag) {
        public Metadata {
            sampleIds = sampleIds == null ? null : sampleIds.clone();
        }

        @Override
        public String[] sampleIds() {
            return sampleIds == null ? null : sampleIds.clone();
        }
    }

    record Block(long rowOffset, String[] rowIds, double[][] values) {
        public int rowCount() {
            return values.length;
        }

        public int columnCount() {
            return values.length == 0 ? 0 : values[0].length;
        }
    }

    interface BlockReader extends AutoCloseable {
        Block readBlock(int maximumRows) throws IOException;

        @Override
        void close() throws IOException;
    }

    Metadata metadata();

    BlockReader open(int[] columnOrder) throws IOException;
}
