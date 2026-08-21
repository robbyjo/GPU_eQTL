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

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

import com.csvreader.CsvReader;

import gov.nih.utils.QDataUtils;

/**
 * Metadata-first, re-readable delimited numeric matrix. The first column is a
 * row identifier and the first row contains sample identifiers.
 */
public final class QDelimitedMatrixSource implements QMatrixRowSource {
    private static final int BUFFER_SIZE = 16 * 1024 * 1024;

    public final class BlockReader implements QMatrixRowSource.BlockReader {
        private final CsvReader reader;
        private final int[] columnOrder;
        private long nextRowOffset;
        private boolean closed;

        private BlockReader(int[] requestedColumnOrder) throws IOException {
            reader = createCsvReader();
            if (!reader.readRecord()) {
                close();
                throw new IOException("Matrix file is empty: " + path);
            }
            columnOrder = normalizeColumnOrder(requestedColumnOrder, metadata.columnCount());
        }

        public QMatrixRowSource.Block readBlock(int maximumRows) throws IOException {
            if (maximumRows <= 0)
                throw new IllegalArgumentException("maximumRows must be positive");
            if (closed)
                throw new IOException("Matrix reader is closed");

            List<String> rowIds = new ArrayList<>(maximumRows);
            List<double[]> rows = new ArrayList<>(maximumRows);
            long blockOffset = nextRowOffset;
            while (rows.size() < maximumRows && reader.readRecord()) {
                String[] tokens = reader.getValues();
                long lineNumber = nextRowOffset + 2;
                validateFieldCount(tokens, lineNumber);
                String rowId = tokens[0].trim();
                if (rowId.isEmpty())
                    throw new IOException("Blank row identifier at line " + lineNumber + " in " + path);

                double[] values = new double[columnOrder.length];
                for (int i = 0; i < columnOrder.length; i++) {
                    String value = tokens[columnOrder[i] + 1].trim();
                    if (isMissingToken(value)) {
                        values[i] = QDataUtils.kUndefinedValue;
                    } else {
                        try {
                            values[i] = Double.parseDouble(value);
                            if (!Double.isFinite(values[i]))
                                throw new NumberFormatException("not finite");
                        } catch (NumberFormatException e) {
                            throw new IOException("Invalid numeric value at line " + lineNumber
                                + ", matrix column " + (columnOrder[i] + 2) + " in " + path, e);
                        }
                    }
                }
                rowIds.add(rowId);
                rows.add(values);
                nextRowOffset++;
            }
            if (rows.isEmpty())
                return null;
            return new QMatrixRowSource.Block(blockOffset, rowIds.toArray(String[]::new),
                rows.toArray(double[][]::new));
        }

        @Override
        public void close() {
            if (!closed) {
                closed = true;
                reader.close();
            }
        }
    }

    private final Path path;
    private final char[] delimiters;
    private final String commentMarkers;
    private final QMatrixRowSource.Metadata metadata;

    public QDelimitedMatrixSource(Path path, char[] delimiters, String commentMarkers) throws IOException {
        this.path = path.toAbsolutePath().normalize();
        this.delimiters = delimiters.clone();
        this.commentMarkers = commentMarkers;
        metadata = scanMetadata();
    }

    @Override
    public QMatrixRowSource.Metadata metadata() {
        return metadata;
    }

    @Override
    public BlockReader open(int[] columnOrder) throws IOException {
        return new BlockReader(columnOrder);
    }

    private QMatrixRowSource.Metadata scanMetadata() throws IOException {
        Set<String> rowIds = new HashSet<>();
        long rowCount = 0;
        CsvReader reader = null;
        try {
            reader = createCsvReader();
            if (!reader.readRecord())
                throw new IOException("Matrix file is empty: " + path);
            String[] header = reader.getValues();
            if (header.length < 2)
                throw new IOException("Matrix header must contain a row-ID column and at least one sample: " + path);
            String[] sampleIds = Arrays.copyOfRange(header, 1, header.length);
            validateUniqueIdentifiers(sampleIds, "sample", path);

            while (reader.readRecord()) {
                String[] tokens = reader.getValues();
                long lineNumber = rowCount + 2;
                if (tokens.length != header.length)
                    throw new IOException("Expected " + header.length + " fields but found " + tokens.length
                        + " at line " + lineNumber + " in " + path);
                String rowId = tokens[0].trim();
                if (rowId.isEmpty())
                    throw new IOException("Blank row identifier at line " + lineNumber + " in " + path);
                if (!rowIds.add(rowId))
                    throw new IOException("Duplicate row identifier '" + rowId + "' in " + path);
                rowCount++;
            }
            if (rowCount == 0)
                throw new IOException("Matrix contains no data rows: " + path);
            return new QMatrixRowSource.Metadata(path, rowCount, sampleIds.length, sampleIds, null);
        } finally {
            if (reader != null)
                reader.close();
        }
    }

    private CsvReader createCsvReader() throws IOException {
        CsvReader reader = new CsvReader(openReader(path));
        reader.setTrimWhitespace(true);
        reader.setDelimiter(delimiters);
        if (commentMarkers != null && !commentMarkers.isEmpty()) {
            reader.setUseComments(true);
            reader.setComment(commentMarkers.charAt(0));
        }
        return reader;
    }

    private void validateFieldCount(String[] tokens, long lineNumber) throws IOException {
        int expected = metadata.columnCount() + 1;
        if (tokens.length != expected)
            throw new IOException("Expected " + expected + " fields but found " + tokens.length
                + " at line " + lineNumber + " in " + path);
    }

    private static int[] normalizeColumnOrder(int[] requested, int columnCount) {
        if (requested == null) {
            int[] identity = new int[columnCount];
            for (int i = 0; i < columnCount; i++)
                identity[i] = i;
            return identity;
        }
        if (requested.length == 0 || requested.length > columnCount)
            throw new IllegalArgumentException("Column selection has " + requested.length
                + " entries; expected between 1 and " + columnCount);
        boolean[] seen = new boolean[columnCount];
        for (int value : requested) {
            if (value < 0 || value >= columnCount || seen[value])
                throw new IllegalArgumentException("Column selection contains an invalid or duplicate index");
            seen[value] = true;
        }
        return requested.clone();
    }

    private static boolean isMissingToken(String value) {
        String normalized = value.trim().toLowerCase(Locale.ROOT);
        return normalized.isEmpty() || normalized.equals("na") || normalized.equals("n/a")
            || normalized.equals("null") || normalized.equals("nan") || normalized.equals(".");
    }

    private static void validateUniqueIdentifiers(String[] identifiers, String kind, Path path) throws IOException {
        Set<String> seen = new HashSet<>();
        for (String raw : identifiers) {
            String identifier = raw == null ? "" : raw.trim();
            if (identifier.isEmpty())
                throw new IOException("Blank " + kind + " identifier in " + path);
            if (!seen.add(identifier))
                throw new IOException("Duplicate " + kind + " identifier '" + identifier + "' in " + path);
        }
    }

    public static Reader openReader(Path path) throws IOException {
        InputStream input = new BufferedInputStream(new FileInputStream(path.toFile()), BUFFER_SIZE);
        String name = path.getFileName().toString().toLowerCase();
        try {
            if (name.endsWith(".gz"))
                input = new GZIPInputStream(input, BUFFER_SIZE);
            else if (name.endsWith(".bz2"))
                input = new BZip2CompressorInputStream(input, true);
            return new BufferedReader(new InputStreamReader(input, Charset.defaultCharset()), BUFFER_SIZE);
        } catch (IOException | RuntimeException e) {
            input.close();
            throw e;
        }
    }
}
