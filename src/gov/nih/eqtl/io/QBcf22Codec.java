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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;
import htsjdk.variant.vcf.VCFHeaderLineTranslator;
import htsjdk.variant.vcf.VCFHeaderVersion;

/**
 * HTSJDK's BCF decoder implements the BCF 2 record layout but rejects the
 * current 2.2 version byte and ignores explicit header IDX values. This
 * adapter extends its version policy and restores the on-disk dictionaries. The
 * matrix reader separately restricts accepted records to biallelic diploid
 * GT/DS values, avoiding the variable-ploidy vector-end case that differs
 * between the old and current specifications.
 */
final class QBcf22Codec extends BCF2Codec {
    private static final int MAX_HEADER_BYTES = 0x08000000; // Same bound as HTSJDK 5.0.0.

    @Override
    public FeatureCodecHeader readHeader(PositionalBufferedStream input) {
        try {
            byte[] prefix = input.readNBytes(9);
            if (prefix.length != 9)
                throw new TribbleException("Truncated BCF header");
            int length = ByteBuffer.wrap(prefix, 5, 4).order(ByteOrder.LITTLE_ENDIAN).getInt();
            if (length <= 0 || length > MAX_HEADER_BYTES)
                throw new TribbleException("Invalid BCF header length: " + length);
            byte[] text = input.readNBytes(length);
            if (text.length != length)
                throw new TribbleException("Truncated BCF header text");

            // Keep the original stream at its original record offset. In particular,
            // Tribble indexes must not see a header length from a rewritten stream.
            byte[] headerBytes = new byte[9 + length];
            System.arraycopy(prefix, 0, headerBytes, 0, 9);
            System.arraycopy(text, 0, headerBytes, 9, length);
            FeatureCodecHeader parsed;
            try (PositionalBufferedStream header = new PositionalBufferedStream(
                    new ByteArrayInputStream(headerBytes))) {
                parsed = super.readHeader(header);
            }
            installExplicitDictionaries(new String(text, StandardCharsets.UTF_8));
            return new FeatureCodecHeader(parsed.getHeaderValue(), input.getPosition());
        } catch (IOException e) {
            throw new TribbleException("I/O error while reading BCF header", e);
        }
    }

    private void installExplicitDictionaries(String text) {
        IndexedNames tags = new IndexedNames("FILTER/INFO/FORMAT");
        IndexedNames contigs = new IndexedNames("contig");
        tags.add("PASS", "0"); // BCF reserves dictionary slot zero, even without a PASS line.
        boolean tagIndices = false, contigIndices = false;
        for (String line : text.split("\n")) {
            int equals = line.indexOf('=');
            if (equals < 0) continue;
            String kind = line.substring(0, equals);
            boolean tag = kind.equals("##FILTER") || kind.equals("##INFO") || kind.equals("##FORMAT");
            if (!tag && !kind.equals("##contig")) continue;
            Map<String, String> fields = VCFHeaderLineTranslator.parseLine(
                VCFHeaderVersion.VCF4_2, line.substring(equals + 1).strip(), java.util.List.of());
            String id = fields.get("ID"), index = fields.get("IDX");
            if (id == null || id.isBlank())
                throw new TribbleException("Missing BCF header ID in " + kind);
            if (tag) { tags.add(id, index); tagIndices |= index != null; }
            else { contigs.add(id, index); contigIndices |= index != null; }
        }
        if (tagIndices) installLookup("dictionary", tags.resolve());
        if (contigIndices) installLookup("contigNames", contigs.resolve());
    }

    /**
     * The pinned HTSJDK codec has private lookup fields and final accessors, with
     * no dictionary-initialization hook. Keep this compatibility bridge in one
     * checked location rather than copying its record/lazy-genotype decoder.
     * An incompatible dependency upgrade fails explicitly; no dense fallback.
     */
    private void installLookup(String name, ArrayList<String> lookup) {
        try {
            Field field = BCF2Codec.class.getDeclaredField(name);
            if (field.getType() != ArrayList.class || !field.trySetAccessible())
                throw new ReflectiveOperationException("Unsupported lookup field " + name);
            field.set(this, lookup);
        } catch (ReflectiveOperationException | SecurityException e) {
            throw new TribbleException("Cannot install explicit BCF IDX dictionary in HTSJDK", e);
        }
    }

    private static final class IndexedNames {
        private final String kind;
        private final Map<String, Integer> indices = new HashMap<>();
        private final Map<Integer, String> names = new HashMap<>();
        private int next;

        IndexedNames(String kind) { this.kind = kind; }

        void add(String id, String text) {
            Integer index = null;
            if (text != null) {
                if (!text.matches("[0-9]+"))
                    throw new TribbleException("Invalid BCF IDX for " + id + ": " + text);
                try { index = Integer.valueOf(text); }
                catch (NumberFormatException e) { throw new TribbleException("BCF IDX exceeds integer range for " + id, e); }
                if (index >= Integer.MAX_VALUE - 1)
                    throw new TribbleException("BCF IDX exceeds supported list range for " + id);
            }
            Integer previous = indices.get(id);
            if (previous != null && index != null && !previous.equals(index))
                throw new TribbleException("Conflicting BCF IDX values for " + id);
            if (previous != null) return;
            // HTSlib registers each first occurrence immediately. Later header
            // lines in another namespace must not move an already assigned ID.
            int assigned = index == null ? next : index;
            if (assigned >= Integer.MAX_VALUE - 1)
                throw new TribbleException("BCF dictionary exceeds supported index range");
            String occupant = names.putIfAbsent(assigned, id);
            if (occupant != null && !occupant.equals(id))
                throw new TribbleException("Duplicate BCF " + kind + " IDX " + assigned);
            indices.put(id, assigned);
            next = Math.max(next, Math.addExact(assigned, 1));
        }

        ArrayList<String> resolve() {
            return new SparseLookup(kind, names, next);
        }
    }

    private static final class SparseLookup extends ArrayList<String> {
        private static final long serialVersionUID = 1L;
        private final String kind;
        private final Map<Integer, String> names;
        private final int length;
        SparseLookup(String kind, Map<Integer, String> names, int length) {
            this.kind = kind; this.names = Map.copyOf(names); this.length = length;
        }
        @Override public int size() { return length; }
        @Override public String get(int index) {
            String name = names.get(index);
            if (name == null)
                throw new TribbleException("Undeclared BCF " + kind + " dictionary index " + index);
            return name;
        }
    }

    @Override
    protected void validateVersionCompatibility(BCFVersion supportedVersion, BCFVersion actualVersion) {
        if (actualVersion.getMajorVersion() != 2
            || (actualVersion.getMinorVersion() != 1 && actualVersion.getMinorVersion() != 2)) {
            throw new TribbleException("GPU eQTL supports BCF 2.1 and 2.2; found " + actualVersion);
        }
    }
}
