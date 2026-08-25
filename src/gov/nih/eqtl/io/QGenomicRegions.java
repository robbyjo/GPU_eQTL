/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;

/** Parsed, normalized genomic intervals with deterministic overlapping set membership. */
public final class QGenomicRegions {
    public enum Coordinates {
        ONE_BASED, BED;

        public static Coordinates parse(String value) {
            if (value == null || value.isBlank())
                return ONE_BASED;
            String normalized = value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
            if (normalized.equals("one-based") || normalized.equals("vcf")
                || normalized.equals("1-based"))
                return ONE_BASED;
            if (normalized.equals("bed") || normalized.equals("zero-based-half-open"))
                return BED;
            throw new IllegalArgumentException(
                "region_coordinates must be one-based or bed");
        }
    }

    public record QueryInterval(String contig, int start, int end) {
        public QueryInterval {
            if (contig == null || contig.isBlank() || start < 1 || end < start)
                throw new IllegalArgumentException("Invalid genomic query interval");
        }
    }

    private record Definition(String setId, String contig, int start, int end, int order) { }

    private final List<Definition> definitions;
    private final List<QueryInterval> queryIntervals;
    private final List<String> setIds;
    private final String signature;

    private QGenomicRegions(List<Definition> definitions, List<QueryInterval> queryIntervals) {
        this.definitions = List.copyOf(definitions);
        this.queryIntervals = List.copyOf(queryIntervals);
        LinkedHashSet<String> ids = new LinkedHashSet<>();
        for (Definition definition : definitions)
            ids.add(definition.setId());
        setIds = List.copyOf(ids);
        signature = digestDefinitions(definitions);
    }

    /** Parse inline one-based regions and an optional region file, then resolve contig aliases. */
    public static QGenomicRegions load(String inline, Path regionFile,
        Coordinates fileCoordinates, List<String> availableContigs) throws IOException {
        List<Definition> parsed = new ArrayList<>();
        int order = 0;
        if (inline != null && !inline.isBlank()) {
            for (String token : inline.split("[;,]")) {
                String value = token.trim();
                if (!value.isEmpty())
                    parsed.add(parseInline(value, order++));
            }
        }
        if (regionFile != null) {
            Path normalized = regionFile.toAbsolutePath().normalize();
            if (!Files.isRegularFile(normalized))
                throw new IOException("Region file does not exist: " + normalized);
            try (BufferedReader reader = Files.newBufferedReader(normalized, StandardCharsets.UTF_8)) {
                String line;
                int lineNumber = 0;
                while ((line = reader.readLine()) != null) {
                    lineNumber++;
                    String trimmed = line.trim();
                    if (trimmed.isEmpty() || trimmed.startsWith("#"))
                        continue;
                    String[] fields = line.split("\\t", -1);
                    if (fields.length != 3 && fields.length != 4)
                        throw new IOException("Region file " + normalized + " line " + lineNumber
                            + " must contain CHROM, START, END or SET_ID, CHROM, START, END");
                    if (isHeader(fields))
                        continue;
                    int offset = fields.length == 4 ? 1 : 0;
                    String setId = fields.length == 4 ? fields[0].trim()
                        : fields[0].trim() + ":" + fields[1].trim() + "-" + fields[2].trim();
                    String contig = fields[offset].trim();
                    int start = parsePosition(fields[offset + 1], normalized + " line " + lineNumber);
                    int end = parsePosition(fields[offset + 2], normalized + " line " + lineNumber);
                    if (fileCoordinates == Coordinates.BED) {
                        if (start < 0 || end <= start)
                            throw new IOException("Invalid BED half-open interval in " + normalized
                                + " line " + lineNumber);
                        start++;
                    } else if (start < 1 || end < start) {
                        throw new IOException("Invalid one-based inclusive interval in " + normalized
                            + " line " + lineNumber);
                    }
                    validateSetId(setId);
                    if (contig.isEmpty())
                        throw new IOException("Blank contig in " + normalized + " line " + lineNumber);
                    parsed.add(new Definition(setId, contig, start, end, order++));
                }
            }
        }
        if (parsed.isEmpty())
            return null;
        if (availableContigs == null || availableContigs.isEmpty())
            throw new IOException("Genomic regions require contig declarations in the variant header or index");

        Map<String, Integer> contigOrder = new LinkedHashMap<>();
        for (String contig : availableContigs)
            if (contig != null && !contig.isBlank())
                contigOrder.putIfAbsent(contig, contigOrder.size());
        List<Definition> resolved = new ArrayList<>(parsed.size());
        for (Definition definition : parsed) {
            String contig = resolveContig(definition.contig(), contigOrder.keySet());
            resolved.add(new Definition(definition.setId(), contig,
                definition.start(), definition.end(), definition.order()));
        }
        List<QueryInterval> queries = mergeQueries(resolved, contigOrder);
        return new QGenomicRegions(resolved, queries);
    }

    public List<QueryInterval> queryIntervals() {
        return queryIntervals;
    }

    public List<String> setIds() {
        return setIds;
    }

    public String signature() {
        return signature;
    }

    /** Semicolon-separated set IDs in definition order; dot when no interval overlaps. */
    public String memberships(VariantContext variant) {
        LinkedHashSet<String> matches = new LinkedHashSet<>();
        for (Definition definition : definitions) {
            if (definition.contig().equals(variant.getContig())
                && variant.getStart() <= definition.end()
                && variant.getEnd() >= definition.start())
                matches.add(definition.setId());
        }
        return matches.isEmpty() ? "." : String.join(";", matches);
    }

    public List<String> emptySetIds(Set<String> observedSetIds) {
        Set<String> observed = observedSetIds == null ? Set.of() : observedSetIds;
        List<String> empty = new ArrayList<>();
        for (String setId : setIds)
            if (!observed.contains(setId))
                empty.add(setId);
        return List.copyOf(empty);
    }

    private static Definition parseInline(String value, int order) {
        String setId = null;
        String location = value;
        int equals = value.indexOf('=');
        if (equals >= 0) {
            setId = value.substring(0, equals).trim();
            location = value.substring(equals + 1).trim();
        }
        int colon = location.lastIndexOf(':');
        int dash = location.indexOf('-', colon + 1);
        if (colon <= 0 || dash <= colon + 1 || dash == location.length() - 1)
            throw new IllegalArgumentException("Invalid region '" + value
                + "'; expected [SET_ID=]CHROM:START-END");
        String contig = location.substring(0, colon).trim();
        int start = parsePosition(location.substring(colon + 1, dash), "region " + value);
        int end = parsePosition(location.substring(dash + 1), "region " + value);
        if (start < 1 || end < start)
            throw new IllegalArgumentException("Invalid one-based inclusive region '" + value + "'");
        if (setId == null || setId.isEmpty())
            setId = contig + ":" + start + "-" + end;
        validateSetId(setId);
        return new Definition(setId, contig, start, end, order);
    }

    private static boolean isHeader(String[] fields) {
        int offset = fields.length == 4 ? 1 : 0;
        return fields[offset].trim().equalsIgnoreCase("chrom")
            && fields[offset + 1].trim().equalsIgnoreCase("start")
            && fields[offset + 2].trim().equalsIgnoreCase("end");
    }

    private static int parsePosition(String value, String context) {
        try {
            long parsed = Long.parseLong(value.trim());
            if (parsed < 0 || parsed > Integer.MAX_VALUE)
                throw new NumberFormatException();
            return (int) parsed;
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid genomic coordinate '" + value
                + "' in " + context, e);
        }
    }

    private static void validateSetId(String setId) {
        if (setId == null || setId.isBlank() || setId.indexOf('\t') >= 0
            || setId.indexOf(';') >= 0 || setId.indexOf(',') >= 0)
            throw new IllegalArgumentException(
                "Region set IDs must be nonblank and must not contain tabs, commas, or semicolons");
    }

    private static String resolveContig(String requested, Set<String> available) throws IOException {
        if (available.contains(requested))
            return requested;
        List<String> candidates = new ArrayList<>();
        String toggled = requested.regionMatches(true, 0, "chr", 0, 3)
            ? requested.substring(3) : "chr" + requested;
        for (String candidate : available) {
            if (candidate.equalsIgnoreCase(requested) || candidate.equalsIgnoreCase(toggled)
                || isMitochondrialAlias(requested, candidate))
                candidates.add(candidate);
        }
        if (candidates.size() == 1)
            return candidates.get(0);
        if (candidates.isEmpty())
            throw new IOException("Region contig '" + requested
                + "' is absent from the variant header/index");
        throw new IOException("Region contig '" + requested + "' is ambiguous: " + candidates);
    }

    private static boolean isMitochondrialAlias(String requested, String candidate) {
        return mitochondrialName(requested) && mitochondrialName(candidate);
    }

    private static boolean mitochondrialName(String value) {
        String normalized = value.toUpperCase(Locale.ROOT);
        if (normalized.startsWith("CHR"))
            normalized = normalized.substring(3);
        return normalized.equals("M") || normalized.equals("MT");
    }

    private static List<QueryInterval> mergeQueries(List<Definition> definitions,
        Map<String, Integer> contigOrder) {
        List<Definition> sorted = new ArrayList<>(definitions);
        sorted.sort(Comparator
            .comparingInt((Definition value) -> contigOrder.get(value.contig()))
            .thenComparingInt(Definition::start)
            .thenComparingInt(Definition::end)
            .thenComparingInt(Definition::order));
        List<QueryInterval> merged = new ArrayList<>();
        for (Definition definition : sorted) {
            if (merged.isEmpty()) {
                merged.add(new QueryInterval(definition.contig(), definition.start(), definition.end()));
                continue;
            }
            QueryInterval previous = merged.get(merged.size() - 1);
            if (previous.contig().equals(definition.contig())
                && (long) definition.start() <= (long) previous.end() + 1) {
                merged.set(merged.size() - 1, new QueryInterval(previous.contig(),
                    previous.start(), Math.max(previous.end(), definition.end())));
            } else {
                merged.add(new QueryInterval(definition.contig(), definition.start(), definition.end()));
            }
        }
        return merged;
    }

    private static String digestDefinitions(List<Definition> definitions) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            digest.update("gpu-eqtl-regions-v1".getBytes(StandardCharsets.UTF_8));
            for (Definition definition : definitions) {
                updateDigest(digest, definition.setId());
                updateDigest(digest, definition.contig());
                updateDigest(digest, Integer.toString(definition.start()));
                updateDigest(digest, Integer.toString(definition.end()));
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
    }

    private static void updateDigest(MessageDigest digest, String value) {
        digest.update((byte) 0);
        digest.update(value.getBytes(StandardCharsets.UTF_8));
    }
}
