/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.settest;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Deterministic explicit variant-to-set membership and effect-allele definitions.
 * Input is tab-separated with required columns SET_ID, VARIANT_ID, REF, ALT, and
 * EFFECT_ALLELE, plus an optional positive WEIGHT column.
 */
public final class QVariantSetTable {
    public record Entry(String setId, String variantId, String ref, String alt,
        String effectAllele, double weight, int sourceLine) {
        public Entry {
            setId = requireIdentifier(setId, "set ID", sourceLine);
            variantId = requireIdentifier(variantId, "variant ID", sourceLine);
            ref = normalizeAllele(ref, "REF", sourceLine);
            alt = normalizeAllele(alt, "ALT", sourceLine);
            effectAllele = normalizeAllele(effectAllele, "effect allele", sourceLine);
            if (ref.equals(alt))
                throw new IllegalArgumentException("REF and ALT are identical at set-definition line "
                    + sourceLine);
            if (!effectAllele.equals(ref) && !effectAllele.equals(alt))
                throw new IllegalArgumentException("Effect allele " + effectAllele
                    + " is neither REF nor ALT at set-definition line " + sourceLine);
            if (!(weight > 0) || !Double.isFinite(weight))
                throw new IllegalArgumentException("Variant weight must be finite and positive at "
                    + "set-definition line " + sourceLine);
        }

        public boolean effectIsAlt() {
            return effectAllele.equals(alt);
        }

        public void requireExactAlleles(String sourceRef, String sourceAlt) {
            String normalizedRef = normalizeAllele(sourceRef, "source REF", sourceLine);
            String normalizedAlt = normalizeAllele(sourceAlt, "source ALT", sourceLine);
            if (!ref.equals(normalizedRef) || !alt.equals(normalizedAlt))
                throw new IllegalArgumentException("Allele mismatch for variant '" + variantId
                    + "' in set '" + setId + "': definition " + ref + "/" + alt
                    + ", source " + normalizedRef + "/" + normalizedAlt
                    + ". Strand flips and implicit allele swaps are not enabled.");
        }
    }

    public record SetDefinition(String id, List<Entry> entries) {
        public SetDefinition {
            entries = List.copyOf(entries);
            if (id == null || id.isBlank() || entries.isEmpty())
                throw new IllegalArgumentException("A variant set must have an ID and entries");
        }
    }

    private static final List<String> REQUIRED_COLUMNS = List.of(
        "SET_ID", "VARIANT_ID", "REF", "ALT", "EFFECT_ALLELE");

    private final Path source;
    private final List<SetDefinition> sets;
    private final String signature;

    private QVariantSetTable(Path source, List<SetDefinition> sets, String signature) {
        this.source = source;
        this.sets = List.copyOf(sets);
        this.signature = signature;
    }

    public static QVariantSetTable load(Path path) throws IOException {
        if (path == null)
            throw new IllegalArgumentException("Variant-set file is required");
        Path normalized = path.toAbsolutePath().normalize();
        if (!Files.isRegularFile(normalized))
            throw new IOException("Variant-set file does not exist: " + normalized);
        LinkedHashMap<String, List<Entry>> bySet = new LinkedHashMap<>();
        Set<String> seenMemberships = new LinkedHashSet<>();
        int[] columns = null;
        boolean weightColumn = false;
        try (BufferedReader reader = Files.newBufferedReader(normalized, StandardCharsets.UTF_8)) {
            String line;
            int lineNumber = 0;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                if (line.isBlank() || line.stripLeading().startsWith("#"))
                    continue;
                String[] fields = line.split("\\t", -1);
                if (columns == null) {
                    columns = parseHeader(fields, normalized, lineNumber);
                    weightColumn = columns[5] >= 0;
                    continue;
                }
                int requiredWidth = 0;
                for (int column : columns)
                    requiredWidth = Math.max(requiredWidth, column + 1);
                if (fields.length < requiredWidth)
                    throw new IOException("Variant-set file " + normalized + " line " + lineNumber
                        + " has " + fields.length + " fields; expected at least " + requiredWidth);
                double weight = 1.0;
                if (weightColumn && !fields[columns[5]].isBlank()) {
                    try {
                        weight = Double.parseDouble(fields[columns[5]].trim());
                    } catch (NumberFormatException e) {
                        throw new IOException("Invalid WEIGHT at " + normalized + " line "
                            + lineNumber, e);
                    }
                }
                Entry entry;
                try {
                    entry = new Entry(fields[columns[0]], fields[columns[1]], fields[columns[2]],
                        fields[columns[3]], fields[columns[4]], weight, lineNumber);
                } catch (IllegalArgumentException e) {
                    throw new IOException("Invalid variant-set row in " + normalized + " line "
                        + lineNumber + ": " + e.getMessage(), e);
                }
                String membership = entry.setId() + "\u0000" + entry.variantId();
                if (!seenMemberships.add(membership))
                    throw new IOException("Duplicate variant '" + entry.variantId() + "' in set '"
                        + entry.setId() + "' at " + normalized + " line " + lineNumber);
                bySet.computeIfAbsent(entry.setId(), ignored -> new ArrayList<>()).add(entry);
            }
        }
        if (columns == null)
            throw new IOException("Variant-set file has no header: " + normalized);
        if (bySet.isEmpty())
            throw new IOException("Variant-set file has no data rows: " + normalized);
        List<SetDefinition> sets = new ArrayList<>(bySet.size());
        for (Map.Entry<String, List<Entry>> entry : bySet.entrySet())
            sets.add(new SetDefinition(entry.getKey(), entry.getValue()));
        return new QVariantSetTable(normalized, sets, signature(sets));
    }

    public Path source() { return source; }
    public List<SetDefinition> sets() { return sets; }
    public String signature() { return signature; }

    private static int[] parseHeader(String[] fields, Path path, int lineNumber)
        throws IOException {
        Map<String, Integer> indices = new LinkedHashMap<>();
        for (int i = 0; i < fields.length; i++) {
            String name = fields[i].trim().toUpperCase(Locale.ROOT);
            if (!name.isEmpty() && indices.putIfAbsent(name, i) != null)
                throw new IOException("Duplicate column '" + name + "' in " + path
                    + " line " + lineNumber);
        }
        int[] columns = new int[6];
        for (int i = 0; i < REQUIRED_COLUMNS.size(); i++) {
            Integer column = indices.get(REQUIRED_COLUMNS.get(i));
            if (column == null)
                throw new IOException("Variant-set header in " + path + " is missing "
                    + REQUIRED_COLUMNS.get(i));
            columns[i] = column;
        }
        columns[5] = indices.getOrDefault("WEIGHT", -1);
        return columns;
    }

    private static String requireIdentifier(String value, String kind, int line) {
        if (value == null || value.isBlank() || value.indexOf('\t') >= 0
            || value.indexOf(';') >= 0 || value.indexOf(',') >= 0)
            throw new IllegalArgumentException(kind + " is blank or contains a reserved delimiter at "
                + "set-definition line " + line);
        return value.trim();
    }

    private static String normalizeAllele(String value, String kind, int line) {
        if (value == null || value.isBlank() || value.indexOf(',') >= 0)
            throw new IllegalArgumentException(kind + " is blank or multiallelic at "
                + "set-definition line " + line);
        return value.trim().toUpperCase(Locale.ROOT);
    }

    private static String signature(List<SetDefinition> sets) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            digest.update("gpu-eqtl-variant-sets-v1".getBytes(StandardCharsets.UTF_8));
            for (SetDefinition set : sets) {
                update(digest, set.id());
                for (Entry entry : set.entries()) {
                    update(digest, entry.variantId());
                    update(digest, entry.ref());
                    update(digest, entry.alt());
                    update(digest, entry.effectAllele());
                    update(digest, Long.toString(Double.doubleToLongBits(entry.weight())));
                }
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
    }

    private static void update(MessageDigest digest, String value) {
        digest.update((byte) 0);
        digest.update(value.getBytes(StandardCharsets.UTF_8));
    }
}
