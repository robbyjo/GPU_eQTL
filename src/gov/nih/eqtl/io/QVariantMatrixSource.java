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

import static gov.nih.utils.QDataUtils.kUndefinedValue;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/** Re-readable VCF/BCF dosage rows with deterministic variant QC. */
public final class QVariantMatrixSource implements QMatrixRowSource {
    private static final int BUFFER_SIZE = 16 * 1024 * 1024;
    private static final double INTEGER_TOLERANCE = 1e-8;

    public enum Format {
        VCF, BCF;

        public static Format parse(String value, Path path) {
            String normalized = value == null ? "auto" : value.trim().toLowerCase(Locale.ROOT);
            if (normalized.equals("auto")) {
                String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
                return name.endsWith(".bcf") ? BCF : VCF;
            }
            if (normalized.equals("vcf") || normalized.equals("vcf.gz") || normalized.equals("vcfgz"))
                return VCF;
            if (normalized.equals("bcf"))
                return BCF;
            throw new IllegalArgumentException("Variant genotype format must be auto, vcf, vcf.gz, or bcf");
        }
    }

    public enum GenotypeField {
        AUTO, DS, GT;

        public static GenotypeField parse(String value) {
            if (value == null || value.isBlank())
                return AUTO;
            try {
                return valueOf(value.trim().toUpperCase(Locale.ROOT));
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("genotype_field must be auto, DS, or GT", e);
            }
        }
    }

    public enum MissingPolicy {
        ERROR, MEAN, ZERO, EXCLUDE_VARIANT, PRESERVE;

        public static MissingPolicy parse(String value) {
            if (value == null || value.isBlank())
                return MEAN;
            String normalized = value.trim().toUpperCase(Locale.ROOT).replace('-', '_');
            if (normalized.equals("EXCLUDE_ROW")) normalized = "EXCLUDE_VARIANT";
            try {
                return valueOf(normalized);
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException(
                    "genotype_missing must be error, mean, zero, or exclude-variant", e);
            }
        }
    }

    public enum MultiallelicPolicy {
        ERROR, EXCLUDE;

        public static MultiallelicPolicy parse(String value) {
            if (value == null || value.isBlank())
                return EXCLUDE;
            try {
                return valueOf(value.trim().toUpperCase(Locale.ROOT));
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("multiallelic must be error or exclude", e);
            }
        }
    }

    public record Options(Format format, GenotypeField genotypeField, MissingPolicy missingPolicy,
        MultiallelicPolicy multiallelicPolicy, double minimumMaf, double minimumMac,
        Path qcOutput) {
        public Options {
            if (format == null || genotypeField == null || missingPolicy == null || multiallelicPolicy == null)
                throw new IllegalArgumentException("Variant reader options must not be null");
            if (!Double.isFinite(minimumMaf) || minimumMaf < 0 || minimumMaf > 0.5)
                throw new IllegalArgumentException("min_maf must be between 0 and 0.5");
            if (!Double.isFinite(minimumMac) || minimumMac < 0)
                throw new IllegalArgumentException("min_mac must not be negative");
            qcOutput = qcOutput == null ? null : qcOutput.toAbsolutePath().normalize();
        }
    }

    public record Summary(long inputRecords, long includedVariants, long monomorphicVariants,
        long singletons, long doubletons, long excludedByMaf, long excludedByMac,
        long excludedForMissingness, long multiallelicRecords) { }

    private interface VariantCursor extends AutoCloseable {
        VCFHeader header();
        boolean hasNext();
        VariantContext next();
        @Override void close() throws IOException;
    }

    private static final class VcfCursor implements VariantCursor {
        private final VCFIterator iterator;

        VcfCursor(Path path) throws IOException {
            iterator = new VCFIteratorBuilder().open(path);
        }

        @Override public VCFHeader header() { return iterator.getHeader(); }
        @Override public boolean hasNext() { return iterator.hasNext(); }
        @Override public VariantContext next() { return iterator.next(); }
        @Override public void close() { iterator.close(); }
    }

    private static final class BcfCursor implements VariantCursor {
        private final QBcf22Codec codec = new QBcf22Codec();
        private final PositionalBufferedStream input;
        private final VCFHeader header;

        BcfCursor(Path path) throws IOException {
            InputStream opened = new BufferedInputStream(Files.newInputStream(path), BUFFER_SIZE);
            try {
                if (IOUtil.isGZIPInputStream(opened))
                    opened = new BufferedInputStream(new GZIPInputStream(opened, BUFFER_SIZE), BUFFER_SIZE);
                input = codec.makeSourceFromStream(opened);
                header = (VCFHeader) codec.readHeader(input).getHeaderValue();
            } catch (IOException | RuntimeException e) {
                opened.close();
                throw e;
            }
        }

        @Override public VCFHeader header() { return header; }
        @Override public boolean hasNext() { return !codec.isDone(input); }

        @Override
        public VariantContext next() {
            if (!hasNext())
                throw new NoSuchElementException();
            return codec.decode(input);
        }

        @Override public void close() { input.close(); }
    }

    private record Evaluation(String rowId, double[] values, int calledSamples, int missingSamples,
        double alternateAlleleCount, double alleleNumber, double eaf, double maf, double mac,
        double hweP, int hweSamples, String classification, boolean included, String exclusionReason) { }

    private final class Reader implements BlockReader {
        private final VariantCursor cursor;
        private final int[] columnOrder;
        private long nextRowOffset;
        private boolean closed;
        private boolean exhausted;

        Reader(int[] requestedColumnOrder) throws IOException {
            cursor = openCursor();
            columnOrder = normalizeColumnOrder(requestedColumnOrder, metadata.columnCount());
        }

        @Override
        public Block readBlock(int maximumRows) throws IOException {
            if (maximumRows <= 0)
                throw new IllegalArgumentException("maximumRows must be positive");
            if (closed)
                throw new IOException("Variant reader is closed");
            if (exhausted)
                return null;

            List<String> rowIds = new ArrayList<>(maximumRows);
            List<double[]> rows = new ArrayList<>(maximumRows);
            long blockOffset = nextRowOffset;
            while (rows.size() < maximumRows && cursor.hasNext()) {
                Evaluation evaluation = evaluate(cursor.next());
                if (!evaluation.included())
                    continue;
                double[] reordered = new double[columnOrder.length];
                for (int i = 0; i < columnOrder.length; i++)
                    reordered[i] = evaluation.values()[columnOrder[i]];
                rowIds.add(evaluation.rowId());
                rows.add(reordered);
                nextRowOffset++;
            }
            if (!cursor.hasNext()) {
                exhausted = true;
                if (nextRowOffset != metadata.rowCount())
                    throw new IOException("Variant source changed while it was being read: expected "
                        + metadata.rowCount() + " accepted rows, found " + nextRowOffset);
            }
            if (rows.isEmpty())
                return null;
            return new Block(blockOffset, rowIds.toArray(String[]::new), rows.toArray(double[][]::new));
        }

        @Override
        public void close() throws IOException {
            if (!closed) {
                closed = true;
                cursor.close();
            }
        }
    }

    private final Path path;
    private final Options options;
    private final GenotypeField selectedField;
    private final String[] sourceSampleIds;
    private final Metadata metadata;
    private final Summary summary;

    public QVariantMatrixSource(Path path, Options options) throws IOException {
        this.path = path.toAbsolutePath().normalize();
        this.options = options;
        if (!Files.isRegularFile(this.path))
            throw new IOException("Variant file does not exist: " + this.path);
        try (VariantCursor cursor = openCursor()) {
            validateSampleIds(cursor.header().getGenotypeSamples());
            sourceSampleIds = cursor.header().getGenotypeSamples().toArray(String[]::new);
            selectedField = selectField(cursor.header(), options.genotypeField());
        }
        ScanResult scan = scanMetadata();
        metadata = scan.metadata();
        summary = scan.summary();
    }

    @Override public Metadata metadata() { return metadata; }
    public Summary summary() { return summary; }
    public GenotypeField selectedField() { return selectedField; }

    @Override
    public BlockReader open(int[] columnOrder) throws IOException {
        return new Reader(columnOrder);
    }

    private record ScanResult(Metadata metadata, Summary summary) { }

    private ScanResult scanMetadata() throws IOException {
        Path qcOutput = options.qcOutput();
        Path temporary = null;
        BufferedWriter writer = null;
        boolean complete = false;
        long inputRecords = 0;
        long included = 0;
        long monomorphic = 0;
        long singletons = 0;
        long doubletons = 0;
        long excludedByMaf = 0;
        long excludedByMac = 0;
        long excludedForMissingness = 0;
        long multiallelic = 0;
        Set<String> rowIds = new HashSet<>();
        String[] samples;
        try {
            if (qcOutput != null) {
                if (path.equals(qcOutput))
                    throw new IOException("Variant QC output must differ from the genotype input: " + path);
                Path parent = qcOutput.getParent();
                if (parent != null)
                    Files.createDirectories(parent);
                temporary = Files.createTempFile(parent, qcOutput.getFileName().toString(), ".partial");
                writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8);
                writer.write("variant_id\tchromosome\tposition\trs_id\tref\talt\tvcf_filter\tgenotype_field"
                    + "\tcalled_samples\tmissing_samples\teffect_allele_count\tallele_number\teaf\tmaf\tmac"
                    + "\thwe_p\thwe_samples\tclassification\tincluded\texclusion_reason\n");
            }
            try (VariantCursor cursor = openCursor()) {
                samples = cursor.header().getGenotypeSamples().toArray(String[]::new);
                while (cursor.hasNext()) {
                    VariantContext variant = cursor.next();
                    inputRecords++;
                    Evaluation evaluation = evaluate(variant);
                    if (evaluation.classification().equals("monomorphic")) monomorphic++;
                    if (evaluation.classification().equals("singleton")) singletons++;
                    if (evaluation.classification().equals("doubleton")) doubletons++;
                    if (evaluation.classification().equals("multiallelic")) multiallelic++;
                    if (evaluation.exclusionReason().contains("maf")) excludedByMaf++;
                    if (evaluation.exclusionReason().contains("mac")) excludedByMac++;
                    if (evaluation.exclusionReason().contains("missing")) excludedForMissingness++;
                    if (!rowIds.add(evaluation.rowId()))
                        throw new IOException("Duplicate canonical variant identifier '"
                            + evaluation.rowId() + "' in " + path);
                    if (evaluation.included())
                        included++;
                    if (writer != null)
                        writeQc(writer, variant, evaluation);
                }
            }
            if (included == 0)
                throw new IOException("No variants remain after QC/filtering in " + path);
            if (writer != null) {
                writer.close();
                writer = null;
                moveAtomically(temporary, qcOutput);
            }
            complete = true;
        } finally {
            if (writer != null)
                writer.close();
            if (!complete && temporary != null)
                Files.deleteIfExists(temporary);
        }

        String signatureTag = String.format(Locale.ROOT,
            "variant-v1;format=%s;field=%s;missing=%s;multiallelic=%s;min-maf=%.17g;min-mac=%.17g",
            options.format(), selectedField, options.missingPolicy(), options.multiallelicPolicy(),
            options.minimumMaf(), options.minimumMac());
        Metadata scanned = new Metadata(path, included, samples.length, samples, signatureTag);
        Summary counts = new Summary(inputRecords, included, monomorphic, singletons, doubletons,
            excludedByMaf, excludedByMac, excludedForMissingness, multiallelic);
        return new ScanResult(scanned, counts);
    }

    private Evaluation evaluate(VariantContext variant) throws IOException {
        if (variant.getNAlleles() != 2) {
            if (options.multiallelicPolicy() == MultiallelicPolicy.ERROR)
                throw new IOException("Multiallelic variant encountered at " + variant.getContig()
                    + ":" + variant.getStart() + "; use --multiallelic exclude");
            return excludedEvaluation(canonicalId(variant), "multiallelic", "multiallelic");
        }

        Allele reference = variant.getReference();
        Allele alternate = variant.getAlternateAllele(0);
        int sampleCount = sourceSampleIds.length;
        double[] values = new double[sampleCount];
        boolean[] missing = new boolean[sampleCount];
        double alternateAlleleCount = 0;
        int calledSamples = 0;
        int missingSamples = 0;
        int homRef = 0;
        int heterozygous = 0;
        int homAlt = 0;
        int hweSamples = 0;
        if (variant.getNSamples() != sampleCount)
            throw new IOException("Variant sample count changed at " + variant.getContig()
                + ":" + variant.getStart());

        for (int sample = 0; sample < sampleCount; sample++) {
            Genotype genotype = variant.getGenotype(sourceSampleIds[sample]);
            Double dosage = selectedField == GenotypeField.DS
                ? dosageFromDs(genotype, variant) : dosageFromGt(genotype, alternate, variant);
            if (dosage == null) {
                missing[sample] = true;
                missingSamples++;
            } else {
                values[sample] = dosage;
                alternateAlleleCount += dosage;
                calledSamples++;
            }

            if (genotype != null && !genotype.isFiltered() && genotype.isCalled()
                && genotype.getPloidy() == 2) {
                int alternateCopies = 0;
                boolean valid = true;
                for (Allele allele : genotype.getAlleles()) {
                    if (allele.equals(reference)) {
                        // Reference copy.
                    } else if (allele.equals(alternate)) {
                        alternateCopies++;
                    } else {
                        valid = false;
                    }
                }
                if (valid) {
                    if (alternateCopies == 0) homRef++;
                    else if (alternateCopies == 1) heterozygous++;
                    else homAlt++;
                    hweSamples++;
                }
            }
        }

        String rowId = canonicalId(variant);
        if (calledSamples == 0) {
            if (options.missingPolicy() == MissingPolicy.ERROR)
                throw new IOException("No called " + selectedField + " values at " + rowId
                    + "; use --genotype-missing exclude-variant after confirming the source field");
            return new Evaluation(rowId, values, 0, missingSamples, 0, 0,
                Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0, "no_calls", false,
                "no_called_genotypes");
        }

        double alleleNumber = 2.0 * calledSamples;
        double eaf = alternateAlleleCount / alleleNumber;
        if (eaf < -INTEGER_TOLERANCE || eaf > 1 + INTEGER_TOLERANCE)
            throw new IOException("Allele dosage is outside [0,2] at " + rowId);
        eaf = Math.max(0, Math.min(1, eaf));
        double maf = Math.min(eaf, 1 - eaf);
        double mac = Math.min(alternateAlleleCount, alleleNumber - alternateAlleleCount);
        String classification = classify(mac);
        double hweP = hweSamples == 0 ? Double.NaN : hardyWeinbergExact(homRef, heterozygous, homAlt);

        List<String> reasons = new ArrayList<>();
        if (classification.equals("monomorphic"))
            reasons.add("monomorphic");
        if (missingSamples > 0) {
            if (options.missingPolicy() == MissingPolicy.ERROR)
                throw new IOException("Missing " + selectedField + " value at " + rowId
                    + "; use --genotype-missing mean or exclude-variant");
            if (options.missingPolicy() == MissingPolicy.EXCLUDE_VARIANT)
                reasons.add("missing_genotype");
        }
        if (maf + INTEGER_TOLERANCE < options.minimumMaf())
            reasons.add("maf_below_minimum");
        if (mac + INTEGER_TOLERANCE < options.minimumMac())
            reasons.add("mac_below_minimum");
        boolean accepted = reasons.isEmpty();
        if (accepted && missingSamples > 0) {
            double meanDosage = options.missingPolicy() == MissingPolicy.ZERO ? 0
                : options.missingPolicy() == MissingPolicy.PRESERVE ? kUndefinedValue : 2 * eaf;
            for (int i = 0; i < values.length; i++)
                if (missing[i])
                    values[i] = meanDosage;
        }
        return new Evaluation(rowId, values, calledSamples, missingSamples, alternateAlleleCount,
            alleleNumber, eaf, maf, mac, hweP, hweSamples, classification, accepted,
            String.join(";", reasons));
    }

    private Double dosageFromDs(Genotype genotype, VariantContext variant) throws IOException {
        if (genotype == null || genotype.isFiltered())
            return null;
        Object raw = genotype.getAnyAttribute("DS");
        if (raw == null)
            return null;
        if (raw instanceof List<?> list) {
            if (list.size() != 1)
                throw new IOException("Expected one DS value at " + canonicalId(variant));
            raw = list.get(0);
        }
        String text = raw.toString().trim();
        if (text.isEmpty() || text.equals("."))
            return null;
        try {
            double value = raw instanceof Number number ? number.doubleValue() : Double.parseDouble(text);
            if (!Double.isFinite(value) || value < -INTEGER_TOLERANCE || value > 2 + INTEGER_TOLERANCE)
                throw new IOException("Invalid DS value '" + text + "' at " + canonicalId(variant));
            return Math.max(0, Math.min(2, value));
        } catch (NumberFormatException e) {
            throw new IOException("Invalid DS value '" + text + "' at " + canonicalId(variant), e);
        }
    }

    private static Double dosageFromGt(Genotype genotype, Allele alternate, VariantContext variant)
        throws IOException {
        if (genotype == null || genotype.isFiltered() || !genotype.isCalled())
            return null;
        if (genotype.getPloidy() != 2)
            throw new IOException("Only diploid GT values are supported at " + canonicalId(variant));
        return (double) genotype.countAllele(alternate);
    }

    private static String classify(double mac) {
        if (mac <= INTEGER_TOLERANCE)
            return "monomorphic";
        if (Math.abs(mac - 1) <= INTEGER_TOLERANCE)
            return "singleton";
        if (Math.abs(mac - 2) <= INTEGER_TOLERANCE)
            return "doubleton";
        return "polymorphic";
    }

    /** Wigginton-style exact HWE test, conditional on the observed rare-allele count. */
    static double hardyWeinbergExact(int observedHomRef, int observedHets, int observedHomAlt) {
        if (observedHomRef < 0 || observedHets < 0 || observedHomAlt < 0)
            throw new IllegalArgumentException("Genotype counts must not be negative");
        int observations = observedHomRef + observedHets + observedHomAlt;
        if (observations == 0)
            return Double.NaN;
        int rareCopies = 2 * Math.min(observedHomRef, observedHomAlt) + observedHets;
        int commonCopies = 2 * observations - rareCopies;
        double[] probabilities = new double[rareCopies + 1];
        int midpoint = (int) ((long) rareCopies * commonCopies / (2L * observations));
        if ((midpoint & 1) != (rareCopies & 1))
            midpoint++;
        probabilities[midpoint] = 1.0;
        double sum = 1.0;

        int rareHom = (rareCopies - midpoint) / 2;
        int commonHom = observations - midpoint - rareHom;
        for (int het = midpoint; het >= 2; het -= 2) {
            probabilities[het - 2] = probabilities[het] * het * (het - 1.0)
                / (4.0 * (rareHom + 1.0) * (commonHom + 1.0));
            sum += probabilities[het - 2];
            rareHom++;
            commonHom++;
        }
        rareHom = (rareCopies - midpoint) / 2;
        commonHom = observations - midpoint - rareHom;
        for (int het = midpoint; het <= rareCopies - 2; het += 2) {
            probabilities[het + 2] = probabilities[het] * 4.0 * rareHom * commonHom
                / ((het + 2.0) * (het + 1.0));
            sum += probabilities[het + 2];
            rareHom--;
            commonHom--;
        }
        for (int i = 0; i < probabilities.length; i++)
            probabilities[i] /= sum;
        double observedProbability = probabilities[observedHets];
        double p = 0;
        for (double probability : probabilities)
            if (probability <= observedProbability + 1e-12)
                p += probability;
        return Math.min(1.0, p);
    }

    private static Evaluation excludedEvaluation(String rowId, String classification, String reason) {
        return new Evaluation(rowId, new double[0], 0, 0, Double.NaN, Double.NaN,
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0, classification, false, reason);
    }

    private static String canonicalId(VariantContext variant) {
        String alternate = variant.getAlternateAlleles().isEmpty() ? "."
            : String.join(",", variant.getAlternateAlleles().stream().map(Allele::getDisplayString).toList());
        return variant.getContig() + ":" + variant.getStart() + ":"
            + variant.getReference().getDisplayString() + ":" + alternate;
    }

    private void writeQc(BufferedWriter writer, VariantContext variant, Evaluation evaluation) throws IOException {
        String rsId = variant.hasID() ? variant.getID() : ".";
        String alternate = variant.getAlternateAlleles().isEmpty() ? "."
            : String.join(",", variant.getAlternateAlleles().stream().map(Allele::getDisplayString).toList());
        String filter = variant.isFiltered() ? String.join(";", variant.getFilters())
            : (variant.filtersWereApplied() ? "PASS" : ".");
        writer.write(tsv(evaluation.rowId())); writer.write('\t');
        writer.write(tsv(variant.getContig())); writer.write('\t');
        writer.write(Long.toString(variant.getStart())); writer.write('\t');
        writer.write(tsv(rsId)); writer.write('\t');
        writer.write(tsv(variant.getReference().getDisplayString())); writer.write('\t');
        writer.write(tsv(alternate)); writer.write('\t');
        writer.write(tsv(filter)); writer.write('\t');
        writer.write(selectedField.name()); writer.write('\t');
        writer.write(Integer.toString(evaluation.calledSamples())); writer.write('\t');
        writer.write(Integer.toString(evaluation.missingSamples())); writer.write('\t');
        writer.write(number(evaluation.alternateAlleleCount())); writer.write('\t');
        writer.write(number(evaluation.alleleNumber())); writer.write('\t');
        writer.write(number(evaluation.eaf())); writer.write('\t');
        writer.write(number(evaluation.maf())); writer.write('\t');
        writer.write(number(evaluation.mac())); writer.write('\t');
        writer.write(number(evaluation.hweP())); writer.write('\t');
        writer.write(Integer.toString(evaluation.hweSamples())); writer.write('\t');
        writer.write(evaluation.classification()); writer.write('\t');
        writer.write(Boolean.toString(evaluation.included())); writer.write('\t');
        writer.write(evaluation.exclusionReason().isEmpty() ? "." : evaluation.exclusionReason());
        writer.write('\n');
    }

    private VariantCursor openCursor() throws IOException {
        return options.format() == Format.BCF ? new BcfCursor(path) : new VcfCursor(path);
    }

    private static GenotypeField selectField(VCFHeader header, GenotypeField requested) {
        if (requested == GenotypeField.AUTO) {
            if (header.getFormatHeaderLine("DS") != null)
                return GenotypeField.DS;
            if (header.getFormatHeaderLine("GT") != null)
                return GenotypeField.GT;
            throw new IllegalArgumentException("Variant header declares neither FORMAT/DS nor FORMAT/GT");
        }
        if (header.getFormatHeaderLine(requested.name()) == null)
            throw new IllegalArgumentException("Variant header does not declare FORMAT/" + requested);
        return requested;
    }

    private static void validateSampleIds(List<String> sampleIds) throws IOException {
        if (sampleIds.isEmpty())
            throw new IOException("Variant file has no genotype samples");
        Set<String> seen = new HashSet<>();
        for (String sampleId : sampleIds) {
            if (sampleId == null || sampleId.isBlank())
                throw new IOException("Variant file contains a blank sample identifier");
            if (!seen.add(sampleId))
                throw new IOException("Variant file contains duplicate sample identifier '" + sampleId + "'");
        }
    }

    private static int[] normalizeColumnOrder(int[] requested, int columnCount) {
        if (requested == null) {
            int[] identity = new int[columnCount];
            for (int i = 0; i < columnCount; i++) identity[i] = i;
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

    private static String number(double value) {
        return Double.isFinite(value) ? Double.toString(value) : "NA";
    }

    private static String tsv(String value) {
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
