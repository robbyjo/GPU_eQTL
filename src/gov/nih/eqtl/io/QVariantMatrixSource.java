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
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/** Re-readable VCF/BCF dosage rows with deterministic variant QC. */
public final class QVariantMatrixSource implements QMatrixRowSource {
    private static final int BUFFER_SIZE = 16 * 1024 * 1024;
    private static final double INTEGER_TOLERANCE = 1e-8;
    private static final long PROGRESS_CHECK_RECORDS = 1_000;
    private static final long PROGRESS_REPORT_RECORDS = 1_000_000;
    private static final long PROGRESS_REPORT_NANOS = 15_000_000_000L;
    private static final int MAX_AUTOMATIC_QC_THREADS = 16;
    private static final int QC_TASKS_PER_THREAD = 2;
    private static final int DECISION_CHUNK_SIZE = 1 << 16;
    private static final int DEFAULT_QC_CHECKPOINT_RECORDS = 10_000;
    private static final String QC_CHECKPOINT_RECORDS_PROPERTY = "eqtl.variant.qc.checkpoint.records";
    private static final String TEST_FAIL_AFTER_PROPERTY = "eqtl.test.variant.qc.fail.after";
    private static final String QC_HEADER = "variant_id\tchromosome\tposition\trs_id\tref\talt\tvcf_filter"
        + "\tgenotype_field\tcalled_samples\tmissing_samples\teffect_allele_count\tallele_number"
        + "\tdosage_sum_squares\teaf\tmaf\tmac\thwe_p\thwe_samples\thom_ref\theterozygous\thom_alt"
        + "\tclassification\tincluded\texclusion_reason"
        + "\tregion_sets\tfrequency_scope\taligned_frequency_reason\n";

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

    public enum FrequencyScope {
        ALIGNED, PATTERN;

        public static FrequencyScope parse(String value) {
            if (value == null || value.isBlank())
                return ALIGNED;
            try {
                return valueOf(value.trim().toUpperCase(Locale.ROOT));
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("frequency_scope must be aligned or pattern", e);
            }
        }
    }

    public record Options(Format format, GenotypeField genotypeField, MissingPolicy missingPolicy,
        MultiallelicPolicy multiallelicPolicy, double minimumMaf, double minimumMac,
        Path qcOutput, int qcThreads, Path qcCheckpoint, Path variantIndex,
        String regions, Path regionsFile, QGenomicRegions.Coordinates regionCoordinates,
        FrequencyScope frequencyScope) {
        public Options(Format format, GenotypeField genotypeField, MissingPolicy missingPolicy,
            MultiallelicPolicy multiallelicPolicy, double minimumMaf, double minimumMac,
            Path qcOutput, int qcThreads, Path qcCheckpoint) {
            this(format, genotypeField, missingPolicy, multiallelicPolicy, minimumMaf, minimumMac,
                qcOutput, qcThreads, qcCheckpoint, null, null, null,
                QGenomicRegions.Coordinates.ONE_BASED, FrequencyScope.ALIGNED);
        }

        public Options(Format format, GenotypeField genotypeField, MissingPolicy missingPolicy,
            MultiallelicPolicy multiallelicPolicy, double minimumMaf, double minimumMac,
            Path qcOutput, int qcThreads) {
            this(format, genotypeField, missingPolicy, multiallelicPolicy, minimumMaf, minimumMac,
                qcOutput, qcThreads, null);
        }

        public Options {
            if (format == null || genotypeField == null || missingPolicy == null || multiallelicPolicy == null)
                throw new IllegalArgumentException("Variant reader options must not be null");
            if (!Double.isFinite(minimumMaf) || minimumMaf < 0 || minimumMaf > 0.5)
                throw new IllegalArgumentException("min_maf must be between 0 and 0.5");
            if (!Double.isFinite(minimumMac) || minimumMac < 0)
                throw new IllegalArgumentException("min_mac must not be negative");
            if (qcThreads < 0)
                throw new IllegalArgumentException("variant_qc_threads must not be negative");
            qcOutput = qcOutput == null ? null : qcOutput.toAbsolutePath().normalize();
            qcCheckpoint = qcCheckpoint == null ? null : qcCheckpoint.toAbsolutePath().normalize();
            variantIndex = variantIndex == null ? null : variantIndex.toAbsolutePath().normalize();
            regionsFile = regionsFile == null ? null : regionsFile.toAbsolutePath().normalize();
            if (regionCoordinates == null || frequencyScope == null)
                throw new IllegalArgumentException("Region coordinates and frequency scope are required");
        }
    }

    public record Summary(long inputRecords, long includedVariants, long monomorphicVariants,
        long singletons, long doubletons, long excludedByMaf, long excludedByMac,
        long excludedForMissingness, long multiallelicRecords, int regionSets,
        List<String> emptyRegionSets) {
        public Summary {
            emptyRegionSets = emptyRegionSets == null ? List.of() : List.copyOf(emptyRegionSets);
        }
    }

    private interface VariantCursor extends AutoCloseable {
        VCFHeader header();
        boolean hasNext();
        VariantContext next() throws IOException;
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

    private static final class IndexedCursor implements VariantCursor {
        private final FeatureReader<VariantContext> reader;
        private final VCFHeader header;
        private final List<QGenomicRegions.QueryInterval> intervals;
        private CloseableTribbleIterator<VariantContext> iterator;
        private final Map<String, Integer> emittedIntervals = new HashMap<>();
        private int nextInterval;
        private int activeInterval = -1;
        private int activeIntervalEnd;
        private String activeContig;
        private VariantContext next;

        IndexedCursor(Path path, Path index, Format format,
            List<QGenomicRegions.QueryInterval> intervals)
            throws IOException {
            try {
                String featurePath = path.toUri().toString();
                String indexPath = index.toUri().toString();
                reader = format == Format.BCF
                    ? AbstractFeatureReader.getFeatureReader(featurePath, indexPath,
                        new QBcf22Codec(), true)
                    : AbstractFeatureReader.getFeatureReader(featurePath, indexPath,
                        new VCFCodec(), true);
                header = (VCFHeader) reader.getHeader();
                this.intervals = List.copyOf(intervals);
            } catch (RuntimeException e) {
                throw new IOException("Cannot open indexed variant source " + path
                    + " with index " + index + ": " + e.getMessage(), e);
            }
            try {
                advance();
            } catch (IOException | RuntimeException e) {
                try {
                    reader.close();
                } catch (IOException closeFailure) {
                    e.addSuppressed(closeFailure);
                }
                if (e instanceof IOException ioException)
                    throw ioException;
                throw new IOException("Cannot start indexed variant query for " + path, e);
            }
        }

        @Override public VCFHeader header() { return header; }

        @Override public boolean hasNext() { return next != null; }

        @Override
        public VariantContext next() throws IOException {
            if (next == null)
                throw new NoSuchElementException();
            VariantContext result = next;
            advance();
            return result;
        }

        private void advance() throws IOException {
            while (true) {
                if (iterator != null && iterator.hasNext()) {
                    VariantContext candidate = iterator.next();
                    String id = canonicalId(candidate);
                    Integer previousInterval = emittedIntervals.get(id);
                    if (previousInterval != null && previousInterval != activeInterval)
                        continue;
                    if (candidate.getEnd() > activeIntervalEnd)
                        emittedIntervals.putIfAbsent(id, activeInterval);
                    next = candidate;
                    return;
                }
                if (iterator != null) {
                    iterator.close();
                    iterator = null;
                }
                if (nextInterval >= intervals.size()) {
                    next = null;
                    return;
                }
                QGenomicRegions.QueryInterval interval = intervals.get(nextInterval++);
                activeInterval = nextInterval - 1;
                activeIntervalEnd = interval.end();
                if (!interval.contig().equals(activeContig)) {
                    emittedIntervals.clear();
                    activeContig = interval.contig();
                }
                try {
                    iterator = reader.query(interval.contig(), interval.start(), interval.end());
                } catch (RuntimeException e) {
                    throw new IOException("Indexed query failed for " + interval.contig() + ":"
                        + interval.start() + "-" + interval.end(), e);
                }
            }
        }

        @Override
        public void close() throws IOException {
            if (iterator != null)
                iterator.close();
            reader.close();
        }
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

    private record Evaluation(String rowId, int calledSamples, int missingSamples,
        double alternateAlleleCount, double alleleNumber, double dosageSumSquares,
        double eaf, double maf, double mac, double hweP, int hweSamples,
        int homRef, int heterozygous, int homAlt,
        String classification, boolean included, String exclusionReason,
        String regionSets, String alignedFrequencyReason) { }

    private record EvaluatedVariant(VariantContext variant, Evaluation evaluation) { }

    /** Compact input-record decisions retained for later sequential source passes. */
    private static final class VariantDecisions {
        private final List<double[]> effectAlleleFrequencyChunks = new ArrayList<>();
        private long size;

        void add(Evaluation evaluation) {
            int indexInChunk = (int) (size % DECISION_CHUNK_SIZE);
            if (indexInChunk == 0)
                effectAlleleFrequencyChunks.add(new double[DECISION_CHUNK_SIZE]);
            effectAlleleFrequencyChunks.get(effectAlleleFrequencyChunks.size() - 1)[indexInChunk]
                = evaluation.included() ? evaluation.eaf() : Double.NaN;
            size++;
        }

        boolean included(long inputRecord) throws IOException {
            return Double.isFinite(effectAlleleFrequency(inputRecord));
        }

        double effectAlleleFrequency(long inputRecord) throws IOException {
            if (inputRecord < 0 || inputRecord >= size)
                throw new IOException("Variant source has more records than the completed QC scan");
            int chunk = Math.toIntExact(inputRecord / DECISION_CHUNK_SIZE);
            int offset = (int) (inputRecord % DECISION_CHUNK_SIZE);
            return effectAlleleFrequencyChunks.get(chunk)[offset];
        }

        long size() {
            return size;
        }
    }

    private static final class ScanState {
        long inputRecords;
        long included;
        long monomorphic;
        long singletons;
        long doubletons;
        long excludedByMaf;
        long excludedByMac;
        long excludedForMissingness;
        long multiallelic;
        String lastRowId;
        final Set<String> rowIds = new HashSet<>();
        final Set<String> observedRegionSets = new LinkedHashSet<>();
        final VariantDecisions decisions = new VariantDecisions();
    }

    private static final class VariantProgress {
        private final String label;
        private final long totalRecords;
        private final long initialProcessedRecords;
        private final long startNanos = System.nanoTime();
        private long lastCheckRecords;
        private long lastReportRecords;
        private long lastReportNanos = startNanos;
        private boolean reported;
        private boolean completed;

        VariantProgress(String label, long totalRecords) {
            this(label, totalRecords, 0);
        }

        VariantProgress(String label, long totalRecords, long initialProcessedRecords) {
            this.label = label;
            this.totalRecords = totalRecords;
            this.initialProcessedRecords = initialProcessedRecords;
            lastCheckRecords = initialProcessedRecords;
            lastReportRecords = initialProcessedRecords;
        }

        void update(long processedRecords, long retainedRecords) {
            if (processedRecords - lastCheckRecords < PROGRESS_CHECK_RECORDS)
                return;
            lastCheckRecords = processedRecords;
            long now = System.nanoTime();
            if (processedRecords - lastReportRecords < PROGRESS_REPORT_RECORDS
                && now - lastReportNanos < PROGRESS_REPORT_NANOS)
                return;
            print(processedRecords, retainedRecords, now, false);
        }

        void complete(long processedRecords, long retainedRecords) {
            if (completed)
                return;
            completed = true;
            long now = System.nanoTime();
            if (reported || now - startNanos >= PROGRESS_REPORT_NANOS)
                print(processedRecords, retainedRecords, now, true);
        }

        private void print(long processedRecords, long retainedRecords, long now, boolean complete) {
            System.out.println(progressMessage(label, processedRecords, totalRecords,
                retainedRecords, now - startNanos, complete,
                Math.max(0, processedRecords - initialProcessedRecords)));
            System.out.flush();
            reported = true;
            lastReportRecords = processedRecords;
            lastReportNanos = now;
        }
    }

    private final class Reader implements BlockReader {
        private final VariantCursor cursor;
        private final int[] columnOrder;
        private final VariantProgress progress;
        private long inputRecords;
        private long nextRowOffset;
        private boolean closed;
        private boolean exhausted;

        Reader(int[] requestedColumnOrder) throws IOException {
            cursor = openCursor();
            columnOrder = normalizeColumnOrder(requestedColumnOrder, metadata.columnCount());
            progress = new VariantProgress("Variant input pass (" + path.getFileName() + ")",
                summary.inputRecords());
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
                VariantContext variant = cursor.next();
                long inputRecord = inputRecords;
                inputRecords++;
                if (decisions.included(inputRecord)) {
                    rowIds.add(canonicalId(variant));
                    rows.add(readSelectedValues(variant, columnOrder,
                        decisions.effectAlleleFrequency(inputRecord)));
                    nextRowOffset++;
                }
                progress.update(inputRecords, nextRowOffset);
            }
            if (!cursor.hasNext()) {
                exhausted = true;
                if (inputRecords != decisions.size())
                    throw new IOException("Variant source changed while it was being read: expected "
                        + decisions.size() + " input records, found " + inputRecords);
                progress.complete(inputRecords, nextRowOffset);
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
    private final int qcThreadCount;
    private final int qcCheckpointRecordBatchSize;
    private final Path qcCheckpointRoot;
    private final Path variantIndex;
    private final QGenomicRegions regions;
    private final List<QGenomicRegions.QueryInterval> wholeSourceIntervals;
    private final long sourceSize;
    private final long sourceModifiedMillis;
    private final long indexSize;
    private final long indexModifiedMillis;
    private Metadata metadata;
    private Summary summary;
    private VariantDecisions decisions;
    private int[] analysisSampleIndices;
    private int analysisSampleCount;
    private boolean scanComplete;
    private Path activeQcCheckpointDirectory;
    private long resumedQcRecords;
    private boolean indexedResumeUsed;

    public QVariantMatrixSource(Path path, Options options) throws IOException {
        this(path, options, false);
    }

    public static QVariantMatrixSource openForAlignment(Path path, Options options) throws IOException {
        return new QVariantMatrixSource(path, options, true);
    }

    private QVariantMatrixSource(Path path, Options options, boolean deferScan) throws IOException {
        this.path = path.toAbsolutePath().normalize();
        this.options = options;
        if (!Files.isRegularFile(this.path))
            throw new IOException("Variant file does not exist: " + this.path);
        sourceSize = Files.size(this.path);
        sourceModifiedMillis = Files.getLastModifiedTime(this.path).toMillis();
        VCFHeader sourceHeader;
        try (VariantCursor cursor = openSequentialCursor()) {
            validateSampleIds(cursor.header().getGenotypeSamples());
            sourceSampleIds = cursor.header().getGenotypeSamples().toArray(String[]::new);
            selectedField = selectField(cursor.header(), options.genotypeField());
            sourceHeader = cursor.header();
        }
        variantIndex = resolveVariantIndex(this.path, options.variantIndex());
        indexSize = variantIndex == null ? -1 : Files.size(variantIndex);
        indexModifiedMillis = variantIndex == null ? -1
            : Files.getLastModifiedTime(variantIndex).toMillis();
        List<String> availableContigs = sourceHeader.getContigLines().stream()
            .map(line -> line.getID()).toList();
        List<QGenomicRegions.QueryInterval> wholeIntervals = new ArrayList<>();
        sourceHeader.getContigLines().forEach(line -> {
            int length = line.getSAMSequenceRecord().getSequenceLength();
            wholeIntervals.add(new QGenomicRegions.QueryInterval(line.getID(), 1,
                length > 0 ? length : Integer.MAX_VALUE));
        });
        wholeSourceIntervals = List.copyOf(wholeIntervals);
        regions = QGenomicRegions.load(options.regions(), options.regionsFile(),
            options.regionCoordinates(), availableContigs);
        if (regions != null && variantIndex == null)
            throw new IOException("Region-limited VCF/BCF access requires a .tbi or .csi index; "
                + "provide --variant-index or place the index beside " + this.path);
        qcThreadCount = recommendQcThreads(options.qcThreads());
        qcCheckpointRecordBatchSize = qcCheckpointBatchSize();
        qcCheckpointRoot = options.qcCheckpoint() != null ? options.qcCheckpoint()
            : options.qcOutput() == null ? null
            : Path.of(options.qcOutput().toString() + ".checkpoint").toAbsolutePath().normalize();
        if (qcCheckpointRoot != null && this.path.equals(qcCheckpointRoot))
            throw new IOException("Variant-QC checkpoint directory must differ from the genotype input: "
                + this.path);
        analysisSampleIndices = identityColumnOrder(sourceSampleIds.length);
        analysisSampleCount = sourceSampleIds.length;
        if (deferScan) {
            metadata = new Metadata(this.path, -1, sourceSampleIds.length,
                sourceSampleIds, "variant-unscanned");
        } else {
            completeScan();
        }
    }

    @Override public Metadata metadata() { return metadata; }
    public Summary summary() {
        ensureScanned();
        return summary;
    }
    public GenotypeField selectedField() { return selectedField; }
    public int analysisSampleCount() { return analysisSampleCount; }
    public int qcThreadCount() { return qcThreadCount; }
    public Path qcCheckpointDirectory() { return activeQcCheckpointDirectory; }
    public long resumedQcRecords() { return resumedQcRecords; }
    public boolean indexedResumeUsed() { return indexedResumeUsed; }
    public Path variantIndex() { return variantIndex; }
    public QGenomicRegions regions() { return regions; }

    public void selectAnalysisSamples(int[] columnOrder) throws IOException {
        if (scanComplete)
            throw new IllegalStateException("Variant QC sample selection is already finalized");
        int[] selected = normalizeColumnOrder(columnOrder, sourceSampleIds.length);
        analysisSampleIndices = selected;
        analysisSampleCount = selected.length;
        completeScan();
    }

    @Override
    public BlockReader open(int[] columnOrder) throws IOException {
        ensureScanned();
        return new Reader(columnOrder);
    }

    private void completeScan() throws IOException {
        ScanResult scan = scanMetadata();
        metadata = scan.metadata();
        summary = scan.summary();
        decisions = scan.decisions();
        scanComplete = true;
    }

    private void ensureScanned() {
        if (!scanComplete)
            throw new IllegalStateException("Variant QC must be finalized after sample alignment before reading rows");
    }

    private record ScanResult(Metadata metadata, Summary summary, VariantDecisions decisions) { }

    private ScanResult scanMetadata() throws IOException {
        try (QVariantQcCheckpoint checkpoint = openQcCheckpoint()) {
            return scanMetadata(checkpoint);
        }
    }

    private ScanResult scanMetadata(QVariantQcCheckpoint checkpoint) throws IOException {
        validateSourceIdentity();
        Path qcOutput = options.qcOutput();
        ScanState state = new ScanState();
        QVariantQcCheckpoint.Snapshot checkpointSnapshot = checkpoint == null
            ? new QVariantQcCheckpoint.Snapshot(0, false)
            : checkpoint.load(line -> acceptEvaluation(parseQcEvaluation(line), state));
        resumedQcRecords = checkpointSnapshot.records();
        if (checkpoint != null)
            activeQcCheckpointDirectory = checkpoint.directory();
        if (checkpointSnapshot.complete()) {
            System.out.println("Reusing completed variant-QC checkpoint with "
                + String.format(Locale.ROOT, "%,d", state.inputRecords) + " records: "
                + checkpoint.directory());
            System.out.flush();
            if (qcOutput != null)
                checkpoint.materialize(qcOutput, QC_HEADER);
            return scanResult(state, sourceSampleIds);
        }

        if (state.inputRecords > 0) {
            System.out.println("Resuming variant QC after "
                + String.format(Locale.ROOT, "%,d", state.inputRecords)
                + " completed records from " + checkpoint.directory());
            if (variantIndex != null && (!wholeSourceIntervals.isEmpty() || regions != null))
                System.out.println("A matching variant index is available; resume will seek to and verify the durable boundary.");
            else
                System.out.println("Sequential VCF/BCF input must decode past completed records; their sample QC and HWE will not be recomputed.");
            System.out.flush();
        }

        VariantProgress progress = new VariantProgress("Variant QC (" + path.getFileName() + ")",
            -1, state.inputRecords);
        ExecutorService executor = qcThreadCount > 1
            ? Executors.newFixedThreadPool(qcThreadCount) : null;
        Deque<Future<EvaluatedVariant>> pending = new ArrayDeque<>();
        int maximumPending = Math.max(1, qcThreadCount * QC_TASKS_PER_THREAD);
        List<String> checkpointBatch = checkpoint == null ? null
            : new ArrayList<>(qcCheckpointRecordBatchSize);
        try {
            try (VariantCursor cursor = openCursorForQcResume(checkpoint, state)) {
                while (cursor.hasNext()) {
                    VariantContext variant = cursor.next();
                    if (executor == null) {
                        acceptQcEvaluation(new EvaluatedVariant(variant, evaluateQc(variant)),
                            state, checkpoint, checkpointBatch, progress);
                    } else {
                        forceGenotypeDecoding(variant);
                        pending.addLast(executor.submit(
                            () -> new EvaluatedVariant(variant, evaluateQc(variant))));
                        if (pending.size() >= maximumPending)
                            acceptQcEvaluation(awaitQcEvaluation(pending.removeFirst()),
                                state, checkpoint, checkpointBatch, progress);
                    }
                }
                while (!pending.isEmpty())
                    acceptQcEvaluation(awaitQcEvaluation(pending.removeFirst()),
                        state, checkpoint, checkpointBatch, progress);
            }
            progress.complete(state.inputRecords, state.included);
            if (checkpoint != null) {
                checkpoint.commit(checkpointBatch);
                checkpointBatch.clear();
                checkpoint.markComplete(state.inputRecords);
                if (qcOutput != null)
                    checkpoint.materialize(qcOutput, QC_HEADER);
            }
        } finally {
            shutdownQcExecutor(executor);
        }
        return scanResult(state, sourceSampleIds);
    }

    private void acceptQcEvaluation(EvaluatedVariant result, ScanState state,
        QVariantQcCheckpoint checkpoint, List<String> checkpointBatch,
        VariantProgress progress) throws IOException {
        Evaluation evaluation = result.evaluation();
        acceptEvaluation(evaluation, state);
        if (checkpoint != null) {
            checkpointBatch.add(formatQc(result.variant(), evaluation));
            if (checkpointBatch.size() >= qcCheckpointRecordBatchSize) {
                checkpoint.commit(checkpointBatch);
                checkpointBatch.clear();
            }
        }
        progress.update(state.inputRecords, state.included);
        maybeFailAfterCheckpoint(state.inputRecords);
    }

    private void acceptEvaluation(Evaluation evaluation, ScanState state) throws IOException {
        state.inputRecords++;
        state.lastRowId = evaluation.rowId();
        if (evaluation.classification().equals("monomorphic")) state.monomorphic++;
        if (evaluation.classification().equals("singleton")) state.singletons++;
        if (evaluation.classification().equals("doubleton")) state.doubletons++;
        if (evaluation.classification().equals("multiallelic")) state.multiallelic++;
        if (evaluation.alignedFrequencyReason().contains("maf")) state.excludedByMaf++;
        if (evaluation.alignedFrequencyReason().contains("mac")) state.excludedByMac++;
        if (evaluation.exclusionReason().contains("missing")) state.excludedForMissingness++;
        if (!evaluation.regionSets().equals("."))
            Collections.addAll(state.observedRegionSets, evaluation.regionSets().split(";"));
        if (!state.rowIds.add(evaluation.rowId()))
            throw new IOException("Duplicate canonical variant identifier '"
                + evaluation.rowId() + "' in " + path);
        if (evaluation.included())
            state.included++;
        state.decisions.add(evaluation);
    }

    private ScanResult scanResult(ScanState state, String[] samples) throws IOException {
        if (state.included == 0)
            throw new IOException("No variants remain after QC/filtering in " + path);
        String signatureTag = String.format(Locale.ROOT,
            "variant-v3;format=%s;field=%s;missing=%s;multiallelic=%s;min-maf=%.17g;min-mac=%.17g;frequency-scope=%s;qc-samples=%d;regions=%s;index=%s:%d:%d",
            options.format(), selectedField, options.missingPolicy(), options.multiallelicPolicy(),
            options.minimumMaf(), options.minimumMac(), options.frequencyScope(), analysisSampleCount,
            regions == null ? "all" : regions.signature(),
            variantIndex == null ? "none" : variantIndex, indexSize, indexModifiedMillis);
        Metadata scanned = new Metadata(path, state.included, samples.length, samples, signatureTag);
        Summary counts = new Summary(state.inputRecords, state.included, state.monomorphic,
            state.singletons, state.doubletons, state.excludedByMaf, state.excludedByMac,
            state.excludedForMissingness, state.multiallelic,
            regions == null ? 0 : regions.setIds().size(),
            regions == null ? List.of() : regions.emptySetIds(state.observedRegionSets));
        return new ScanResult(scanned, counts, state.decisions);
    }

    private QVariantQcCheckpoint openQcCheckpoint() throws IOException {
        if (options.qcOutput() != null && path.equals(options.qcOutput()))
            throw new IOException("Variant QC output must differ from the genotype input: " + path);
        if (qcCheckpointRoot == null)
            return null;
        String signature = qcCheckpointSignature();
        String description = String.format(Locale.ROOT,
            "input=%s; format=%s; field=%s; missing=%s; multiallelic=%s; min_maf=%.17g; min_mac=%.17g; frequency_scope=%s; aligned_samples=%d; regions=%s",
            path, options.format(), selectedField, options.missingPolicy(),
            options.multiallelicPolicy(), options.minimumMaf(), options.minimumMac(),
            options.frequencyScope(), analysisSampleCount,
            regions == null ? "all" : regions.signature());
        return new QVariantQcCheckpoint(qcCheckpointRoot, signature, description);
    }

    private String qcCheckpointSignature() throws IOException {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is unavailable", e);
        }
        updateDigest(digest, "gpu-eqtl-variant-qc-checkpoint-v3");
        updateDigest(digest, path.toString());
        updateDigest(digest, Long.toString(sourceSize));
        updateDigest(digest, Long.toString(sourceModifiedMillis));
        updateDigest(digest, options.format().name());
        updateDigest(digest, selectedField.name());
        updateDigest(digest, options.missingPolicy().name());
        updateDigest(digest, options.multiallelicPolicy().name());
        updateDigest(digest, Double.toHexString(options.minimumMaf()));
        updateDigest(digest, Double.toHexString(options.minimumMac()));
        updateDigest(digest, options.frequencyScope().name());
        updateDigest(digest, regions == null ? "all-regions" : regions.signature());
        if (variantIndex != null) {
            updateDigest(digest, variantIndex.toString());
            updateDigest(digest, Long.toString(indexSize));
            updateDigest(digest, Long.toString(indexModifiedMillis));
        } else {
            updateDigest(digest, "no-index");
        }
        updateDigest(digest, Integer.toString(sourceSampleIds.length));
        for (String sampleId : sourceSampleIds)
            updateDigest(digest, sampleId);
        updateDigest(digest, Integer.toString(analysisSampleIndices.length));
        for (int index : analysisSampleIndices) {
            updateDigest(digest, Integer.toString(index));
            updateDigest(digest, sourceSampleIds[index]);
        }
        return HexFormat.of().formatHex(digest.digest());
    }

    private static void updateDigest(MessageDigest digest, String value) {
        digest.update(value.getBytes(StandardCharsets.UTF_8));
        digest.update((byte) 0);
    }

    private void validateCheckpointPrefix(VariantCursor cursor, QVariantQcCheckpoint checkpoint,
        long completedRecords) throws IOException {
        if (checkpoint == null || completedRecords == 0)
            return;
        long[] validated = {0};
        checkpoint.forEachLine(line -> {
            if (!cursor.hasNext())
                throw new IOException("Variant input ended while validating resumed checkpoint after "
                    + validated[0] + " records: " + path);
            VariantContext variant = cursor.next();
            String expected = qcRowId(line);
            String observed = canonicalId(variant);
            if (!expected.equals(observed))
                throw new IOException("Variant input no longer matches resumed checkpoint at record "
                    + (validated[0] + 1) + ": expected " + expected + ", found " + observed);
            validated[0]++;
        });
        if (validated[0] != completedRecords)
            throw new IOException("Variant-QC checkpoint loaded " + completedRecords
                + " records but prefix validation found " + validated[0]);
    }

    private VariantCursor openCursorForQcResume(QVariantQcCheckpoint checkpoint,
        ScanState state) throws IOException {
        if (checkpoint == null || state.inputRecords == 0)
            return openCursor();
        List<QGenomicRegions.QueryInterval> intervals = regions == null
            ? wholeSourceIntervals : regions.queryIntervals();
        if (variantIndex != null && !intervals.isEmpty() && state.lastRowId != null) {
            validateSourceIdentity();
            VariantCursor indexed = indexedCursorAfter(state.lastRowId, intervals);
            indexedResumeUsed = true;
            System.out.println("Indexed QC resume sought directly to durable boundary "
                + state.lastRowId + " using " + variantIndex);
            System.out.flush();
            return indexed;
        }
        VariantCursor sequential = openCursor();
        validateCheckpointPrefix(sequential, checkpoint, state.inputRecords);
        return sequential;
    }

    private VariantCursor indexedCursorAfter(String boundaryId,
        List<QGenomicRegions.QueryInterval> intervals) throws IOException {
        Boundary boundary = Boundary.parse(boundaryId);
        List<QGenomicRegions.QueryInterval> remaining = new ArrayList<>();
        boolean foundContig = false;
        for (QGenomicRegions.QueryInterval interval : intervals) {
            if (!foundContig) {
                if (!interval.contig().equals(boundary.contig()) || interval.end() < boundary.position())
                    continue;
                foundContig = true;
                remaining.add(new QGenomicRegions.QueryInterval(interval.contig(),
                    Math.max(interval.start(), boundary.position()), interval.end()));
            } else {
                remaining.add(interval);
            }
        }
        if (remaining.isEmpty())
            throw new IOException("Indexed QC resume boundary is outside the selected source regions: "
                + boundaryId);
        IndexedCursor cursor = new IndexedCursor(path, variantIndex, options.format(), remaining);
        try {
            while (cursor.hasNext()) {
                VariantContext variant = cursor.next();
                if (canonicalId(variant).equals(boundaryId))
                    return cursor;
            }
        } catch (IOException | RuntimeException e) {
            cursor.close();
            throw e;
        }
        cursor.close();
        throw new IOException("Indexed QC resume could not find durable boundary " + boundaryId
            + " in " + path);
    }

    private record Boundary(String contig, int position) {
        static Boundary parse(String canonicalId) throws IOException {
            int first = canonicalId.indexOf(':');
            int second = first < 0 ? -1 : canonicalId.indexOf(':', first + 1);
            if (first <= 0 || second <= first + 1)
                throw new IOException("Malformed canonical variant boundary " + canonicalId);
            try {
                return new Boundary(canonicalId.substring(0, first),
                    Integer.parseInt(canonicalId.substring(first + 1, second)));
            } catch (NumberFormatException e) {
                throw new IOException("Malformed canonical variant boundary " + canonicalId, e);
            }
        }
    }

    private static String qcRowId(String line) throws IOException {
        int tab = line.indexOf('\t');
        if (tab <= 0)
            throw new IOException("Malformed variant-QC checkpoint record");
        return line.substring(0, tab);
    }

    private static Evaluation parseQcEvaluation(String line) throws IOException {
        String[] fields = line.split("\\t", -1);
        if (fields.length != 27)
            throw new IOException("Malformed variant-QC checkpoint record: expected 27 fields, found "
                + fields.length);
        try {
            boolean included;
            if (fields[22].equals("true")) included = true;
            else if (fields[22].equals("false")) included = false;
            else throw new IOException("Malformed included flag in variant-QC checkpoint: " + fields[22]);
            return new Evaluation(fields[0], Integer.parseInt(fields[8]), Integer.parseInt(fields[9]),
                checkpointNumber(fields[10]), checkpointNumber(fields[11]), checkpointNumber(fields[12]),
                checkpointNumber(fields[13]), checkpointNumber(fields[14]), checkpointNumber(fields[15]),
                checkpointNumber(fields[16]), Integer.parseInt(fields[17]), Integer.parseInt(fields[18]),
                Integer.parseInt(fields[19]), Integer.parseInt(fields[20]), fields[21], included,
                fields[23].equals(".") ? "" : fields[23], fields[24],
                fields[26].equals(".") ? "" : fields[26]);
        } catch (NumberFormatException e) {
            throw new IOException("Malformed numeric field in variant-QC checkpoint record " + fields[0], e);
        }
    }

    private static double checkpointNumber(String value) {
        return value.equals("NA") ? Double.NaN : Double.parseDouble(value);
    }

    private static int qcCheckpointBatchSize() {
        String configured = System.getProperty(QC_CHECKPOINT_RECORDS_PROPERTY);
        if (configured == null || configured.isBlank())
            return DEFAULT_QC_CHECKPOINT_RECORDS;
        try {
            int records = Integer.parseInt(configured);
            if (records <= 0)
                throw new IllegalArgumentException(QC_CHECKPOINT_RECORDS_PROPERTY + " must be positive");
            return records;
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(QC_CHECKPOINT_RECORDS_PROPERTY
                + " must be a positive integer", e);
        }
    }

    private static void maybeFailAfterCheckpoint(long records) throws IOException {
        String configured = System.getProperty(TEST_FAIL_AFTER_PROPERTY);
        if (configured == null || configured.isBlank())
            return;
        try {
            if (records == Long.parseLong(configured))
                throw new IOException("Injected variant-QC interruption after " + records + " records");
        } catch (NumberFormatException e) {
            throw new IOException("Invalid test interruption record count", e);
        }
    }

    private static EvaluatedVariant awaitQcEvaluation(Future<EvaluatedVariant> future)
        throws IOException {
        try {
            return future.get();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted while evaluating variant QC", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof IOException ioException)
                throw ioException;
            if (cause instanceof RuntimeException runtimeException)
                throw runtimeException;
            if (cause instanceof Error error)
                throw error;
            throw new IOException("Variant QC worker failed", cause);
        }
    }

    private static void forceGenotypeDecoding(VariantContext variant) {
        if (variant.getNAlleles() == 2 && variant.getNSamples() > 0)
            variant.getGenotype(0);
    }

    private static void shutdownQcExecutor(ExecutorService executor) {
        if (executor == null)
            return;
        executor.shutdownNow();
        try {
            if (!executor.awaitTermination(1, TimeUnit.MINUTES))
                System.err.println("WARNING: Variant QC workers did not terminate within one minute.");
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    private double[] readSelectedValues(VariantContext variant, int[] columnOrder,
        double effectAlleleFrequency) throws IOException {
        if (variant.getNAlleles() != 2)
            throw new IOException("Variant source changed after QC at " + variant.getContig()
                + ":" + variant.getStart());
        if (variant.getNSamples() != sourceSampleIds.length)
            throw new IOException("Variant sample count changed at " + variant.getContig()
                + ":" + variant.getStart());
        Allele alternate = variant.getAlternateAllele(0);
        double[] values = new double[columnOrder.length];
        for (int outputColumn = 0; outputColumn < columnOrder.length; outputColumn++) {
            int sourceColumn = columnOrder[outputColumn];
            Genotype genotype = variant.getGenotype(sourceSampleIds[sourceColumn]);
            Double dosage = selectedField == GenotypeField.DS
                ? dosageFromDs(genotype, variant) : dosageFromGt(genotype, alternate, variant);
            if (dosage != null) {
                values[outputColumn] = dosage;
            } else if (options.missingPolicy() == MissingPolicy.ZERO) {
                values[outputColumn] = 0;
            } else if (options.missingPolicy() == MissingPolicy.PRESERVE) {
                values[outputColumn] = kUndefinedValue;
            } else {
                values[outputColumn] = 2 * effectAlleleFrequency;
            }
        }
        return values;
    }

    private Evaluation evaluateQc(VariantContext variant) throws IOException {
        if (variant.getNAlleles() != 2) {
            if (options.multiallelicPolicy() == MultiallelicPolicy.ERROR)
                throw new IOException("Multiallelic variant encountered at " + variant.getContig()
                    + ":" + variant.getStart() + "; use --multiallelic exclude");
            return excludedEvaluation(variant, "multiallelic", "multiallelic");
        }

        Allele reference = variant.getReference();
        Allele alternate = variant.getAlternateAllele(0);
        int sampleCount = sourceSampleIds.length;
        double alternateAlleleCount = 0;
        double dosageSumSquares = 0;
        int calledSamples = 0;
        int missingSamples = 0;
        int homRef = 0;
        int heterozygous = 0;
        int homAlt = 0;
        int hweSamples = 0;
        if (variant.getNSamples() != sampleCount)
            throw new IOException("Variant sample count changed at " + variant.getContig()
                + ":" + variant.getStart());

        for (int sample : analysisSampleIndices) {
            Genotype genotype = variant.getGenotype(sourceSampleIds[sample]);
            Double dosage = selectedField == GenotypeField.DS
                ? dosageFromDs(genotype, variant) : dosageFromGt(genotype, alternate, variant);
            if (dosage == null) {
                missingSamples++;
            } else {
                alternateAlleleCount += dosage;
                dosageSumSquares += dosage * dosage;
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
            return new Evaluation(rowId, 0, missingSamples, 0, 0,
                0, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0,
                0, 0, 0, "no_calls", false,
                "no_called_genotypes", regionMemberships(variant), "");
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
        List<String> frequencyReasons = new ArrayList<>();
        if (maf + INTEGER_TOLERANCE < options.minimumMaf())
            frequencyReasons.add("maf_below_minimum");
        if (mac + INTEGER_TOLERANCE < options.minimumMac())
            frequencyReasons.add("mac_below_minimum");
        if (options.frequencyScope() == FrequencyScope.ALIGNED)
            reasons.addAll(frequencyReasons);
        boolean accepted = reasons.isEmpty();
        return new Evaluation(rowId, calledSamples, missingSamples, alternateAlleleCount,
            alleleNumber, dosageSumSquares, eaf, maf, mac, hweP, hweSamples,
            homRef, heterozygous, homAlt, classification, accepted,
            String.join(";", reasons), regionMemberships(variant),
            String.join(";", frequencyReasons));
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

    static String progressMessage(String label, long processedRecords, long totalRecords,
        long retainedRecords, long elapsedNanos, boolean complete) {
        return progressMessage(label, processedRecords, totalRecords, retainedRecords,
            elapsedNanos, complete, processedRecords);
    }

    private static String progressMessage(String label, long processedRecords, long totalRecords,
        long retainedRecords, long elapsedNanos, boolean complete, long timedRecords) {
        double elapsedSeconds = Math.max(0, elapsedNanos) / 1_000_000_000.0;
        double recordsPerSecond = elapsedSeconds > 0 ? timedRecords / elapsedSeconds : 0;
        StringBuilder message = new StringBuilder(label)
            .append(complete ? " complete: " : " progress: ");
        if (totalRecords > 0) {
            double percent = Math.min(100, 100.0 * processedRecords / totalRecords);
            message.append(String.format(Locale.ROOT, "%,d / %,d records (%.1f%%)",
                processedRecords, totalRecords, percent));
        } else {
            message.append(String.format(Locale.ROOT, "%,d records scanned", processedRecords));
        }
        message.append(String.format(Locale.ROOT, "; retained=%,d; elapsed=%s; rate=%,.1f records/s",
            retainedRecords, duration(elapsedSeconds), recordsPerSecond));
        if (!complete && totalRecords > processedRecords && recordsPerSecond > 0) {
            double etaSeconds = (totalRecords - processedRecords) / recordsPerSecond;
            message.append("; ETA=").append(duration(etaSeconds));
        }
        return message.toString();
    }

    private static String duration(double seconds) {
        if (seconds < 60)
            return String.format(Locale.ROOT, "%.1f s", seconds);
        if (seconds < 3600)
            return String.format(Locale.ROOT, "%.1f min", seconds / 60);
        return String.format(Locale.ROOT, "%.2f h", seconds / 3600);
    }

    private Evaluation excludedEvaluation(VariantContext variant, String classification, String reason) {
        return new Evaluation(canonicalId(variant), 0, 0, Double.NaN, Double.NaN,
            Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN, 0,
            0, 0, 0, classification, false, reason,
            regionMemberships(variant), "");
    }

    private static String canonicalId(VariantContext variant) {
        String alternate = variant.getAlternateAlleles().isEmpty() ? "."
            : String.join(",", variant.getAlternateAlleles().stream().map(Allele::getDisplayString).toList());
        return variant.getContig() + ":" + variant.getStart() + ":"
            + variant.getReference().getDisplayString() + ":" + alternate;
    }

    private String formatQc(VariantContext variant, Evaluation evaluation) {
        String rsId = variant.hasID() ? variant.getID() : ".";
        String alternate = variant.getAlternateAlleles().isEmpty() ? "."
            : String.join(",", variant.getAlternateAlleles().stream().map(Allele::getDisplayString).toList());
        String filter = variant.isFiltered() ? String.join(";", variant.getFilters())
            : (variant.filtersWereApplied() ? "PASS" : ".");
        return String.join("\t",
            tsv(evaluation.rowId()),
            tsv(variant.getContig()),
            Long.toString(variant.getStart()),
            tsv(rsId),
            tsv(variant.getReference().getDisplayString()),
            tsv(alternate),
            tsv(filter),
            selectedField.name(),
            Integer.toString(evaluation.calledSamples()),
            Integer.toString(evaluation.missingSamples()),
            number(evaluation.alternateAlleleCount()),
            number(evaluation.alleleNumber()),
            number(evaluation.dosageSumSquares()),
            number(evaluation.eaf()),
            number(evaluation.maf()),
            number(evaluation.mac()),
            number(evaluation.hweP()),
            Integer.toString(evaluation.hweSamples()),
            Integer.toString(evaluation.homRef()),
            Integer.toString(evaluation.heterozygous()),
            Integer.toString(evaluation.homAlt()),
            evaluation.classification(),
            Boolean.toString(evaluation.included()),
            evaluation.exclusionReason().isEmpty() ? "." : evaluation.exclusionReason(),
            evaluation.regionSets(),
            options.frequencyScope().name().toLowerCase(Locale.ROOT),
            evaluation.alignedFrequencyReason().isEmpty() ? "." : evaluation.alignedFrequencyReason());
    }

    private String regionMemberships(VariantContext variant) {
        return regions == null ? "." : regions.memberships(variant);
    }

    private VariantCursor openCursor() throws IOException {
        validateSourceIdentity();
        if (regions != null)
            return new IndexedCursor(path, variantIndex, options.format(), regions.queryIntervals());
        return openSequentialCursor();
    }

    private VariantCursor openSequentialCursor() throws IOException {
        return options.format() == Format.BCF ? new BcfCursor(path) : new VcfCursor(path);
    }

    private void validateSourceIdentity() throws IOException {
        if (Files.size(path) != sourceSize
            || Files.getLastModifiedTime(path).toMillis() != sourceModifiedMillis)
            throw new IOException("Variant source changed during analysis: " + path);
        if (variantIndex != null && (Files.size(variantIndex) != indexSize
            || Files.getLastModifiedTime(variantIndex).toMillis() != indexModifiedMillis))
            throw new IOException("Variant index changed during analysis: " + variantIndex);
    }

    private static Path resolveVariantIndex(Path variantPath, Path configured) throws IOException {
        if (configured != null) {
            if (!Files.isRegularFile(configured))
                throw new IOException("Variant index does not exist: " + configured);
            return configured;
        }
        for (String suffix : List.of(".csi", ".tbi", ".idx")) {
            Path candidate = Path.of(variantPath.toString() + suffix).toAbsolutePath().normalize();
            if (Files.isRegularFile(candidate))
                return candidate;
        }
        return null;
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
        if (requested == null)
            return identityColumnOrder(columnCount);
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

    private static int[] identityColumnOrder(int columnCount) {
        int[] identity = new int[columnCount];
        for (int i = 0; i < columnCount; i++)
            identity[i] = i;
        return identity;
    }

    static int recommendQcThreads(int configuredThreads) {
        if (configuredThreads < 0)
            throw new IllegalArgumentException("variant_qc_threads must not be negative");
        if (configuredThreads > 0)
            return configuredThreads;
        int processors = Runtime.getRuntime().availableProcessors();
        if (processors <= 1)
            return 1;
        return Math.min(MAX_AUTOMATIC_QC_THREADS, Math.max(2, processors - 1));
    }

    private static String number(double value) {
        return Double.isFinite(value) ? Double.toString(value) : "NA";
    }

    private static String tsv(String value) {
        return value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }

}
