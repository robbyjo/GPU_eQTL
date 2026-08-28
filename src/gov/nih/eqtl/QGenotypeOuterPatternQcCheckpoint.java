/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Properties;
import java.util.zip.CheckedInputStream;
import java.util.zip.CheckedOutputStream;
import java.util.zip.CRC32;

/** Durable source-ordered genotype-block aggregates for exact trait-pattern variant QC. */
final class QGenotypeOuterPatternQcCheckpoint {
    private static final String FORMAT = "gpu-eqtl-genotype-outer-pattern-qc-v1";
    private static final String MANIFEST = "manifest.properties";
    private static final int MAGIC = 0x51475051; // QGPQ
    private static final int VERSION = 1;

    record BlockCounts(long[] inputVariants, long[] includedVariants,
        long[] monomorphicVariants, long[] belowMinimumMaf, long[] belowMinimumMac,
        long[] noCallVariants, long[] missingGenotypes) {
        BlockCounts(int patternCount) {
            this(new long[patternCount], new long[patternCount], new long[patternCount],
                new long[patternCount], new long[patternCount], new long[patternCount],
                new long[patternCount]);
        }

        BlockCounts {
            int patterns = inputVariants.length;
            if (includedVariants.length != patterns || monomorphicVariants.length != patterns
                || belowMinimumMaf.length != patterns || belowMinimumMac.length != patterns
                || noCallVariants.length != patterns || missingGenotypes.length != patterns)
                throw new IllegalArgumentException("Pattern-QC count arrays have different lengths");
        }

        int patternCount() { return inputVariants.length; }
    }

    private final Path directory;
    private final int totalBlocks;
    private final int patternCount;
    private final boolean keepAfterSuccess;
    private final boolean reusedAtStart;

    private QGenotypeOuterPatternQcCheckpoint(Path directory, int totalBlocks,
        int patternCount, boolean keepAfterSuccess, boolean reusedAtStart) {
        this.directory = directory;
        this.totalBlocks = totalBlocks;
        this.patternCount = patternCount;
        this.keepAfterSuccess = keepAfterSuccess;
        this.reusedAtStart = reusedAtStart;
    }

    static QGenotypeOuterPatternQcCheckpoint open(Path directory, String signature,
        int totalBlocks, int patternCount, boolean resume, boolean keepAfterSuccess)
        throws IOException {
        Path normalized = directory.toAbsolutePath().normalize();
        Path manifest = normalized.resolve(MANIFEST);
        if (Files.exists(normalized)) {
            if (!resume)
                throw new IOException("Pattern-QC checkpoint directory already exists; use --resume "
                    + "or choose another directory: " + normalized);
            Properties properties = load(manifest);
            if (!FORMAT.equals(properties.getProperty("format"))
                || !signature.equals(properties.getProperty("signature"))
                || totalBlocks != Integer.parseInt(properties.getProperty("total_blocks", "-1"))
                || patternCount != Integer.parseInt(properties.getProperty("pattern_count", "-1")))
                throw new IOException("Pattern-QC checkpoint does not match this analysis: "
                    + normalized);
            System.out.println("Resuming pattern-QC checkpoint: " + normalized);
        } else {
            Files.createDirectories(normalized);
            Properties properties = new Properties();
            properties.setProperty("format", FORMAT);
            properties.setProperty("signature", signature);
            properties.setProperty("total_blocks", Integer.toString(totalBlocks));
            properties.setProperty("pattern_count", Integer.toString(patternCount));
            Path temporary = Files.createTempFile(normalized, "manifest", ".partial");
            try (var output = new BufferedOutputStream(Files.newOutputStream(temporary))) {
                properties.store(output, "GPU eQTL genotype-outer pattern-QC checkpoint");
            }
            moveAtomically(temporary, manifest);
            System.out.println("Pattern-QC checkpoint directory: " + normalized);
        }
        QGenotypeOuterPatternQcCheckpoint checkpoint =
            new QGenotypeOuterPatternQcCheckpoint(normalized, totalBlocks, patternCount,
                keepAfterSuccess, false);
        boolean allComplete = totalBlocks > 0;
        for (int block = 0; block < totalBlocks; block++)
            allComplete &= checkpoint.isComplete(block);
        return new QGenotypeOuterPatternQcCheckpoint(normalized, totalBlocks, patternCount,
            keepAfterSuccess, allComplete);
    }

    boolean isComplete(int blockNumber) throws IOException {
        validateBlock(blockNumber);
        Path part = part(blockNumber);
        if (!Files.isRegularFile(part)) return false;
        readBlock(blockNumber);
        return true;
    }

    void writeBlock(int blockNumber, BlockCounts counts) throws IOException {
        validateBlock(blockNumber);
        if (counts.patternCount() != patternCount)
            throw new IOException("Pattern-QC block has " + counts.patternCount()
                + " patterns; expected " + patternCount);
        Path partial = partial(blockNumber);
        Files.deleteIfExists(partial);
        boolean complete = false;
        CheckedOutputStream checked = new CheckedOutputStream(
            new BufferedOutputStream(Files.newOutputStream(partial)), new CRC32());
        try (DataOutputStream output = new DataOutputStream(checked)) {
            output.writeInt(MAGIC);
            output.writeInt(VERSION);
            output.writeInt(blockNumber);
            output.writeInt(patternCount);
            for (int pattern = 0; pattern < patternCount; pattern++) {
                output.writeLong(counts.inputVariants()[pattern]);
                output.writeLong(counts.includedVariants()[pattern]);
                output.writeLong(counts.monomorphicVariants()[pattern]);
                output.writeLong(counts.belowMinimumMaf()[pattern]);
                output.writeLong(counts.belowMinimumMac()[pattern]);
                output.writeLong(counts.noCallVariants()[pattern]);
                output.writeLong(counts.missingGenotypes()[pattern]);
            }
            output.writeLong(checked.getChecksum().getValue());
            output.flush();
            complete = true;
        } finally {
            if (complete) moveAtomically(partial, part(blockNumber));
            else Files.deleteIfExists(partial);
        }
    }

    void assemble(Path outputPath, QTraitPatternModelSet models, double minimumMaf,
        double minimumMac) throws IOException {
        BlockCounts totals = new BlockCounts(patternCount);
        for (int block = 0; block < totalBlocks; block++)
            add(totals, readBlock(block));
        Path output = outputPath.toAbsolutePath().normalize();
        Path parent = output.getParent();
        if (parent != null) Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, output.getFileName().toString(), ".partial");
        boolean complete = false;
        try (BufferedWriter writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {
            writer.write("pattern_id\ttrait_rows\tobserved_samples\tinput_variants"
                + "\tincluded_variants\tmonomorphic_variants\tbelow_min_maf\tbelow_min_mac"
                + "\tno_call_variants\tmissing_genotypes\tfrequency_scope\tstatistics_cache\treused");
            writer.newLine();
            for (QTraitPatternModelSet.Model model : models.models()) {
                if (!model.estimable()) continue;
                int id = model.id;
                writer.write(id + "\t" + model.traitRows + "\t" + model.observed.length
                    + "\t" + totals.inputVariants()[id]
                    + "\t" + totals.includedVariants()[id]
                    + "\t" + totals.monomorphicVariants()[id]
                    + "\t" + totals.belowMinimumMaf()[id]
                    + "\t" + totals.belowMinimumMac()[id]
                    + "\t" + totals.noCallVariants()[id]
                    + "\t" + totals.missingGenotypes()[id]
                    + "\taligned\t" + directory + "\t" + reusedAtStart);
                writer.newLine();
            }
            complete = true;
        } finally {
            if (complete) moveAtomically(temporary, output);
            else Files.deleteIfExists(temporary);
        }
        System.out.println("Genotype-outer pattern variant QC: " + output
            + " (candidate pattern MAF >= " + minimumMaf + ", MAC >= " + minimumMac
            + " counts are audited but not applied)");
    }

    void finishSuccess() throws IOException {
        if (!keepAfterSuccess) cleanCompletedFiles();
    }

    Path directory() { return directory; }

    private BlockCounts readBlock(int blockNumber) throws IOException {
        validateBlock(blockNumber);
        Path path = part(blockNumber);
        if (!Files.isRegularFile(path))
            throw new IOException("Incomplete pattern-QC checkpoint block " + blockNumber);
        BlockCounts counts = new BlockCounts(patternCount);
        CheckedInputStream checked = new CheckedInputStream(
            new BufferedInputStream(Files.newInputStream(path)), new CRC32());
        try (DataInputStream input = new DataInputStream(checked)) {
            if (input.readInt() != MAGIC || input.readInt() != VERSION
                || input.readInt() != blockNumber || input.readInt() != patternCount)
                throw new IOException("Invalid pattern-QC checkpoint block " + path);
            for (int pattern = 0; pattern < patternCount; pattern++) {
                counts.inputVariants()[pattern] = readNonnegative(input, path);
                counts.includedVariants()[pattern] = readNonnegative(input, path);
                counts.monomorphicVariants()[pattern] = readNonnegative(input, path);
                counts.belowMinimumMaf()[pattern] = readNonnegative(input, path);
                counts.belowMinimumMac()[pattern] = readNonnegative(input, path);
                counts.noCallVariants()[pattern] = readNonnegative(input, path);
                counts.missingGenotypes()[pattern] = readNonnegative(input, path);
            }
            long observedChecksum = checked.getChecksum().getValue();
            long expectedChecksum = input.readLong();
            if (observedChecksum != expectedChecksum)
                throw new IOException("Pattern-QC checkpoint checksum failure in " + path);
            if (input.read() != -1)
                throw new IOException("Pattern-QC checkpoint has trailing data: " + path);
        } catch (EOFException e) {
            throw new IOException("Truncated pattern-QC checkpoint " + path, e);
        }
        return counts;
    }

    private static long readNonnegative(DataInputStream input, Path path) throws IOException {
        long value = input.readLong();
        if (value < 0) throw new IOException("Negative aggregate in pattern-QC checkpoint " + path);
        return value;
    }

    private static void add(BlockCounts target, BlockCounts source) {
        for (int pattern = 0; pattern < target.patternCount(); pattern++) {
            target.inputVariants()[pattern] = Math.addExact(target.inputVariants()[pattern], source.inputVariants()[pattern]);
            target.includedVariants()[pattern] = Math.addExact(target.includedVariants()[pattern], source.includedVariants()[pattern]);
            target.monomorphicVariants()[pattern] = Math.addExact(target.monomorphicVariants()[pattern], source.monomorphicVariants()[pattern]);
            target.belowMinimumMaf()[pattern] = Math.addExact(target.belowMinimumMaf()[pattern], source.belowMinimumMaf()[pattern]);
            target.belowMinimumMac()[pattern] = Math.addExact(target.belowMinimumMac()[pattern], source.belowMinimumMac()[pattern]);
            target.noCallVariants()[pattern] = Math.addExact(target.noCallVariants()[pattern], source.noCallVariants()[pattern]);
            target.missingGenotypes()[pattern] = Math.addExact(target.missingGenotypes()[pattern], source.missingGenotypes()[pattern]);
        }
    }

    private void cleanCompletedFiles() throws IOException {
        for (int block = 0; block < totalBlocks; block++) {
            Files.deleteIfExists(part(block));
            Files.deleteIfExists(partial(block));
        }
        Files.deleteIfExists(directory.resolve(MANIFEST));
        Files.deleteIfExists(directory);
    }

    private Path part(int blockNumber) {
        return directory.resolve(String.format("block-%08d.qc", blockNumber));
    }

    private Path partial(int blockNumber) {
        return directory.resolve(String.format("block-%08d.partial", blockNumber));
    }

    private void validateBlock(int blockNumber) {
        if (blockNumber < 0 || blockNumber >= totalBlocks)
            throw new IllegalArgumentException("Invalid pattern-QC checkpoint block " + blockNumber);
    }

    private static Properties load(Path manifest) throws IOException {
        Properties properties = new Properties();
        try (var input = new BufferedInputStream(Files.newInputStream(manifest))) {
            properties.load(input);
        }
        return properties;
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE,
                StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
