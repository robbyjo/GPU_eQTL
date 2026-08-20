/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.Properties;

/** Atomic genotype-block result parts used for restart and deterministic assembly. */
public final class QAnalysisCheckpoint {
    @FunctionalInterface
    public interface BlockWriter {
        void write(Writer writer) throws Exception;
    }

    private static final String FORMAT = "gpu-eqtl-checkpoint-v1";
    private static final String MANIFEST = "manifest.properties";

    private final Path directory;
    private final String signature;
    private final int totalBlocks;
    private final boolean keepAfterSuccess;

    private QAnalysisCheckpoint(Path directory, String signature, int totalBlocks,
        boolean keepAfterSuccess) {
        this.directory = directory;
        this.signature = signature;
        this.totalBlocks = totalBlocks;
        this.keepAfterSuccess = keepAfterSuccess;
    }

    public static QAnalysisCheckpoint open(Path directory, String signature, int totalBlocks,
        boolean resume, boolean keepAfterSuccess) throws IOException {
        Path normalized = directory.toAbsolutePath().normalize();
        Path manifest = normalized.resolve(MANIFEST);
        if (Files.exists(normalized)) {
            if (!resume)
                throw new IOException("Checkpoint directory already exists; use --resume or choose another directory: " + normalized);
            Properties properties = load(manifest);
            if (!FORMAT.equals(properties.getProperty("format"))
                || !signature.equals(properties.getProperty("signature"))
                || totalBlocks != Integer.parseInt(properties.getProperty("total_blocks", "-1")))
                throw new IOException("Checkpoint does not match this analysis: " + normalized);
            System.out.println("Resuming checkpoint: " + normalized);
        } else {
            Files.createDirectories(normalized);
            Properties properties = new Properties();
            properties.setProperty("format", FORMAT);
            properties.setProperty("signature", signature);
            properties.setProperty("total_blocks", Integer.toString(totalBlocks));
            Path temporary = Files.createTempFile(normalized, "manifest", ".partial");
            try (OutputStream output = new BufferedOutputStream(Files.newOutputStream(temporary))) {
                properties.store(output, "GPU eQTL checkpoint");
            }
            moveAtomically(temporary, manifest);
            System.out.println("Checkpoint directory: " + normalized);
        }
        return new QAnalysisCheckpoint(normalized, signature, totalBlocks, keepAfterSuccess);
    }

    public boolean isComplete(int blockNumber) {
        validateBlock(blockNumber);
        return Files.isRegularFile(part(blockNumber));
    }

    public int completedCount() {
        int completed = 0;
        for (int block = 0; block < totalBlocks; block++)
            if (isComplete(block))
                completed++;
        return completed;
    }

    public void writeBlock(int blockNumber, BlockWriter action) throws Exception {
        validateBlock(blockNumber);
        Path partial = partial(blockNumber);
        Files.deleteIfExists(partial);
        boolean complete = false;
        try (Writer writer = new BufferedWriter(Files.newBufferedWriter(partial,
            StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE), 1024 * 1024)) {
            action.write(writer);
            writer.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(partial, part(blockNumber));
        }
    }

    public void assemble(Path outputPath, String header) throws IOException {
        for (int block = 0; block < totalBlocks; block++)
            if (!isComplete(block))
                throw new IOException("Cannot assemble output; genotype block " + block + " is incomplete");
        Path output = outputPath.toAbsolutePath().normalize();
        Path parent = output.getParent();
        if (parent != null)
            Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, output.getFileName().toString(), ".partial");
        boolean complete = false;
        try (OutputStream destination = new BufferedOutputStream(Files.newOutputStream(temporary), 16 * 1024 * 1024)) {
            destination.write(header.getBytes(StandardCharsets.UTF_8));
            destination.write(System.lineSeparator().getBytes(StandardCharsets.UTF_8));
            for (int block = 0; block < totalBlocks; block++) {
                try (BufferedInputStream source = new BufferedInputStream(Files.newInputStream(part(block)), 16 * 1024 * 1024)) {
                    source.transferTo(destination);
                }
            }
            complete = true;
        } finally {
            if (complete)
                moveAtomically(temporary, output);
            else
                Files.deleteIfExists(temporary);
        }
        if (!keepAfterSuccess)
            cleanCompletedFiles();
    }

    public Path directory() {
        return directory;
    }

    public String signature() {
        return signature;
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
        return directory.resolve(String.format("block-%08d.part", blockNumber));
    }

    private Path partial(int blockNumber) {
        return directory.resolve(String.format("block-%08d.partial", blockNumber));
    }

    private void validateBlock(int blockNumber) {
        if (blockNumber < 0 || blockNumber >= totalBlocks)
            throw new IllegalArgumentException("Invalid checkpoint block " + blockNumber);
    }

    private static Properties load(Path manifest) throws IOException {
        Properties properties = new Properties();
        try (BufferedInputStream input = new BufferedInputStream(Files.newInputStream(manifest))) {
            properties.load(input);
        }
        return properties;
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
