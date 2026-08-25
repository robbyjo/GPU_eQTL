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
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.Properties;

/** Durable, ordered result groups for exact trait-missingness-pattern analysis. */
final class QTraitPatternCheckpoint {
    private static final String FORMAT = "gpu-eqtl-trait-pattern-checkpoint-v1";
    private static final String MANIFEST = "manifest.properties";
    private static final int BUFFER_BYTES = 16 * 1024 * 1024;

    private final Path directory;
    private final String signature;
    private final int totalPatterns;
    private final boolean keepAfterSuccess;

    private QTraitPatternCheckpoint(Path directory, String signature, int totalPatterns,
        boolean keepAfterSuccess) {
        this.directory = directory;
        this.signature = signature;
        this.totalPatterns = totalPatterns;
        this.keepAfterSuccess = keepAfterSuccess;
    }

    static QTraitPatternCheckpoint open(Path directory, String signature, int totalPatterns,
        boolean resume, boolean keepAfterSuccess) throws IOException {
        if (totalPatterns <= 0)
            throw new IllegalArgumentException("Trait-pattern checkpoint needs at least one pattern");
        Path normalized = directory.toAbsolutePath().normalize();
        Path manifest = normalized.resolve(MANIFEST);
        if (Files.exists(normalized)) {
            if (!Files.isDirectory(normalized))
                throw new IOException("Trait-pattern checkpoint is not a directory: " + normalized);
            if (!resume)
                throw new IOException("Trait-pattern checkpoint already exists; use --resume or choose another directory: "
                    + normalized);
            Properties properties = load(manifest);
            if (!FORMAT.equals(properties.getProperty("format"))
                || !signature.equals(properties.getProperty("signature"))
                || totalPatterns != Integer.parseInt(properties.getProperty("total_patterns", "-1")))
                throw new IOException("Trait-pattern checkpoint does not match this analysis: " + normalized);
            System.out.println("Resuming trait-pattern checkpoint: " + normalized);
        } else {
            Files.createDirectories(normalized);
            Properties properties = new Properties();
            properties.setProperty("format", FORMAT);
            properties.setProperty("signature", signature);
            properties.setProperty("total_patterns", Integer.toString(totalPatterns));
            Path temporary = normalized.resolve(MANIFEST + ".partial");
            Files.deleteIfExists(temporary);
            try (OutputStream output = new BufferedOutputStream(Files.newOutputStream(temporary,
                StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE))) {
                properties.store(output, "GPU eQTL trait-pattern checkpoint");
            }
            moveAtomically(temporary, manifest);
            System.out.println("Trait-pattern checkpoint directory: " + normalized);
        }
        QTraitPatternCheckpoint checkpoint = new QTraitPatternCheckpoint(normalized, signature,
            totalPatterns, keepAfterSuccess);
        checkpoint.deleteKnownPartials();
        return checkpoint;
    }

    boolean isResultComplete(int patternIndex) {
        validatePattern(patternIndex);
        return Files.isRegularFile(resultPart(patternIndex));
    }

    boolean isQcComplete(int patternIndex) {
        validatePattern(patternIndex);
        return Files.isRegularFile(qcPart(patternIndex));
    }

    int completedResultCount() {
        int completed = 0;
        for (int pattern = 0; pattern < totalPatterns; pattern++)
            if (isResultComplete(pattern)) completed++;
        return completed;
    }

    Path groupOutput(int patternIndex) {
        validatePattern(patternIndex);
        return directory.resolve(String.format("pattern-%08d.assembled.csv", patternIndex));
    }

    Path blockCheckpointDirectory(int patternIndex) {
        validatePattern(patternIndex);
        return directory.resolve(String.format("pattern-%08d.blocks", patternIndex));
    }

    void commitResultFromOutput(int patternIndex) throws IOException {
        validatePattern(patternIndex);
        Path source = groupOutput(patternIndex);
        if (!Files.isRegularFile(source))
            throw new IOException("Completed trait-pattern group output is missing: " + source);
        Path partial = resultPartial(patternIndex);
        Files.deleteIfExists(partial);
        boolean complete = false;
        try (BufferedReader reader = Files.newBufferedReader(source, StandardCharsets.UTF_8);
             BufferedWriter writer = new BufferedWriter(Files.newBufferedWriter(partial,
                 StandardCharsets.UTF_8, StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE),
                 1024 * 1024)) {
            if (reader.readLine() == null)
                throw new IOException("Trait-pattern group output has no header: " + source);
            String line;
            while ((line = reader.readLine()) != null) {
                writer.write(line);
                writer.newLine();
            }
            writer.flush();
            complete = true;
        } finally {
            if (complete)
                moveAtomically(partial, resultPart(patternIndex));
            else
                Files.deleteIfExists(partial);
        }
        Files.deleteIfExists(source);
    }

    void commitEmptyResult(int patternIndex) throws IOException {
        validatePattern(patternIndex);
        writeTextPart(resultPartial(patternIndex), resultPart(patternIndex), "");
    }

    void commitQcLine(int patternIndex, String line) throws IOException {
        validatePattern(patternIndex);
        writeTextPart(qcPartial(patternIndex), qcPart(patternIndex), line + System.lineSeparator());
    }

    void assembleResults(Path output, String header) throws IOException {
        assemble(output, header, false);
    }

    void assembleQc(Path output, String header) throws IOException {
        assemble(output, header, true);
    }

    void finishSuccess() throws IOException {
        if (keepAfterSuccess)
            return;
        for (int pattern = 0; pattern < totalPatterns; pattern++) {
            Files.deleteIfExists(resultPart(pattern));
            Files.deleteIfExists(qcPart(pattern));
            Files.deleteIfExists(resultPartial(pattern));
            Files.deleteIfExists(qcPartial(pattern));
            Files.deleteIfExists(groupOutput(pattern));
        }
        Files.deleteIfExists(directory.resolve(MANIFEST + ".partial"));
        Files.deleteIfExists(directory.resolve(MANIFEST));
        Files.delete(directory);
    }

    Path directory() {
        return directory;
    }

    String signature() {
        return signature;
    }

    private void assemble(Path outputPath, String header, boolean qc) throws IOException {
        for (int pattern = 0; pattern < totalPatterns; pattern++) {
            Path part = qc ? qcPart(pattern) : resultPart(pattern);
            if (!Files.isRegularFile(part))
                throw new IOException("Cannot assemble trait-pattern " + (qc ? "QC" : "results")
                    + "; pattern " + pattern + " is incomplete");
        }
        Path output = outputPath.toAbsolutePath().normalize();
        Path parent = output.getParent();
        if (parent != null)
            Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, output.getFileName().toString(), ".partial");
        boolean complete = false;
        try (OutputStream destination = new BufferedOutputStream(
            Files.newOutputStream(temporary), BUFFER_BYTES)) {
            destination.write(header.getBytes(StandardCharsets.UTF_8));
            destination.write(System.lineSeparator().getBytes(StandardCharsets.UTF_8));
            for (int pattern = 0; pattern < totalPatterns; pattern++) {
                Path part = qc ? qcPart(pattern) : resultPart(pattern);
                try (BufferedInputStream source = new BufferedInputStream(
                    Files.newInputStream(part), BUFFER_BYTES)) {
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
    }

    private void deleteKnownPartials() throws IOException {
        Files.deleteIfExists(directory.resolve(MANIFEST + ".partial"));
        for (int pattern = 0; pattern < totalPatterns; pattern++) {
            Files.deleteIfExists(resultPartial(pattern));
            Files.deleteIfExists(qcPartial(pattern));
        }
    }

    private static void writeTextPart(Path partial, Path target, String text) throws IOException {
        Files.deleteIfExists(partial);
        boolean complete = false;
        try {
            Files.writeString(partial, text, StandardCharsets.UTF_8,
                StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE);
            complete = true;
        } finally {
            if (complete)
                moveAtomically(partial, target);
            else
                Files.deleteIfExists(partial);
        }
    }

    private Path resultPart(int patternIndex) {
        return directory.resolve(String.format("pattern-%08d.results.part", patternIndex));
    }

    private Path resultPartial(int patternIndex) {
        return directory.resolve(String.format("pattern-%08d.results.partial", patternIndex));
    }

    private Path qcPart(int patternIndex) {
        return directory.resolve(String.format("pattern-%08d.qc.part", patternIndex));
    }

    private Path qcPartial(int patternIndex) {
        return directory.resolve(String.format("pattern-%08d.qc.partial", patternIndex));
    }

    private void validatePattern(int patternIndex) {
        if (patternIndex < 0 || patternIndex >= totalPatterns)
            throw new IllegalArgumentException("Invalid trait-pattern checkpoint index " + patternIndex);
    }

    private static Properties load(Path path) throws IOException {
        if (!Files.isRegularFile(path))
            throw new IOException("Trait-pattern checkpoint manifest is missing: " + path);
        Properties properties = new Properties();
        try (BufferedInputStream input = new BufferedInputStream(Files.newInputStream(path))) {
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
