/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Properties;

/** Atomic, source-ordered persistence for long variant-QC scans. */
final class QVariantQcCheckpoint implements AutoCloseable {
    private static final String FORMAT_VERSION = "1";
    private static final String MANIFEST_NAME = "manifest.properties";

    @FunctionalInterface
    interface LineConsumer {
        void accept(String line) throws IOException;
    }

    record Snapshot(long records, boolean complete) { }

    private final Path directory;
    private final Path manifestPath;
    private final String signature;
    private final String description;
    private final FileChannel lockChannel;
    private final FileLock lock;
    private boolean complete;
    private long completedRecords;
    private int nextPart;

    QVariantQcCheckpoint(Path root, String signature, String description) throws IOException {
        if (root == null || signature == null || signature.isBlank())
            throw new IllegalArgumentException("Variant-QC checkpoint root and signature are required");
        this.signature = signature;
        this.description = description == null ? "" : description;
        Path normalizedRoot = root.toAbsolutePath().normalize();
        Files.createDirectories(normalizedRoot);
        directory = normalizedRoot.resolve(signature);
        Files.createDirectories(directory);
        manifestPath = directory.resolve(MANIFEST_NAME);
        lockChannel = FileChannel.open(directory.resolve(".lock"),
            StandardOpenOption.CREATE, StandardOpenOption.WRITE);
        FileLock acquired;
        try {
            acquired = lockChannel.tryLock();
        } catch (OverlappingFileLockException e) {
            acquired = null;
        }
        if (acquired == null) {
            lockChannel.close();
            throw new IOException("Variant-QC checkpoint is already in use by another process: "
                + directory);
        }
        lock = acquired;
        try {
            initializeOrValidateManifest();
        } catch (IOException | RuntimeException e) {
            try {
                close();
            } catch (IOException closeFailure) {
                e.addSuppressed(closeFailure);
            }
            throw e;
        }
    }

    Path directory() {
        return directory;
    }

    Snapshot load(LineConsumer consumer) throws IOException {
        List<Path> parts = partFiles();
        long records = 0;
        for (int index = 0; index < parts.size(); index++) {
            Path part = parts.get(index);
            String expected = partName(index);
            if (!part.getFileName().toString().equals(expected))
                throw new IOException("Variant-QC checkpoint has a missing or out-of-order part: expected "
                    + expected + " in " + directory);
            try (BufferedReader reader = Files.newBufferedReader(part, StandardCharsets.UTF_8)) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.isEmpty())
                        throw new IOException("Variant-QC checkpoint contains a blank record in " + part);
                    consumer.accept(line);
                    records++;
                }
            }
        }
        nextPart = parts.size();
        completedRecords = records;
        if (complete && records != manifestRecordCount())
            throw new IOException("Completed variant-QC checkpoint expected " + manifestRecordCount()
                + " records but contains " + records + " in " + directory);
        return new Snapshot(records, complete);
    }

    void commit(List<String> lines) throws IOException {
        if (lines == null || lines.isEmpty())
            return;
        if (complete)
            throw new IOException("Cannot append to completed variant-QC checkpoint " + directory);
        Path target = directory.resolve(partName(nextPart));
        if (Files.exists(target))
            throw new IOException("Variant-QC checkpoint part already exists: " + target);
        Path temporary = Files.createTempFile(directory, "part-", ".partial");
        boolean moved = false;
        try (BufferedWriter writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {
            for (String line : lines) {
                if (line == null || line.indexOf('\n') >= 0 || line.indexOf('\r') >= 0)
                    throw new IOException("Variant-QC checkpoint records must be single non-null lines");
                writer.write(line);
                writer.write('\n');
            }
        }
        try {
            moveAtomically(temporary, target);
            moved = true;
            nextPart++;
            completedRecords += lines.size();
        } finally {
            if (!moved)
                Files.deleteIfExists(temporary);
        }
    }

    void markComplete(long records) throws IOException {
        if (records != completedRecords)
            throw new IOException("Cannot complete variant-QC checkpoint: state has " + records
                + " records but durable parts contain " + completedRecords);
        complete = true;
        writeManifest(true, records);
    }

    void materialize(Path output, String header) throws IOException {
        if (!complete)
            throw new IOException("Cannot write QC output from incomplete checkpoint " + directory);
        Path normalized = output.toAbsolutePath().normalize();
        Path parent = normalized.getParent();
        if (parent != null)
            Files.createDirectories(parent);
        Path temporary = Files.createTempFile(parent, normalized.getFileName().toString(), ".partial");
        boolean moved = false;
        try (BufferedWriter writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {
            writer.write(header);
            if (!header.endsWith("\n"))
                writer.write('\n');
            for (Path part : partFiles()) {
                try (BufferedReader reader = Files.newBufferedReader(part, StandardCharsets.UTF_8)) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        writer.write(line);
                        writer.write('\n');
                    }
                }
            }
        }
        try {
            moveAtomically(temporary, normalized);
            moved = true;
        } finally {
            if (!moved)
                Files.deleteIfExists(temporary);
        }
    }

    void forEachLine(LineConsumer consumer) throws IOException {
        for (Path part : partFiles()) {
            try (BufferedReader reader = Files.newBufferedReader(part, StandardCharsets.UTF_8)) {
                String line;
                while ((line = reader.readLine()) != null)
                    consumer.accept(line);
            }
        }
    }

    @Override
    public void close() throws IOException {
        try {
            if (lock.isValid())
                lock.release();
        } finally {
            lockChannel.close();
        }
    }

    private void initializeOrValidateManifest() throws IOException {
        if (!Files.exists(manifestPath)) {
            try (var entries = Files.list(directory)) {
                if (entries.anyMatch(path -> {
                    String name = path.getFileName().toString();
                    return !name.equals(".lock") && !name.endsWith(".partial");
                }))
                    throw new IOException("Variant-QC checkpoint directory has no manifest but is not empty: "
                        + directory);
            }
            complete = false;
            completedRecords = 0;
            writeManifest(false, 0);
            return;
        }
        Properties properties = readManifest();
        if (!FORMAT_VERSION.equals(properties.getProperty("format")))
            throw new IOException("Unsupported variant-QC checkpoint format in " + manifestPath);
        if (!signature.equals(properties.getProperty("signature")))
            throw new IOException("Variant-QC checkpoint signature mismatch in " + manifestPath);
        complete = Boolean.parseBoolean(properties.getProperty("complete", "false"));
    }

    private Properties readManifest() throws IOException {
        Properties properties = new Properties();
        try (Reader reader = Files.newBufferedReader(manifestPath, StandardCharsets.UTF_8)) {
            properties.load(reader);
        }
        return properties;
    }

    private long manifestRecordCount() throws IOException {
        String value = readManifest().getProperty("records", "0");
        try {
            return Long.parseLong(value);
        } catch (NumberFormatException e) {
            throw new IOException("Invalid record count in " + manifestPath, e);
        }
    }

    private void writeManifest(boolean isComplete, long records) throws IOException {
        Properties properties = new Properties();
        properties.setProperty("format", FORMAT_VERSION);
        properties.setProperty("signature", signature);
        properties.setProperty("complete", Boolean.toString(isComplete));
        properties.setProperty("records", Long.toString(records));
        properties.setProperty("description", description);
        Path temporary = Files.createTempFile(directory, "manifest-", ".partial");
        boolean moved = false;
        try (Writer writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {
            properties.store(writer, "GPU eQTL resumable variant-QC checkpoint");
        }
        try {
            moveAtomically(temporary, manifestPath);
            moved = true;
        } finally {
            if (!moved)
                Files.deleteIfExists(temporary);
        }
    }

    private List<Path> partFiles() throws IOException {
        List<Path> parts = new ArrayList<>();
        try (var entries = Files.list(directory)) {
            entries.filter(path -> path.getFileName().toString().matches("part-[0-9]{8}\\.tsv"))
                .sorted(Comparator.comparing(path -> path.getFileName().toString()))
                .forEach(parts::add);
        }
        return parts;
    }

    private static String partName(int index) {
        return String.format("part-%08d.tsv", index);
    }

    private static void moveAtomically(Path source, Path target) throws IOException {
        try {
            Files.move(source, target, StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
        } catch (AtomicMoveNotSupportedException e) {
            Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
        }
    }
}
