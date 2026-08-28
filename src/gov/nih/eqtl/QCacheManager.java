/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.channels.OverlappingFileLockException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.FileTime;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Read-only cache inventory and conservative, recoverable stale-artifact pruning. */
final class QCacheManager {
    record Options(Path root, Path report, boolean prune, boolean apply, int olderThanDays) { }
    record Summary(long files, long bytes, long candidates, long candidateBytes,
        long moved, long movedBytes, long locked) { }
    private record Entry(Path path, String type, long bytes, FileTime modified,
        String sha256, String integrity, String action) { }

    private static final DateTimeFormatter TRASH_STAMP = DateTimeFormatter
        .ofPattern("uuuuMMdd-HHmmss'Z'").withZone(ZoneOffset.UTC);

    private QCacheManager() { }

    static Summary run(Options options) throws IOException {
        if (options == null || options.root() == null)
            throw new IllegalArgumentException("--cache-dir is required for cache operations");
        Path root = options.root().toAbsolutePath().normalize();
        if (!Files.isDirectory(root, LinkOption.NOFOLLOW_LINKS))
            throw new IllegalArgumentException("Cache directory does not exist: " + root);
        if (options.apply() && !options.prune())
            throw new IllegalArgumentException("--apply-cache-prune requires --prune-cache");
        if (options.prune() && options.olderThanDays() < 1)
            throw new IllegalArgumentException("--prune-cache requires --prune-cache-older-than-days of at least 1");
        Instant cutoff = options.prune() ? Instant.now().minusSeconds(
            Math.multiplyExact((long) options.olderThanDays(), 24L * 60 * 60)) : Instant.MIN;
        Path trash = root.resolve(".trash").resolve(TRASH_STAMP.format(Instant.now())).normalize();
        List<Path> files;
        try (var stream = Files.walk(root)) {
            files = stream.filter(path -> Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS))
                .filter(path -> !path.startsWith(root.resolve(".trash")))
                .sorted(Comparator.comparing(path -> root.relativize(path).toString()))
                .toList();
        }
        List<Entry> entries = new ArrayList<>();
        long bytes = 0, candidates = 0, candidateBytes = 0, moved = 0, movedBytes = 0, locked = 0;
        for (Path path : files) {
            long size = Files.size(path);
            bytes = Math.addExact(bytes, size);
            FileTime modified = Files.getLastModifiedTime(path, LinkOption.NOFOLLOW_LINKS);
            String type = classify(path);
            String digest = "";
            String integrity = "not-cache-artifact";
            if (!type.equals("other") && !type.equals("lock")) {
                try {
                    digest = sha256(path);
                    integrity = size == 0 ? "empty" : "readable-sha256";
                } catch (IOException error) {
                    integrity = "unreadable:" + oneLine(error.getMessage());
                }
            }
            String action = "keep";
            boolean stale = options.prune() && isPrunable(type)
                && !modified.toInstant().isAfter(cutoff);
            if (stale) {
                candidates++;
                candidateBytes = Math.addExact(candidateBytes, size);
                if (hasActiveLock(root, path, files)) {
                    action = "skip-active-lock";
                    locked++;
                } else if (!integrity.equals("readable-sha256")) {
                    action = "skip-integrity-" + integrity;
                } else if (!options.apply()) {
                    action = "would-move-to-trash";
                } else {
                    Path relative = root.relativize(path);
                    Path destination = trash.resolve(relative).normalize();
                    requireInside(root, path);
                    requireInside(trash, destination);
                    if (destination.getParent() != null) Files.createDirectories(destination.getParent());
                    Files.move(path, destination, StandardCopyOption.ATOMIC_MOVE);
                    action = "moved-to-trash:" + root.relativize(destination);
                    moved++;
                    movedBytes = Math.addExact(movedBytes, size);
                }
            }
            entries.add(new Entry(root.relativize(path), type, size, modified,
                digest, integrity, action));
        }
        Path report = options.report() == null ? root.resolve("cache-inventory.tsv")
            : options.report().toAbsolutePath().normalize();
        writeReport(report, entries);
        Map<String, long[]> totals = new LinkedHashMap<>();
        for (Entry entry : entries) {
            long[] total = totals.computeIfAbsent(entry.type(), ignored -> new long[2]);
            total[0]++; total[1] = Math.addExact(total[1], entry.bytes());
        }
        System.out.println("Cache root: " + root);
        for (Map.Entry<String, long[]> total : totals.entrySet())
            System.out.println("  " + total.getKey() + ": " + total.getValue()[0]
                + " file(s), " + humanBytes(total.getValue()[1]));
        System.out.println("Cache inventory: " + report);
        if (options.prune())
            System.out.println("Stale prune: candidates=" + candidates + " ("
                + humanBytes(candidateBytes) + "), active-lock skips=" + locked
                + (options.apply() ? ", moved=" + moved + " (" + humanBytes(movedBytes)
                    + ") to " + trash : ", dry-run; use --apply-cache-prune to move candidates"));
        return new Summary(files.size(), bytes, candidates, candidateBytes, moved, movedBytes, locked);
    }

    private static void writeReport(Path report, List<Entry> entries) throws IOException {
        if (report.getParent() != null) Files.createDirectories(report.getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(report, StandardCharsets.UTF_8)) {
            writer.write("type\tpath\tbytes\tmodified_utc\tsha256\tintegrity\taction\n");
            for (Entry entry : entries) {
                writer.write(entry.type()); writer.write('\t');
                writer.write(entry.path().toString()); writer.write('\t');
                writer.write(Long.toString(entry.bytes())); writer.write('\t');
                writer.write(entry.modified().toInstant().toString()); writer.write('\t');
                writer.write(entry.sha256()); writer.write('\t');
                writer.write(entry.integrity()); writer.write('\t');
                writer.write(entry.action()); writer.newLine();
            }
        }
    }

    private static String classify(Path path) {
        String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        if (name.endsWith(".lock") || name.equals("lock")) return "lock";
        if (name.endsWith(".qcache")) return "prepared-matrix";
        if (name.endsWith(".qraw") || name.endsWith(".qraw.idx")) return "aligned-raw";
        if (name.endsWith(".qmiss")) return "missingness";
        if (name.endsWith(".qpvs")) return "pattern-variant-statistics";
        if (name.endsWith(".partial") || name.endsWith(".part")) return "partial";
        return "other";
    }

    private static boolean isPrunable(String type) {
        return !(type.equals("other") || type.equals("lock"));
    }

    private static boolean hasActiveLock(Path root, Path candidate, List<Path> files) {
        Path sibling = candidate.resolveSibling(candidate.getFileName() + ".lock");
		if (Files.exists(sibling, LinkOption.NOFOLLOW_LINKS) && isActiveLock(sibling)) return true;
        for (Path lock : files) {
            if (!classify(lock).equals("lock")) continue;
            Path parent = lock.getParent();
			String lockName = lock.getFileName().toString();
			if (parent != null && (lockName.equals(".lock") || lockName.equals("lock"))) {
				if (candidate.startsWith(parent) && isActiveLock(lock)) return true;
			} else if (lockName.toLowerCase(Locale.ROOT).endsWith(".lock")) {
				Path protectedPath = lock.resolveSibling(
					lockName.substring(0, lockName.length() - ".lock".length()));
				if (candidate.equals(protectedPath) && isActiveLock(lock)) return true;
			}
        }
        return false;
    }

	private static boolean isActiveLock(Path path) {
		try (FileChannel channel = FileChannel.open(path, StandardOpenOption.WRITE)) {
			try (FileLock lock = channel.tryLock()) {
				return lock == null;
			}
		} catch (OverlappingFileLockException | IOException error) {
			return true;
		}
	}

    private static String sha256(Path path) throws IOException {
        MessageDigest digest;
        try {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException error) {
            throw new IllegalStateException(error);
        }
        byte[] buffer = new byte[1024 * 1024];
        try (InputStream input = new BufferedInputStream(Files.newInputStream(path))) {
            for (int count; (count = input.read(buffer)) >= 0; )
                if (count > 0) digest.update(buffer, 0, count);
        }
        return HexFormat.of().formatHex(digest.digest());
    }

    private static void requireInside(Path root, Path path) {
        if (!path.toAbsolutePath().normalize().startsWith(root.toAbsolutePath().normalize()))
            throw new IllegalArgumentException("Cache operation escaped its validated root: " + path);
    }

    private static String humanBytes(long bytes) {
        if (bytes < 1024) return bytes + " B";
        return String.format(Locale.ROOT, "%.2f MiB", bytes / (1024.0 * 1024));
    }

    private static String oneLine(String value) {
        return value == null ? "unknown" : value.replace('\t', ' ').replace('\r', ' ').replace('\n', ' ');
    }
}
