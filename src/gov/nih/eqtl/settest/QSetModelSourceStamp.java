/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.io.IOException;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.security.*;
import java.util.*;

/** Portable model-input identity, with process-local metadata checks between commits. */
final class QSetModelSourceStamp {
    private record Entry(BasicFileAttributes attributes,String digest) { }
    private final Map<Path,Entry> cache=new HashMap<>();

    String capture(Path... files) throws IOException {
        StringBuilder result=new StringBuilder();
        for(Path file:files) {
            Path path=file.toAbsolutePath().normalize();
            BasicFileAttributes before=Files.readAttributes(path,BasicFileAttributes.class);
            Entry entry=cache.get(path);
            // A file key is useful only within this process: its device number
            // can differ across cluster nodes mounting the same shared files.
            if(entry==null || before.fileKey()==null || !same(entry.attributes,before)) {
                MessageDigest digest;
                try { digest=MessageDigest.getInstance("SHA-256"); }
                catch(NoSuchAlgorithmException e) { throw new IllegalStateException(e); }
                try(var input=Files.newInputStream(path)) {
                    byte[] buffer=new byte[65536];int count;
                    while((count=input.read(buffer))!=-1)digest.update(buffer,0,count);
                }
                BasicFileAttributes after=Files.readAttributes(path,BasicFileAttributes.class);
                if(!same(before,after))throw new IOException("Model input changed while fingerprinting: "+path);
                entry=new Entry(after,HexFormat.of().formatHex(digest.digest()));cache.put(path,entry);
            }
            result.append(path).append(':').append(before.size()).append(':').append(before.lastModifiedTime())
                .append(":sha256=").append(entry.digest).append('\n');
        }
        return result.toString();
    }

    private static boolean same(BasicFileAttributes a,BasicFileAttributes b) {
        return a.size()==b.size() && a.lastModifiedTime().equals(b.lastModifiedTime())
            && a.creationTime().equals(b.creationTime()) && Objects.equals(a.fileKey(),b.fileKey());
    }
}
