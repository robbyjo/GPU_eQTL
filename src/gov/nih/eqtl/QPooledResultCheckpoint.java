/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.security.*;
import java.util.Arrays;

/** New-mode trait parts have model/index-bound SHA-256 records in one bounded index file. */
final class QPooledResultCheckpoint implements AutoCloseable {
    private final QAnalysisCheckpoint parts;
    private final FileChannel digests;
    private final Path index,retainedDirectory;
    private final int count,headerLength;
    private final String signature;
    private QPooledResultCheckpoint(QAnalysisCheckpoint parts,FileChannel digests,Path index,int count,int headerLength,String signature,Path retainedDirectory){this.parts=parts;this.digests=digests;this.index=index;this.count=count;this.headerLength=headerLength;this.signature=signature;this.retainedDirectory=retainedDirectory;}
    static QPooledResultCheckpoint open(Path directory,String signature,int count,boolean resume)throws IOException {
        Path normalized=directory.toAbsolutePath().normalize();
        return open(normalized,signature,count,resume,normalized.resolveSibling(normalized.getFileName()+".uncommitted"));
    }
    static QPooledResultCheckpoint open(Path directory,String signature,int count,boolean resume,Path retainedDirectory)throws IOException {
        boolean existingDirectory=Files.exists(directory);
        retainedDirectory=retainedDirectory.toAbsolutePath().normalize();
        if(retainedDirectory.startsWith(directory.toAbsolutePath().normalize()))throw new IllegalArgumentException("Retained parts must be outside the cleaned checkpoint");
        QAnalysisCheckpoint parts=QAnalysisCheckpoint.open(directory,signature,count,resume,true);
        ByteArrayOutputStream bytes=new ByteArrayOutputStream();try(DataOutputStream out=new DataOutputStream(bytes)){out.writeInt(0x51504331);out.writeUTF(signature);out.writeInt(count);}
        byte[] header=bytes.toByteArray();Path index=parts.directory().resolve("trait-digests.bin");boolean existing=Files.exists(index);
        // A missing index cannot distinguish uncommitted parts from committed
        // parts whose integrity records were lost. Never recreate that knowledge.
        if(existingDirectory&&!existing)throw new IOException("Pooled checkpoint digest index is missing: "+index);
        FileChannel channel=FileChannel.open(index,StandardOpenOption.CREATE,StandardOpenOption.READ,StandardOpenOption.WRITE);
        try {
            long size=Math.addExact(header.length,Math.multiplyExact((long)count,32));
            if(existing) {
                if(channel.size()!=size)throw new IOException("Pooled checkpoint digest index is truncated or has extra data");
                ByteBuffer observed=ByteBuffer.allocate(header.length);readFully(channel,observed,0);
                if(!Arrays.equals(header,observed.array()))throw new IOException("Pooled checkpoint digest model mismatch");
            } else {
                ByteBuffer initial=ByteBuffer.wrap(header);long position=0;
                while(initial.hasRemaining())position+=channel.write(initial,position);
                if(size>header.length)channel.write(ByteBuffer.wrap(new byte[]{0}),size-1);
                channel.force(true);
            }
            return new QPooledResultCheckpoint(parts,channel,index,count,header.length,signature,retainedDirectory);
        }catch(Throwable error){channel.close();throw error;}
    }
    boolean isComplete(int block)throws IOException {
        requireBlock(block);byte[] expected=record(block);boolean empty=true;for(byte b:expected)empty&=b==0;
        if(empty)return false;
        if(!parts.isComplete(block))throw new IOException("Committed pooled checkpoint part is missing: "+block);
        if(!MessageDigest.isEqual(expected,digest(block,part(block))))throw new IOException("Pooled checkpoint part checksum mismatch: "+block);
        return true;
    }
    void writeBlock(int block,QAnalysisCheckpoint.BlockWriter action)throws Exception {
        if(isComplete(block))throw new IOException("Refuse to overwrite committed pooled trait "+block);
        retainUncommitted(part(block));retainUncommitted(parts.directory().resolve(String.format("block-%08d.partial",block)));
        parts.writeBlock(block,action);
        ByteBuffer value=ByteBuffer.wrap(digest(block,part(block)));long position=headerLength+32L*block;
        while(value.hasRemaining())position+=digests.write(value,position);
    }
    @FunctionalInterface interface BatchWriter {void write(Writer[] writers)throws Exception;}
    /** Open one bounded group without recursively nesting one callback per trait. */
    void writeBlocks(int[] blocks,BatchWriter action)throws Exception {
        if(blocks.length<1||blocks.length>256)throw new IllegalArgumentException("Expected 1..256 checkpoint parts");
        var unique=new java.util.HashSet<Integer>();
        for(int block:blocks){requireBlock(block);if(!unique.add(block))throw new IllegalArgumentException("Duplicate checkpoint part");}
        Writer[] writers=new Writer[blocks.length];Path[] temporary=new Path[blocks.length];
        try {
            for(int i=0;i<blocks.length;i++) {
                int block=blocks[i];
                if(isComplete(block)){writers[i]=Writer.nullWriter();continue;}
                Path path=parts.directory().resolve(String.format("block-%08d.partial",block));
                retainUncommitted(part(block));retainUncommitted(path);Files.deleteIfExists(path);
                temporary[i]=path;writers[i]=new BufferedWriter(Files.newBufferedWriter(path,StandardCharsets.UTF_8,StandardOpenOption.CREATE_NEW,StandardOpenOption.WRITE),1024*1024);
            }
            action.write(writers);
            // Close every writer successfully before committing any new part.
            for(int i=0;i<writers.length;i++){writers[i].close();writers[i]=null;}
            for(int i=0;i<blocks.length;i++)if(temporary[i]!=null) {
                int block=blocks[i];
                try{Files.move(temporary[i],part(block),StandardCopyOption.ATOMIC_MOVE,StandardCopyOption.REPLACE_EXISTING);}
                catch(AtomicMoveNotSupportedException e){Files.move(temporary[i],part(block),StandardCopyOption.REPLACE_EXISTING);}
                ByteBuffer value=ByteBuffer.wrap(digest(block,part(block)));long position=headerLength+32L*block;
                while(value.hasRemaining())position+=digests.write(value,position);
            }
        } finally {
            IOException failure=null;
            for(var writer:writers)if(writer!=null)try{writer.close();}catch(IOException e){if(failure==null)failure=e;else failure.addSuppressed(e);}
            if(failure!=null)throw failure;
        }
    }
    void assemble(Path output,String header)throws IOException {
        for(int block=0;block<count;block++)if(!isComplete(block))throw new IOException("Uncommitted pooled trait "+block);
        parts.assemble(output,header);
    }
    void removeCompleted()throws IOException {
        for(int block=0;block<count;block++)if(!isComplete(block))throw new IOException("Cannot clean incomplete pooled checkpoint");
        digests.force(true);digests.close();Files.delete(index);parts.removeCompleted();
    }
    public void close()throws IOException {if(digests.isOpen()){digests.force(false);digests.close();}}
    private void retainUncommitted(Path path)throws IOException {
        if(!Files.isRegularFile(path))return;
        Files.createDirectories(retainedDirectory);
        String hash=hexFile(path);Path target=retainedDirectory.resolve(path.getFileName()+"."+hash);
        if(Files.exists(target)){if(!hexFile(target).equals(hash))throw new IOException("Retained uncommitted part changed");}
        else Files.copy(path,target);
    }
    private byte[] record(int block)throws IOException {ByteBuffer value=ByteBuffer.allocate(32);readFully(digests,value,headerLength+32L*block);return value.array();}
    private byte[] digest(int block,Path path)throws IOException {
        try {
            MessageDigest digest=MessageDigest.getInstance("SHA-256");digest.update((signature+":"+block+":").getBytes(StandardCharsets.UTF_8));
            try(InputStream input=Files.newInputStream(path)){byte[] buffer=new byte[65536];int read;while((read=input.read(buffer))!=-1)digest.update(buffer,0,read);}
            return digest.digest();
        }catch(NoSuchAlgorithmException impossible){throw new IllegalStateException(impossible);}
    }
    private static String hexFile(Path path)throws IOException {
        try(InputStream input=Files.newInputStream(path)){MessageDigest digest=MessageDigest.getInstance("SHA-256");byte[] b=new byte[65536];int n;while((n=input.read(b))!=-1)digest.update(b,0,n);return java.util.HexFormat.of().formatHex(digest.digest());}
        catch(NoSuchAlgorithmException impossible){throw new IllegalStateException(impossible);}
    }
    private static void readFully(FileChannel channel,ByteBuffer buffer,long position)throws IOException {while(buffer.hasRemaining()){int n=channel.read(buffer,position);if(n<0)throw new EOFException("Truncated pooled digest index");position+=n;}}
    private void requireBlock(int block){if(block<0||block>=count)throw new IllegalArgumentException("Invalid pooled trait index");}
    private Path part(int block){return parts.directory().resolve(String.format("block-%08d.part",block));}
}
