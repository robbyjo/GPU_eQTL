/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.file.*;
import java.security.*;
import java.util.Arrays;

/** Immutable, checksummed per-trait nulls shared safely across chromosome processes. */
final class QPooledNullCache {
    private static final int MAGIC=0x51504e32;
    private QPooledNullCache(){ }
    static synchronized QPooledPedigreeModel.Fit fit(QPooledPedigreeModel.Geometry geometry,double[] response,
        String trait,Path directory,boolean rebuild)throws Exception {
        String key=geometry.nullCacheKey(trait,response);Path parent=directory.resolve(key.substring(0,2));Files.createDirectories(parent);
        Path file=parent.resolve(key+".qnull"),lockPath=parent.resolve(key+".lock");
        // Keep lock inodes stable: unlinking a lock after release races a third process.
        try(FileChannel channel=FileChannel.open(lockPath,StandardOpenOption.CREATE,StandardOpenOption.WRITE);var lock=channel.lock()) {
            if(Files.exists(file)) {
                var previous=read(file,key,geometry);
                if(!rebuild)return previous;
                // Explicit rebuild recomputes but never overwrites the immutable verified fit.
                return geometry.fit(response);
            }
            var fit=geometry.fit(response);byte[] encoded=encode(key,fit);
            Path temporary=Files.createTempFile(parent,key,".partial");
            try {
                Files.write(temporary,encoded);
                try{Files.move(temporary,file,StandardCopyOption.ATOMIC_MOVE);}catch(AtomicMoveNotSupportedException error){Files.move(temporary,file);}
            }finally{Files.deleteIfExists(temporary);}
            return fit;
        }
    }
    private static byte[] encode(String key,QPooledPedigreeModel.Fit fit)throws Exception {
        ByteArrayOutputStream bytes=new ByteArrayOutputStream();
        try(DataOutputStream out=new DataOutputStream(bytes)) {
            out.writeInt(MAGIC);out.writeUTF(key);out.writeInt(fit.geometry.n);out.writeInt(fit.geometry.p);out.writeInt(fit.variance.length);
            for(double value:fit.variance)out.writeDouble(value);for(double value:fit.weights)out.writeDouble(value);for(double value:fit.py)out.writeDouble(value);
            for(double[] row:fit.inverse)for(double value:row)out.writeDouble(value);
            out.writeDouble(fit.logLik);out.writeDouble(fit.maximumScore);out.writeInt(fit.iterations);
        }
        byte[] payload=bytes.toByteArray(),checksum=MessageDigest.getInstance("SHA-256").digest(payload);
        byte[] result=Arrays.copyOf(payload,payload.length+checksum.length);System.arraycopy(checksum,0,result,payload.length,checksum.length);return result;
    }
    private static QPooledPedigreeModel.Fit read(Path file,String key,QPooledPedigreeModel.Geometry geometry)throws Exception {
        long maximum=256+8L*(geometry.names.size()+2L*geometry.n+(long)geometry.p*geometry.p);
        long size=Files.size(file);if(size<128||size>maximum)throw new IOException("Invalid pooled null-cache length: "+file);
        byte[] encoded=Files.readAllBytes(file),payload=Arrays.copyOf(encoded,encoded.length-32),checksum=Arrays.copyOfRange(encoded,encoded.length-32,encoded.length);
        if(!MessageDigest.isEqual(checksum,MessageDigest.getInstance("SHA-256").digest(payload)))throw new IOException("Pooled null-cache checksum mismatch: "+file);
        try(DataInputStream in=new DataInputStream(new ByteArrayInputStream(payload))) {
            if(in.readInt()!=MAGIC||!in.readUTF().equals(key)||in.readInt()!=geometry.n||in.readInt()!=geometry.p||in.readInt()!=geometry.names.size())throw new IOException("Pooled null-cache model mismatch");
            double[] variance=values(in,geometry.names.size()),weights=values(in,geometry.n),py=values(in,geometry.n);double[][] inverse=new double[geometry.p][];
            for(int i=0;i<inverse.length;i++)inverse[i]=values(in,geometry.p);
            double logLik=in.readDouble(),score=in.readDouble();int iterations=in.readInt();
            if(in.available()!=0||!Double.isFinite(logLik)||!Double.isFinite(score)||score<0||score>QPooledPedigreeModel.SCORE_TOLERANCE||iterations<1)
                throw new IOException("Invalid pooled null-cache convergence record");
            for(double value:variance)if(value<0)throw new IOException("Negative cached variance");
            for(double value:weights)if(value<=0)throw new IOException("Nonpositive cached covariance weight");
            return new QPooledPedigreeModel.Fit(geometry,variance,weights,py,inverse,logLik,score,iterations,true);
        }
    }
    private static double[] values(DataInputStream input,int n)throws IOException {
        double[] values=new double[n];for(int i=0;i<n;i++){values[i]=input.readDouble();if(!Double.isFinite(values[i]))throw new IOException("Nonfinite pooled cached state");}return values;
    }
}
