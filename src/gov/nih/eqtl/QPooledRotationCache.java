/* Copyright 2026 Roby Joehanes; GNU GPL version 3.
 * Cache exact original centering and pedigree rotation across trait batches;
 * model fitting, score products and inference remain the released implementation.
 */
package gov.nih.eqtl;

import gov.nih.eqtl.io.*;
import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.*;

final class QPooledRotationCache {
    static QRawMatrixCache open(Path directory, QRawMatrixCache raw, int[] mask,
            QPooledPedigreeModel.Geometry geometry, int block) throws Exception {
        String stamp = "eqtl-centered-rotation-v1;" + QRawMatrixCache.signature(raw,mask)
            + ";" + geometry.nullCacheKey("rotation-cache-only",new double[geometry.n]);
        String key=HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256")
            .digest(stamp.getBytes(StandardCharsets.UTF_8)));
        Files.createDirectories(directory);
        String[] all=raw.metadata().sampleIds();
        String[] ids=Arrays.stream(mask).mapToObj(i->all[i]).toArray(String[]::new);
        QMatrixRowSource source=new QMatrixRowSource() {
            public Metadata metadata(){return new Metadata(raw.path(),raw.metadata().rowCount(),mask.length,ids,"centered-rotation-v1;"+key);}
            public BlockReader open(int[] columns)throws IOException {
                if(columns!=null){if(columns.length!=mask.length)throw new IOException("Rotation-cache column count changed");for(int i=0;i<columns.length;i++)if(columns[i]!=i)throw new IOException("Rotation-cache columns must be identity");}
                BlockReader reader=raw.open(mask);
                return new BlockReader() {
                    public Block readBlock(int rows)throws IOException {
                        Block tile=reader.readBlock(rows);if(tile==null)return null;
                        double[][] centered=new double[tile.rowCount()][];
                        for(int a=0;a<centered.length;a++) {
                            double[] g=tile.values()[a].clone();
                            for(double v:g)if(!Double.isFinite(v)||QMissingnessScan.isMissing(v))throw new IOException("Rotation cache requires complete finite genotypes");
                            double mean=0;for(double v:g)mean+=v;mean/=g.length;
                            for(int i=0;i<g.length;i++)g[i]-=mean;
                            centered[a]=g;
                        }
                        return new Block(tile.rowOffset(),tile.rowIds(),geometry.rotateRows(centered));
                    }
                    public void close()throws IOException{reader.close();}
                };
            }
        };
        try(FileChannel channel=FileChannel.open(directory.resolve(key+".lock"),StandardOpenOption.CREATE,StandardOpenOption.WRITE);var lock=channel.lock()) {
            // Existing corruption is fatal; it must never trigger a silent rebuild.
            QRawMatrixCache previous=QRawMatrixCache.openIfPresent(directory,key,source);
            if(previous!=null)return previous;
            return QRawMatrixCache.openOrBuild(directory,key,source,null,block,false);
        }
    }
}

