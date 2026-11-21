/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.io.*;
import java.nio.file.*;
import java.nio.channels.*;
import java.util.*;
import gov.nih.eqtl.io.*;
import gov.nih.gpu.*;

/** Opt-in pooled score scan; bounded SNP tiles, raw indexed cache and atomic trait parts. */
final class QPooledScoreRunner {
    private static final String ALGORITHM="pooled-pedigree-block-v3-normalized-null-integrity";
    private static final String HEADER="SNP,Probe,N,DF,Effect,SE,Statistic,log10P,Test,Status,VarianceComponents,REMLLogLik,Iterations,MaxProjectedScore,Rank";
    private static final String AUDIT="Trait\tN\tRank\tCohortN\tTestedVariants\tDegenerateVariants\tStatus\tVarianceComponents\tREMLLogLik\tMaxProjectedScore\tIdentityReduction";
    static void validate(QeQTLAnalysisConfig config,QPooledPedigreeModel model,QMatrixRowSource genotypes,int[] genotypeColumns,
        QMatrixRowSource traits,int[] traitColumns)throws Exception {
        CisIndex cis=null;
        if(config.getTraitPositionsFilename()!=null) {
            List<String> ids=new ArrayList<>();
            try(var reader=genotypes.open(new int[]{genotypeColumns[0]})) {QMatrixRowSource.Block tile;while((tile=reader.readBlock(4096))!=null)Collections.addAll(ids,tile.rowIds());}
            cis=new CisIndex(Path.of(config.getTraitPositionsFilename()),ids.toArray(String[]::new),Long.parseLong(config.get("mixed_cis_window")));
        }
        try(var reader=traits.open(traitColumns)) {
            QMatrixRowSource.Block row;while((row=reader.readBlock(1))!=null) {
                String id=row.rowIds()[0];requireTsv(id);if(cis!=null)cis.rows(id);
                int[] mask=observedMask(row.values()[0],id);
                try{model.geometry(mask);}catch(QPooledPedigreeModel.Unestimable error){if(!config.getSkipUnestimableTraits())throw error;}
            }
        }
        System.out.println("Pooled alignment, pedigree, observed-mask rank/DF, covariance identifiability and cis metadata validated; no fit or scan.");
    }
    private static int[] observedMask(double[] complete,String id) {
        for(double v:complete)if(Double.isInfinite(v))throw new IllegalArgumentException("Infinite trait value: "+id);
        return java.util.stream.IntStream.range(0,complete.length).filter(i->!QMissingnessScan.isMissing(complete[i])).toArray();
    }
    static void run(QeQTLAnalysisConfig config,QPooledPedigreeModel model,QMatrixRowSource genotypes,int[] genotypeColumns,
        QMatrixRowSource traits,int[] traitColumns,GpuContext[] contexts)throws Exception {
        String stamp=QMixedModelRunner.sourceStamp(config);
        String signature=QMixedModelRunner.signature(config,genotypeColumns,traitColumns,ALGORITHM+";"+model.engine()+";"+QMixedModelRunner.execution(contexts));
        Path out=Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
        Path checkpointDir=config.getCheckpointDirectory()==null?Path.of(out+".pooled-checkpoint"):Path.of(config.getCheckpointDirectory());
        checkpointDir=checkpointDir.toAbsolutePath().normalize();
        Path ownership=checkpointDir.resolveSibling(checkpointDir.getFileName()+".owner.lock");Files.createDirectories(ownership.getParent());
        try(FileChannel owner=FileChannel.open(ownership,StandardOpenOption.CREATE,StandardOpenOption.WRITE);FileLock owned=acquire(owner)) {
        int traitCount=Math.toIntExact(traits.metadata().rowCount());int block=config.getGenotypeBlockRows()>0?config.getGenotypeBlockRows():256;
        Path cacheDir=config.getCacheDirectory()==null?Path.of(out+".pooled-cache"):Path.of(config.getCacheDirectory());
        Path nullCache=config.getMixedNullCacheDirectory()==null?cacheDir.resolve("pooled-nulls"):Path.of(config.getMixedNullCacheDirectory());
        try(QPooledResultCheckpoint results=QPooledResultCheckpoint.open(checkpointDir,signature,traitCount,config.getResume());
            QPooledResultCheckpoint audits=QPooledResultCheckpoint.open(checkpointDir.resolve("audit"),"pooled-audit-v2;"+signature,traitCount,config.getResume())) {
        QMixedModelRunner.requireUnchanged(config,stamp);
        try(QRawMatrixCache raw=QRawMatrixCache.openOrBuild(cacheDir.resolve("aligned-raw"),QRawMatrixCache.signature(genotypes,genotypeColumns),genotypes,genotypeColumns,block,config.getRebuildCache());
            GpuContextPool pool=new GpuContextPool(contexts);var traitReader=traits.open(traitColumns)) {
            compile(pool);
            CisIndex cis=config.getTraitPositionsFilename()==null?null:new CisIndex(Path.of(config.getTraitPositionsFilename()),raw.indexedRowIds(),Long.parseLong(config.get("mixed_cis_window")));
            for(int t=0;t<traitCount;t++) {
                var row=traitReader.readBlock(1);if(row==null||row.rowOffset()!=t||row.rowCount()!=1)throw new IOException("Pooled trait source changed");
                if(results.isComplete(t)&&audits.isComplete(t))continue;
                QMixedModelRunner.requireUnchanged(config,stamp);
                String id=row.rowIds()[0];requireTsv(id);double[] complete=row.values()[0];
                int[] mask=observedMask(complete,id);
                int[] cisRows=cis==null?null:cis.rows(id);
                double[] y=Arrays.stream(mask).mapToDouble(i->complete[i]).toArray();
                QPooledPedigreeModel.Geometry geometry;
                QPooledPedigreeModel.Fit fit;
                try{geometry=model.geometry(mask);fit=QPooledNullCache.fit(geometry,y,id,nullCache,config.getRebuildCache());}catch(QPooledPedigreeModel.Unestimable error) {
                    if(!config.getSkipUnestimableTraits())throw error;
                    if(!results.isComplete(t))results.writeBlock(t,writer->{});
                    if(!audits.isComplete(t))audits.writeBlock(t,writer->writer.write(id+"\t"+mask.length+"\tNA\t"+model.cohortCounts(mask)+"\t0\t0\tskipped:"+error.getMessage()+"\tNA\tNA\tNA\tNA\n"));
                    continue;
                }
                long[] counts=new long[2];
                if(!results.isComplete(t))results.writeBlock(t,writer->{
                    if(cisRows!=null) {
                        for(int start=0;start<cisRows.length;start+=block) {
                            int[] selected=Arrays.copyOfRange(cisRows,start,Math.min(cisRows.length,start+block));
                            scan(raw.readRows(selected,mask),geometry,fit,pool,config,id,writer,counts);
                        }
                    } else try(var reader=raw.open(mask)) {QMatrixRowSource.Block tile;while((tile=reader.readBlock(block))!=null)scan(tile,geometry,fit,pool,config,id,writer,counts);}
                    QMixedModelRunner.requireUnchanged(config,stamp);
                });
                if(!audits.isComplete(t)) {
                    // An interrupted paired audit can be reconstructed, but counts need a fresh scan.
                    if(counts[0]==0&&results.isComplete(t)) {
                        if(cisRows!=null)for(int start=0;start<cisRows.length;start+=block)scan(raw.readRows(Arrays.copyOfRange(cisRows,start,Math.min(cisRows.length,start+block)),mask),geometry,fit,pool,config,id,Writer.nullWriter(),counts);
                        else try(var reader=raw.open(mask)){QMatrixRowSource.Block tile;while((tile=reader.readBlock(block))!=null)scan(tile,geometry,fit,pool,config,id,Writer.nullWriter(),counts);}
                    }
                    audits.writeBlock(t,writer->writer.write(id+"\t"+geometry.n+"\t"+geometry.p+"\t"+model.cohortCounts(mask)+"\t"+counts[0]+"\t"+counts[1]+"\tok\t"+fit.variances()+"\t"+fit.logLik+"\t"+fit.maximumScore+"\t"+geometry.combinedIdentity+"\n"));
                }
                System.out.println("Pooled score trait "+id+": N="+geometry.n+", rank="+geometry.p+", projected score="+fit.maximumScore+", geometry builds="+model.geometryBuilds()+", null="+(fit.reused?"cache":"computed"));
            }
            if(traitReader.readBlock(1)!=null)throw new IOException("Pooled trait source grew");
        }
        if(!signature.equals(QMixedModelRunner.signature(config,genotypeColumns,traitColumns,ALGORITHM+";"+model.engine()+";"+QMixedModelRunner.execution(contexts))))throw new IOException("Pooled inputs changed");
        // Assemble audit first because the parent checkpoint may be cleaned after result assembly.
        audits.assemble(Path.of(out+".pooled.tsv"),AUDIT);results.assemble(out,HEADER);
        if(!config.getKeepCheckpoints()){audits.removeCompleted();results.removeCompleted();}
        }
        }
    }
    private static FileLock acquire(FileChannel owner)throws IOException {
        try{FileLock lock=owner.tryLock();if(lock==null)throw new IOException("Another process owns this pooled checkpoint");return lock;}
        catch(OverlappingFileLockException error){throw new IOException("Another pooled scan in this JVM owns this checkpoint",error);}
    }
    private static void scan(QMatrixRowSource.Block tile,QPooledPedigreeModel.Geometry geometry,QPooledPedigreeModel.Fit fit,GpuContextPool pool,
        QeQTLAnalysisConfig config,String trait,Writer writer,long[] counts)throws Exception {
        double[][] rotated=new double[tile.rowCount()][];
        for(int a=0;a<rotated.length;a++) {
            double[] g=tile.values()[a].clone();double sum=0;int called=0;
            for(double v:g){if(Double.isInfinite(v))throw new IllegalArgumentException("Infinite genotype");if(!QMissingnessScan.isMissing(v)){sum+=v;called++;}}
            if(called<g.length&&config.getPredictorMissingPolicy()==QMissingValuePolicy.ERROR)throw new IllegalArgumentException("Missing genotype in pooled observed mask");
            double fill=config.getPredictorMissingPolicy()==QMissingValuePolicy.ZERO||called==0?0:sum/called;
            for(int i=0;i<g.length;i++)if(QMissingnessScan.isMissing(g[i]))g[i]=fill;
            // Centering by a constant improves collinearity cancellation; intercept is in the design.
            double mean=0;for(double v:g)mean+=v;mean/=g.length;for(int i=0;i<g.length;i++)g[i]-=mean;
            rotated[a]=g;
        }
        rotated=geometry.rotateRows(rotated);
        double[][] products=products(rotated,fit,pool);
        for(int a=0;a<rotated.length;a++) {
            counts[0]++;var result=fit.score(rotated[a],products[a]);if(!result.status().equals("ok"))counts[1]++;
            if(config.getThresholdType().equals("pval")&&result.log10P()>Math.log10(config.getThresholdValue()))continue;
            writer.write(csv(tile.rowIds()[a])+","+csv(trait)+","+geometry.n+",1,"+result.beta()+","+result.se()+","+result.statistic()+","+result.log10P()+",pooled-pedigree-block-score,"+result.status()+","+csv(fit.variances())+","+fit.logLik+","+fit.iterations+","+fit.maximumScore+","+geometry.p+"\n");
        }
    }
    static double[][] products(double[][] g,QPooledPedigreeModel.Fit fit,GpuContextPool pool)throws Exception {
        if(pool==null||pool.getNumContexts()==0){double[][] result=new double[g.length][];for(int a=0;a<g.length;a++)result[a]=fit.products(g[a]);return result;}
        int n=fit.geometry.n,p=fit.geometry.p+1,rows=round(g.length,16),cols=round(p,16),shared=round(n,64);
        double[] left=new double[Math.multiplyExact(rows,shared)],right=new double[Math.multiplyExact(shared,cols)];
        for(int a=0;a<g.length;a++)System.arraycopy(g[a],0,left,a*shared,n);
        double[][] scores=fit.scoreColumns();for(int i=0;i<n;i++)System.arraycopy(scores[i],0,right,i*cols,p);
        GpuContext context=pool.reserveContext();double[] packed;
        try{packed=context.executeDoubleKernel(left,right,Math.multiplyExact(rows,cols),17L*64*8,shared,cols,new long[]{cols,rows},new long[]{16,16});}finally{pool.releaseContext(context);}
        double[][] result=new double[g.length][p];for(int a=0;a<g.length;a++)for(int j=0;j<p;j++){double v=packed[a*cols+j];if(!Double.isFinite(v))throw new ArithmeticException("Nonfinite pooled GPU score product");result[a][j]=v;}return result;
    }
    static void compile(GpuContextPool pool)throws Exception {
        String kernel="#define BLOCK_SIZE 16\n#define DATATYPE double\n#define N_MIN_1 1\n#if defined(cl_khr_fp64)\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n#elif defined(cl_amd_fp64)\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n"+QeQTLAnalysis.eqtlReal;
        for(GpuContext context:pool.getAllContexts())context.compileKernel(kernel,"eqtlReal",GpuPrecision.FP64);
    }
    private static int round(int n,int b){return Math.multiplyExact(Math.addExact(n,b-1)/b,b);}
    private static String csv(String s){return "\""+s.replace("\"","\"\"")+"\"";}
    private static void requireTsv(String id){if(id.contains("\t")||id.contains("\n")||id.contains("\r"))throw new IllegalArgumentException("Trait ID cannot be represented in pooled audit TSV");}
    static final class CisIndex {
        private record Position(String chromosome,long position) { }
        private record Variant(long position,int row) { }
        private final Map<String,Position> traits=new HashMap<>();
        private final Map<String,Variant[]> variants=new HashMap<>();private final long radius;
        CisIndex(Path positions,String[] ids,long radius)throws IOException {
            if(radius<0)throw new IllegalArgumentException("Cis radius must be nonnegative");
            this.radius=radius;List<String> lines=Files.readAllLines(positions);
            if(lines.isEmpty()||!lines.get(0).equals("ID\tCHROM\tPOS"))throw new IOException("Trait positions require ID/CHROM/POS TSV header");
            for(String line:lines.subList(1,lines.size())){String[] v=line.split("\t",-1);if(v.length!=3||v[0].isBlank()||traits.putIfAbsent(v[0],new Position(chrom(v[1]),positive(v[2])))!=null)throw new IOException("Invalid or duplicate trait position");}
            Map<String,List<Variant>> byChromosome=new HashMap<>();
            for(int i=0;i<ids.length;i++){String[] v=ids[i].split(":",3);if(v.length<2)throw new IOException("Cis score mode requires CHROM:POS-prefixed variant IDs");byChromosome.computeIfAbsent(chrom(v[0]),unused->new ArrayList<>()).add(new Variant(positive(v[1]),i));}
            for(var entry:byChromosome.entrySet()){entry.getValue().sort(Comparator.comparingLong(Variant::position).thenComparingInt(Variant::row));variants.put(entry.getKey(),entry.getValue().toArray(Variant[]::new));}
        }
        int[] rows(String id)throws IOException {
            Position p=traits.get(id);if(p==null)throw new IOException("No exact versioned trait position for "+id);
            Variant[] rows=variants.getOrDefault(p.chromosome,new Variant[0]);long lo=Math.max(1,p.position-Math.min(radius,p.position-1)),hi=radius>Long.MAX_VALUE-p.position?Long.MAX_VALUE:p.position+radius;
            int a=lower(rows,lo),b=hi==Long.MAX_VALUE?rows.length:lower(rows,hi+1);
            return Arrays.stream(rows,a,b).mapToInt(Variant::row).sorted().toArray();
        }
        private static int lower(Variant[] a,long value){int lo=0,hi=a.length;while(lo<hi){int m=(lo+hi)>>>1;if(a[m].position<value)lo=m+1;else hi=m;}return lo;}
        private static long positive(String v)throws IOException{try{long value=Long.parseLong(v);if(value>0)return value;}catch(NumberFormatException ignored){}throw new IOException("Position must be a positive integer");}
        private static String chrom(String v)throws IOException{String s=v.startsWith("chr")?v.substring(3):v;if(s.isBlank())throw new IOException("Blank chromosome");return s;}
    }
}
