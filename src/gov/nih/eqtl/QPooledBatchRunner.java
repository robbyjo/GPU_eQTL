/* Copyright 2026 Roby Joehanes; GNU GPL version 3. Opt-in bounded FP64 pedigree batching. */
package gov.nih.eqtl;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import gov.nih.eqtl.io.*;
import gov.nih.gpu.*;

final class QPooledBatchRunner {
    private static final class Trait {
        final int index;final String id;final int[] mask;final QPooledPedigreeModel.Fit fit;final long[] counts=new long[2];Writer writer;
        Trait(int index,String id,int[] mask,QPooledPedigreeModel.Fit fit){this.index=index;this.id=id;this.mask=mask;this.fit=fit;}
    }
    static void run(QeQTLAnalysisConfig config,QPooledPedigreeModel model,QRawMatrixCache raw,QMatrixRowSource.BlockReader reader,
            GpuContextPool pool,QPooledResultCheckpoint results,QPooledResultCheckpoint audits,Path cacheDir,Path nullCache,String stamp,int traitCount,int block) throws Exception {
        int limit=config.getMixedTraitBatchSize();
        if(limit<1||limit>256)throw new IllegalArgumentException("mixed-trait-batch-size must be 1..256");
        if(config.getTraitPositionsFilename()!=null||config.getPredictorMissingPolicy()!=QMissingValuePolicy.ERROR)throw new IllegalArgumentException("Pedigree batching requires genome-wide complete predictors");
        List<Trait> batch=new ArrayList<>();
        for(int t=0;t<traitCount;t++) {
            var row=reader.readBlock(1);if(row==null||row.rowOffset()!=t||row.rowCount()!=1)throw new IOException("Pooled trait source changed");
            if(results.isComplete(t)&&audits.isComplete(t))continue;
            QMixedModelRunner.requireUnchanged(config,stamp);String id=row.rowIds()[0];double[] complete=row.values()[0];
            if(id.contains("\t")||id.contains("\n")||id.contains("\r"))throw new IllegalArgumentException("Invalid audit trait ID");
            for(double v:complete)if(Double.isInfinite(v))throw new IllegalArgumentException("Infinite trait value");
            int[] mask=java.util.stream.IntStream.range(0,complete.length).filter(i->!QMissingnessScan.isMissing(complete[i])).toArray();
            if(!batch.isEmpty()&&(batch.size()>=limit||!Arrays.equals(batch.get(0).mask,mask))){flush(batch,config,model,raw,pool,results,audits,cacheDir,stamp,block);batch.clear();}
            QPooledPedigreeModel.Fit fit;
            try{var geometry=model.geometry(mask);double[] y=Arrays.stream(mask).mapToDouble(i->complete[i]).toArray();fit=QPooledNullCache.fit(geometry,y,id,nullCache,config.getRebuildCache());}
            catch(QPooledPedigreeModel.Unestimable error){
                if(!config.getSkipUnestimableTraits())throw error;
                if(!results.isComplete(t))results.writeBlock(t,w->{});
                if(!audits.isComplete(t))audits.writeBlock(t,w->w.write(id+"\t"+mask.length+"\tNA\t"+model.cohortCounts(mask)+"\t0\t0\tskipped:"+error.getMessage()+"\tNA\tNA\tNA\tNA\n"));
                continue;
            }
            batch.add(new Trait(t,id,mask,fit));
        }
        if(!batch.isEmpty())flush(batch,config,model,raw,pool,results,audits,cacheDir,stamp,block);
        if(reader.readBlock(1)!=null)throw new IOException("Pooled trait source grew");
    }
    private static void flush(List<Trait> batch,QeQTLAnalysisConfig config,QPooledPedigreeModel model,QRawMatrixCache raw,
            GpuContextPool pool,QPooledResultCheckpoint results,QPooledResultCheckpoint audits,Path cacheDir,String stamp,int block) throws Exception {
        var first=batch.get(0);var geometry=first.fit.geometry;var fits=batch.stream().map(t->t.fit).toArray(QPooledPedigreeModel.Fit[]::new);
        var plan=new QPooledBatchProducts(fits,pool);
        double logThreshold=config.getThresholdType().equals("pval")?Math.log10(config.getThresholdValue()):Double.POSITIVE_INFINITY;
        double critical=config.getThresholdType().equals("pval")?jdistlib.ChiSquare.quantile(config.getThresholdValue(),1,false,false):Double.NEGATIVE_INFINITY;
        try(var rotated=QPooledRotationCache.open(cacheDir.resolve("centered-rotated-v1"),raw,first.mask,geometry,block)) {
            results.writeBlocks(batch.stream().mapToInt(t->t.index).toArray(),writers->{
                for(int t=0;t<batch.size();t++)batch.get(t).writer=writers[t];
                try(var source=rotated.open(null)) {
                    QMatrixRowSource.Block tile;
                    while((tile=source.readBlock(block))!=null) {
                        var products=plan.calculate(tile.values());
                        for(int a=0;a<tile.rowCount();a++)for(int t=0;t<batch.size();t++) {
                            var trait=batch.get(t);trait.counts[0]++;var result=plan.score(tile.values()[a],products,a,t,critical);
                            if(result==null)continue;
                            if(!result.status().equals("ok"))trait.counts[1]++;
                            if(result.log10P()>logThreshold)continue;
                            var fit=trait.fit;
                            trait.writer.write(csv(tile.rowIds()[a])+","+csv(trait.id)+","+geometry.n+",1,"+result.beta()+","+result.se()+","+result.statistic()+","+result.log10P()+",pooled-pedigree-block-score,"+result.status()+","+csv(fit.variances())+","+fit.logLik+","+fit.iterations+","+fit.maximumScore+","+geometry.p+"\n");
                        }
                    }
                }
                QMixedModelRunner.requireUnchanged(config,stamp);
            });
        }
        for(var trait:batch) {
            var fit=trait.fit;
            if(!audits.isComplete(trait.index))audits.writeBlock(trait.index,w->w.write(trait.id+"\t"+geometry.n+"\t"+geometry.p+"\t"+model.cohortCounts(trait.mask)+"\t"+trait.counts[0]+"\t"+trait.counts[1]+"\tok\t"+fit.variances()+"\t"+fit.logLik+"\t"+fit.maximumScore+"\t"+geometry.combinedIdentity+"\n"));
        }
        System.out.println("Pedigree batch completed: traits="+batch.size()+", variants="+first.counts[0]+", N="+geometry.n+", rank="+geometry.p);
    }
    private static String csv(String s){return "\""+s.replace("\"","\"\"")+"\"";}
}

