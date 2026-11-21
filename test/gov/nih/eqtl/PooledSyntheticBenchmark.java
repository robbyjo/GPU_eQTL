/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
/* Synthetic engineering check only; no controlled cohort observations. */
package gov.nih.eqtl;
import java.nio.file.*;
import java.util.*;
import gov.nih.gpu.*;
import gov.nih.gpu.cuda.CudaGpuBackend;
public final class PooledSyntheticBenchmark {
 static double seconds(long start){return (System.nanoTime()-start)/1e9;}
 public static void main(String[] args)throws Exception {
  Path root=Path.of(args[0]);if(Files.exists(root))throw new IllegalArgumentException("Use a fresh benchmark output directory");Files.createDirectories(root);int n=8955,p=55,nf=4746;Random random=new Random(20261120);
  double[][] x=new double[n][p];double[] y=new double[n],genetic=new double[n];String[] ids=new String[n],groups=new String[n];StringBuilder ped=new StringBuilder("ID\tSIRE\tDAM\n");
  for(int first=0;first<nf;){int size=first==0?512:Math.min(4,nf-first);genetic[first]=random.nextGaussian();if(size>1)genetic[first+1]=random.nextGaussian();
   for(int j=0;j<size;j++){int i=first+j;ped.append("synthetic").append(i).append('\t').append(j<2?"0":"synthetic"+first).append('\t').append(j<2?"0":"synthetic"+(first+1)).append('\n');if(j>=2)genetic[i]=.5*(genetic[first]+genetic[first+1])+Math.sqrt(.5)*random.nextGaussian();}first+=size;
  }
  Path pedigree=root.resolve("synthetic-pedigree.tsv");Files.writeString(pedigree,ped);
  for(int i=0;i<n;i++){ids[i]="synthetic"+i;groups[i]=i<nf?"FHS":i<nf+2005?"WHI":"CARDIA";x[i][0]=1;x[i][1]=groups[i].equals("WHI")?1:0;x[i][2]=groups[i].equals("CARDIA")?1:0;for(int j=3;j<p;j++)x[i][j]=random.nextGaussian();y[i]=3+.2*x[i][3]+(i<nf?Math.sqrt(.6)*genetic[i]:0)+random.nextGaussian()*(i<nf?Math.sqrt(.7):i<nf+2005?1.2:.9);}
  long start=System.nanoTime();var model=new QPooledPedigreeModel(x,groups,"FHS",ids,QMixedModel.readPedigree(pedigree),1000);double graph=seconds(start);
  int[] mask=java.util.stream.IntStream.range(0,n).toArray();start=System.nanoTime();var geometry=model.geometry(mask);double geometryTime=seconds(start);
  start=System.nanoTime();var fit=QPooledNullCache.fit(geometry,y,"synthetic-full",root.resolve("nulls"),false);double nullTime=seconds(start);
  start=System.nanoTime();var cached=QPooledNullCache.fit(geometry,y,"synthetic-full",root.resolve("nulls"),false);double reuseTime=seconds(start);if(!cached.reused||model.fitExecutions()!=1)throw new AssertionError("Null fit not reused");
  double[][] g=new double[17][n];for(double[] row:g)for(int i=0;i<n;i++)row[i]=(random.nextBoolean()?1:0)+(random.nextBoolean()?1:0);g=geometry.rotateRows(g);
  double productTime,worst=0;var devices=new GpuRuntime(new CudaGpuBackend()).getGpuDevices(true,true);if(devices.isEmpty())throw new AssertionError("CUDA device absent");
  String device=devices.get(0).getName();
  try(var pool=new GpuContextPool(new GpuContext[]{devices.get(0).openContext()})){QPooledScoreRunner.compile(pool);start=System.nanoTime();double[][] products=QPooledScoreRunner.products(g,fit,pool);productTime=seconds(start);
   for(int a=0;a<g.length;a++){var reference=fit.score(g[a],fit.products(g[a]));var actual=fit.score(g[a],products[a]);double error=Math.abs(reference.log10P()-actual.log10P());worst=Math.max(worst,error);if(error>1e-8)throw new AssertionError("GPU/scalar difference "+error);}}
  int[] partial=java.util.stream.IntStream.range(0,n).filter(i->i%17!=0).toArray();start=System.nanoTime();var partialGeometry=model.geometry(partial);double partialGeometryTime=seconds(start);double[] yy=Arrays.stream(partial).mapToDouble(i->y[i]).toArray();start=System.nanoTime();var missing=partialGeometry.fit(yy);double partialFitTime=seconds(start);
  System.out.printf(Locale.ROOT,"{\"synthetic_only\":true,\"samples\":%d,\"design_rank\":%d,\"largest_pedigree_block\":%d,\"pedigree_blocks\":%d,\"graph_seconds\":%.9f,\"geometry_seconds\":%.9f,\"null_fit_seconds\":%.9f,\"null_cache_hit_seconds\":%.9f,\"GPU_17_row_products_seconds\":%.9f,\"maximum_scalar_GPU_log10P_error\":%.12g,\"partial_N\":%d,\"partial_geometry_seconds\":%.9f,\"partial_null_fit_seconds\":%.9f,\"full_projected_score\":%.12g,\"partial_projected_score\":%.12g,\"engine\":\"%s\",\"device\":\"%s\",\"passed\":true}%n",n,geometry.p,model.largestPedigreeBlock(),model.pedigreeBlockCount(),graph,geometryTime,nullTime,reuseTime,productTime,worst,partial.length,partialGeometryTime,partialFitTime,fit.maximumScore,missing.maximumScore,model.engine(),device);
 }
}
