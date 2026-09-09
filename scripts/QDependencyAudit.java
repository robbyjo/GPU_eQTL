/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package org.jlinalg.settest;

import java.nio.file.*;
import java.util.*;
import gov.nih.jama.QRDecomposition;
import gov.nih.jama.CholeskyDecomposition;
import gov.nih.jama.Matrix;
import org.jlinalg.ols.*;
import org.jlinalg.compute.BackendPolicy;
import jdistlib.accelerator.CpuComputeBackend;
import org.jlinalg.reml.*;

/** Separate audit classpath only: does not add either candidate library to production. */
public final class QDependencyAudit {
    private static volatile double sink;
    public static void main(String[] args) throws Exception {
        double maxOld=0,maxNew=0,maxBetween=0,maxRelative=0;
        var tails=Files.readAllLines(Path.of("test/resources/audit-reference/t-tail.csv"));
        for(String line:tails.subList(1,tails.size())) {
            String[] v=line.split(","); double df=Double.parseDouble(v[0]),t=Double.parseDouble(v[1]),r=Double.parseDouble(v[2]);
            maxOld=Math.max(maxOld,Math.abs(net.sourceforge.jdistlib.T.cumulative(t,df,false,true)-r));
            maxNew=Math.max(maxNew,Math.abs(jdistlib.T.cumulative(t,df,false,true)-r));
            maxBetween=Math.max(maxBetween,Math.abs(net.sourceforge.jdistlib.T.cumulative(t,df,false,true)-jdistlib.T.cumulative(t,df,false,true)));
            maxRelative=Math.max(maxRelative,Math.abs(jdistlib.T.cumulative(t,df,false,true)-r)/Math.max(1,Math.abs(r)));
        }
        System.out.println("T_TAIL_R cases="+(tails.size()-1)+" old_max_abs="+maxOld+" new_max_abs="+maxNew);
        System.out.println("T_TAIL between_versions_max_abs="+maxBetween+" relative_log_error="+maxRelative);
        time("old-jdistlib-t-100000", ()->{for(int i=0;i<100000;i++)sink+=net.sourceforge.jdistlib.T.cumulative(i*.0001,120,false,true);});
        time("new-jdistlib-t-100000", ()->{for(int i=0;i<100000;i++)sink+=jdistlib.T.cumulative(i*.0001,120,false,true);});
        Random rng=new Random(4815); int n=512,p=8;
        double[][] x=new double[n][p]; double[] flat=new double[n*p], y=new double[n];
        for(int i=0;i<n;i++){ y[i]=rng.nextGaussian(); for(int j=0;j<p;j++)flat[i*p+j]=x[i][j]=j==0?1:rng.nextGaussian(); }
        CpuComputeBackend cpu=new CpuComputeBackend();
        double[] rhs=y.clone();
        double[] current=new QRDecomposition(x).solve(new Matrix(rhs,n)).getColumnPackedCopy();
        double[] candidate=cpu.dgeqp3(flat,n,p).solveLeastSquares(rhs);
        double error=0; for(int j=0;j<p;j++)error=Math.max(error,Math.abs(current[j]-candidate[j]));
        System.out.println("QR coefficient_max_abs="+error);
        time("current-qr-factor-solve-512x8",()->sink+=new QRDecomposition(x).solve(new Matrix(rhs,n)).get(0,0));
        time("candidate-qr-factor-solve-512x8",()->sink+=cpu.dgeqp3(flat,n,p).solveLeastSquares(rhs)[0]);
        time("jlinalg-full-ols-fit-512x8",()->sink+=Ols.fit(rhs,x,OlsOptions.defaults(),BackendPolicy.CPU).coefficients()[0]);
        int d=128; double[][] v=new double[d][d]; double[] vf=new double[d*d], b=new double[d];
        for(int i=0;i<d;i++){b[i]=rng.nextGaussian();for(int j=0;j<d;j++)vf[i*d+j]=v[i][j]=Math.pow(.4,Math.abs(i-j));}
        time("current-cholesky-factor-solve-128",()->sink+=new CholeskyDecomposition(new Matrix(v)).solve(new Matrix(b,d)).get(0,0));
        time("candidate-cholesky-factor-solve-128",()->sink+=cpu.dpotrf(vf,d).solve(b)[0]);
        // Package placement permits evaluating the published internal numerical tail,
        // without proposing its package-private API as a production integration point.
        int success=0,failures=0;double maxMixtureError=0;
        var mixtures=Files.readAllLines(Path.of("test/resources/audit-reference/mixtures.csv"));
        for(String line:mixtures.subList(1,mixtures.size())) {
            String[] row=line.split(","); double[] lambda=Arrays.stream(row[0].split(";")).mapToDouble(Double::parseDouble).toArray();
            try {
                var tail=QuadraticFormDistribution.survival(Double.parseDouble(row[1]),lambda);
                maxMixtureError=Math.max(maxMixtureError,Math.abs(tail.pValue()-Double.parseDouble(row[2])));success++;
            } catch(IllegalArgumentException ex) { failures++;System.out.println("JLIN_MIXTURE_UNRESOLVED "+row[0]+" q="+row[1]+" "+ex.getMessage()); }
        }
        System.out.println("JLIN_MIXTURE_R success="+success+" failures="+failures+" max_abs="+maxMixtureError);
        var mixedData=Files.readAllLines(Path.of("test/resources/audit-reference/mixed-data.csv"));
        int nm=mixedData.size()-1;double[][] xm=new double[nm][2];double[] ym=new double[nm];List<String> groups=new ArrayList<>();
        for(int i=0;i<nm;i++){String[] row=mixedData.get(i+1).split(",");groups.add(row[0]);xm[i][0]=1;xm[i][1]=Double.parseDouble(row[1]);ym[i]=Double.parseDouble(row[2]);}
        var components=List.of(VarianceComponent.randomIntercept("group",groups),VarianceComponent.identity("residual",nm));
        var mixedFit=Reml.fit(ym,xm,components,RemlOptions.defaults(),BackendPolicy.CPU);
        System.out.println("JLIN_REML converged="+mixedFit.converged()+" coefficients="+Arrays.toString(mixedFit.fixedEffects())+" variances="+Arrays.toString(mixedFit.varianceComponents()));
        double[] expected=Arrays.stream(Files.readAllLines(Path.of("test/resources/audit-reference/mixed-reference.csv")).get(1).split(",")).mapToDouble(Double::parseDouble).toArray();
        double remlError=0;
        for(int i=0;i<2;i++){remlError=Math.max(remlError,Math.abs(expected[i]-mixedFit.fixedEffects()[i]));remlError=Math.max(remlError,Math.abs(expected[2+i]-mixedFit.varianceComponents()[i]));}
        if(!mixedFit.converged()||remlError>2e-5)throw new AssertionError("REML differs from nlme: "+remlError);
        System.out.println("JLIN_REML_R max_abs_coefficient_or_variance="+remlError);
        time("jlinalg-reml-n40-g10",()->sink+=Reml.fit(ym,xm,components,RemlOptions.defaults(),BackendPolicy.CPU).fixedEffects()[0]);
        System.out.println("checksum="+sink);
    }
    private static void time(String name,Runnable job) {
        for(int i=0;i<100;i++)job.run();
        long[] times=new long[11];
        for(int i=0;i<times.length;i++){long t=System.nanoTime();for(int batch=0;batch<10;batch++)job.run();times[i]=(System.nanoTime()-t)/10;}
        Arrays.sort(times);System.out.printf(Locale.ROOT,"DEPENDENCY %s median_ms=%.6f min_ms=%.6f max_ms=%.6f%n",name,times[5]/1e6,times[0]/1e6,times[10]/1e6);
    }
}
