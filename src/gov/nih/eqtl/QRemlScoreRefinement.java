/*
 * Copyright (C) 2026 JLinAlg contributors
 * Adapted for GPU eQTL by Roby Joehanes: score-checked REML plateau refinement.
 * SPDX-License-Identifier: GPL-2.0-or-later
 */
package gov.nih.eqtl;

import java.util.Arrays;
import java.util.List;
import jdistlib.accelerator.ComputeBackend;
import jdistlib.accelerator.MatrixTranspose;
import org.jlinalg.compute.BackendContext;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.internal.MatrixOps;
import org.jlinalg.reml.*;

/**
 * Refines a stalled null in stable coordinates using JLinAlg's REML score and
 * expected-information equations. A strict likelihood comparison can reject
 * the final Fisher step when its improvement rounds to zero. Only within a
 * 16-ulp likelihood plateau may a smaller projected score authorize a step;
 * success still requires the original projected-score tolerance. Neither the
 * statistical objective nor the variance bounds/tolerance is changed.
 */
final class QRemlScoreRefinement {
    private QRemlScoreRefinement() { }

    static RemlResult fit(double[] y,double[][] design,List<VarianceComponent> components,RemlOptions options) {
        int n=y.length,p=design[0].length,k=components.size();
        double[] x=MatrixOps.rowMajor(design,n),logs=options.initialVariances();
        if(logs==null || logs.length!=k)throw new IllegalArgumentException("REML refinement requires initial variances");
        for(int i=0;i<k;i++)logs[i]=Math.log(logs[i]);
        double lower=Math.log(options.minimumVariance()),upper=Math.log(options.maximumVariance());
        double[][] bases=components.stream().map(VarianceComponent::covariance).toArray(double[][]::new);
        try(var context=BackendContext.select(BackendPolicy.CPU)) {
            ComputeBackend backend=context.backend();
            Evaluation current=evaluate(y,x,n,p,bases,logs,backend);
            for(int iteration=1;iteration<=options.maximumIterations();iteration++) {
                if(scoreMaximum(current.score,logs,lower,upper)<=options.scoreTolerance()) {
                    System.out.println("REML null plateau resolved; maximum projected log-variance score="
                        +scoreMaximum(current.score,logs,lower,upper));
                    double[] se=new double[p],df=new double[p];Arrays.fill(df,n-p-1.0);
                    for(int j=0;j<p;j++)se[j]=Math.sqrt(current.betaCov[j*p+j]);
                    var association=AssociationStatistics.studentT(current.beta,se,df,
                        DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
                    return RemlResult.fromCoefficientSpace(components.stream().map(VarianceComponent::name).toList(),
                        Arrays.stream(logs).map(Math::exp).toArray(),association,current.betaCov,
                        MatrixOps.multiply(backend,x,n,p,current.beta),current.residual,current.likelihood,
                        VarianceEstimation.REML,n,p,iteration,true,context.provenance());
                }
                double[] step=informationStep(current.information,current.score,logs,lower,upper,backend);
                Search next=search(y,x,n,p,bases,logs,current,step,lower,upper,options,backend);
                if(next==null) {
                    double[] gradient=current.score.clone();
                    for(int j=0;j<k;j++)if(blocked(gradient[j],logs[j],lower,upper))gradient[j]=0;
                    next=search(y,x,n,p,bases,logs,current,gradient,lower,upper,options,backend);
                }
                if(next==null)throw new IllegalArgumentException("REML score refinement could not improve likelihood or projected score");
                logs=next.logs;current=next.evaluation;
            }
        }
        throw new IllegalArgumentException("REML score refinement exhausted the remaining iteration budget");
    }

    private static Search search(double[] y,double[] x,int n,int p,double[][] bases,double[] logs,
        Evaluation current,double[] step,double lower,double upper,RemlOptions options,ComputeBackend backend) {
        double largest=Arrays.stream(step).map(Math::abs).max().orElse(0);
        if(largest>options.maximumLogVarianceStep())
            for(int j=0;j<step.length;j++)step[j]*=options.maximumLogVarianceStep()/largest;
        double oldScore=scoreMaximum(current.score,logs,lower,upper);
        double plateau=16*Math.ulp(Math.max(1,Math.abs(current.likelihood)));
        for(int attempt=0;attempt<24;attempt++) {
            double scale=Math.scalb(1.0,-attempt);double[] candidate=logs.clone();
            for(int j=0;j<candidate.length;j++)candidate[j]=Math.max(lower,Math.min(upper,logs[j]+scale*step[j]));
            if(Arrays.equals(candidate,logs))continue;
            Evaluation next;
            try{next=evaluate(y,x,n,p,bases,candidate,backend);}
            catch(IllegalArgumentException | IllegalStateException invalid){continue;}
            if(next.likelihood>current.likelihood || (next.likelihood>=current.likelihood-plateau
                    && scoreMaximum(next.score,candidate,lower,upper)<oldScore))return new Search(candidate,next);
        }
        return null;
    }

    private static double[] informationStep(double[] information,double[] score,double[] logs,double lower,double upper,
        ComputeBackend backend) {
        int k=score.length;
        // A component whose score points out of an active bound is fixed for
        // this step. Its coupled Fisher direction must not move free components.
        int[] free=java.util.stream.IntStream.range(0,k).filter(j->!blocked(score[j],logs[j],lower,upper)).toArray();
        int dimension=free.length;double maximum=0;
        double[] reduced=new double[dimension*dimension],rhs=new double[dimension];
        for(int a=0;a<dimension;a++) {
            rhs[a]=score[free[a]];
            for(int b=0;b<dimension;b++)reduced[a*dimension+b]=information[free[a]*k+free[b]];
            maximum=Math.max(maximum,Math.abs(reduced[a*dimension+a]));
        }
        double ridge=Math.max(1e-12,maximum*1e-10);
        for(int attempt=0;attempt<10;attempt++,ridge*=10) {
            double[] regularized=reduced.clone();
            for(int j=0;j<dimension;j++)regularized[j*dimension+j]+=ridge;
            try {
                double[] direction=backend.dpotrf(regularized,dimension).solve(rhs),step=new double[k];
                for(int j=0;j<dimension;j++)step[free[j]]=direction[j];
                return step;
            }
            catch(IllegalArgumentException | IllegalStateException invalid){ /* Try a damped search direction. */ }
        }
        throw new IllegalArgumentException("REML refinement information is not numerically identifiable");
    }

    private static Evaluation evaluate(double[] y,double[] x,int n,int p,double[][] bases,double[] logs,
        ComputeBackend backend) {
        int k=bases.length;double[] variances=Arrays.stream(logs).map(Math::exp).toArray(),covariance=new double[n*n];
        for(int j=0;j<k;j++)for(int i=0;i<covariance.length;i++)covariance[i]+=variances[j]*bases[j][i];
        var factor=backend.dpotrf(covariance,n);
        double[] inverseX=factor.solve(x,p);
        double[] fixedInformation=MatrixOps.transposeMultiply(backend,x,n,p,inverseX,p);symmetrize(fixedInformation,p);
        var fixedFactor=backend.dpotrf(fixedInformation,p);
        double[] rhs=new double[p],inverseY=factor.solve(y);
        backend.dgemv(MatrixTranspose.TRANSPOSE,n,p,1,x,inverseY,0,rhs);
        double[] beta=fixedFactor.solve(rhs),residual=MatrixOps.subtract(y,MatrixOps.multiply(backend,x,n,p,beta));
        double[] py=factor.solve(residual),betaCov=fixedFactor.solve(MatrixOps.identity(p),p);
        double likelihood=-.5*((n-p)*Math.log(2*Math.PI)+factor.logDeterminant()+fixedFactor.logDeterminant()
            +backend.ddot(n,residual,0,1,py,0,1));
        double[] inverse=factor.solve(MatrixOps.identity(n),n),middle=MatrixOps.multiply(backend,inverseX,n,p,betaCov,p);
        double[] correction=new double[n*n];
        backend.dgemm(MatrixTranspose.NONE,MatrixTranspose.TRANSPOSE,n,n,p,1,middle,inverseX,0,correction);
        double[] projection=MatrixOps.subtract(inverse,correction);symmetrize(projection,n);
        double[] score=new double[k],information=new double[k*k];double[][] products=new double[k][];
        for(int j=0;j<k;j++) {
            products[j]=MatrixOps.multiply(backend,projection,n,n,bases[j],n);
            double trace=0;for(int i=0;i<n;i++)trace+=products[j][i*n+i];
            double[] kp=MatrixOps.multiply(backend,bases[j],n,n,py);
            score[j]=.5*variances[j]*(backend.ddot(n,py,0,1,kp,0,1)-trace);
            if(!Double.isFinite(score[j]))throw new IllegalArgumentException("Nonfinite REML refinement score");
        }
        for(int a=0;a<k;a++)for(int b=0;b<=a;b++) {
            double trace=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++)trace+=products[a][i*n+j]*products[b][j*n+i];
            information[a*k+b]=information[b*k+a]=.5*variances[a]*variances[b]*trace;
        }
        if(!Double.isFinite(likelihood))throw new IllegalArgumentException("Nonfinite REML refinement likelihood");
        return new Evaluation(beta,betaCov,residual,score,information,likelihood);
    }

    private static void symmetrize(double[] matrix,int n) {
        for(int i=0;i<n;i++)for(int j=0;j<i;j++)matrix[i*n+j]=matrix[j*n+i]=(matrix[i*n+j]+matrix[j*n+i])/2;
    }
    private static double scoreMaximum(double[] score,double[] logs,double lower,double upper) {
        double maximum=0;
        for(int j=0;j<score.length;j++) {
            if(blocked(score[j],logs[j],lower,upper))continue;
            maximum=Math.max(maximum,Math.abs(score[j]));
        }
        return maximum;
    }
    private static boolean blocked(double score,double log,double lower,double upper) {
        return (log<=lower+1e-12 && score<0)||(log>=upper-1e-12 && score>0);
    }
    private record Search(double[] logs,Evaluation evaluation) { }
    private record Evaluation(double[] beta,double[] betaCov,double[] residual,double[] score,double[] information,double likelihood) { }
}
