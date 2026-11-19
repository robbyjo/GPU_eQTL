/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.util.*;
import gov.nih.gpu.SpectralLinearAlgebra;
import gov.nih.jama.QRDecomposition;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.compute.BackendProvenance;
import org.jlinalg.inference.AssociationStatistics;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.reml.RemlResult;
import org.jlinalg.reml.VarianceEstimation;

/** Two-component Gaussian REML in exact fixed-effect error contrasts. No low-rank approximation. */
final class QSpectralReml {
    record Fitted(RemlResult fit,double[] factor,int columns,double[] response,double lower,double upper,
                  double maximumProjectedScore,boolean combinedIdentity) { }
    private record Evaluation(double t,double objective,double genetic,double residual,double scoreG,double scoreE,double derivative) { }
    private final int n,p,d;
    private final double[][] design,q;
    private final double[] qFlat,k,h,eigen,transform;
    private final double logR;
    private final boolean identity;
    private final String backendIdentity;
    private int evaluations;

    QSpectralReml(double[][] design,double[] covariance) {
        this.design=design;n=design.length;p=design[0].length;d=n-p;
        if(covariance.length!=Math.multiplyExact(n,n))throw new IllegalArgumentException("Spectral covariance dimensions differ");
        k=covariance.clone();boolean same=true;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++) {
            double value=k[i*n+j];
            if(!Double.isFinite(value)||Math.abs(value-k[j*n+i])>1e-12*Math.max(1,Math.abs(value)))
                throw new IllegalArgumentException("Spectral pedigree covariance must be finite and symmetric");
            same&=value==(i==j?1:0);
        }
        identity=same;
        QRDecomposition qr=new QRDecomposition(design);
        if(!qr.isFullRank()||d<=1)throw new IllegalArgumentException("Spectral null requires full rank and positive association DF");
        q=qr.getQArray();qFlat=flatten(q);double[][] r=qr.getRArray();
        double determinant=0;for(int i=0;i<p;i++)determinant+=Math.log(Math.abs(r[i][i]));logR=determinant;
        double[][] inverse=new double[p][p];transform=new double[p*p];
        for(int column=0;column<p;column++)for(int i=p-1;i>=0;i--) {
            double value=i==column?1:0;for(int j=i+1;j<p;j++)value-=r[i][j]*inverse[j][column];
            inverse[i][column]=value/r[i][i];
        }
        int[] pivot=qr.getPivot();for(int i=0;i<p;i++)System.arraycopy(inverse[i],0,transform,pivot[i]*p,p);
        try(var operations=SpectralLinearAlgebra.open()) {
            backendIdentity=operations.identity();
            // Cholesky cheaply verifies the usual positive-definite pedigree. A singular
            // PSD kernel uses a full eigen check rather than being rejected or regularized.
            if(!identity)try{operations.solvePositive(k,n,new double[n],1);}catch(ArithmeticException|IllegalArgumentException error) {
                double[] spectrum=operations.eigen(k,n).values();checkSpectrum(spectrum);
            }
            double[] b=operations.complement(q);
            if(identity){h=b;eigen=new double[d];Arrays.fill(eigen,1);}
            else {
                double[] kb=operations.multiply(k,false,b,false,n,d,n);
                double[] kernel=operations.multiply(b,true,kb,false,d,d,n);symmetrize(kernel,d);
                var decomposition=operations.eigen(kernel,d);eigen=decomposition.values();checkSpectrum(eigen);
                h=operations.multiply(b,false,decomposition.vectors(),false,n,d,d);
            }
        }
    }
    private static void checkSpectrum(double[] values) {
        double scale=1;for(double value:values){if(!Double.isFinite(value))throw new ArithmeticException("Nonfinite spectral eigenvalue");scale=Math.max(scale,Math.abs(value));}
        for(int i=0;i<values.length;i++) {
            if(values[i]<-1e-10*scale)throw new IllegalArgumentException("Spectral pedigree covariance is not positive semidefinite");
            // Only roundoff-negative eigenvalues of a PSD input are clamped; positive
            // eigenvalues are never truncated, merged or approximated.
            if(values[i]<0)values[i]=0;
        }
    }
    private static double[] flatten(double[][] matrix) {
        int columns=matrix[0].length;double[] values=new double[Math.multiplyExact(matrix.length,columns)];
        for(int i=0;i<matrix.length;i++)System.arraycopy(matrix[i],0,values,i*columns,columns);return values;
    }
    private static void symmetrize(double[] a,int n) {
        for(int i=0;i<n;i++)for(int j=0;j<i;j++)a[i*n+j]=a[j*n+i]=.5*a[i*n+j]+.5*a[j*n+i];
    }
    Fitted fit(double[] y) {
        if(y.length!=n)throw new IllegalArgumentException("Spectral response dimensions differ");
        double originalNorm=0;for(double value:y){if(!Double.isFinite(value))throw new IllegalArgumentException("Nonfinite spectral response");originalNorm=Math.hypot(originalNorm,value);}
        boolean intercept=design[0][0]!=0;for(int i=1;i<n;i++)intercept&=design[i][0]==design[0][0];
        double anchor=intercept?y[0]:0;double[] centered=y.clone();for(int i=0;i<n;i++)centered[i]-=anchor;
        double[] residual=centered.clone();
        for(int a=0;a<p;a++){double dot=0;for(int i=0;i<n;i++)dot+=q[i][a]*centered[i];for(int i=0;i<n;i++)residual[i]-=q[i][a]*dot;}
        double norm=0;for(double value:residual)norm=Math.hypot(norm,value);
        if(!(norm>1e-12*originalNorm))throw new IllegalArgumentException("Spectral trait has no residual variation after fixed effects");
        double sd=norm/Math.sqrt(d),pilot=sd*sd;
        if(!Double.isFinite(pilot)||pilot<Double.MIN_NORMAL)throw new IllegalArgumentException("Spectral variance is outside the supported FP64 range; rescale the trait");
        double lo=Math.max(Double.MIN_NORMAL,pilot*1e-10),hi=Math.min(Double.MAX_VALUE/n,pilot*1e10);
        double[] z=new double[d],z2=new double[d];
        for(int j=0;j<d;j++){double dot=0;for(int i=0;i<n;i++)dot+=h[i*d+j]*residual[i];z[j]=dot;z2[j]=dot*dot;}
        double[] theta,v=new double[d];double maximumScore,objective;evaluations=0;
        if(identity) {
            double total=0;for(double value:z2)total+=value;total=Math.max(2*lo,Math.min(2*hi,total/d));
            theta=new double[]{total};Arrays.fill(v,total);double score=0;for(double value:z2)score+=.5*(value/total-1);
            maximumScore=Math.abs(score);evaluations=1;
        } else {
            Evaluation optimum=optimize(z2,lo,hi);
            theta=new double[]{optimum.genetic,optimum.residual};
            for(int j=0;j<d;j++)v[j]=theta[0]*eigen[j]+theta[1];
            double sg=optimum.scoreG,se=optimum.scoreE;
            if(theta[0]<=lo*(1+1e-8)&&sg<0||theta[0]>=hi*(1-1e-8)&&sg>0)sg=0;
            if(theta[1]<=lo*(1+1e-8)&&se<0||theta[1]>=hi*(1-1e-8)&&se>0)se=0;
            maximumScore=Math.max(Math.abs(sg),Math.abs(se));
        }
        if(!Double.isFinite(maximumScore)||maximumScore>1e-6)
            throw new ArithmeticException("Spectral REML failed its projected-score check: "+maximumScore);
        objective=d*Math.log(2*Math.PI);for(int j=0;j<d;j++)objective+=Math.log(v[j])+z2[j]/v[j];
        double logLikelihood=-.5*objective-logR;
        double[] response=new double[n],factor=h.clone();
        for(int i=0;i<n;i++)for(int j=0;j<d;j++){response[i]+=h[i*d+j]*(z[j]/v[j]);factor[i*d+j]/=Math.sqrt(v[j]);}
        // Supply genuine original-coordinate coefficient diagnostics for the public
        // RemlResult adapter; no fabricated beta, covariance or residual entries.
        RemlResult fit;
        try(var operations=SpectralLinearAlgebra.open()) {
            double[] covariance=k.clone();
            for(int i=0;i<covariance.length;i++)covariance[i]*=identity?0:theta[0];
            for(int i=0;i<n;i++)covariance[i*n+i]+=identity?theta[0]:theta[1];
            double[] rhs=new double[n*(p+1)];
            for(int i=0;i<n;i++){System.arraycopy(q[i],0,rhs,i*(p+1),p);rhs[i*(p+1)+p]=centered[i];}
            double[] solved=operations.solvePositive(covariance,n,rhs,p+1),vq=new double[n*p],vy=new double[n];
            for(int i=0;i<n;i++){System.arraycopy(solved,i*(p+1),vq,i*p,p);vy[i]=solved[i*(p+1)+p];}
            double[] information=operations.multiply(qFlat,true,vq,false,p,p,n);symmetrize(information,p);
            double[] unit=new double[p*p];for(int i=0;i<p;i++)unit[i*p+i]=1;
            double[] covarianceQ=operations.solvePositive(information,p,unit,p),scoreQ=new double[p];
            for(int a=0;a<p;a++)for(int i=0;i<n;i++)scoreQ[a]+=q[i][a]*vy[i];
            double[] betaQ=operations.multiply(covarianceQ,false,scoreQ,false,p,1,p);
            double[] beta=operations.multiply(transform,false,betaQ,false,p,1,p);
            if(intercept)beta[0]+=anchor/design[0][0];
            double[] middle=operations.multiply(transform,false,covarianceQ,false,p,p,p);
            double[] betaCov=operations.multiply(middle,false,transform,true,p,p,p),se=new double[p];
            for(int a=0;a<p;a++)se[a]=Math.sqrt(Math.max(0,betaCov[a*p+a]));
            double[] fitted=new double[n],errors=new double[n];
            for(int i=0;i<n;i++){double value=0;for(int a=0;a<p;a++)value+=q[i][a]*betaQ[a];fitted[i]=anchor+value;errors[i]=centered[i]-value;}
            var statistics=AssociationStatistics.studentT(beta,se,d,DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION);
            fit=RemlResult.fromCoefficientSpace(identity?List.of("total-genetic-plus-residual"):List.of("additive-genetic","residual"),
                theta,statistics,betaCov,fitted,errors,logLikelihood,VarianceEstimation.REML,n,p,evaluations,true,
                new BackendProvenance(BackendPolicy.CPU,backendIdentity,"Spectral Gaussian error contrasts",!backendIdentity.startsWith("portable"),false));
        }
        return new Fitted(fit,factor,d,response,lo,hi,maximumScore,identity);
    }
    private Evaluation evaluate(double t,double[] z2,double lo,double hi) {
        if(++evaluations>65536)throw new ArithmeticException("Spectral REML exceeded 65536 profile evaluations");
        double ratio=Math.exp(t),sum=0;
        for(int j=0;j<d;j++){double shape=1+ratio*eigen[j];if(!(shape>0))return new Evaluation(t,Double.POSITIVE_INFINITY,0,0,0,0,Double.NaN);sum+=z2[j]/shape;}
        double lower=Math.max(lo,lo/ratio),upper=Math.min(hi,hi/ratio),raw=sum/d;
        if(lower>upper*(1+1e-12))return new Evaluation(t,Double.POSITIVE_INFINITY,0,0,0,0,Double.NaN);
        double e=Math.max(lower,Math.min(upper,raw)),g=e*ratio,objective=0,sg=0,se=0;
        for(int j=0;j<d;j++) {
            double v=e*(1+ratio*eigen[j]),difference=z2[j]/v-1;
            objective+=.5*(Math.log(v)+z2[j]/v);
            sg+=.5*(g*eigen[j]/v)*difference;se+=.5*(e/v)*difference;
        }
        double derivative=-sg;
        if((raw<lower||raw>upper)&&(e==lo/ratio||e==hi/ratio))derivative=se;
        return new Evaluation(t,objective,g,e,sg,se,derivative);
    }
    private Evaluation optimize(double[] z2,double lo,double hi) {
        double bound=Math.log(hi)-Math.log(lo),step=2*bound/1024;
        Evaluation[] grid=new Evaluation[1025];List<Evaluation> candidates=new ArrayList<>();
        for(int i=0;i<grid.length;i++)grid[i]=evaluate(-bound+i*step,z2,lo,hi);
        candidates.add(grid[0]);candidates.add(grid[1024]);
        for(int i=1;i<1024;i++)if(grid[i].objective<=grid[i-1].objective&&grid[i].objective<=grid[i+1].objective) {
            double a=grid[i-1].t,b=grid[i+1].t,ratio=(Math.sqrt(5)-1)/2;
            double left=b-ratio*(b-a),right=a+ratio*(b-a);
            Evaluation fl=evaluate(left,z2,lo,hi),fr=evaluate(right,z2,lo,hi);
            for(int iteration=0;iteration<100&&b-a>1e-12;iteration++) {
                if(fl.objective<fr.objective){b=right;right=left;fr=fl;left=b-ratio*(b-a);fl=evaluate(left,z2,lo,hi);}
                else{a=left;left=right;fl=fr;right=a+ratio*(b-a);fr=evaluate(right,z2,lo,hi);}
            }
            candidates.add(fl.objective<fr.objective?fl:fr);
            a=grid[i-1].t;b=grid[i+1].t;
            Evaluation fa=grid[i-1],fb=grid[i+1];
            if(fa.derivative*fb.derivative<0) {
                for(int iteration=0;iteration<100&&b-a>1e-13;iteration++) {
                    double mid=.5*a+.5*b;Evaluation fm=evaluate(mid,z2,lo,hi);
                    if(fa.derivative*fm.derivative<=0){b=mid;fb=fm;}else{a=mid;fa=fm;}
                }
                candidates.add(evaluate(.5*a+.5*b,z2,lo,hi));
            }
        }
        double minimum=candidates.stream().mapToDouble(Evaluation::objective).min().orElseThrow();
        return candidates.stream().filter(e->e.objective-minimum<=1e-10)
            .min(Comparator.comparingDouble(e->Math.abs(e.derivative))).orElseThrow();
    }
}
