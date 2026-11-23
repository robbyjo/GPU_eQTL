/* Copyright 2026 Roby Joehanes; GNU GPL version 3.
 * Bounded FP64 score products across traits with the same pedigree geometry.
 */
package gov.nih.eqtl;
import java.util.*;
import gov.nih.gpu.*;

final class QPooledBatchProducts {
    private final QPooledPedigreeModel.Fit[] fits;
    private final GpuContextPool pool;
    private final int n,p,shared,columns,normColumns;
    private final double[] right,normRight;
    record Products(double[] score,double[] raw,boolean[] scalarNorm,int scoreStride,int normStride) { }
    QPooledBatchProducts(QPooledPedigreeModel.Fit[] fits,GpuContextPool pool) {
        if(fits.length==0||fits.length>256)throw new IllegalArgumentException("Batch requires 1..256 null fits");
        this.fits=fits.clone();this.pool=pool;n=fits[0].geometry.n;p=fits[0].geometry.p+1;
        for(var fit:fits)if(fit.geometry!=fits[0].geometry)throw new IllegalArgumentException("Batch requires identical observed-mask geometry");
        shared=round(n,64);columns=round(Math.multiplyExact(p,fits.length),16);normColumns=round(fits.length,16);
        checkMemory(0);
        right=new double[Math.multiplyExact(shared,columns)];normRight=new double[Math.multiplyExact(shared,normColumns)];
        for(int t=0;t<fits.length;t++) {
            double[][] values=fits[t].scoreColumns();
            for(int i=0;i<n;i++){System.arraycopy(values[i],0,right,i*columns+t*p,p);normRight[i*normColumns+t]=fits[t].weights[i];}
        }
    }
    private void checkMemory(int rows) {
        long left=Math.multiplyExact((long)rows,shared),rhs=Math.multiplyExact((long)shared,columns),normRhs=Math.multiplyExact((long)shared,normColumns);
        long scoreOut=Math.multiplyExact((long)rows,columns),normOut=Math.multiplyExact((long)rows,normColumns);
        long bytes=Math.multiplyExact(8L,Math.addExact(Math.addExact(2*left,rhs+normRhs),scoreOut+normOut));
        Runtime rt=Runtime.getRuntime();long available=rt.maxMemory()-(rt.totalMemory()-rt.freeMemory());
        if(bytes>available/2)throw new IllegalArgumentException("Pedigree batch exceeds bounded JVM heap budget; reduce batch or tile size");
        if(pool!=null)for(var c:pool.getAllContexts()) {
            var d=c.getDevice();long total=d.getGlobalMemoryBytes(),allocation=d.getMaxAllocationBytes();
            long largest=8*Math.max(Math.max(left,rhs),Math.max(normRhs,Math.max(scoreOut,normOut)));
            if((allocation>0&&largest>allocation)||(total>0&&bytes>total*3/4))throw new IllegalArgumentException("Pedigree batch exceeds device memory budget");
        }
    }
    Products calculate(double[][] g) throws Exception {
        int rows=round(g.length,16);checkMemory(rows);double[] left=new double[Math.multiplyExact(rows,shared)],squared=new double[left.length];boolean[] scalar=new boolean[g.length];
        for(int a=0;a<g.length;a++) {
            if(g[a].length!=n)throw new IllegalArgumentException("Genotype mask differs from batch");
            for(int i=0;i<n;i++){double v=g[a][i];if(!Double.isFinite(v))throw new ArithmeticException("Nonfinite rotated genotype");left[a*shared+i]=v;double vv=v*v;if(!Double.isFinite(vv)||(v!=0&&vv==0)){scalar[a]=true;vv=0;}squared[a*shared+i]=vv;}
        }
        double[] products,norms;
        if(pool==null||pool.getNumContexts()==0) {
            products=new double[Math.multiplyExact(rows,columns)];norms=new double[Math.multiplyExact(rows,normColumns)];
            for(int a=0;a<g.length;a++)for(int t=0;t<fits.length;t++){double[] q=fits[t].products(g[a]);System.arraycopy(q,0,products,a*columns+t*p,p);double raw=0;for(int i=0;i<n;i++){double v=g[a][i]*Math.sqrt(fits[t].weights[i]);raw+=v*v;}norms[a*normColumns+t]=raw;}
        } else {
            GpuContext context=pool.reserveContext();
            try {
                products=context.executeDoubleKernel(left,right,Math.multiplyExact(rows,columns),17L*64*8,shared,columns,new long[]{columns,rows},new long[]{16,16});
                norms=context.executeDoubleKernel(squared,normRight,Math.multiplyExact(rows,normColumns),17L*64*8,shared,normColumns,new long[]{normColumns,rows},new long[]{16,16});
            } finally {pool.releaseContext(context);}
        }
        return new Products(products,norms,scalar,columns,normColumns);
    }
    QPooledPedigreeModel.Association score(double[] g,Products q,int row,int trait,double criticalStatistic) {
        var fit=fits[trait];int offset=row*q.scoreStride+trait*p;double raw=q.raw[row*q.normStride+trait];double correction=0;
        for(int j=0;j<p;j++)if(!Double.isFinite(q.score[offset+j]))throw new ArithmeticException("Nonfinite batched pedigree score product");
        for(int j=1;j<p;j++)correction+=q.score[offset+j]*q.score[offset+j];
        double info=raw-correction,u=q.score[offset];
        // Preserve the original explicit residual/reprojection calculation at
        // variance bounds and collinearity, with a wider guard around its .01 cutoff.
        if(q.scalarNorm[row]||!Double.isFinite(raw)||raw<0||!Double.isFinite(info)||info<=.02*raw) {
            return fit.score(g,Arrays.copyOfRange(q.score,offset,offset+p));
        }
        double z=u/Math.sqrt(info),stat=z*z,beta=(u/info)*fit.responseScale,se=fit.responseScale/Math.sqrt(info);
        if(!Double.isFinite(stat)||!Double.isFinite(beta)||!Double.isFinite(se)||se<=0)throw new ArithmeticException("Nonfinite batched score statistic/effect/SE");
        // Avoid tail evaluation only safely below the output threshold. The
        // caller still applies the original log-p comparison to returned rows.
        if(stat<criticalStatistic-1e-10*Math.max(1,criticalStatistic))return null;
        return new QPooledPedigreeModel.Association(beta,se,stat,jdistlib.ChiSquare.cumulative(stat,1,false,true)/Math.log(10),"ok");
    }
    private static int round(int n,int b){return Math.multiplyExact(Math.addExact(n,b-1)/b,b);}
}

