/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.gpu.cpu;

import java.util.concurrent.locks.ReentrantLock;
import gov.nih.gpu.SpectralLinearAlgebra;
import gov.nih.jama.Matrix;
import gov.nih.jama.QRDecomposition;
import jdistlib.accelerator.CpuComputeBackend;
import org.bytedeco.javacpp.Loader;
import static org.bytedeco.openblas.global.openblas_full.*;

/** LAPACK geometry, with a deterministic portable reference selected by cpu.blas=java. */
public final class CpuSpectralLinearAlgebra implements SpectralLinearAlgebra {
    private static final ReentrantLock LOCK = new ReentrantLock();
    private final boolean nativeBlas;
    private final int previousThreads;
    private boolean closed;

    public CpuSpectralLinearAlgebra() {
        LOCK.lock();
        boolean loaded = false; int previous = 0;
        try {
            String requested = System.getProperty("eqtl.cpu.blas", "auto");
            if (!requested.equals("java")) {
                try {
                    Loader.load(org.bytedeco.openblas.global.openblas_full.class);
                    previous = openblas_get_num_threads();
                    int threads = Integer.getInteger("eqtl.cpu.threads", 1);
                    if (threads < 1) throw new IllegalArgumentException("CPU thread count must be positive");
                    openblas_set_num_threads(threads);
                    loaded = true;
                } catch (LinkageError | RuntimeException error) {
                    if (!requested.equals("auto")) throw error;
                }
            }
        } catch (Throwable error) {
            LOCK.unlock(); throw error;
        }
        nativeBlas = loaded; previousThreads = previous;
    }
    public String identity() { return nativeBlas ? "openblas-full-0.3.34-lapack-fp64" : "portable-spectral-fp64"; }
    private static void require(int code, String operation) {
        if (code != 0) throw new ArithmeticException(operation+" failed: LAPACK info="+code);
    }
    public double[] complement(double[][] q) {
        int n=q.length,p=q[0].length,d=n-p;
        double[] b=new double[Math.multiplyExact(n,d)];
        if (nativeBlas) {
            double[] full=new double[Math.multiplyExact(n,n)],tau=new double[p];
            for(int i=0;i<n;i++)System.arraycopy(q[i],0,full,i*n,p);
            require(LAPACKE_dgeqrf(101,n,p,full,n,tau),"Spectral QR");
            require(LAPACKE_dorgqr(101,n,n,p,full,n,tau),"Spectral complete Q");
            for(int i=0;i<n;i++)System.arraycopy(full,i*n+p,b,i*d,d);
        } else {
            double[][] packed=new QRDecomposition(q).getQRArray();
            for(int j=0;j<d;j++)b[(p+j)*d+j]=1;
            for(int a=p-1;a>=0;a--)if(packed[a][a]!=0)for(int j=0;j<d;j++) {
                double dot=0;for(int i=a;i<n;i++)dot+=packed[i][a]*b[i*d+j];
                dot=-dot/packed[a][a];for(int i=a;i<n;i++)b[i*d+j]+=dot*packed[i][a];
            }
        }
        return b;
    }
    public Eigen eigen(double[] matrix,int n) {
        double[] values=new double[n],vectors=matrix.clone();
        if(nativeBlas)require(LAPACKE_dsyevd(101,(byte)'V',(byte)'L',n,vectors,n,values),"Spectral symmetric eigenproblem");
        else {
            double[][] copy=new double[n][n];for(int i=0;i<n;i++)System.arraycopy(matrix,i*n,copy[i],0,n);
            var eig=new Matrix(copy).eig();values=eig.getRealEigenvalues();double[][] v=eig.getV().getArray();
            for(int i=0;i<n;i++)System.arraycopy(v[i],0,vectors,i*n,n);
        }
        return new Eigen(values,vectors);
    }
    public QrProjection projectResidual(double[] matrix,double[] response,int n,int p) {
        if(n<p||p<1||matrix.length!=Math.multiplyExact(n,p)||response.length!=n)
            throw new IllegalArgumentException("Spectral QR dimensions differ");
        double[] residual=response.clone(),orthogonal=new double[matrix.length],upper=new double[p*p];
        if(nativeBlas) {
            double[] factor=matrix.clone(),tau=new double[p];
            require(LAPACKE_dgeqrf(101,n,p,factor,p,tau),"Weighted thin QR");
            for(int i=0;i<p;i++)for(int j=i;j<p;j++)upper[i*p+j]=factor[i*p+j];
            require(LAPACKE_dormqr(101,(byte)'L',(byte)'T',n,1,p,factor,p,tau,residual,1),"Weighted QR response");
            java.util.Arrays.fill(residual,0,p,0);
            require(LAPACKE_dormqr(101,(byte)'L',(byte)'N',n,1,p,factor,p,tau,residual,1),"Weighted QR residual");
            require(LAPACKE_dorgqr(101,n,p,p,factor,p,tau),"Weighted thin Q");
            orthogonal=factor;
        } else {
            double[][] a=new double[n][p];for(int i=0;i<n;i++)System.arraycopy(matrix,i*p,a[i],0,p);
            var qr=new QRDecomposition(a,0);
            if(qr.getRank()!=p)throw new ArithmeticException("Weighted design lost rank");
            double[][] packed=qr.getQRArray(),q=qr.getQArray(),r=qr.getRArray();
            for(int k=0;k<p;k++)reflect(packed,residual,k,n);
            java.util.Arrays.fill(residual,0,p,0);
            for(int k=p-1;k>=0;k--)reflect(packed,residual,k,n);
            for(int i=0;i<n;i++)System.arraycopy(q[i],0,orthogonal,i*p,p);
            for(int i=0;i<p;i++)System.arraycopy(r[i],0,upper,i*p,p);
        }
        for(int i=0;i<p;i++)if(!Double.isFinite(upper[i*p+i])||upper[i*p+i]==0)
            throw new ArithmeticException("Weighted design lost rank");
        return new QrProjection(residual,orthogonal,upper);
    }
    private static void reflect(double[][] packed,double[] vector,int k,int n) {
        double sum=0;for(int i=k;i<n;i++)sum+=packed[i][k]*vector[i];
        sum=-sum/packed[k][k];for(int i=k;i<n;i++)vector[i]+=sum*packed[i][k];
    }
    public double[] multiply(double[] a,boolean ta,double[] b,boolean tb,int m,int n,int k) {
        if(a.length!=Math.multiplyExact(m,k)||b.length!=Math.multiplyExact(k,n))throw new IllegalArgumentException("Spectral product dimensions differ");
        double[] c=new double[Math.multiplyExact(m,n)];
        if(nativeBlas)cblas_dgemm(101,ta?112:111,tb?112:111,m,n,k,1,a,ta?m:k,b,tb?k:n,0,c,n);
        else for(int i=0;i<m;i++)for(int j=0;j<n;j++) {
            double sum=0;for(int h=0;h<k;h++)sum+=a[ta?h*m+i:i*k+h]*b[tb?j*k+h:h*n+j];c[i*n+j]=sum;
        }
        return c;
    }
    public double[] solvePositive(double[] a,int n,double[] b,int columns) {
        if(a.length!=Math.multiplyExact(n,n)||b.length!=Math.multiplyExact(n,columns))throw new IllegalArgumentException("Spectral solve dimensions differ");
        if(nativeBlas) {
            double[] factor=a.clone(),answer=b.clone();
            require(LAPACKE_dpotrf(101,(byte)'L',n,factor,n),"Spectral covariance Cholesky");
            require(LAPACKE_dpotrs(101,(byte)'L',n,columns,factor,n,answer,columns),"Spectral covariance solve");
            return answer;
        }
        try(var backend=new CpuComputeBackend()){return backend.dpotrf(a,n).solve(b,columns);}
    }
    public void close() {
        if(!closed){closed=true;try{if(nativeBlas)openblas_set_num_threads(previousThreads);}finally{LOCK.unlock();}}
    }
}
