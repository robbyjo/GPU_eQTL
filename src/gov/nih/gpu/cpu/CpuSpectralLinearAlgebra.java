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
