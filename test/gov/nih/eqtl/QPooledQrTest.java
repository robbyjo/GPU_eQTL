/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import gov.nih.gpu.SpectralLinearAlgebra;
import org.junit.jupiter.api.Test;

class QPooledQrTest {
    @Test void portableProjectionPreservesAnalyticResidualAtExtremeCohortWeights() {
        check("java");
    }
    @Test void nativeProjectionPreservesAnalyticResidualAtExtremeCohortWeights() {
        check("openblas");
    }
    private static void check(String engine) {
        String previous=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try(var algebra=SpectralLinearAlgebra.open()) {
            int n=12,p=3;double[] a=new double[n*p],y=new double[n],expected=new double[n];
            // The noise is orthogonal to both the intercept and linear trend in
            // every cohort. Weighting makes two design columns almost parallel.
            double[] noise={1,-1,-1,1},scale={1e5,1,1e-3};
            for(int i=0;i<n;i++) {
                int group=i/4,j=i%4;double s=scale[group];
                a[i*p]=s;a[i*p+1]=group==0?s:0;a[i*p+2]=s*j;
                expected[i]=group==0?0:s*noise[j];
                y[i]=2*a[i*p]-.5*a[i*p+1]+.25*a[i*p+2]+expected[i];
            }
            double[] originalA=a.clone(),originalY=y.clone();
            var result=algebra.projectResidual(a,y,n,p);
            assertArrayEquals(originalA,a,0);assertArrayEquals(originalY,y,0);
            assertArrayEquals(expected,result.residual(),1e-9);
            double[] q=result.orthonormalColumns(),r=result.upper();
            for(int j=0;j<p;j++)for(int k=0;k<p;k++) {
                double dot=0;for(int i=0;i<n;i++)dot+=q[i*p+j]*q[i*p+k];
                assertEquals(j==k?1:0,dot,2e-14);
            }
            for(int i=0;i<n;i++)for(int j=0;j<p;j++) {
                double reconstructed=0;for(int k=0;k<p;k++)reconstructed+=q[i*p+k]*r[k*p+j];
                assertEquals(a[i*p+j],reconstructed,2e-10);
            }
            assertThrows(IllegalArgumentException.class,()->algebra.projectResidual(new double[2],new double[1],1,2));
        } finally {if(previous==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",previous);}
    }
}
