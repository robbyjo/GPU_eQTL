/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import org.junit.jupiter.api.Test;

class QSpectralScaleTest {
    @Test void traitUnitsDoNotChangeTheProfileOptimumOrConvergence() {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas","openblas");
        try {
            int n=512;double[][] x=new double[n][3];double[] k=new double[n*n],y=new double[n];
            Random random=new Random(20261119);
            for(int group=0;group<n/4;group++) {
                double family=random.nextGaussian();
                for(int member=0;member<4;member++) {
                    int i=group*4+member;x[i]=new double[]{1,(i-128.)/90,Math.sin(.15*i)};
                    y[i]=.3+.2*x[i][1]+.1*x[i][2]+.25*family+.4*random.nextGaussian();
                    for(int j=0;j<n;j++)k[i*n+j]=i==j?1:i/4==j/4?.5:0;
                }
            }
            var geometry=new QSpectralReml(x,k);var reference=geometry.fit(y);
            for(double scale:new double[]{1e-100,1e-50,1e50,1e100}) {
                double[] response=Arrays.stream(y).map(v->v*scale).toArray();
                var actual=geometry.fit(response);
                assertTrue(actual.maximumProjectedScore()<1e-8,"scale="+scale+" score="+actual.maximumProjectedScore());
                double[] a=actual.fit().varianceComponents(),b=reference.fit().varianceComponents();
                for(int i=0;i<a.length;i++)assertEquals(b[i],a[i]/scale/scale,Math.max(Math.abs(b[i])*2e-7,Arrays.stream(b).sum()*1e-8),"scale="+scale);
                assertEquals(reference.fit().restrictedLogLikelihood(),actual.fit().restrictedLogLikelihood()+(n-3)*Math.log(scale),1e-7);
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
}
