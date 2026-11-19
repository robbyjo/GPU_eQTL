/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.util.*;
import gov.nih.gpu.SpectralLinearAlgebra;
import org.jlinalg.reml.VarianceComponent;
import org.junit.jupiter.api.Test;

class QSpectralRemlTest {
    @Test void resolvedGeometryEngineBindsSignatureAndCannotChangeDuringFit() {
        String old=System.getProperty("eqtl.cpu.blas");
        try {
            Data data=data(false);System.setProperty("eqtl.cpu.blas","java");
            var portable=model(data,true);String portableSignature=portable.signature();
            System.setProperty("eqtl.cpu.blas","openblas");
            var nativeModel=model(data,true);
            assertNotEquals(portableSignature,nativeModel.signature());
            assertThrows(IllegalStateException.class,()->portable.fitSetNull(data.y));
            System.setProperty("eqtl.cpu.blas","auto");
            // On this bundled-native test host auto resolves to the same actual engine.
            assertEquals(nativeModel.signature(),model(data,true).signature());
            System.setProperty("eqtl.cpu.blas","java");
            var fitted=portable.fitSetNull(data.y);assertTrue(fitted.usesSpectralSetNull());
            System.setProperty("eqtl.cpu.blas","openblas");
            assertThrows(IllegalStateException.class,()->portable.fitSetNull(data.y));
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    record Data(double[][] x,double[] k,double[] y) { }
    private Data data(boolean identity) {
        int n=56;double[][] x=new double[n][3];double[] k=new double[n*n],y=new double[n];Random random=new Random(191119);
        for(int i=0;i<n;i++) {
            x[i][0]=1;x[i][1]=(i-27.)/17;x[i][2]=Math.cos(.31*i);
            for(int j=0;j<n;j++)k[i*n+j]=i==j?1:!identity&&i/4==j/4?.5:0;
        }
        for(int group=0;group<n/4;group++) {
            double shared=random.nextGaussian();
            for(int j=0;j<4;j++){int i=group*4+j;y[i]=.2+.04*x[i][1]+.1*x[i][2]+.25*shared+.2*random.nextGaussian();}
        }
        return new Data(x,k,y);
    }
    private QMixedModel model(Data data,boolean spectral) {
        var model=new QMixedModel(data.x,List.of(new VarianceComponent("additive-genetic",data.y.length,data.k),
            VarianceComponent.identity("residual",data.y.length)),1000,"residual");
        return spectral?model.withSpectralSetNull():model;
    }
    private void compare(boolean identity,String engine) {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            Data data=data(identity);var fast=model(data,true).fitSetNull(data.y);var reference=model(data,false).fitSetNull(data.y);
            assertTrue(fast.usesSpectralSetNull());assertTrue(fast.fit().converged());
            assertEquals(identity,fast.combinedIdentityVariance());assertTrue(fast.spectralMaximumProjectedScore()<=1e-6);
            assertEquals(reference.fit().restrictedLogLikelihood(),fast.fit().restrictedLogLikelihood(),1e-6);
            assertArrayEquals(reference.fit().fixedEffects(),fast.fit().fixedEffects(),2e-6);
            assertArrayEquals(reference.fit().fixedEffectCovariance(),fast.fit().fixedEffectCovariance(),2e-6);
            assertArrayEquals(reference.projection(),fast.projection(),2e-5);
            assertArrayEquals(reference.projectedResponse(),fast.projectedResponse(),2e-5);
            double[][] factor=fast.spectralFactor();double[] projection=fast.projection();int n=data.y.length;
            for(int i=0;i<n;i++)for(int j=0;j<n;j++) {
                double value=0;for(int a=0;a<factor[0].length;a++)value+=factor[i][a]*factor[j][a];
                assertEquals(projection[i*n+j],value,Math.abs(value)*1e-12+1e-10);
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    @Test void portableFamilyFitMatchesIndependentDenseReml(){compare(false,"java");}
    @Test void portableIdentityFitPreservesCombinedVariance(){compare(true,"java");}
    @Test void nativeFamilyFitMatchesIndependentDenseReml(){compare(false,"openblas");}
    @Test void nativeIdentityFitPreservesCombinedVariance(){compare(true,"openblas");}
    @Test void geometryCacheUsesTheExactOrderedObservedMaskAndFitsEachResponse() {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas","java");
        try {
            Data data=data(false);var full=model(data,true);int[] mask=java.util.stream.IntStream.range(0,48).toArray();
            double[] y=Arrays.copyOf(data.y,48);var first=full.fitSetNull(y,mask);
            double[] changed=y.clone();for(int i=0;i<changed.length;i++)changed[i]=.3*y[i]*y[i]+.03*Math.sin(i);
            var second=full.fitSetNull(changed,mask);var independent=full.subset(mask).fitSetNull(changed);
            assertArrayEquals(independent.projection(),second.projection(),1e-10);
            assertFalse(Arrays.equals(first.fit().varianceComponents(),second.fit().varianceComponents()));
            int[] reverse=mask.clone();double[] reversedY=y.clone();
            for(int i=0;i<mask.length;i++){reverse[i]=mask[mask.length-1-i];reversedY[i]=y[y.length-1-i];}
            var reversed=full.fitSetNull(reversedY,reverse);
            assertEquals(first.fit().restrictedLogLikelihood(),reversed.fit().restrictedLogLikelihood(),1e-9);
            assertNotEquals(full.signature(),model(data,false).signature());
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    @Test void spectralGeometryRejectsInvalidCovarianceAndConstantTraits() {
        String old=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas","java");
        try {
            Data data=data(false);double[] wrong=data.k.clone();wrong[0]=-1;
            assertThrows(IllegalArgumentException.class,()->new QSpectralReml(data.x,wrong));
            double[] constant=new double[data.y.length];Arrays.fill(constant,1);
            assertThrows(IllegalArgumentException.class,()->model(data,true).fitSetNull(constant));
            int[] duplicate={0,1,1,3,4,5};assertThrows(IllegalArgumentException.class,()->model(data,true).fitSetNull(new double[6],duplicate));
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    @Test void nativeAndPortableGeometryAgreeInRowMajorLayout() {
        double[][] q={{.5, .5},{.5,-.5},{.5,.5},{.5,-.5}};
        String old=System.getProperty("eqtl.cpu.blas");
        try {
            for(String engine:List.of("java","openblas")) {
                System.setProperty("eqtl.cpu.blas",engine);
                try(var operations=SpectralLinearAlgebra.open()) {
                    double[] b=operations.complement(q);double[] gram=operations.multiply(b,true,b,false,2,2,4);
                    assertArrayEquals(new double[]{1,0,0,1},gram,1e-12);
                    var eigen=operations.eigen(new double[]{2,1,1,2},2);
                    assertArrayEquals(new double[]{1,3},eigen.values(),1e-12);
                    double[] weighted=eigen.vectors().clone();for(int i=0;i<2;i++)for(int j=0;j<2;j++)weighted[i*2+j]*=eigen.values()[j];
                    assertArrayEquals(new double[]{2,1,1,2},operations.multiply(weighted,false,eigen.vectors(),true,2,2,2),1e-12);
                }
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
}
