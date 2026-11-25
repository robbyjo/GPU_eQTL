/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.*;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.concurrent.*;
import gov.nih.eqtl.QAnalysisCheckpoint;
import gov.nih.gpu.SpectralLinearAlgebra;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QSetNumericsPerformanceTest {
    @Test void nativeSkatOmnibusVectorsMatchPortableAtMultipleScales() {
        String old=System.getProperty("eqtl.cpu.blas");
        double[][] covariance={{1.5,.5,0},{.5,1.5,0},{0,0,0}};
        double[] scores={1,.3,0},rho={0,.5,1};
        var reference=new QSetScoreDistribution(covariance).omnibus(scores,rho,false,100,7,1e-6,1048576);
        try {
            for(String engine:new String[]{"java","openblas"}) {
                System.setProperty("eqtl.cpu.blas",engine);
                try(var algebra=SpectralLinearAlgebra.open()) {
                    for(double scale:new double[]{1e-100,1,1e100}) {
                        double[][] c=Arrays.stream(covariance).map(row->Arrays.stream(row).map(v->v*scale).toArray()).toArray(double[][]::new);
                        double[] u=Arrays.stream(scores).map(v->v*Math.sqrt(scale)).toArray();
                        var actual=new QSetScoreDistribution(c,algebra,true).omnibus(u,rho,false,100,7,1e-6,1048576);
                        assertEquals(reference.p(),actual.p(),reference.p()*2e-10);
                        assertEquals(reference.minimumP(),actual.minimumP(),reference.minimumP()*2e-10);
                    }
                }
            }
        } finally {if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old);}
    }
    @Test void interruptedProductionTailsStopWithoutReturningAProbability() {
        var distribution=new QSetScoreDistribution(new double[][]{{1,0},{0,2}});
        try {
            Thread.currentThread().interrupt();
            assertThrows(CancellationException.class,()->distribution.omnibus(new double[]{1,.3},new double[]{0,1},true,1,7,1e-5,1048576));
            assertThrows(CancellationException.class,()->distribution.omnibus(new double[]{1,.3},new double[]{0,1},false,1,7,1e-5,1048576));
            assertThrows(CancellationException.class,()->QQuadraticFormDistribution.survival(1.25,new double[]{1,2}));
            assertThrows(CancellationException.class,()->QSkatOmnibusIntegration.integrate(new double[]{1,2},new double[]{1,1},new double[]{0,1},new double[]{3,3},.1,1e-5,1048576,7));
        } finally {Thread.interrupted();}
    }
    @Test void failureCancelsARealSimulationWithoutHangingCleanup(@TempDir Path root)throws Exception {
        Path log=root.resolve("simulation-cancellation.log");
        var child=new ProcessBuilder(Path.of(System.getProperty("java.home"),"bin","java").toString(),"-Xmx128m","-cp",System.getProperty("java.class.path"),SimulationChild.class.getName())
            .redirectErrorStream(true).redirectOutput(log.toFile()).start();
        try {
            assertTrue(child.waitFor(15,TimeUnit.SECONDS),"Real simulation cancellation timed out");
            assertEquals(0,child.exitValue(),java.nio.file.Files.readString(log));
        } finally {if(child.isAlive()){child.destroyForcibly();child.waitFor();}}
    }
    public static final class SimulationChild {
        public static void main(String[] args)throws Exception {
            CountDownLatch started=new CountDownLatch(1),stopped=new CountDownLatch(1);
            var distribution=new QSetScoreDistribution(new double[][]{{1,0},{0,2}});
            try(var pipeline=new QSetResultPipeline(2)) {
                StringWriter writer=new StringWriter();
                pipeline.submit(writer,()->{assertTrue(started.await(5,TimeUnit.SECONDS));Thread.sleep(200);throw new ArithmeticException("synthetic failure");});
                pipeline.submit(writer,()->{started.countDown();try {distribution.omnibus(new double[]{1,.3},new double[]{0,1},true,Integer.MAX_VALUE,7,1e-5,1048576);return "unexpected";}finally{stopped.countDown();}});
                assertThrows(ArithmeticException.class,()->pipeline.drain(writer));
            }
            assertTrue(stopped.await(1,TimeUnit.SECONDS));
        }
    }
    @Test void nativeValuesOnlyAndVectorsMatchScalarOnRotatedSingularKernels() {
        String old=System.getProperty("eqtl.cpu.blas");
        try {
            for (String engine:new String[]{"java","openblas"}) {
                System.setProperty("eqtl.cpu.blas",engine);
                try (var algebra=SpectralLinearAlgebra.open()) {
                    int n=33;
                    double[][] c=new double[n][n];double[] u=new double[n];
                    // Known orthonormal DCT basis; zero, weak and repeated modes.
                    for (int k=0;k<n;k++) {
                        double lambda=k<4?0:k==4?1e-9:(k%7+1)/7.;
                        for (int i=0;i<n;i++) for (int j=0;j<n;j++)
                            c[i][j]+=lambda*Math.sqrt((k==0?1.:2.)/n)*Math.cos(Math.PI*(i+.5)*k/n)
                                *Math.sqrt((k==0?1.:2.)/n)*Math.cos(Math.PI*(j+.5)*k/n);
                    }
                    for (int i=0;i<n;i++) u[i]=Math.sin(i*.37)*3;
                    var scalar=new QSetScoreDistribution(c).skat(u);
                    for (double scale:new double[]{1e-100,1,1e100}) {
                        double[][] scaled=Arrays.stream(c).map(row->Arrays.stream(row).map(x->x*scale).toArray()).toArray(double[][]::new);
                        double[] scores=Arrays.stream(u).map(x->x*Math.sqrt(scale)).toArray();
                        var only=new QSetScoreDistribution(scaled,algebra,false);
                        var vectors=new QSetScoreDistribution(scaled,algebra,true);
                        assertEquals(scalar.p(),only.skat(scores).p(),scalar.p()*2e-10,engine);
                        assertEquals(scalar.p(),vectors.skat(scores).p(),scalar.p()*2e-10,engine);
                        assertThrows(IllegalStateException.class,()->only.omnibus(scores,new double[]{0,1},true,10,1,1e-6,1024));
                    }
                    double[] flat={2,1,1,2},before=flat.clone();
                    var values=algebra.eigen(flat,2,false);
                    assertArrayEquals(new double[]{1,3},values.values(),1e-14);
                    assertNull(values.vectors());assertArrayEquals(before,flat);
                    assertThrows(ArithmeticException.class,()->new QSetScoreDistribution(new double[][]{{1,2},{2,1}},algebra,false));
                }
            }
        } finally { if(old==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",old); }
    }
    @Test void costlySeriesHandsOffWithoutChangingTheTailLaw() {
        // Independent closed-form convolution: epsilon*chi-square(1)+chi-square(2).
        for (double epsilon:new double[]{.001,.005,.01}) for (double q:new double[]{.1,10,100,1000}) {
            double expected=Math.exp(-q/2)/Math.sqrt(1-epsilon)
                *jdistlib.ChiSquare.cumulative((1-epsilon)*q/epsilon,1,true,false)
                +jdistlib.ChiSquare.cumulative(q/epsilon,1,false,false);
            var tail=QQuadraticFormDistribution.survival(q,new double[]{epsilon,1,1});
            assertEquals(expected,tail.pValue(),expected*3e-11);
            assertTrue(tail.method().startsWith("tilted-characteristic-inversion"));
        }
        assertEquals("positive-gamma-series",QQuadraticFormDistribution.survival(3.25,new double[]{.5,2,7}).method());
    }
    @Test void earlyInversionFailureRetainsFullSeriesAndCancellationIsFatal() {
        double[] lambda={.01,.02,.05,.1,1};
        var recovered=QQuadraticFormDistribution.survival(100,lambda,(q,spectrum,lower)->{
            throw new ArithmeticException("Injected inversion work-limit failure");
        });
        assertEquals("positive-gamma-series",recovered.method());
        assertEquals(1.6748563178173308e-23,recovered.pValue(),recovered.pValue()*3e-11);
        assertThrows(CancellationException.class,()->QQuadraticFormDistribution.survival(100,lambda,
            (q,spectrum,lower)->{throw new CancellationException("Injected cancellation");}));
        var failed=assertThrows(IllegalArgumentException.class,()->QQuadraticFormDistribution.survival(100,new double[]{1e-12,.1,1},
            (q,spectrum,lower)->{throw new ArithmeticException("Injected inversion failure");}));
        assertEquals(1,failed.getSuppressed().length);
    }
    @Test void pendingTailsOverlapAndWriteOnlyInInputOrder() throws Exception {
        CountDownLatch second=new CountDownLatch(1);StringWriter writer=new StringWriter();
        try (var pipeline=new QSetResultPipeline(2)) {
            pipeline.submit(writer,()->{assertTrue(second.await(5,TimeUnit.SECONDS));return "first\n";});
            pipeline.submit(writer,()->{second.countDown();return "second\n";});
            pipeline.submit(writer,()->"third\n");pipeline.drain(writer);
        }
        assertEquals("first\nsecond\nthird\n",writer.toString());
        assertEquals(8,QSetResultPipeline.workers(0,64,32,1L<<30));
        assertEquals(32,QSetResultPipeline.workers(64,64,32,1L<<30));
        assertEquals(1,QSetResultPipeline.workers(64,64,64,32L<<20));
        assertEquals(1,QSetResultPipeline.workers(64,1,64,1L<<30));
    }
    @Test void failureCancelsOutstandingTailAndCannotCommitCheckpoint(@TempDir Path root) throws Exception {
        var checkpoint=QAnalysisCheckpoint.open(root.resolve("checkpoint"),"pipeline-test",1,false,true);
        CountDownLatch started=new CountDownLatch(1),stopped=new CountDownLatch(1);
        try (var pipeline=new QSetResultPipeline(2)) {
            assertThrows(ArithmeticException.class,()->checkpoint.writeBlock(0,writer->{
                pipeline.submit(writer,()->{assertTrue(started.await(5,TimeUnit.SECONDS));throw new ArithmeticException("tail failed");});
                pipeline.submit(writer,()->{started.countDown();try{new CountDownLatch(1).await();return "unreachable";}finally{stopped.countDown();}});
                pipeline.drain(writer);
            }));
            assertFalse(checkpoint.isComplete(0));
        }
        assertTrue(stopped.await(5,TimeUnit.SECONDS));
        assertFalse(checkpoint.isComplete(0));
    }
}
