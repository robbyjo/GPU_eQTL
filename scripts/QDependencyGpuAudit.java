/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
import java.util.*;
import jdistlib.accelerator.*;

/** Same transfer-inclusive GEMM shapes/data as GPU_eQTL's GpuBackendBenchmark. */
public final class QDependencyGpuAudit {
    private static volatile double sink;
    public static void main(String[] args) {
        for(Compute requested:new Compute[]{Compute.CUDA,Compute.OPENCL}) {
            try(var selection=ComputeBackends.select(requested)) {
                ComputeBackend backend=selection.backend();
                System.out.println("CANDIDATE_BACKEND "+backend.deviceInfo());
                for(int[] shape:new int[][]{{128,2048},{512,1024},{2048,512}}) {
                    int rows=shape[0],width=shape[1];
                    double[] a=new double[width*rows],b=new double[rows*width],c=new double[width*width];
                    for(int i=0;i<a.length;i++)a[i]=((i*17L)%101-50)/51.0;
                    for(int i=0;i<b.length;i++)b[i]=((i*31L)%97-48)/49.0;
                    Runnable job=()->backend.dgemm(MatrixTranspose.NONE,MatrixTranspose.NONE,width,width,rows,1.0/(rows-1),a,b,0,c);
                    for(int i=0;i<3;i++)job.run();long[] times=new long[9];
                    for(int i=0;i<times.length;i++){long start=System.nanoTime();job.run();times[i]=System.nanoTime()-start;sink+=c[i];}
                    double maxError=0;
                    for(int i:new int[]{0,1,width/2,width-1})for(int j:new int[]{0,2,width/2,width-1}) {
                        double expected=0;for(int k=0;k<rows;k++)expected+=a[i*rows+k]*b[k*width+j];
                        maxError=Math.max(maxError,Math.abs(expected/(rows-1)-c[i*width+j]));
                    }
                    if(maxError>1e-11)throw new AssertionError("GEMM reference mismatch "+maxError);
                    Arrays.sort(times);System.out.printf(Locale.ROOT,"CANDIDATE_GEMM %s rows=%d width=%d median_ms=%.6f max_reference_error=%.3g%n",requested,rows,width,times[4]/1e6,maxError);
                }
            } catch(IllegalStateException | LinkageError ex) { System.out.println("CANDIDATE_UNAVAILABLE "+requested+" "+ex); }
        }
    }
}
