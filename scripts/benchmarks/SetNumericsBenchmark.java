/* Copyright 2026 Roby Joehanes; GNU GPL version 3. Synthetic public benchmark. */
package gov.nih.eqtl.settest;
import java.util.*;
import gov.nih.gpu.SpectralLinearAlgebra;

public final class SetNumericsBenchmark {
    public static void main(String[] args) throws Exception {
        String mode=args[0];int repeats=Integer.parseInt(args[1]);
        System.out.printf(Locale.ROOT,"{\"kind\":\"environment\",\"mode\":\"%s\",\"java\":\"%s\",\"cores\":%d,\"threads\":\"%s\",\"heap\":%d}%n",
            mode,System.getProperty("java.version"),Runtime.getRuntime().availableProcessors(),System.getProperty("eqtl.cpu.threads"),Runtime.getRuntime().maxMemory());
        if (mode.equals("tail")) {
            QQuadraticFormDistribution.survival(3.25,new double[]{.5,2,7});
            for (double q:new double[]{10,100,500}) {
                long[] times=new long[repeats];double p=0;String method="";
                for (int j=0;j<repeats;j++) {
                    // Alternate spectra to exclude the thread-local coefficient-cache advantage.
                    QQuadraticFormDistribution.survival(3.25,new double[]{.5,2,7});
                    long start=System.nanoTime();var value=QQuadraticFormDistribution.survival(q,new double[]{.01,.02,.05,.1,1});
                    times[j]=System.nanoTime()-start;p=value.pValue();method=value.method();
                }
                Arrays.sort(times);
                System.out.printf(Locale.ROOT,"{\"kind\":\"tail\",\"q\":%.17g,\"p\":%.17g,\"median_ms\":%.9f,\"method\":\"%s\"}%n",q,p,times[times.length/2]/1e6,method);
            }
            return;
        }
        boolean nativeEigen=!mode.equals("eigen-java");
        try (SpectralLinearAlgebra algebra=nativeEigen?SpectralLinearAlgebra.open():null) {
            for (int n:new int[]{64,256,512,1024}) {
                double[][] c=new double[n][n];double[] u=new double[n];Arrays.fill(u,Math.sqrt(2));
                for(int i=0;i<n;i++){c[i][i]=2;if(i>0)c[i][i-1]=c[i-1][i]=.3;}
                long[] times=new long[repeats];double p=0;
                for (int j=0;j<repeats;j++) {
                    long start=System.nanoTime();
                    var distribution=nativeEigen?new QSetScoreDistribution(c,algebra,mode.equals("eigen-vectors")):new QSetScoreDistribution(c);
                    times[j]=System.nanoTime()-start;p=distribution.skat(u).p();
                }
                Arrays.sort(times);
                System.out.printf(Locale.ROOT,"{\"kind\":\"eigen\",\"n\":%d,\"p\":%.17g,\"median_ms\":%.9f}%n",n,p,times[times.length/2]/1e6);
            }
        }
    }
}
