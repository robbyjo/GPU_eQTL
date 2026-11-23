/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledBoundaryTest {
    @TempDir Path temporary;

    @Test void zeroResidualCohortPreservesOtherCohortsAnalyticRemlVariances() throws Exception {
        checkBoundary("java");
    }
    @Test void nativeQrPreservesTheSameAnalyticBoundarySolution() throws Exception {
        checkBoundary("openblas");
    }
    private void checkBoundary(String engine) throws Exception {
        String previous=System.getProperty("eqtl.cpu.blas");System.setProperty("eqtl.cpu.blas",engine);
        try {
            int perGroup=24,n=3*perGroup;
            double[][] x=new double[n][3];double[] y=new double[n];String[] group=new String[n],ids=new String[n];
            StringBuilder pedigree=new StringBuilder("ID\tSIRE\tDAM\n");
            for(int i=0;i<n;i++) {
                int g=i/perGroup,j=i%perGroup;group[i]=new String[]{"FHS","WHI","CARDIA"}[g];ids[i]="synthetic-"+i;
                x[i][0]=1;x[i][1]=g==1?1:0;x[i][2]=g==2?1:0;
                y[i]=g==1?0:Math.sin(j*.71)*(g==0?1:.3)+Math.cos(j*.23)*.1;
                if(g==0)pedigree.append(ids[i]).append("\t0\t0\n");
            }
            Path path=temporary.resolve("pedigree.tsv");Files.writeString(path,pedigree);
            var model=new QPooledPedigreeModel(x,group,"FHS",ids,QMixedModel.readPedigree(path),1000);
            var geometry=model.geometry(java.util.stream.IntStream.range(0,n).toArray());
            var fit=geometry.fit(y);
            assertTrue(fit.maximumScore<=1e-6);assertTrue(fit.iterations<=1000);
            for(int g:new int[]{0,2}) {
                double mean=0,sum=0;for(int j=0;j<perGroup;j++)mean+=y[g*perGroup+j]/perGroup;
                for(int j=0;j<perGroup;j++)sum+=Math.pow(y[g*perGroup+j]-mean,2);
                int component=geometry.names.indexOf(g==0?"total:FHS":"residual:CARDIA");
                assertEquals(sum/(perGroup-1),fit.variance[component],2e-7);
            }
            assertTrue(fit.variance[geometry.names.indexOf("residual:WHI")]<1e-8);
        } finally {if(previous==null)System.clearProperty("eqtl.cpu.blas");else System.setProperty("eqtl.cpu.blas",previous);}
    }
}
