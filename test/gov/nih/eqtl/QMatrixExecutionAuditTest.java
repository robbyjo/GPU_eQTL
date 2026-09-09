/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import gov.nih.utils.matrix.*;

class QMatrixExecutionAuditTest {
    @Test
    void workerExceptionsDoNotReturnPartialScientificResults() {
        double[][] x=new double[32][1024], y=new double[1024][128];
        y[17]=new double[127];
        assertThrows(IllegalStateException.class,()->QMatrixUtils.parallelMatrixMultiplication(x,y,null,1,32,128,EMultiplicationMode.XY));
    }
    @Test
    void threadedProductMatchesExplicitScalarForEveryCell() {
        int n=48,k=512,m=48;double[][] x=new double[n][k],y=new double[k][m];
        for(int i=0;i<n;i++)for(int j=0;j<k;j++)x[i][j]=Math.sin(i+j*.3);
        for(int i=0;i<k;i++)for(int j=0;j<m;j++)y[i][j]=Math.cos(i*.7-j);
        double[][] actual=QMatrixUtils.parallelMatrixMultiplication(x,y,null,1,n,m,EMultiplicationMode.XY);
        for(int i=0;i<n;i++)for(int j=0;j<m;j++){
            double expected=0;for(int z=0;z<k;z++)expected+=x[i][z]*y[z][j];
            assertEquals(expected,actual[i][j],0);
        }
    }
}
