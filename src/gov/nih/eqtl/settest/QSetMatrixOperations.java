/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import gov.nih.eqtl.QeQTLAnalysis;
import gov.nih.gpu.*;
import java.util.*;

/** Tiled FP64 row products; scalar implementation is the independent execution oracle. */
public final class QSetMatrixOperations implements AutoCloseable {
    private final GpuContextPool pool;
    private final boolean cpuProjection;
    private final Map<String,Long> calls = new LinkedHashMap<>();
    private static final int TILE = 128;

    public QSetMatrixOperations(GpuContext[] contexts) {
        this(contexts,false);
    }
    public QSetMatrixOperations(GpuContext[] contexts,boolean cpuProjection) {
        this.cpuProjection = cpuProjection;
        pool = contexts == null || contexts.length == 0 ? null : GpuContextPool.borrowed(contexts);
        if (pool != null) {
            String source = "#define BLOCK_SIZE 16\n#define DATATYPE double\n#define N_MIN_1 1\n"
                + "#if defined(cl_khr_fp64)\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
                + "#elif defined(cl_amd_fp64)\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n"
                + QeQTLAnalysis.eqtlReal;
            for (GpuContext context : contexts) context.compileKernel(source, "eqtlReal", GpuPrecision.FP64);
        }
    }

    /** A B' with bounded device tiles, including partial final rows and columns. */
    public double[][] cross(double[][] a, double[][] b, String operation) {
        int k = a[0].length;
        for (double[] row : a) validateRow(row,k);
        for (double[] row : b) validateRow(row,k);
        double[][] result = new double[a.length][b.length];
        if (pool == null) {
            for (int i = 0; i < a.length; i++) for (int j = 0; j < b.length; j++)
                for (int h = 0; h < k; h++) result[i][j] += a[i][h]*b[j][h];
        } else {
            int width = round(k,64);
            List<GpuDevice> devices = pool.getAllContexts().stream().map(GpuContext::getDevice).toList();
            int sizingRows = Math.multiplyExact(TILE*4,devices.size());
            long availableHeap = Runtime.getRuntime().maxMemory()-Runtime.getRuntime().totalMemory()+Runtime.getRuntime().freeMemory();
            int recommended = devices.stream().allMatch(d -> d.getBackendName().equalsIgnoreCase("cpu"))
                ? GpuTuning.recommendCpuBlockSize(Math.max(2,k),sizingRows,sizingRows,GpuPrecision.FP64,Math.max(1,availableHeap)).blockSize()
                : GpuTuning.recommendBlockSize(devices,Math.max(2,k),sizingRows,sizingRows,GpuPrecision.FP64).blockSize();
            int tile = Math.min(TILE,recommended);
            for (int i0 = 0; i0 < a.length; i0 += tile) for (int j0 = 0; j0 < b.length; j0 += tile) {
                int nr = Math.min(tile,a.length-i0), nc = Math.min(tile,b.length-j0);
                int rows = round(nr,16), columns = round(nc,16);
                double[] left = new double[Math.multiplyExact(rows,width)];
                double[] right = new double[Math.multiplyExact(width,columns)];
                for (int i = 0; i < nr; i++) System.arraycopy(a[i0+i],0,left,i*width,k);
                for (int j = 0; j < nc; j++) for (int h = 0; h < k; h++) right[h*columns+j] = b[j0+j][h];
                GpuContext context = pool.reserveContext();
                try {
                    double[] product = context.executeDoubleKernel(left,right,rows*columns,
                        17L*64*8,width,columns,new long[]{columns,rows},new long[]{16,16});
                    for (int i = 0; i < nr; i++) System.arraycopy(product,i*columns,result[i0+i],j0,nc);
                    calls.merge(context.getDevice().getBackendName()+":"+operation,1L,Long::sum);
                } finally { pool.releaseContext(context); }
            }
        }
        for (double[] row : result) validateRow(row,b.length);
        return result;
    }

    public double[][] residualize(double[][] rows, double[][] q) {
        if (cpuProjection && pool != null) {
            try (var scalar = new QSetMatrixOperations(null)) { return scalar.residualize(rows,q); }
        }
        double[][] qt = transpose(q);
        double[][] coefficients = cross(rows,qt,"fixed-projection-coefficients");
        double[][] fitted = cross(coefficients,q,"fixed-projection-fitted");
        double[][] residual = new double[rows.length][rows[0].length];
        for (int i = 0; i < rows.length; i++) {
            double original = 0, remaining = 0;
            for (int j = 0; j < rows[i].length; j++) {
                residual[i][j] = rows[i][j]-fitted[i][j];
                original = Math.hypot(original,rows[i][j]); remaining = Math.hypot(remaining,residual[i][j]);
            }
            if (remaining <= 1e-12*original) Arrays.fill(residual[i],0);
        }
        return residual;
    }

    public String executionSignature() {
        StringBuilder result = new StringBuilder("set-products-fp64-v2;cpu-projection=").append(cpuProjection).append(";blas=")
            .append(System.getProperty("eqtl.cpu.blas","auto"))
            .append(";threads=").append(System.getProperty("eqtl.cpu.threads","auto"));
        if (pool == null) return result.append(";scalar").toString();
        for (GpuContext context : pool.getAllContexts()) {
            GpuDevice d = context.getDevice();
            result.append(';').append(d.getBackendName()).append(':').append(d.getName())
                .append(':').append(d.getDriverVersion()).append(':').append(d.getComputeApiVersion());
        }
        return result.toString();
    }
    public Map<String,Long> executedProducts() { return Map.copyOf(calls); }
    public static double[][] transpose(double[][] matrix) {
        double[][] result = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) for (int j = 0; j < matrix[i].length; j++) result[j][i] = matrix[i][j];
        return result;
    }
    private static void validateRow(double[] row,int size) {
        if (row == null || row.length != size) throw new IllegalArgumentException("Ragged set matrix");
        for (double value : row) if (!Double.isFinite(value)) throw new ArithmeticException("Non-finite FP64 set matrix product/input");
    }
    private static int round(int n,int multiple) { return Math.multiplyExact(Math.addExact(n,multiple-1)/multiple,multiple); }
    @Override public void close() { if (pool != null) pool.close(); }
}
