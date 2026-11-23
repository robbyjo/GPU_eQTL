/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.gpu;

/** CPU geometry for explicitly selected Gaussian spectral null fits. Row-major FP64. */
public interface SpectralLinearAlgebra extends AutoCloseable {
    record Eigen(double[] values, double[] vectors) { }
    record QrProjection(double[] residual, double[] orthonormalColumns, double[] upper) { }
    static SpectralLinearAlgebra open() {
        return new gov.nih.gpu.cpu.CpuSpectralLinearAlgebra();
    }
    String identity();
    double[] complement(double[][] orthonormalColumns);
    Eigen eigen(double[] symmetric, int dimension);
    /** Thin full-rank QR and (I-QQ')response; returns row-major Q and R without
     * forming A'A or a full Q. Input arrays are not modified. */
    QrProjection projectResidual(double[] matrix, double[] response, int rows, int columns);
    double[] multiply(double[] left, boolean transposeLeft, double[] right,
        boolean transposeRight, int rows, int columns, int shared);
    double[] solvePositive(double[] matrix, int dimension, double[] right, int columns);
    @Override void close();
}
