/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.gpu;

/** CPU geometry for explicitly selected Gaussian spectral null fits. Row-major FP64. */
public interface SpectralLinearAlgebra extends AutoCloseable {
    record Eigen(double[] values, double[] vectors) { }
    static SpectralLinearAlgebra open() {
        return new gov.nih.gpu.cpu.CpuSpectralLinearAlgebra();
    }
    String identity();
    double[] complement(double[][] orthonormalColumns);
    Eigen eigen(double[] symmetric, int dimension);
    double[] multiply(double[] left, boolean transposeLeft, double[] right,
        boolean transposeRight, int rows, int columns, int shared);
    double[] solvePositive(double[] matrix, int dimension, double[] right, int columns);
    @Override void close();
}
