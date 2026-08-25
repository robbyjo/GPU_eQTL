/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.eqtl;

import java.util.List;
import java.util.stream.Collectors;

import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;
import gov.nih.gpu.GpuExecutionMetrics;
import gov.nih.gpu.GpuPrecision;

/** Applies one validated orthonormal covariate projection across exclusive compute contexts. */
final class QGpuResidualizer implements QeQTLPreprocessor.Residualizer, AutoCloseable {
    private final GpuContextPool contextPool;
    private final List<GpuContext> contexts;
    private final double[][] covariateQ;
    private final int sampleCount;
    private final int covariateRank;
    private final GpuPrecision precision;
    private final double[] flattenedQ;
    private final float[] flattenedQ32;
    private final QeQTLProfiler profiler;
    private final String cacheSignatureTag;
    private boolean closed;

    QGpuResidualizer(GpuContext[] contexts, double[][] covariateQ,
        GpuPrecision precision, QeQTLProfiler profiler) {
        if (contexts == null || contexts.length == 0)
            throw new IllegalArgumentException("At least one compute context is required for accelerated residualization");
        if (covariateQ == null || covariateQ.length == 0 || covariateQ[0].length == 0)
            throw new IllegalArgumentException("A non-empty covariate Q matrix is required");
        if (precision == null)
            throw new IllegalArgumentException("precision must not be null");
        this.contextPool = GpuContextPool.borrowed(contexts);
        this.contexts = contextPool.getAllContexts();
        this.covariateQ = covariateQ;
        this.sampleCount = covariateQ.length;
        this.covariateRank = covariateQ[0].length;
        this.precision = precision;
        this.profiler = profiler == null ? new QeQTLProfiler(false) : profiler;
		for (GpuContext context : this.contexts)
			context.setProfilingEnabled(this.profiler.isEnabled());
        for (int sample = 0; sample < sampleCount; sample++) {
            if (covariateQ[sample].length != covariateRank)
                throw new IllegalArgumentException("Covariate Q matrix is not rectangular");
            for (double value : covariateQ[sample])
                if (!Double.isFinite(value))
                    throw new IllegalArgumentException("Covariate Q matrix contains a non-finite value");
        }
        int projectionElements = Math.multiplyExact(sampleCount, covariateRank);
        if (precision == GpuPrecision.FP32) {
            flattenedQ = null;
            flattenedQ32 = new float[projectionElements];
            for (int sample = 0; sample < sampleCount; sample++)
                for (int column = 0; column < covariateRank; column++)
                    flattenedQ32[sample * covariateRank + column] = (float) covariateQ[sample][column];
        } else {
            flattenedQ = new double[projectionElements];
            flattenedQ32 = null;
            for (int sample = 0; sample < sampleCount; sample++)
                System.arraycopy(covariateQ[sample], 0, flattenedQ,
                    sample * covariateRank, covariateRank);
        }
        String backends = this.contexts.stream()
            .map(context -> context.getDevice().getBackendName())
            .collect(Collectors.joining("+"));
        cacheSignatureTag = "gpu-projection-v1-" + precision.optionName() + "-" + backends;
    }

    @Override
    public double[][] residualize(double[][] values, double[][] projection, String matrixName) {
        ensureOpen();
        if (projection != covariateQ)
            throw new IllegalArgumentException("Backend residualizer received a different covariate projection");
        if (values == null || values.length == 0)
            throw new IllegalArgumentException(matrixName + " block must not be empty");
        int rowCount = values.length;
        int elements = Math.multiplyExact(rowCount, sampleCount);
        long waitStarted = profiler.start();
        GpuContext context = contextPool.reserveContext();
        profiler.record(QeQTLProfiler.Phase.GPU_CONTEXT_WAIT, waitStarted, 1, 0);
        try {
            double[][] result = new double[rowCount][sampleCount];
            if (precision == GpuPrecision.FP32) {
                float[] flattened = new float[elements];
                for (int row = 0; row < rowCount; row++) {
                    validateRow(values[row], matrixName, row);
                    for (int sample = 0; sample < sampleCount; sample++)
                        flattened[row * sampleCount + sample] = (float) values[row][sample];
                }
                float[] residuals = context.residualizeFloatRows(flattened, flattenedQ32,
                    rowCount, sampleCount, covariateRank);
                for (int row = 0; row < rowCount; row++)
                    for (int sample = 0; sample < sampleCount; sample++)
                        result[row][sample] = residuals[row * sampleCount + sample];
            } else {
                double[] flattened = new double[elements];
                for (int row = 0; row < rowCount; row++) {
                    validateRow(values[row], matrixName, row);
                    System.arraycopy(values[row], 0, flattened, row * sampleCount, sampleCount);
                }
                double[] residuals = context.residualizeDoubleRows(flattened, flattenedQ,
                    rowCount, sampleCount, covariateRank);
                for (int row = 0; row < rowCount; row++)
                    System.arraycopy(residuals, row * sampleCount, result[row], 0, sampleCount);
            }
            recordMetrics(context.getLastExecutionMetrics());
            return result;
        } finally {
            contextPool.releaseContext(context);
        }
    }

    private void validateRow(double[] row, String matrixName, int rowNumber) {
        if (row == null || row.length != sampleCount)
            throw new IllegalArgumentException(matrixName + " row " + rowNumber + " has the wrong sample count");
    }

    private void recordMetrics(GpuExecutionMetrics metrics) {
		profiler.recordElapsed(QeQTLProfiler.Phase.GPU_RESIDUALIZATION_BUFFER_SETUP,
			metrics.bufferSetupNanoseconds(), 1, 0);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_RESIDUALIZATION_UPLOAD,
            metrics.uploadNanoseconds(), 1, metrics.uploadedBytes());
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_RESIDUALIZATION_COMPUTE,
            metrics.computeNanoseconds(), 1, 0);
        profiler.recordElapsed(QeQTLProfiler.Phase.GPU_RESIDUALIZATION_DOWNLOAD,
            metrics.downloadNanoseconds(), 1, metrics.downloadedBytes());
    }

    @Override
    public int concurrency() {
        return contexts.size();
    }

	@Override
	public int estimatedHostBytesPerValue() {
		// Raw/prepared doubles plus flattened input, backend output, and converted residual rows.
		return precision == GpuPrecision.FP32 ? 24 : 32;
	}

    @Override
    public String cacheSignatureTag() {
        return cacheSignatureTag;
    }

    private void ensureOpen() {
        if (closed)
            throw new IllegalStateException("Backend residualizer is closed");
    }

    @Override
    public void close() {
        if (closed)
            return;
        RuntimeException failure = null;
        for (GpuContext context : contexts) {
            try {
                context.releaseResidualizationResources();
            } catch (RuntimeException e) {
                if (failure == null)
                    failure = e;
                else
                    failure.addSuppressed(e);
            }
        }
        contextPool.close();
        closed = true;
        if (failure != null)
            throw failure;
    }
}
