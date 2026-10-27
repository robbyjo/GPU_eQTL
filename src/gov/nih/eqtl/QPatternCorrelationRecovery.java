/* Copyright 2026 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import gov.nih.eqtl.io.QMissingnessScan;
import gov.nih.gpu.GpuPrecision;
import gov.nih.jama.QRDecomposition;

/** Bounded scalar QR verification of a failed exact-pattern correlation. */
final class QPatternCorrelationRecovery {
    // This is a verification window, never a correlation-clamping tolerance.
    private static final double VERIFICATION_WINDOW = 1e-6;

    record Result(double correlation, double genotypeSd, double traitScale) { }

    private QPatternCorrelationRecovery() { }

    static boolean canVerify(double correlation) {
        return Double.isFinite(correlation)
            && Math.abs(correlation) > 1 + QeQTLStatistics.FP64_CORRELATION_TOLERANCE
            && Math.abs(correlation) <= 1 + VERIFICATION_WINDOW;
    }

    static Result verify(double correlation, double genotypeSd, double[] rawGenotype,
        double[] preparedTrait, double replacement, int[] observed, double[][] design) {
        if (!canVerify(correlation))
            throw new IllegalArgumentException("Correlation is outside the scalar verification window: "
                + correlation);
        double[][] selectedDesign = new double[observed.length][];
        double[] genotype = new double[observed.length];
        double[] trait = new double[observed.length];
        for (int i = 0; i < observed.length; i++) {
            int sample = observed[i];
            selectedDesign[i] = design[sample].clone();
            double value = rawGenotype[sample];
            genotype[i] = QMissingnessScan.isMissing(value) ? replacement : value;
            trait[i] = preparedTrait[sample];
        }
        QRDecomposition qr = new QRDecomposition(selectedDesign);
        if (!qr.isFullRank())
            throw new IllegalArgumentException("Scalar verification design is rank deficient");
        double[][] q = qr.getQ().getArray();
        double verifiedSd = residualizeAndStandardize(genotype, q);
        double traitScale = residualizeAndStandardize(trait, q);
        double verified = QeQTLStatistics.validateCorrelation(
            QeQTLPreprocessor.correlation(genotype, trait), GpuPrecision.FP64);
        if (Math.abs(verified - correlation) > VERIFICATION_WINDOW
            || !Double.isFinite(genotypeSd) || !(genotypeSd > 0)
            || Math.abs(verifiedSd / genotypeSd - 1) > VERIFICATION_WINDOW
            || Math.abs(traitScale - 1) > VERIFICATION_WINDOW)
            throw new IllegalArgumentException("Scalar QR verification disagrees with pattern statistics: "
                + "correlation=" + correlation + ", verified=" + verified
                + ", genotype SD=" + genotypeSd + ", verified SD=" + verifiedSd
                + ", trait scale=" + traitScale);
        return new Result(verified, verifiedSd, traitScale);
    }

    private static double residualizeAndStandardize(double[] values, double[][] q) {
        // Exact-pattern designs contain an intercept. Remove the origin first so a
        // large raw mean cannot contaminate a small residual through cancellation.
        center(values);
        for (int pass = 0; pass < 2; pass++) {
            double[] coefficients = new double[q[0].length];
            for (int sample = 0; sample < values.length; sample++)
                for (int column = 0; column < coefficients.length; column++)
                    coefficients[column] += q[sample][column] * values[sample];
            for (int sample = 0; sample < values.length; sample++)
                for (int column = 0; column < coefficients.length; column++)
                    values[sample] -= q[sample][column] * coefficients[column];
        }
        center(values);
        double norm = 0;
        for (double value : values) norm = Math.hypot(norm, value);
        double sd = norm / Math.sqrt(values.length - 1.0);
        if (!(sd > 0) || !Double.isFinite(sd))
            throw new IllegalArgumentException("Zero or invalid variance during scalar QR verification");
        for (int i = 0; i < values.length; i++) values[i] /= sd;
        return sd;
    }

    private static void center(double[] values) {
        double origin = values[0];
        double sum = 0;
        for (int i = 0; i < values.length; i++) {
            values[i] -= origin;
            sum += values[i];
        }
        double mean = sum / values.length;
        for (int i = 0; i < values.length; i++) values[i] -= mean;
    }
}
