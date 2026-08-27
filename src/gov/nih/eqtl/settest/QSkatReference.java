/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import gov.nih.jama.EigenvalueDecomposition;
import gov.nih.jama.Matrix;
import net.sourceforge.jdistlib.ChiSquare;
import gov.nih.utils.matrix.EMultiplicationMode;
import gov.nih.utils.matrix.QMatrixUtils;

/** Deterministic FP64 continuous-trait SKAT and simulation-adjusted SKAT-O reference. */
public final class QSkatReference {
    private static final double IMHOF_ABSOLUTE_TOLERANCE = 1e-9;
    private static final int IMHOF_INTERVALS_PER_BLOCK = 128;
    private static final int IMHOF_MAX_BLOCKS = 512;
    private static final int IMHOF_QUIET_BLOCKS = 6;

    public record MixturePValue(double pValue, String method, double degreesOfFreedom,
        double scale) { }
    public record Result(double statistic, double pValue, double log10P,
        double[] eigenvalues, String pValueMethod) {
        public Result { eigenvalues = eigenvalues.clone(); }
        @Override public double[] eigenvalues() { return eigenvalues.clone(); }
    }
    public record Component(double rho, Result result) { }
    public record OmnibusResult(List<Component> components, double minimumComponentP,
        double adjustedP, int simulations, long seed, String adjustmentMethod) {
        public OmnibusResult { components = List.copyOf(components); }
    }

    private QSkatReference() { }

    /** Rows are already covariate-residualized variants; weights multiply kernel columns. */
    public static Result calculate(double[][] residualVariants, double[] weights,
        double[] residualTrait, int nullResidualDf) {
        validate(residualVariants, weights, residualTrait, nullResidualDf);
        double sigmaSquared = dot(residualTrait, residualTrait) / nullResidualDf;
        if (!(sigmaSquared > 0) || !Double.isFinite(sigmaSquared))
            throw new IllegalArgumentException("SKAT null residual variance is invalid");
        double[][] weighted = weighted(residualVariants, weights);
        double statistic = scoreStatistic(weighted, residualTrait) / sigmaSquared;
        double[] eigenvalues = eigenvalues(gram(weighted));
        MixturePValue mixture = mixturePValue(statistic, eigenvalues);
        return new Result(statistic, mixture.pValue(), Math.log10(mixture.pValue()),
            eigenvalues, mixture.method());
    }

    /**
     * Deterministic parametric-null adjustment over the supplied rho grid. Simulations use
     * the same Gaussian vector for every component, preserving their null correlation.
     */
    public static OmnibusResult calculateOmnibus(double[][] residualVariants, double[] weights,
        double[] residualTrait, int nullResidualDf, double[] rhoGrid, int simulations, long seed) {
        validate(residualVariants, weights, residualTrait, nullResidualDf);
        if (rhoGrid == null || rhoGrid.length == 0 || simulations < 1)
            throw new IllegalArgumentException("SKAT-O requires a rho grid and positive simulation count");
        double[][] base = weighted(residualVariants, weights);
        double[] burden = new double[residualTrait.length];
        for (double[] variant : base)
            for (int sample = 0; sample < burden.length; sample++) burden[sample] += variant[sample];
        double sigmaSquared = dot(residualTrait, residualTrait) / nullResidualDf;
        List<FactorState> factors = new ArrayList<>();
        List<Component> components = new ArrayList<>();
        double minimum = 1;
        for (double rho : rhoGrid) {
            if (!Double.isFinite(rho) || rho < 0 || rho > 1)
                throw new IllegalArgumentException("SKAT-O rho values must be within [0,1]");
            double[][] factor = factor(base, burden, rho);
            double[] eigenvalues = eigenvalues(gram(factor));
            factors.add(new FactorState(factor, eigenvalues));
            Result result = calculateFactor(factor, eigenvalues, residualTrait, sigmaSquared);
            components.add(new Component(rho, result));
            minimum = Math.min(minimum, result.pValue());
        }
        double[] criticalStatistics = new double[factors.size()];
        for (int component = 0; component < factors.size(); component++)
            criticalStatistics[component] = criticalStatistic(
                factors.get(component).eigenvalues(), minimum);
        Random random = new Random(seed);
        int atLeastAsExtreme = 0;
        double[] gaussian = new double[residualTrait.length];
        for (int simulation = 0; simulation < simulations; simulation++) {
            for (int sample = 0; sample < gaussian.length; sample++) gaussian[sample] = random.nextGaussian();
            boolean atLeastAsExtremeSimulation = false;
            for (int component = 0; component < factors.size(); component++) {
                FactorState state = factors.get(component);
                double q = 0;
                for (double[] column : state.factor()) {
                    double score = dot(column, gaussian);
                    q += score * score;
                }
                if (q >= criticalStatistics[component]) {
                    atLeastAsExtremeSimulation = true;
                    break;
                }
            }
            if (atLeastAsExtremeSimulation) atLeastAsExtreme++;
        }
        double adjusted = (atLeastAsExtreme + 1.0) / (simulations + 1.0);
        return new OmnibusResult(components, minimum, adjusted, simulations, seed,
            "correlated-parametric-null");
    }

    public static MixturePValue mixturePValue(double statistic, double[] eigenvalues) {
        if (!Double.isFinite(statistic) || statistic < 0 || eigenvalues == null)
            throw new IllegalArgumentException("Invalid quadratic-form inputs");
        double first = 0, second = 0;
        int positive = 0;
        double firstPositive = Double.NaN;
        boolean equalPositive = true;
        for (double value : eigenvalues) if (value > 1e-12 && Double.isFinite(value)) {
            if (positive == 0) firstPositive = value;
            else if (Math.abs(value - firstPositive) > 1e-12 * Math.max(value, firstPositive))
                equalPositive = false;
            positive++;
            first += value; second += value * value;
        }
        if (!(first > 0) || !(second > 0))
            throw new IllegalArgumentException("Quadratic-form kernel has no positive eigenvalue");
        double df = first * first / second;
        double scale = second / first;
        if (positive == 1) {
            double p = ChiSquare.cumulative(statistic / scale, 1, false, false);
            return new MixturePValue(clampProbability(p), "exact-scaled-chi-square", 1, scale);
        }
        if (equalPositive) {
            double p = ChiSquare.cumulative(statistic / firstPositive,
                positive, false, false);
            return new MixturePValue(clampProbability(p),
                "exact-equal-eigenvalue-chi-square", positive, firstPositive);
        }
        Double imhof = imhofSurvival(statistic, eigenvalues);
        if (imhof != null)
            return new MixturePValue(clampProbability(imhof), "imhof-converged", df, scale);
        double p = ChiSquare.cumulative(statistic / scale, df, false, false);
        return new MixturePValue(clampProbability(p), "satterthwaite-fallback", df, scale);
    }

    private static Result calculateFactor(double[][] factor, double[] eigenvalues,
        double[] trait, double variance) {
        double q = 0;
        for (double[] column : factor) { double score = dot(column, trait); q += score * score / variance; }
        MixturePValue p = mixturePValue(q, eigenvalues);
        return new Result(q, p.pValue(), Math.log10(p.pValue()), eigenvalues, p.method());
    }

    /** Imhof inversion in fixed deterministic blocks; null means the convergence test failed. */
    private static Double imhofSurvival(double statistic, double[] eigenvalues) {
        double maximum = Arrays.stream(eigenvalues).max().orElse(0);
        if (!(maximum > 0)) return null;
        double blockWidth = 1.0 / maximum;
        double integral = 0;
        int quiet = 0;
        for (int block = 0; block < IMHOF_MAX_BLOCKS; block++) {
            double start = block * blockWidth;
            double contribution = simpson(start, start + blockWidth,
                IMHOF_INTERVALS_PER_BLOCK, statistic, eigenvalues);
            if (!Double.isFinite(contribution)) return null;
            integral += contribution;
            if (Math.abs(contribution) < IMHOF_ABSOLUTE_TOLERANCE) quiet++;
            else quiet = 0;
            if (quiet >= IMHOF_QUIET_BLOCKS) {
                double survival = 0.5 + integral / Math.PI;
                return survival >= -1e-7 && survival <= 1 + 1e-7 ? survival : null;
            }
        }
        return null;
    }

    private static double simpson(double start, double end, int intervals,
        double statistic, double[] eigenvalues) {
        double step = (end - start) / intervals;
        double total = imhofIntegrand(start, statistic, eigenvalues)
            + imhofIntegrand(end, statistic, eigenvalues);
        for (int i = 1; i < intervals; i++)
            total += (i % 2 == 0 ? 2 : 4)
                * imhofIntegrand(start + i * step, statistic, eigenvalues);
        return total * step / 3;
    }

    private static double imhofIntegrand(double t, double statistic,
        double[] eigenvalues) {
        if (t == 0) {
            double sum = 0;
            for (double eigenvalue : eigenvalues)
                if (eigenvalue > 1e-12 && Double.isFinite(eigenvalue)) sum += eigenvalue;
            return sum - statistic;
        }
        double angle = -statistic * t;
        double logMagnitudeDenominator = 0;
        for (double eigenvalue : eigenvalues) {
            if (!(eigenvalue > 1e-12) || !Double.isFinite(eigenvalue)) continue;
            double twice = 2 * eigenvalue * t;
            angle += 0.5 * Math.atan(twice);
            logMagnitudeDenominator += 0.25 * Math.log1p(twice * twice);
        }
        return Math.sin(angle) * Math.exp(-logMagnitudeDenominator) / t;
    }

    private static double clampProbability(double p) {
        return Math.max(Double.MIN_VALUE, Math.min(1, p));
    }

    private static double criticalStatistic(double[] eigenvalues, double survivalProbability) {
        if (survivalProbability >= 1) return 0;
        double high = Arrays.stream(eigenvalues).sum();
        for (int expansion = 0; expansion < 128
            && mixturePValue(high, eigenvalues).pValue() > survivalProbability; expansion++)
            high *= 2;
        double low = 0;
        for (int iteration = 0; iteration < 64; iteration++) {
            double middle = (low + high) * 0.5;
            if (mixturePValue(middle, eigenvalues).pValue() > survivalProbability) low = middle;
            else high = middle;
        }
        return high;
    }

    private static double[][] factor(double[][] base, double[] burden, double rho) {
        List<double[]> columns = new ArrayList<>();
        if (rho < 1) for (double[] variant : base) columns.add(scale(variant, Math.sqrt(1 - rho)));
        if (rho > 0) columns.add(scale(burden, Math.sqrt(rho)));
        return columns.toArray(double[][]::new);
    }
    private static double[][] weighted(double[][] rows, double[] weights) {
        double[][] result = new double[rows.length][];
        for (int i = 0; i < rows.length; i++) result[i] = scale(rows[i], weights[i]);
        return result;
    }
    private static double[] scale(double[] values, double scale) {
        double[] result = values.clone(); for (int i = 0; i < result.length; i++) result[i] *= scale; return result;
    }
    private static double[][] gram(double[][] rows) {
        return QMatrixUtils.parallelMatrixMultiplication(rows, QMatrixUtils.transpose(rows),
            null, 1, rows.length, rows.length, EMultiplicationMode.XY);
    }
    private static double scoreStatistic(double[][] rows, double[] trait) {
        double[][] column = new double[trait.length][1];
        for (int sample = 0; sample < trait.length; sample++) column[sample][0] = trait[sample];
        double[][] scores = QMatrixUtils.parallelMatrixMultiplication(rows, column, null, 1,
            rows.length, 1, EMultiplicationMode.XY);
        double statistic = 0;
        for (double[] score : scores) statistic += score[0] * score[0];
        return statistic;
    }
    private static double[] eigenvalues(double[][] matrix) {
        EigenvalueDecomposition decomposition = new EigenvalueDecomposition(new Matrix(matrix));
        return Arrays.stream(decomposition.getRealEigenvalues()).filter(value -> value > 1e-12).sorted().toArray();
    }
    private static double dot(double[] left, double[] right) {
        double value = 0; for (int i = 0; i < left.length; i++) value += left[i] * right[i]; return value;
    }
    private static void validate(double[][] variants, double[] weights, double[] trait, int df) {
        if (variants == null || variants.length == 0 || weights == null || variants.length != weights.length
            || trait == null || trait.length < 2 || df <= 0)
            throw new IllegalArgumentException("Invalid SKAT reference dimensions");
        for (int i = 0; i < variants.length; i++) {
            if (variants[i] == null || variants[i].length != trait.length || !(weights[i] > 0)
                || !Double.isFinite(weights[i])) throw new IllegalArgumentException("Invalid SKAT variant row or weight");
            for (double value : variants[i]) if (!Double.isFinite(value)) throw new IllegalArgumentException("Non-finite SKAT genotype");
        }
        for (double value : trait) if (!Double.isFinite(value)) throw new IllegalArgumentException("Non-finite SKAT trait");
    }

    private record FactorState(double[][] factor, double[] eigenvalues) { }
}
