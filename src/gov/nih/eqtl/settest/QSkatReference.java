/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import gov.nih.jama.EigenvalueDecomposition;
import gov.nih.jama.Matrix;
import net.sourceforge.jdistlib.ChiSquare;

/** Deterministic FP64 continuous-trait SKAT and simulation-adjusted SKAT-O reference. */
public final class QSkatReference {
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

    /** Immutable kernel geometry reused across traits sharing a complete-case null design. */
    public static final class PreparedKernel {
        private final double[][] weighted;
        private final double[] spectrum;
        private PreparedKernel(double[][] variants, double[] weights) {
            validate(variants, weights, new double[variants[0].length], 1);
            weighted = weighted(variants, weights);
            spectrum = eigenvalues(gram(weighted));
        }
        public Result calculate(double[] trait, int df) {
            if (trait == null || trait.length != weighted[0].length || df <= 0)
                throw new IllegalArgumentException("Invalid SKAT trait dimensions or DF");
            for (double value : trait) if (!Double.isFinite(value))
                throw new IllegalArgumentException("Non-finite SKAT trait");
            double variance = dot(trait, trait) / df;
            if (!(variance > 0) || !Double.isFinite(variance))
                throw new IllegalArgumentException("SKAT null residual variance is invalid");
            return calculateFactor(weighted, spectrum, trait, variance);
        }
    }

    public static PreparedKernel prepare(double[][] variants, double[] weights) {
        if (variants == null || variants.length == 0 || variants[0] == null)
            throw new IllegalArgumentException("SKAT variants are required");
        return new PreparedKernel(variants, weights);
    }

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
        return prepareOmnibus(residualVariants, weights, rhoGrid)
            .calculate(residualTrait, nullResidualDf, simulations, seed);
    }

    public static PreparedOmnibus prepareOmnibus(double[][] variants, double[] weights, double[] rhoGrid) {
        if (variants == null || variants.length == 0 || variants[0] == null)
            throw new IllegalArgumentException("SKAT variants are required");
        validate(variants, weights, new double[variants[0].length], 1);
        return new PreparedOmnibus(variants, weights, rhoGrid);
    }

    public static final class PreparedOmnibus {
        private final List<FactorState> factors = new ArrayList<>();
        private final double[] rhoGrid;
        private final double[][] base;
        private final double[] burden;
        private PreparedOmnibus(double[][] variants, double[] weights, double[] rhoGrid) {
            if (rhoGrid == null || rhoGrid.length == 0)
                throw new IllegalArgumentException("SKAT-O requires a rho grid");
            this.rhoGrid = rhoGrid.clone();
            base = weighted(variants, weights);
            burden = new double[variants[0].length];
            for (double[] variant : base)
                for (int sample = 0; sample < burden.length; sample++) burden[sample] += variant[sample];
            for (double rho : rhoGrid) {
                if (!Double.isFinite(rho) || rho < 0 || rho > 1)
                    throw new IllegalArgumentException("SKAT-O rho values must be within [0,1]");
                double[][] factor = factor(base, burden, rho);
                factors.add(new FactorState(factor, rho == 1 && dot(burden, burden) == 0
                    ? new double[0] : eigenvalues(gram(factor))));
            }
        }

        public OmnibusResult calculate(double[] residualTrait, int nullResidualDf, int simulations, long seed) {
            if (residualTrait == null || residualTrait.length != burden.length || nullResidualDf <= 0 || simulations < 1)
                throw new IllegalArgumentException("Invalid SKAT-O trait dimensions, DF, or simulations");
            for (double value : residualTrait) if (!Double.isFinite(value))
                throw new IllegalArgumentException("Non-finite SKAT-O trait");
            double sigmaSquared = dot(residualTrait, residualTrait) / nullResidualDf;
            if (!(sigmaSquared > 0) || !Double.isFinite(sigmaSquared))
                throw new IllegalArgumentException("SKAT null residual variance is invalid");
            List<Component> components = new ArrayList<>();
            double minimum = 1;
            for (int component = 0; component < rhoGrid.length; component++) {
                FactorState state = factors.get(component);
                Result result = calculateFactor(state.factor(), state.eigenvalues(), residualTrait, sigmaSquared);
                components.add(new Component(rhoGrid[component], result));
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
                double independentQ = scoreStatistic(base, gaussian);
                double burdenScore = dot(burden, gaussian);
                boolean atLeastAsExtremeSimulation = false;
                for (int component = 0; component < factors.size(); component++) {
                    double rho = rhoGrid[component];
                    double q = (1 - rho) * independentQ + rho * burdenScore * burdenScore;
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
    }

    public static MixturePValue mixturePValue(double statistic, double[] eigenvalues) {
        var tail = QQuadraticFormDistribution.survival(statistic, eigenvalues);
        return new MixturePValue(tail.pValue(), tail.method(), tail.degreesOfFreedom(), tail.scale());
    }

    private static Result calculateFactor(double[][] factor, double[] eigenvalues,
        double[] trait, double variance) {
        if (eigenvalues.length == 0)
            return new Result(0, 1, 0, eigenvalues, "degenerate-zero-kernel");
        double q = 0;
        for (double[] column : factor) { double score = dot(column, trait); q += score * score / variance; }
        MixturePValue p = mixturePValue(q, eigenvalues);
        return new Result(q, p.pValue(), Math.log10(p.pValue()), eigenvalues, p.method());
    }

    private static double criticalStatistic(double[] eigenvalues, double probability) {
        if (eigenvalues.length == 0) return Double.POSITIVE_INFINITY;
        return QQuadraticFormDistribution.critical(eigenvalues, probability);
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
        double[][] result = new double[rows.length][rows.length];
        for (int i = 0; i < rows.length; i++)
            for (int j = 0; j <= i; j++) result[i][j] = result[j][i] = dot(rows[i], rows[j]);
        return result;
    }
    private static double scoreStatistic(double[][] rows, double[] trait) {
        double statistic = 0;
        for (double[] row : rows) { double score = dot(row, trait); statistic += score * score; }
        return statistic;
    }
    private static double[] eigenvalues(double[][] matrix) {
        double maximum = 0;
        for (double[] row : matrix) for (double value : row) {
            if (!Double.isFinite(value)) throw new IllegalArgumentException("Non-finite SKAT kernel");
            maximum = Math.max(maximum, Math.abs(value));
        }
        if (!(maximum > 0)) throw new IllegalArgumentException("Zero SKAT kernel");
        double[][] normalized = new double[matrix.length][matrix.length];
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix.length; j++) normalized[i][j] = matrix[i][j] / maximum;
        EigenvalueDecomposition decomposition = new EigenvalueDecomposition(new Matrix(normalized));
        double[] spectrum = decomposition.getRealEigenvalues();
        double largest = Arrays.stream(spectrum).max().orElse(0);
        double tolerance = largest * Math.ulp(1.0) * matrix.length * 16;
        for (double value : spectrum) if (value < -tolerance)
            throw new IllegalArgumentException("SKAT kernel is not positive semidefinite");
        final double multiplier = maximum;
        return Arrays.stream(spectrum).filter(value -> value > tolerance)
            .map(value -> value * multiplier).sorted().toArray();
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
