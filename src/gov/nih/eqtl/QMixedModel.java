/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.jama.QRDecomposition;
import org.jlinalg.compute.BackendPolicy;
import org.jlinalg.inference.DegreesOfFreedomMethod;
import org.jlinalg.pedigree.Pedigree;
import org.jlinalg.pedigree.PedigreeIndividual;
import org.jlinalg.reml.*;
import jdistlib.accelerator.CpuComputeBackend;

/** Gaussian REML adapter. JLinAlg owns fitting; GPU eQTL owns data and inference policy. */
public final class QMixedModel {
    private final double[][] design;
    private final List<VarianceComponent> components;
    private final RemlOptions options;
    private final String dfMethod;

    public record Association(double beta, double standardError, double statistic,
        double df, double log10P, String test, String status, RemlResult fit) { }

    public QMixedModel(double[][] design, List<VarianceComponent> components,
        int maximumIterations, String dfMethod) {
        if (design == null || design.length < 3 || design[0].length < 1
            || design.length <= design[0].length + 1)
            throw new IllegalArgumentException("Mixed model requires positive association residual DF");
        QRDecomposition qr = new QRDecomposition(design);
        if (!qr.isFullRank()) throw new IllegalArgumentException("Mixed-model design is rank deficient");
        this.design = Arrays.stream(design).map(double[]::clone).toArray(double[][]::new);
        this.components = List.copyOf(components);
        this.dfMethod = dfMethod;
        DegreesOfFreedomMethod method = switch (dfMethod) {
            case "residual" -> DegreesOfFreedomMethod.RESIDUAL_APPROXIMATION;
            case "satterthwaite" -> DegreesOfFreedomMethod.SATTERTHWAITE;
            case "kenward-roger" -> DegreesOfFreedomMethod.KENWARD_ROGER;
            default -> throw new IllegalArgumentException("--reml-df must be residual, satterthwaite, or kenward-roger");
        };
        options = RemlOptions.builder().maximumIterations(maximumIterations)
            .degreesOfFreedomMethod(method).build();
    }

    public static QMixedModel create(QeQTLAnalysisConfig config, double[][] design,
        QCovariateTable covariates, boolean[] selected, String[] sampleIds) throws Exception {
        int n = sampleIds.length;
        // Several dense n*n workspaces are needed by the fitter; fail before allocation.
        long estimate = Math.multiplyExact(Math.multiplyExact((long)n, n),
            8L * (24 + (config.getRandomCovariates() == null ? 0 : config.getRandomCovariates().length)));
        if (estimate > Runtime.getRuntime().maxMemory() * .65)
            throw new IllegalArgumentException("Dense REML workspace estimate " + estimate
                + " bytes exceeds 65% of JVM heap; increase -Xmx or reduce the cohort");
        List<VarianceComponent> components = new ArrayList<>();
        if (config.getPedigreeFilename() != null) {
            Pedigree pedigree = readPedigree(Path.of(config.getPedigreeFilename()));
            String[] individuals = config.get("pedigree_id_column") == null ? sampleIds
                : covariates.columnValues(config.get("pedigree_id_column"), selected);
            if (config.get("pedigree_id_column") == null && config.getGenotypeIdStripPrefix() != null) {
                individuals = individuals.clone();
                String prefix = config.getGenotypeIdStripPrefix();
                for (int i = 0; i < individuals.length; i++)
                    if (individuals[i].startsWith(prefix)) individuals[i] = individuals[i].substring(prefix.length());
            }
            double[] relationship = new double[Math.multiplyExact(n, n)];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
                relationship[i*n+j] = pedigree.relationship(individuals[i], individuals[j]);
            components.add(new VarianceComponent("additive-genetic", n, relationship));
        }
        if (config.getRandomCovariates() != null)
            for (String column : config.getRandomCovariates())
                components.add(VarianceComponent.randomIntercept("group:" + column,
                    Arrays.asList(covariates.columnValues(column, selected))));
        components.add(VarianceComponent.identity("residual", n));
        if (design == null) { design = new double[n][1]; for (double[] row : design) row[0] = 1; }
        return new QMixedModel(design, components, config.getRemlMaximumIterations(), config.getRemlDf());
    }

    /** Headered TSV: ID, SIRE, DAM; 0 or . denotes an unknown parent. IDs are literal. */
    public static Pedigree readPedigree(Path path) throws Exception {
        List<PedigreeIndividual> rows = new ArrayList<>();
        try (var reader = Files.newBufferedReader(path)) {
            String line = reader.readLine();
            if (!"ID\tSIRE\tDAM".equals(line))
                throw new IllegalArgumentException("Pedigree header must be ID<TAB>SIRE<TAB>DAM");
            int number = 1;
            while ((line = reader.readLine()) != null) {
                number++;
                String[] fields = line.split("\t", -1);
                if (fields.length != 3 || Arrays.stream(fields).anyMatch(String::isBlank))
                    throw new IllegalArgumentException("Malformed pedigree row " + number);
                for (String value : fields) if (!value.equals(value.trim()))
                    throw new IllegalArgumentException("Whitespace around pedigree ID at row " + number);
                if (fields[0].equals("0") || fields[0].equals("."))
                    throw new IllegalArgumentException("Pedigree ID cannot use an unknown-parent sentinel");
                rows.add(new PedigreeIndividual(fields[0], parent(fields[1]), parent(fields[2])));
            }
        }
        return Pedigree.of(rows);
    }
    private static String parent(String value) { return value.equals("0") || value.equals(".") ? null : value; }

    public NullState fitNull(double[] y) {
        if (y == null || y.length != design.length)
            throw new IllegalArgumentException("Mixed-model response length differs from design");
        for (double value : y) if (!Double.isFinite(value))
            throw new IllegalArgumentException("Mixed-model response contains missing/non-finite values");
        double[][] q = new QRDecomposition(design).getQ().getArray();
        double[] residual = y.clone();
        for (int column = 0; column < q[0].length; column++) {
            double product = 0; for (int i = 0; i < y.length; i++) product += q[i][column]*y[i];
            for (int i = 0; i < y.length; i++) residual[i] -= q[i][column]*product;
        }
        double residualNorm = 0, originalNorm = 0;
        for (int i = 0; i < y.length; i++) { residualNorm = Math.hypot(residualNorm,residual[i]); originalNorm = Math.hypot(originalNorm,y[i]); }
        if (!(residualNorm > 1e-12*originalNorm))
            throw new IllegalArgumentException("Mixed-model trait has no residual variation after fixed effects");
        double pilotSd = residualNorm/Math.sqrt(y.length-design[0].length), pilotVariance = pilotSd*pilotSd;
        if (!Double.isFinite(pilotVariance) || pilotVariance < Double.MIN_NORMAL)
            throw new IllegalArgumentException("Mixed-model trait variance is outside the supported FP64 range; rescale the trait");
        double minimumVariance = Math.max(Double.MIN_NORMAL,pilotVariance*1e-10);
        double maximumVariance = Math.min(Double.MAX_VALUE/y.length,pilotVariance*1e10);
        RemlResult fit = Reml.fit(y, design, components,
            RemlOptions.builder().maximumIterations(options.maximumIterations())
                .varianceBounds(minimumVariance,maximumVariance).build(), BackendPolicy.CPU);
        requireConverged(fit);
        int n = y.length, p = design[0].length;
        double[] covariance = new double[n*n], variance = fit.varianceComponents();
        for (int k = 0; k < components.size(); k++) {
            double[] basis = components.get(k).covariance();
            for (int i = 0; i < covariance.length; i++) covariance[i] += variance[k]*basis[i];
        }
        try (CpuComputeBackend backend = new CpuComputeBackend()) {
            var factor = backend.dpotrf(covariance, n);
            double[] identity = new double[n*n]; for (int i = 0; i < n; i++) identity[i*n+i] = 1;
            double[] projection = factor.solve(identity, n), x = new double[n*p];
            for (int i = 0; i < n; i++) System.arraycopy(design[i], 0, x, i*p, p);
            double[] inverseX = factor.solve(x, p), betaCov = fit.fixedEffectCovariance();
            double[] middle = new double[n*p];
            for (int i = 0; i < n; i++) for (int a = 0; a < p; a++)
                for (int b = 0; b < p; b++) middle[i*p+a] += inverseX[i*p+b]*betaCov[b*p+a];
            for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
                double correction = 0;
                for (int a = 0; a < p; a++) correction += middle[i*p+a]*inverseX[j*p+a];
                projection[i*n+j] = projection[j*n+i] = projection[i*n+j]-correction;
            }
            return new NullState(fit, projection, factor.solve(fit.residuals()), design, minimumVariance, maximumVariance);
        }
    }

    public Association refit(double[] y, double[] g, NullState state) {
        double[] centered = state.center(g);
        if (!state.estimable(centered, state.project(centered))) return degenerate(state.fit);
        int n = y.length, p = design[0].length;
        double[][] full = new double[n][p+1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(design[i], 0, full[i], 0, p); full[i][p] = centered[i];
        }
        var refitOptions = RemlOptions.builder().maximumIterations(options.maximumIterations())
            .degreesOfFreedomMethod(options.degreesOfFreedomMethod())
            .varianceBounds(state.minimumVariance,state.maximumVariance)
            .initialVariances(state.fit.varianceComponents()).build();
        RemlResult fit = Reml.fit(y, full, components, refitOptions, BackendPolicy.CPU);
        requireConverged(fit);
        double beta = fit.fixedEffects()[p], se = fit.standardErrors()[p];
        double df = dfMethod.equals("residual") ? n-p-1 : fit.degreesOfFreedom()[p];
        double t = beta/se;
        double logP = (Math.log(2) + jdistlib.T.cumulative(Math.abs(t), df, false, true))/Math.log(10);
        return new Association(beta, se, t, df, logP, "reml-wald-"+dfMethod, "ok", fit);
    }

    static void requireConverged(RemlResult fit) {
        if (!fit.converged()) throw new IllegalArgumentException("REML did not converge: " + fit.convergenceMessage());
    }
    static Association degenerate(RemlResult fit) {
        return new Association(Double.NaN, Double.NaN, 0, 1, 0,
            "unestimable", "monomorphic-or-collinear", fit);
    }
    public static final class NullState {
        private final RemlResult fit;
        private final double[] projection, projectedResponse;
        private final double norm;
        private final double minimumVariance, maximumVariance;
        private final double[][] q;
        NullState(RemlResult fit, double[] projection, double[] response, double[][] design, double minimumVariance, double maximumVariance) {
            this.fit = fit; this.projection = projection; this.projectedResponse = response;
            this.minimumVariance = minimumVariance; this.maximumVariance = maximumVariance;
            q = new QRDecomposition(design).getQ().getArray();
            double maximum = 0; int n = response.length;
            for (int i = 0; i < n; i++) { double sum = 0;
                for (int j = 0; j < n; j++) sum += Math.abs(projection[i*n+j]); maximum = Math.max(maximum, sum); }
            norm = maximum;
        }
        public RemlResult fit() { return fit; }
        public double[] projection() { return projection.clone(); }
        public double[] projectedResponse() { return projectedResponse.clone(); }
        /** Remove fixed-effect span before products to avoid offset cancellation. */
        public double[] center(double[] g) {
            double[] result = g.clone();
            for (double value : g) if (!Double.isFinite(value))
                throw new IllegalArgumentException("Mixed model received an unfilled predictor");
            for (int a = 0; a < q[0].length; a++) {
                double product = 0; for (int i = 0; i < g.length; i++) product += q[i][a]*g[i];
                for (int i = 0; i < g.length; i++) result[i] -= q[i][a]*product;
            }
            double original = 0, residual = 0;
            for (int i = 0; i < g.length; i++) { original = Math.hypot(original,g[i]); residual = Math.hypot(residual,result[i]); }
            if (!Double.isFinite(original) || !Double.isFinite(residual))
                throw new IllegalArgumentException("Mixed-model predictor magnitude exceeds the FP64 range");
            if (residual <= 1e-12*original) Arrays.fill(result, 0);
            return result;
        }
        public double[] project(double[] g) {
            int n = g.length; double[] result = new double[n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) result[i] += projection[i*n+j]*g[j];
            return result;
        }
        boolean estimable(double[] g, double[] pg) {
            double information = dot(g,pg), tolerance = 1e-12*norm*dot(g,g);
            if (!Double.isFinite(information) || !Double.isFinite(tolerance) || information < -tolerance)
                throw new IllegalArgumentException("Invalid mixed-model score information");
            return information > tolerance;
        }
        public Association score(double[] centered, double[] projected) {
            if (!estimable(centered, projected)) return degenerate(fit);
            double information = dot(centered,projected), score = dot(centered,projectedResponse);
            double statistic = score*score/information;
            double logP = jdistlib.ChiSquare.cumulative(statistic, 1, false, true)/Math.log(10);
            return new Association(score/information, 1/Math.sqrt(information), statistic, 1,
                logP, "null-reml-score-chi-square", "ok", fit);
        }
    }
    private static double dot(double[] a, double[] b) { double sum = 0; for (int i = 0; i < a.length; i++) sum += a[i]*b[i]; return sum; }
}
