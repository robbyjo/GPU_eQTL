/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.util.*;
import gov.nih.jama.*;

/** Score-space SKAT/SKAT-O: U~N(0,C), including fitted related-sample covariance. */
public final class QSetScoreDistribution {
    public record Result(double statistic,double p,double minimumP,String components,String method) { }
    private final double[] lambda, burden;
    private final int dimension;

    public QSetScoreDistribution(double[][] covariance) {
        if (covariance == null || covariance.length == 0) throw new IllegalArgumentException("Empty score covariance");
        int m = covariance.length;
        dimension = m;
        for (double[] row : covariance) if (row == null || row.length != m) throw new IllegalArgumentException("Score covariance must be square");
        double scale = 0;
        for (double[] row : covariance) for (double value : row) {
            if (!Double.isFinite(value)) throw new ArithmeticException("Non-finite score covariance");
            scale = Math.max(scale,Math.abs(value));
        }
        if (!(scale > 0)) throw new IllegalArgumentException("Zero score covariance");
        double[][] normalized = new double[m][m];
        for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
            normalized[i][j] = .5*(covariance[i][j]/scale+covariance[j][i]/scale);
        var eig = new EigenvalueDecomposition(new Matrix(normalized));
        double[] all = eig.getRealEigenvalues(); double[][] v = eig.getV().getArray();
        double tolerance = Arrays.stream(all).max().orElse(0)*Math.ulp(1.0)*m*16;
        for (double value : all) if (value < -tolerance)
            throw new ArithmeticException("Set score covariance is not positive semidefinite");
        int rank = (int)Arrays.stream(all).filter(x -> x > tolerance).count();
        if (rank == 0) throw new IllegalArgumentException("Zero projected set kernel");
        lambda = new double[rank]; burden = new double[rank];
        for (int j = 0,k = 0; j < m; j++) if (all[j] > tolerance) {
            lambda[k] = all[j]*scale;
            for (int i = 0; i < m; i++) {
                burden[k] += v[i][j]*Math.sqrt(lambda[k]);
            }
            k++;
        }
    }
    public Result skat(double[] scores) {
        validateScores(scores);
        double q = 0; for (double u : scores) q += u*u;
        var tail = QQuadraticFormDistribution.survival(q,lambda);
        return new Result(q,tail.pValue(),Double.NaN,"",tail.method());
    }
    public Result omnibus(double[] scores,double[] rho,boolean simulation,int simulations,long seed,
        double tolerance,int maxEvaluations) {
        validateScores(scores);
        double skat = 0, b = 0;
        for (double u : scores) { skat += u*u; b += u; }
        List<double[]> spectra = new ArrayList<>();
        double minimum = 1;
        StringJoiner components = new StringJoiner(";");
        for (double r : rho) {
            double[] spectrum = spectrum(r); spectra.add(spectrum);
            double q = (1-r)*skat+r*b*b;
            double p = spectrum.length == 0 ? 1 : QQuadraticFormDistribution.survival(q,spectrum).pValue();
            minimum = Math.min(minimum,p); components.add(r+":"+p);
        }
        double[] critical = new double[rho.length];
        for (int i = 0; i < rho.length; i++) critical[i] = spectra.get(i).length == 0
            ? Double.POSITIVE_INFINITY : QQuadraticFormDistribution.critical(spectra.get(i),minimum);
        if (!simulation) {
            var result = QSkatOmnibusIntegration.integrate(lambda,burden,rho,critical,minimum,tolerance,maxEvaluations,seed);
            return new Result(Double.NaN,result.probability(),minimum,components.toString(),result.method()
                + ";estimated-error="+result.error()+";evaluations="+result.evaluations()+";seed="+seed);
        }
        Random random = new Random(seed); int extreme = 0;
        for (int iteration = 0; iteration < simulations; iteration++) {
            double qs = 0, qb = 0;
            for (int j = 0; j < lambda.length; j++) {
                double z = random.nextGaussian(); qs += lambda[j]*z*z; qb += burden[j]*z;
            }
            for (int i = 0; i < rho.length; i++) if ((1-rho[i])*qs+rho[i]*qb*qb >= critical[i]) { extreme++; break; }
        }
        return new Result(Double.NaN,(extreme+1.0)/(simulations+1.0),minimum,components.toString(),
            "correlated-parametric-null;simulations="+simulations+";seed="+seed);
    }
    private void validateScores(double[] scores) {
        if (scores == null || scores.length != dimension) throw new IllegalArgumentException("Score dimensions differ from covariance");
        for (double value : scores) if (!Double.isFinite(value)) throw new ArithmeticException("Non-finite set score");
    }
    private double[] spectrum(double rho) {
        if (rho == 0) return lambda;
        double b = 0; for (double x : burden) b += x*x;
        if (rho == 1) return b > 1e-28*Arrays.stream(lambda).sum() ? new double[]{b} : new double[0];
        double[][] matrix = new double[lambda.length][lambda.length];
        double scale = Arrays.stream(lambda).max().orElseThrow();
        for (int i = 0; i < lambda.length; i++) for (int j = 0; j < lambda.length; j++)
            matrix[i][j] = ((i == j ? (1-rho)*lambda[i] : 0)+rho*burden[i]*burden[j])/scale;
        return Arrays.stream(new EigenvalueDecomposition(new Matrix(matrix)).getRealEigenvalues())
            .filter(x -> x > 0).map(x -> x*scale).toArray();
    }
}
