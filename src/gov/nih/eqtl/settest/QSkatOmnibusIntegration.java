/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.util.*;
import java.util.function.DoubleUnaryOperator;
import jdistlib.ChiSquare;
import jdistlib.Normal;

/**
 * Joint Gaussian quadratic-form tail, integrating radius analytically.
 * Z=R D, R^2~chi-square(rank), D uniform on the sphere, independently.
 * The union across rho is R^2 >= min_rho critical[rho]/q_rho(D).
 * All tails are evaluated directly (no subtraction from one or event-count floor).
 */
final class QSkatOmnibusIntegration {
    record Result(double probability, double error, int evaluations, String method) { }
    private QSkatOmnibusIntegration() { }

    static Result integrate(double[] lambda, double[] burden, double[] rho, double[] critical,
        double minimumP, double relativeTolerance, int maximumEvaluations, long seed) {
        if (!(minimumP > 0) || minimumP < 1e-280)
            throw new ArithmeticException("SKAT-O tail is outside the supported integration range (p >= 1e-280)");
        if (minimumP == 1) return new Result(1,0,0,"joint-radial-exact");
        if (lambda.length == 1) return new Result(minimumP,0,0,"joint-radial-rank-one-exact");
        double smallest = Arrays.stream(lambda).min().orElseThrow(), largest = Arrays.stream(lambda).max().orElseThrow();
        if (largest-smallest <= 1e-13*largest)
            return isotropic(lambda,burden,rho,critical,minimumP,relativeTolerance,maximumEvaluations);
        if (lambda.length == 2) {
            int[] count = {0};
            DoubleUnaryOperator f = angle -> {
                if (++count[0] > maximumEvaluations)
                    throw new ArithmeticException("SKAT-O angular integration exceeded maximum evaluations");
                return radial(new double[]{Math.cos(angle),Math.sin(angle)},lambda,burden,rho,critical)/minimumP;
            };
            double value = 0, error = 0;
            // Seed multiple panels so a narrow peak cannot hide inside one coarse panel.
            for (int i = 0; i < 64; i++) {
                double a = Math.PI*i/64, b = Math.PI*(i+1)/64, m = (a+b)/2;
                double fa = f.applyAsDouble(a), fm = f.applyAsDouble(m), fb = f.applyAsDouble(b);
                double[] integral = adaptive(f,a,b,fa,fm,fb,(b-a)*(fa+4*fm+fb)/6,
                    relativeTolerance*Math.PI/256,28);
                value += integral[0]; error += integral[1];
            }
            return checked(value*minimumP/Math.PI,error*minimumP/Math.PI,count[0],
                "joint-radial-adaptive-angular",minimumP,rho.length,relativeTolerance);
        }
        // Independently shifted Halton rules integrate directions, not rare-event indicators.
        // Replicate uncertainty is an estimate, not an interval-arithmetic guarantee.
        final int replicates = 16;
        int rank = lambda.length;
        int[] primes = primes(rank);
        Random random = new Random(seed);
        double[][] shifts = new double[replicates][rank];
        for (double[] shift : shifts) for (int j = 0; j < rank; j++) shift[j] = random.nextDouble();
        double[] sums = new double[replicates], direction = new double[rank];
        double previous = Double.NaN;
        int completed = 0, successive = 0;
        for (int perReplicate = 256; (long)perReplicate*replicates <= maximumEvaluations; perReplicate *= 2) {
            for (int r = 0; r < replicates; r++) for (int index = completed+1; index <= perReplicate; index++) {
                double norm = 0;
                for (int j = 0; j < rank; j++) {
                    double u = radicalInverse(index,primes[j])+shifts[r][j]; u -= Math.floor(u);
                    direction[j] = Normal.quantile(Math.max(1e-15,Math.min(1-1e-15,u)),0,1,true,false);
                    norm += direction[j]*direction[j];
                }
                norm = Math.sqrt(norm);
                for (int j = 0; j < rank; j++) direction[j] /= norm;
                sums[r] += radial(direction,lambda,burden,rho,critical)/minimumP;
            }
            double mean = Arrays.stream(sums).sum()/replicates/perReplicate, variance = 0;
            for (double sum : sums) variance += Math.pow(sum/perReplicate-mean,2);
            double error = 4*Math.sqrt(variance/(replicates*(replicates-1)));
            boolean stable = error <= relativeTolerance*mean
                && Double.isFinite(previous) && Math.abs(mean-previous) <= relativeTolerance*mean;
            successive = stable ? successive+1 : 0;
            if (successive >= 2)
                return checked(mean*minimumP,error*minimumP,perReplicate*replicates,
                    "joint-radial-randomized-quadrature",minimumP,rho.length,relativeTolerance);
            previous = mean; completed = perReplicate;
            if (perReplicate > maximumEvaluations/replicates/2) break;
        }
        throw new ArithmeticException("SKAT-O radial quadrature did not converge within " + maximumEvaluations
            + " evaluations at relative tolerance " + relativeTolerance + "; no p-value substituted");
    }

    /** Equal nonzero eigenvalues: condition on the burden direction; remaining radius is chi-square(rank-1). */
    private static Result isotropic(double[] lambda,double[] burden,double[] rho,double[] critical,
        double minimum,double tolerance,int limit) {
        double scale = Arrays.stream(lambda).average().orElseThrow(), b2 = 0;
        for (double b : burden) b2 += b*b;
        double[] coefficient = new double[rho.length]; double bound = Double.POSITIVE_INFINITY;
        for (int i = 0; i < rho.length; i++) {
            coefficient[i] = (1-rho[i])*scale+rho[i]*b2;
            if (coefficient[i] > 0) bound = Math.min(bound,Math.sqrt(critical[i]/coefficient[i]));
        }
        final int[] count = {0};
        DoubleUnaryOperator f = z -> {
            if (++count[0] > limit) throw new ArithmeticException("SKAT-O conditional integration exceeded maximum evaluations");
            double threshold = Double.POSITIVE_INFINITY;
            for (int i = 0; i < rho.length; i++) if (rho[i] < 1)
                threshold = Math.min(threshold,(critical[i]-coefficient[i]*z*z)/((1-rho[i])*scale));
            return 2*Normal.density(z,0,1,false)*ChiSquare.cumulative(Math.max(0,threshold),lambda.length-1,false,false)/minimum;
        };
        double value = 0, error = 0;
        for (int i = 0; i < 64; i++) {
            double a = bound*i/64, b = bound*(i+1)/64, m = (a+b)/2;
            double fa = f.applyAsDouble(a), fm = f.applyAsDouble(m), fb = f.applyAsDouble(b);
            double[] part = adaptive(f,a,b,fa,fm,fb,(b-a)*(fa+4*fm+fb)/6,tolerance/256,28);
            value += part[0]; error += part[1];
        }
        return checked(value*minimum+2*Normal.cumulative(bound,0,1,false,false),error*minimum,count[0],
            "joint-conditional-isotropic-integral",minimum,rho.length,tolerance);
    }

    private static Result checked(double p,double error,int evaluations,String method,
        double minimum,int components,double tolerance) {
        if (!Double.isFinite(p) || !(p > 0) || error > tolerance*p
            || p < minimum*(1-tolerance) || p > Math.min(1,components*minimum)*(1+tolerance))
            throw new ArithmeticException("SKAT-O integration failed error/union probability bounds");
        return new Result(Math.max(minimum,Math.min(1,p)),error,evaluations,method);
    }
    private static double radial(double[] direction,double[] lambda,double[] burden,double[] rho,double[] critical) {
        double skat = 0, b = 0;
        for (int j = 0; j < direction.length; j++) {
            skat += lambda[j]*direction[j]*direction[j]; b += burden[j]*direction[j];
        }
        double boundary = Double.POSITIVE_INFINITY;
        for (int j = 0; j < rho.length; j++) {
            double q = (1-rho[j])*skat+rho[j]*b*b;
            if (q > 0) boundary = Math.min(boundary,critical[j]/q);
        }
        return ChiSquare.cumulative(boundary,direction.length,false,false);
    }
    private static double[] adaptive(DoubleUnaryOperator f,double a,double b,double fa,double fm,double fb,
        double whole,double tolerance,int depth) {
        double m = (a+b)/2, fl = f.applyAsDouble((a+m)/2), fr = f.applyAsDouble((m+b)/2);
        double left = (m-a)*(fa+4*fl+fm)/6, right = (b-m)*(fm+4*fr+fb)/6;
        double delta = left+right-whole;
        if (Math.abs(delta) <= 15*tolerance) return new double[]{left+right+delta/15,Math.abs(delta)/15};
        if (depth == 0) throw new ArithmeticException("SKAT-O angular quadrature did not converge");
        double[] x = adaptive(f,a,m,fa,fl,fm,left,tolerance/2,depth-1);
        double[] y = adaptive(f,m,b,fm,fr,fb,right,tolerance/2,depth-1);
        return new double[]{x[0]+y[0],x[1]+y[1]};
    }
    private static int[] primes(int count) {
        int[] result = new int[count]; int found = 0;
        for (int x = 2; found < count; x++) {
            boolean prime = true;
            for (int y = 2; y <= x/y; y++) if (x%y == 0) { prime = false; break; }
            if (prime) result[found++] = x;
        }
        return result;
    }
    private static double radicalInverse(int index,int base) {
        double value = 0, place = 1.0/base;
        while (index > 0) { value += (index%base)*place; index /= base; place /= base; }
        return value;
    }
}
