/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.util.Arrays;

/**
 * Numerical inversion of the same weighted chi-square law, on a tilted
 * Bromwich contour. Tilting is a change of integration contour, not a
 * saddlepoint approximation. Extreme tails are integrated directly.
 */
final class QTiltedQuadraticForm {
    private static final int ORDER = 12;
    private static final int MAXIMUM_EVALUATIONS = 262144;
    private static final double RELATIVE_TOLERANCE = 2e-11;
    private static final double[] X = {.9914553711208126,.9491079123427585,
        .8648644233597691,.7415311855993944,.5860872354676911,
        .4058451513773972,.2077849550078985};
    private static final double[] WK = {.022935322010529225,.06309209262997855,
        .10479001032225018,.14065325971552592,.1690047266392679,
        .1903505780647854,.20443294007529889,.20948214108472783};
    private static final double[] WG = {.1294849661688697,.27970539148927667,
        .38183005050511894,.4179591836734694};

    record Result(double probability, double estimatedRelativeError) { }
    private QTiltedQuadraticForm() { }

    static Result probability(double q, double[] lambda, boolean lower) {
        if (q == 0) return new Result(lower ? 0 : 1, 0);
        if (!Double.isFinite(q)) return new Result(lower ? 1 : 0, 0);
        boolean upperContour = q >= Arrays.stream(lambda).sum();
        double left = -lambda.length / 2.0, right = Math.nextDown(q / 2);
        for (int i = 0; i < 100; i++) {
            double middle = left + (right-left)/2, derivative = 0;
            for (double value : lambda) derivative += value / (q-2*value*middle);
            if (derivative < 1) left = middle; else right = middle;
        }
        double tau = left + (right-left)/2;
        // Keep the contour away from the pole at zero, including q=E[Q].
        tau = upperContour ? Math.max(.25,tau) : Math.min(-.25,tau);
        double[] b = new double[lambda.length];
        double logScale = -tau;
        for (int i = 0; i < b.length; i++) {
            double denominator = q-2*lambda[i]*tau;
            if (!(denominator > 0)) throw new ArithmeticException("Invalid quadratic-form contour");
            b[i] = 2*lambda[i]/denominator;
            double ratio = -2*lambda[i]*tau/q;
            logScale -= .5*(Math.abs(ratio) < .5 ? Math.log1p(ratio)
                : Math.log(denominator)-Math.log(q));
        }
        // Chernoff's bound on the selected tail is exp(logScale).
        if (logScale < Math.log(Double.MIN_VALUE))
            return new Result(upperContour == lower ? 1 : 0, 0);
        Integrand f = new Integrand(tau,b);
        double tolerance = 2e-13/(1+Math.abs(tau));
        for (int attempt = 0; attempt < 4; attempt++, tolerance /= 8) {
            double end = 8, remainder;
            while ((remainder = remainder(end,b)) > tolerance/8) {
                end *= 2;
                if (end > 32768) throw new ArithmeticException("Quadratic-form inversion truncation limit exceeded");
            }
            int panels = (int)Math.ceil(end/(Math.PI/2));
            double integral = 0, error = remainder;
            for (int panel = 0; panel < panels; panel++) {
                double[] part = integrate(f,end*panel/panels,end*(panel+1)/panels,
                    tolerance/panels,20);
                integral += part[0]; error += part[1];
            }
            integral += tailCorrection(end,tau,b);
            double positive = upperContour ? integral : -integral;
            if (!(positive > 0) || !Double.isFinite(positive))
                throw new ArithmeticException("Quadratic-form inversion returned a nonpositive integral");
            if (error > RELATIVE_TOLERANCE*positive) continue;
            double small = Math.exp(logScale+Math.log(positive/Math.PI));
            if (!Double.isFinite(small) || small > 1+RELATIVE_TOLERANCE)
                throw new ArithmeticException("Quadratic-form inversion failed probability bounds");
            small = Math.min(1,small);
            double relative = error/positive;
            if (upperContour == lower) relative *= small/Math.max(Double.MIN_VALUE,1-small);
            if (relative > RELATIVE_TOLERANCE) continue;
            return new Result(upperContour == lower ? 1-small : small,relative);
        }
        throw new ArithmeticException("Quadratic-form inversion failed its relative error check");
    }

    private static final class Integrand {
        final double tau;
        final double[] b;
        int evaluations;
        double roundingScale;
        Integrand(double tau,double[] b) { this.tau=tau; this.b=b; }
        double value(double u) {
            if (++evaluations > MAXIMUM_EVALUATIONS)
                throw new ArithmeticException("Quadratic-form inversion exceeded maximum evaluations");
            double angle=Math.atan2(u,tau);
            double logAmplitude=0,phase=-u-angle,absolutePhase=0,amplitudeCorrection=0,phaseCorrection=0;
            for (double value : b) {
                double x = value*u;
                double term=-.25*Math.log1p(x*x)-amplitudeCorrection;
                double next=logAmplitude+term;amplitudeCorrection=(next-logAmplitude)-term;logAmplitude=next;
                double addition=.5*Math.atan(x);absolutePhase+=Math.abs(addition);
                term=addition-phaseCorrection;next=phase+term;phaseCorrection=(next-phase)-term;phase=next;
            }
            double envelope=Math.exp(logAmplitude)/Math.hypot(tau,u);
            roundingScale=envelope*(1+Math.abs(u)+Math.abs(angle)+absolutePhase+Math.abs(logAmplitude));
            return envelope*Math.cos(phase);
        }
    }

    private static double[] integrate(Integrand f,double a,double b,double tolerance,int depth) {
        double middle=(a+b)/2, half=(b-a)/2, center=f.value(middle);
        double kronrod=WK[7]*center,gauss=WG[3]*center,rounding=WK[7]*f.roundingScale;
        for (int j=0;j<7;j++) {
            double first=f.value(middle-half*X[j]),firstRounding=f.roundingScale;
            double second=f.value(middle+half*X[j]);
            kronrod += WK[j]*(first+second);rounding += WK[j]*(firstRounding+f.roundingScale);
            if (j%2==1) gauss += WG[(j-1)/2]*(first+second);
        }
        double quadratureError=Math.abs(kronrod-gauss)*half;
        // Use the complex amplitude/phase scale, not |Re(f)|: at a zero of
        // cosine, phase rounding can dominate an arbitrarily small real part.
        double roundoff=8*Math.ulp(1.0)*rounding*half;
        double error=Math.max(quadratureError,roundoff);
        // Refining a panel below its rounding floor cannot improve it. Return
        // that floor to the global relative-error check instead of recursing.
        if (quadratureError <= Math.max(tolerance,roundoff)) return new double[]{kronrod*half,error};
        if (depth==0) throw new ArithmeticException("Quadratic-form inversion quadrature did not converge");
        double[] first=integrate(f,a,middle,tolerance/2,depth-1);
        double[] second=integrate(f,middle,b,tolerance/2,depth-1);
        return new double[]{first[0]+second[0],first[1]+second[1]};
    }

    /**
     * For H(u)=prod(1-i*b_j*u)^(-1/2)/(tau+i*u),
     * |H^(m)(u)| <= (rank/2+1)_m * |phi(u)| / u^(m+1).
     * After m integrations by parts the remaining integral is therefore at
     * most (rank/2+1)_m * |phi(U)| / (m*U^m). This is a truncation bound;
     * Gauss-Kronrod and FP64 rounding errors are estimated separately.
     *
     * A second bound retains r_j=b_j*U/hypot(1,b_j*U). For u>=U, the
     * absolute k-th logarithmic derivative of H is bounded by
     * (k-1)!*(1+sum(r_j^k)/2)/U^k. Exponentiating these positive Taylor
     * coefficients bounds |H^(m)(u)| by m!*c_m*|phi(u)|/(U^m*u).
     * Also -d(log|phi(u)|)/d(log u)>=nu=sum(r_j^2)/2, so integrating
     * |phi(u)|/u gives at most |phi(U)|/nu. Both bounds apply to the
     * same unmodified spectrum; taking their minimum is safe and avoids
     * treating every weak eigenvalue as a full 1/u derivative at U.
     */
    private static double remainder(double end,double[] b) {
        double log = -Math.log(ORDER)-ORDER*Math.log(end);
        for (int k=0;k<ORDER;k++) log += Math.log(b.length/2.0+1+k);
        for (double value:b) log -= .25*Math.log1p(value*end*value*end);
        double rankBound = Math.exp(log);
        // A nominal rank can be much larger than the number of eigenvalues
        // that affect this endpoint. Bound the logarithmic derivatives using
        // their actual denominators, rather than replacing every one by 1/u.
        double[] logarithm = new double[ORDER+1];
        Arrays.fill(logarithm,1); // the pole 1/(tau+i*u) is bounded by 1/u
        double logAmplitude=0, decay=0;
        for (double value:b) {
            double x=value*end, ratio=Math.abs(x)/Math.hypot(1,x), power=ratio;
            logAmplitude -= .25*Math.log1p(x*x);
            decay += .5*ratio*ratio;
            for (int k=1;k<=ORDER;k++,power*=ratio) logarithm[k] += .5*power;
        }
        if (!(decay > 0)) return rankBound;
        double[] coefficients = new double[ORDER+1];
        coefficients[0]=1;
        double logFactorial=0;
        for (int k=1;k<=ORDER;k++) {
            for (int j=1;j<=k;j++) coefficients[k] += logarithm[j]*coefficients[k-j]/k;
            logFactorial += Math.log(k);
        }
        double spectrumBound=Math.exp(logAmplitude+logFactorial+Math.log(coefficients[ORDER])
            -ORDER*Math.log(end)-Math.log(decay));
        return Math.min(rankBound,spectrumBound);
    }

    /** Endpoint terms of integration by parts, using Taylor coefficients of H. */
    private static double tailCorrection(double u,double tau,double[] b) {
        double[][] logarithm=new double[ORDER][2];
        // -log(tau+i*u), then -.5*sum(log(1-i*b_j*u)).
        addLog(logarithm,u/(tau*tau+u*u),tau/(tau*tau+u*u),1);
        double logAmplitude=0, phase=0;
        for (double value:b) {
            double x=value*u;
            logAmplitude -= .25*Math.log1p(x*x); phase += .5*Math.atan(x);
            addLog(logarithm,value*x/(1+x*x),-value/(1+x*x),.5);
        }
        double amplitude=Math.exp(logAmplitude),denominator=tau*tau+u*u;
        double[][] coefficients=new double[ORDER][2];
        coefficients[0][0]=amplitude*(tau*Math.cos(phase)+u*Math.sin(phase))/denominator;
        coefficients[0][1]=amplitude*(tau*Math.sin(phase)-u*Math.cos(phase))/denominator;
        for (int k=1;k<ORDER;k++) for (int j=1;j<=k;j++) {
            double[] a=logarithm[j],c=coefficients[k-j];
            coefficients[k][0] += j*(a[0]*c[0]-a[1]*c[1])/k;
            coefficients[k][1] += j*(a[0]*c[1]+a[1]*c[0])/k;
        }
        double real=0,imaginary=0,factorial=1,fr=0,fi=-1;
        for (int k=0;k<ORDER;k++) {
            if (k>0) factorial *= k;
            real += factorial*(fr*coefficients[k][0]-fi*coefficients[k][1]);
            imaginary += factorial*(fr*coefficients[k][1]+fi*coefficients[k][0]);
            double next=fi; fi=-fr; fr=next;
        }
        return Math.cos(u)*real+Math.sin(u)*imaginary;
    }

    private static void addLog(double[][] coefficients,double real,double imaginary,double weight) {
        double pr=1,pi=0;
        for (int k=1;k<coefficients.length;k++) {
            double next=pr*real-pi*imaginary;
            pi=pr*imaginary+pi*real; pr=next;
            double sign=k%2==0?1:-1;
            coefficients[k][0] += weight*sign*pr/k;
            coefficients[k][1] += weight*sign*pi/k;
        }
    }
}
