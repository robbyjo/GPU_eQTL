/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import java.util.*;
import gov.nih.eqtl.QMixedModel;
import gov.nih.jama.QRDecomposition;
import static gov.nih.utils.QDataUtils.kUndefinedValue;

/** One response and exact ordered observed-sample mask; no cross-trait fitted-null reuse. */
public final class QSetTraitModel {
    public record Scores(double[] u,double[][] covariance) { }
    private final int[] samples;
    private final double[][] q, projection;
    private final double[] response;
    private final QMixedModel.NullState mixed;
    private final double sigma;
    private final int nullDf;

    public QSetTraitModel(double[] y,double[][] design,QMixedModel fullModel) {
        samples = observed(y);
        double[][] x = subsetDesign(design,samples);
        QRDecomposition qr = new QRDecomposition(x);
        if (!qr.isFullRank() || samples.length <= qr.getRank()+1)
            throw new IllegalArgumentException("Set trait mask has rank-deficient design or non-positive association DF");
        q = qr.getQ().getArray(); nullDf = samples.length-qr.getRank();
        double[] values = new double[samples.length];
        for (int i = 0; i < samples.length; i++) values[i] = y[samples[i]];
        if (fullModel != null) {
            mixed = fullModel.subset(samples).fitNull(values);
            response = mixed.projectedResponse(); sigma = 1;
            double[] p = mixed.projection(); projection = new double[samples.length][samples.length];
            for (int i = 0; i < samples.length; i++) System.arraycopy(p,i*samples.length,projection[i],0,samples.length);
        } else {
            mixed = null; projection = null; response = values.clone();
            for (int a = 0; a < q[0].length; a++) {
                double dot = 0; for (int i = 0; i < values.length; i++) dot += q[i][a]*values[i];
                for (int i = 0; i < values.length; i++) response[i] -= q[i][a]*dot;
            }
            double norm = 0, original = 0;
            for (int i = 0; i < values.length; i++) { norm = Math.hypot(norm,response[i]); original = Math.hypot(original,values[i]); }
            if (!(norm > 1e-12*original)) throw new IllegalArgumentException("Set trait has no residual variation");
            sigma = norm/Math.sqrt(nullDf);
            for (int i = 0; i < response.length; i++) response[i] /= sigma;
        }
    }
    public Scores scores(double[][] weightedRows,QSetMatrixOperations operations,boolean burden) {
        double[][] rows = weightedRows;
        if (burden) {
            double[] ones = new double[rows.length]; Arrays.fill(ones,1);
            rows = operations.cross(new double[][]{ones},QSetMatrixOperations.transpose(rows),"burden-aggregation");
        }
        if (mixed != null) {
            double[][] centered = new double[rows.length][];
            for (int i = 0; i < rows.length; i++) centered[i] = mixed.center(rows[i]);
            rows = centered;
        } else rows = operations.residualize(rows,q);
        double[][] projected = mixed == null ? rows : operations.cross(rows,projection,"fitted-P-times-genotypes");
        double[][] products = operations.cross(rows,new double[][]{response},"genotype-scores");
        double[] u = new double[rows.length];
        for (int i = 0; i < u.length; i++) u[i] = products[i][0];
        double[][] covariance = operations.cross(rows,projected,"score-covariance");
        for (int i = 0; i < covariance.length; i++) for (int j = 0; j < i; j++)
            covariance[i][j] = covariance[j][i] = .5*covariance[i][j]+.5*covariance[j][i];
        return new Scores(u,covariance);
    }
    public int[] samples() { return samples.clone(); }
    public int sampleCount() { return samples.length; }
    public int nullDf() { return nullDf; }
    public boolean related() { return mixed != null; }
    public double sigma() { return sigma; }
    public String diagnostics() {
        if (mixed == null) return "fixed;null-df="+nullDf;
        StringJoiner variances = new StringJoiner("|");
        double[] values = mixed.fit().varianceComponents();
        for (int i = 0; i < values.length; i++) variances.add(mixed.fit().componentNames().get(i)+"="+values[i]);
        return "null-reml-score;variance="+variances
            +";logLik="+mixed.fit().restrictedLogLikelihood()+";iterations="+mixed.fit().iterations()
            +";converged="+mixed.fit().converged()+";singular="+mixed.fit().isSingular(1e-8);
    }
    public static int[] observed(double[] y) {
        return java.util.stream.IntStream.range(0,y.length).filter(i -> {
            if (Double.isNaN(y[i]) || y[i] == kUndefinedValue) return false;
            if (!Double.isFinite(y[i])) throw new IllegalArgumentException("Infinite set trait value");
            return true;
        }).toArray();
    }
    public static double[][] subsetDesign(double[][] design,int[] samples) {
        if (samples.length == 0) throw new IllegalArgumentException("Set trait has no observed samples");
        double[][] x = new double[samples.length][];
        for (int i = 0; i < samples.length; i++) x[i] = design[samples[i]].clone();
        return x;
    }
}
