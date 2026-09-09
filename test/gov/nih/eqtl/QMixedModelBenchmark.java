/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import java.util.*;

/** Manual timing of the two deliberately different mixed-model inference modes. */
public final class QMixedModelBenchmark {
    private static volatile double sink;
    public static void main(String[] args) throws Exception {
        var data = QMixedModelTest.data(); var model = QMixedModelTest.model("combined","residual");
        double[][] variants = new double[32][];
        for (int g = 0; g < variants.length; g++) {
            variants[g] = data.g().clone();
            for (int i = 0; i < data.g().length; i++) variants[g][i] += .1*Math.sin((g+1)*(i+1));
        }
        for (boolean refit : new boolean[]{false,true}) {
            long[] times = new long[5];
            for (int run = -2; run < times.length; run++) {
                long start = System.nanoTime(); var state = model.fitNull(data.y());
                for (double[] variant : variants) {
                    if (refit) sink += model.refit(data.y(),variant,state).log10P();
                    else { var g = state.center(variant); sink += state.score(g,state.project(g)).log10P(); }
                }
                if (run >= 0) times[run] = System.nanoTime()-start;
            }
            Arrays.sort(times);
            System.out.printf(Locale.ROOT,"MIXED_BENCH mode=%s n=48 predictors=32 null_fit_included=true median_ms=%.6f%n",
                refit?"per-predictor-reml":"shared-null-score",times[2]/1e6);
        }
        System.out.println("checksum="+sink);
    }
}
