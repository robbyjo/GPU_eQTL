/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.util.Arrays;
import java.util.Random;
import gov.nih.eqtl.io.QMatrixRowSource.Block;
import gov.nih.eqtl.settest.QSkatReference;
import gov.nih.eqtl.settest.*;
import java.nio.file.*;
import java.util.ArrayList;
import gov.nih.utils.matrix.QMatrixUtils;
import gov.nih.utils.matrix.EMultiplicationMode;

/** Manual median timings; identical harness can run against saved baseline classes. */
public final class QModuleAuditBenchmark {
    private static volatile double sink;
    public static void main(String[] args) throws Exception {
        Random random = new Random(20260909);
        int n = 512, m = 16;
        double[][] variants = new double[m][n], design = new double[n][3];
        double[] trait = new double[n], weights = new double[m];
        Arrays.fill(weights, 1);
        for (int j = 0; j < n; j++) {
            trait[j] = random.nextGaussian();
            design[j][0] = 1; design[j][1] = random.nextGaussian(); design[j][2] = random.nextGaussian();
            for (int i = 0; i < m; i++) variants[i][j] = random.nextGaussian();
        }
        double[][] q = new gov.nih.jama.QRDecomposition(design).getQ().getArray();
        String[] ids = new String[m]; Arrays.fill(ids, "benchmark");
        time("fixed-preparation-16x512-p3", () -> sink += QeQTLPreprocessor.prepare(
            new Block(0, ids, copy(variants)), q, "benchmark").standardDeviations()[0]);
        time("burden-score-product-16x512x16", () -> sink += QMatrixUtils.parallelMatrixMultiplication(
            variants, QMatrixUtils.transpose(variants), null, 1, m, m, EMultiplicationMode.XY)[0][0]);
        time("skat-16x512", () -> sink += QSkatReference.calculate(variants, weights, trait, n-3).pValue());
        time("skat-o-16x512-rho3-b99", () -> sink += QSkatReference.calculateOmnibus(
            variants, weights, trait, n-3, new double[]{0,0.5,1}, 99, 813).adjustedP());
        ArrayList<QBurdenReference.Variant> genotype = new ArrayList<>();
        StringBuilder definition = new StringBuilder("SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\tWEIGHT\n");
        for (int i=0;i<m;i++) {
            double[] dosage=new double[n]; for(int j=0;j<n;j++)dosage[j]=random.nextInt(3);
            genotype.add(new QBurdenReference.Variant("v"+i,"A","G",dosage));
            definition.append("set\tv"+i+"\tA\tG\tG\t1\n");
        }
        double[][] traits=new double[8][n];String[] traitIds=new String[8];
        for(int i=0;i<traits.length;i++){traitIds[i]="t"+i;for(int j=0;j<n;j++)traits[i][j]=random.nextGaussian();}
        var nullModel=QSetTestNullModel.create(traitIds,traits,design);
        Path path=Path.of("target/module-benchmark-sets.tsv");Files.writeString(path,definition);
        var sets=QVariantSetTable.load(path);
        var policy=QSetTestPolicy.unfiltered(QMissingValuePolicy.ERROR);
        time("burden-production-16variants-8traits-n512",()->sink+=QBurdenReference.analyzeBatched(sets,genotype,nullModel,policy).results().get(0).tStatistic());
        for(var method:new QSetTestMethod[]{QSetTestMethod.SKAT,QSetTestMethod.SKAT_O}) {
            var options=new QKernelSetReference.Options(method,new double[]{0,.5,1},99,813);
            time(method.optionName()+"-production-16variants-8traits-n512",()->sink+=QKernelSetReference.analyze(sets,genotype,nullModel,policy,options).results().get(0).pValue());
            for(var result:QKernelSetReference.analyze(sets,genotype,nullModel,policy,options).results())
                System.out.println("RESULT "+method.optionName()+" "+result.setId()+" "+result.traitId()+" "+result.statistic()+" "+result.pValue()+" "+result.log10P());
        }
        System.out.println("checksum=" + sink);
    }
    private static double[][] copy(double[][] rows) {
        return Arrays.stream(rows).map(double[]::clone).toArray(double[][]::new);
    }
    private static void time(String name, Runnable job) {
        for (int i=0; i<3; i++) job.run();
        long[] times = new long[7];
        for (int i=0; i<times.length; i++) { long start=System.nanoTime(); job.run(); times[i]=System.nanoTime()-start; }
        Arrays.sort(times);
        System.out.printf(java.util.Locale.ROOT, "AUDIT %s median_ms=%.6f min_ms=%.6f max_ms=%.6f%n", name,
            times[3]/1e6, times[0]/1e6, times[6]/1e6);
    }
}
