/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import java.nio.file.*;
import java.util.*;
import gov.nih.eqtl.io.QCovariateTable;

/** Matched baseline/candidate benchmark for independent repeated-subject blocks. */
public final class QCohortAuditBenchmark {
    private static volatile double sink;
    public static void main(String[] args) throws Exception {
        Path root=Path.of("target/cohort-benchmark");Files.createDirectories(root);
        int n=512;String[] ids=new String[n];boolean[] retained=new boolean[n];Arrays.fill(retained,true);
        StringBuilder table=new StringBuilder("id,age,subject,cohort\n");
        for(int i=0;i<n;i++){ids[i]="s"+i;table.append(ids[i]+","+(i%53)+",p"+(i/2)+",A\n");}
        Files.writeString(root.resolve("cov.csv"),table);
        Files.writeString(root.resolve("model.tsv"),"COHORT\tFIXED_COVARIATES\tSUBJECT_COLUMN\tREPEATED_POLICY\tWITHIN_SUBJECT_CORRELATION\nA\tage\tsubject\tcompound-symmetry\t0.25\n");
        var definitions=QCohortModel.Definitions.load(root.resolve("model.tsv"));
        var covariates=QCovariateTable.load(root.resolve("cov.csv"),new char[]{','},"#");
        time("compound-symmetry-prepare-512",()->sink+=QCohortModel.create(definitions,covariates,"cohort",retained,ids,null,null).design().length);
        var model=QCohortModel.create(definitions,covariates,"cohort",retained,ids,null,null);
        double[][] q=new gov.nih.jama.QRDecomposition(model.design()).getQ().getArray();
        double[][] rows=new double[16][n];Random random=new Random(812);
        for(double[] row:rows)for(int i=0;i<n;i++)row[i]=random.nextGaussian();
        time("compound-symmetry-project-16x512",()->{
            double[] row=model.wrap(null).residualize(rows,q,"benchmark")[0];
            sink+=QCohortModel.stackedDot(row,row);
        });
        System.out.println("checksum="+sink);
    }
    private interface Job {void run() throws Exception;}
    private static void time(String name,Job job) throws Exception {
        for(int i=0;i<3;i++)job.run();long[] times=new long[7];
        for(int i=0;i<times.length;i++){long t=System.nanoTime();job.run();times[i]=System.nanoTime()-t;}
        Arrays.sort(times);System.out.printf(Locale.ROOT,"AUDIT %s median_ms=%.6f min_ms=%.6f max_ms=%.6f%n",name,times[3]/1e6,times[0]/1e6,times[6]/1e6);
    }
}
