/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QMissingValuePolicy;
import gov.nih.eqtl.io.QDelimitedMatrixSource;
import gov.nih.eqtl.settest.QSetTestPolicy.FailurePolicy;

class QSetTestRunnerTest {
    @Test
    void boundedBurdenOutputAndResumeAreDeterministic(@TempDir Path directory)
        throws Exception {
        Path variants = directory.resolve("variants.csv");
        Files.writeString(variants, ",S1,S2,S3,S4,S5,S6,S7,S8\n"
            + "1:100:A:G,0,0,0,1,0,0,1,0\n"
            + "1:110:C:T,0,0,1,0,0,0,,0\n"
            + "1:120:G:A,0,1,0,1,2,1,2,1\n"
            + "1:130:A:C,0,0,0,0,0,0,0,0\n");
        Path traits = directory.resolve("traits.csv");
        Files.writeString(traits, ",S1,S2,S3,S4,S5,S6,S7,S8\n"
            + "trait1,1.2,1.0,1.4,2.1,1.8,2.3,2.8,2.5\n"
            + "trait2,3.0,2.8,2.7,2.4,2.5,2.1,1.9,1.8\n");
        Path sets = directory.resolve("sets.tsv");
        Files.writeString(sets, "SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\tWEIGHT\n"
            + "geneA\t1:100:A:G\tA\tG\tG\t1\n"
            + "geneA\t1:110:C:T\tC\tT\tC\t0.5\n"
            + "geneB\t1:110:C:T\tC\tT\tC\t0.5\n"
            + "geneB\t1:120:G:A\tG\tA\tA\t2\n"
            + "empty\t1:999:A:G\tA\tG\tG\t1\n"
            + "mono\t1:130:A:C\tA\tC\tC\t1\n");
        double[][] design = design();
        QSetTestPolicy policy = new QSetTestPolicy(0, 0.2, 0,
            Double.POSITIVE_INFINITY, QMissingValuePolicy.MEAN,
            FailurePolicy.SKIP, FailurePolicy.SKIP);
        Path output = directory.resolve("burden.csv");
        Path audit = directory.resolve("burden.sets.tsv");
        Path checkpoint = directory.resolve("checkpoint");
        QSetTestRunner.Options first = new QSetTestRunner.Options(QSetTestMethod.BURDEN,
            sets, output, audit, checkpoint, 1, policy, false, true, "none", 0,
            1, new double[] {0, 0.5, 1}, 199, 17);
        QSetTestRunner.run(new QDelimitedMatrixSource(variants, cCommonDelimiter, "#"),
            identity(8), new QDelimitedMatrixSource(traits, cCommonDelimiter, "#"),
            identity(8), design, first);
        byte[] expected = Files.readAllBytes(output);
        List<String> lines = Files.readAllLines(output);
        assertEquals(5, lines.size());
        assertArrayEquals(new String[] {"geneA/trait1", "geneA/trait2",
            "geneB/trait1", "geneB/trait2"}, lines.subList(1, lines.size()).stream()
                .map(line -> {
                    String[] fields = line.split(",", -1);
                    return fields[0] + "/" + fields[1];
                }).toArray(String[]::new));
        assertTrue(Files.readString(audit).contains("empty\t1\t1\t0\t0\tskipped-empty"));
        assertTrue(Files.readString(audit).contains("mono\t1\t0\t0\t1\tskipped-monomorphic"));

        QSetTestRunner.Options resumed = new QSetTestRunner.Options(QSetTestMethod.BURDEN,
            sets, output, audit, checkpoint, 1, policy, true, true, "none", 0,
            1, new double[] {0, 0.5, 1}, 199, 17);
        QSetTestRunner.run(new QDelimitedMatrixSource(variants, cCommonDelimiter, "#"),
            identity(8), new QDelimitedMatrixSource(traits, cCommonDelimiter, "#"),
            identity(8), design(), resumed);
        assertArrayEquals(expected, Files.readAllBytes(output));
    }

    @Test
    void productionSkatAndSkatOUseTheSameBoundedSchema(@TempDir Path directory)
        throws Exception {
        Path variants = directory.resolve("variants.csv");
        Files.writeString(variants, ",S1,S2,S3,S4,S5,S6,S7,S8\n"
            + "1:100:A:G,0,0,0,1,0,0,1,0\n"
            + "1:110:C:T,0,0,1,0,0,0,,0\n");
        Path traits = directory.resolve("traits.csv");
        Files.writeString(traits, ",S1,S2,S3,S4,S5,S6,S7,S8\n"
            + "trait1,1.2,1.0,1.4,2.1,1.8,2.3,2.8,2.5\n");
        Path sets = directory.resolve("sets.tsv");
        Files.writeString(sets, "SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\tWEIGHT\n"
            + "geneA\t1:100:A:G\tA\tG\tG\t1\n"
            + "geneA\t1:110:C:T\tC\tT\tC\t0.5\n");
        QSetTestPolicy policy = QSetTestPolicy.unfiltered(QMissingValuePolicy.MEAN);
        for (QSetTestMethod method : List.of(QSetTestMethod.SKAT, QSetTestMethod.SKAT_O)) {
            Path output = directory.resolve(method.optionName() + ".csv");
            QSetTestRunner.Options options = new QSetTestRunner.Options(method, sets, output,
                directory.resolve(method.optionName() + ".audit.tsv"),
                directory.resolve(method.optionName() + ".checkpoint"), 1, policy,
                false, false, "none", 0, 1, new double[] {0, 0.5, 1}, 99, 19);
            QSetTestRunner.run(new QDelimitedMatrixSource(variants, cCommonDelimiter, "#"),
                identity(8), new QDelimitedMatrixSource(traits, cCommonDelimiter, "#"),
                identity(8), design(), options);
            List<String> lines = Files.readAllLines(output);
            assertEquals(2, lines.size());
            assertEquals(method.optionName(), lines.get(1).split(",", -1)[2]);
            assertTrue(Double.parseDouble(lines.get(1).split(",", -1)[10]) > 0);
            if (method == QSetTestMethod.SKAT_O)
                assertTrue(lines.get(1).contains("correlated-parametric-null"));
        }
    }

    private static double[][] design() {
        double[] age = {20, 22, 25, 28, 30, 35, 38, 42};
        double[][] design = new double[age.length][2];
        for (int i = 0; i < age.length; i++) {
            design[i][0] = 1;
            design[i][1] = age[i];
        }
        return design;
    }

    private static int[] identity(int size) {
        int[] order = new int[size];
        for (int i = 0; i < size; i++) order[i] = i;
        return order;
    }
}
