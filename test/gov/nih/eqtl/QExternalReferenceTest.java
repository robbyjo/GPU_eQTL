/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.eqtl.io.QMatrixRowSource.Block;
import gov.nih.eqtl.settest.QSkatReference;
import gov.nih.eqtl.settest.*;

class QExternalReferenceTest {
    private static final Path ROOT = Path.of("test/resources/audit-reference");
    @Test
    void productionBurdenAndSkatSingletonMatchR(@TempDir Path temp) throws Exception {
        var data=Files.readAllLines(ROOT.resolve("data.csv"));int n=data.size()-1;
        double[][] design=new double[n][2];double[] g=new double[n],y=new double[n];
        for(int i=0;i<n;i++) {String[] row=data.get(i+1).split(",");design[i][0]=1;design[i][1]=Double.parseDouble(row[1]);g[i]=Double.parseDouble(row[2]);y[i]=Double.parseDouble(row[3]);}
        var nullModel=QSetTestNullModel.create(new String[]{"trait"},new double[][]{y},design);
        Path path=temp.resolve("sets.tsv");Files.writeString(path,"SET_ID\tVARIANT_ID\tREF\tALT\tEFFECT_ALLELE\tWEIGHT\nset\tv\tA\tG\tG\t1\n");
        var sets=QVariantSetTable.load(path);
        var variants=java.util.List.of(new QBurdenReference.Variant("v","A","G",g));
        var policy=QSetTestPolicy.unfiltered(QMissingValuePolicy.ERROR);
        var burden=QBurdenReference.analyzeBatched(sets,variants,nullModel,policy).results().get(0);
        String[] r=Files.readAllLines(ROOT.resolve("gls.csv")).get(1).split(",");
        assertEquals(Double.parseDouble(r[1]),burden.effect(),2e-12);
        assertEquals(Double.parseDouble(r[2]),burden.rSquared(),2e-12);
        assertEquals(Double.parseDouble(r[3]),burden.tStatistic(),2e-11);
        assertEquals(Double.parseDouble(r[5]),burden.log10P(),2e-11);
        var skat=QKernelSetReference.analyze(sets,variants,nullModel,policy,
            new QKernelSetReference.Options(QSetTestMethod.SKAT,null,0,0)).results().get(0);
        String[] sr=Files.readAllLines(ROOT.resolve("skat-singleton.csv")).get(1).split(",");
        assertEquals(Double.parseDouble(sr[0]),skat.statistic(),2e-11);
        assertEquals(Double.parseDouble(sr[2]),skat.pValue(),2e-12);
    }
    @Test
    void fixedRepeatedAndSuppliedPedigreeCovarianceMatchR(@TempDir Path temp) throws Exception {
        var data = Files.readAllLines(ROOT.resolve("data.csv"));
        int n = data.size()-1;
        String[] ids = new String[n]; boolean[] selected = new boolean[n]; Arrays.fill(selected,true);
        double[] g=new double[n], y=new double[n];
        StringBuilder cov = new StringBuilder("id,age,subject,cohort\n");
        for (int i=0;i<n;i++) {
            String[] row=data.get(i+1).split(","); ids[i]=row[0];
            g[i]=Double.parseDouble(row[2]); y[i]=Double.parseDouble(row[3]);
            cov.append(row[0]+","+row[1]+","+row[4]+",A\n");
        }
        Files.writeString(temp.resolve("cov.csv"),cov);
        for (String line : Files.readAllLines(ROOT.resolve("gls.csv")).subList(1,4)) {
            String[] expected=line.split(","); String name=expected[0];
            String header="COHORT\tFIXED_COVARIATES\tCOVARIANCE_FILE\tSUBJECT_COLUMN\tREPEATED_POLICY\tWITHIN_SUBJECT_CORRELATION\n";
            // Repeat reference also checks the O(n) compound-symmetry operator.
            String definition=name.equals("repeated") ? "A\tage\t\tsubject\tcompound-symmetry\t0.35\n"
                : "A\tage\t"+ROOT.resolve(name+".tsv").toAbsolutePath()+"\t\tindependent\t\n";
            Files.writeString(temp.resolve("model.tsv"),header+definition);
            QCohortModel model=QCohortModel.create(QCohortModel.Definitions.load(temp.resolve("model.tsv")),
                QCovariateTable.load(temp.resolve("cov.csv"),new char[]{','},"#"),"cohort",selected,ids,null,null);
            double[][] q=new gov.nih.jama.QRDecomposition(model.design()).getQ().getArray();
            var pg=QeQTLPreprocessor.prepare(new Block(0,new String[]{"g"},new double[][]{g.clone()}),q,"g",model.wrap(null));
            var py=QeQTLPreprocessor.prepare(new Block(0,new String[]{"y"},new double[][]{y.clone()}),q,"y",model.wrap(null));
            var actual=QeQTLStatistics.calculate(QeQTLPreprocessor.correlation(pg.values()[0],py.values()[0]),
                py.standardDeviations()[0],pg.standardDeviations()[0],Integer.parseInt(expected[4]),0);
            assertEquals(Double.parseDouble(expected[1]),actual.effect(),2e-12,name);
            assertEquals(Double.parseDouble(expected[2]),actual.rSquared(),2e-12,name);
            assertEquals(Double.parseDouble(expected[3]),actual.tStatistic(),2e-11,name);
            assertEquals(Double.parseDouble(expected[5]),actual.log10P(),2e-11,name);
        }
    }

    @Test
    void tTailsMatchR() throws Exception {
        for (String line:Files.readAllLines(ROOT.resolve("t-tail.csv")).subList(1,55)) {
            String[] row=line.split(","); double expected=Double.parseDouble(row[2]);
            assertEquals(expected,net.sourceforge.jdistlib.T.cumulative(Double.parseDouble(row[1]),
                Double.parseDouble(row[0]),false,true),Math.max(1e-12,Math.abs(expected)*2e-13),line);
        }
    }

    @Test
    void reportedConvergedMixturesMatchIndependentDavies() throws Exception {
        int exact=0, converged=0, fallback=0;
        double maxError=0, fallbackError=0;
        for (String line:Files.readAllLines(ROOT.resolve("mixtures.csv")).subList(1,31)) {
            String[] row=line.split(","); double[] lambda=Arrays.stream(row[0].split(";")).mapToDouble(Double::parseDouble).toArray();
            var result=QSkatReference.mixturePValue(Double.parseDouble(row[1]),lambda);
            double error=Math.abs(result.pValue()-Double.parseDouble(row[2]));
            assertFalse(result.method().contains("fallback"));
            if (result.method().contains("fallback")) { fallback++; fallbackError=Math.max(error,fallbackError); }
            else {
                if(result.method().contains("exact")) exact++; else converged++;
                maxError=Math.max(error,maxError);
                assertTrue(error<2e-7,"False convergence: "+line+" actual="+result);
            }
        }
        System.out.println("QUADRATIC_AUDIT exact="+exact+" converged="+converged+" fallback="+fallback
            +" max_absolute_error="+maxError+" fallback_max_absolute_error="+fallbackError);
    }
}
