/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;

class QMixedModelTest {
    static final Path ROOT = Path.of("test/resources/mixed-reference");
    record Data(double[][] x, double[] g, double[] y, String[] ids) { }
    static Data data() throws Exception {
        var lines = Files.readAllLines(ROOT.resolve("data.csv")); int n = lines.size()-1;
        double[][] x = new double[n][2]; double[] g = new double[n], y = new double[n]; String[] ids = new String[n];
        for (int i = 0; i < n; i++) { String[] v = lines.get(i+1).split(",");
            ids[i] = v[0]; x[i][0] = 1; x[i][1] = Double.parseDouble(v[2]);
            g[i] = Double.parseDouble(v[3]); y[i] = Double.parseDouble(v[4]); }
        return new Data(x,g,y,ids);
    }
    static QMixedModel model(String type, String df) throws Exception {
        var config = new QeQTLAnalysisConfig(new LinkedHashMap<>());
        config.set("covariate_file",ROOT.resolve("data.csv").toString());
        config.set("reml_df",df);
        if (!type.equals("group")) config.set("pedigree_file",ROOT.resolve("pedigree.tsv").toString());
        if (!type.equals("pedigree")) config.set("covariate_random","group");
        Data d = data();
        return QMixedModel.create(config,d.x(),QCovariateTable.load(ROOT.resolve("data.csv"),new char[]{','},"#"),null,d.ids());
    }
    @Test void fittedPedigreeAndGroupingMatchIndependentROptimizer() throws Exception {
        Data d = data(); double max = 0;
        for (String row : Files.readAllLines(ROOT.resolve("association.csv")).subList(1,7)) {
            String[] v = row.split(","); QMixedModel model = model(v[0],"residual");
            var state = model.fitNull(d.y()); double[] g = state.center(d.g());
            var actual = v[1].equals("score") ? state.score(g,state.project(g)) : model.refit(d.y(),d.g(),state);
            double[] observed = {actual.beta(),actual.standardError(),actual.statistic(),actual.df(),actual.log10P(),actual.fit().restrictedLogLikelihood()};
            for (int i = 0; i < observed.length; i++) {
                double error = Math.abs(observed[i]-Double.parseDouble(v[i+2])); max = Math.max(max,error);
                assertEquals(Double.parseDouble(v[i+2]),observed[i],3e-6,row+" field="+i);
            }
            double[] variances = Arrays.stream(v[8].split(";")).mapToDouble(Double::parseDouble).toArray();
            assertArrayEquals(variances,actual.fit().varianceComponents(),3e-6,row);
            assertTrue(actual.fit().converged());
        }
        System.out.println("MIXED_R_REFERENCE six fits; max absolute statistic/log-likelihood error="+max);
    }
    @Test void pedigreeValidationAndRelationshipsAreExplicit(@TempDir Path root) throws Exception {
        var pedigree = QMixedModel.readPedigree(ROOT.resolve("pedigree.tsv"));
        assertEquals(.5,pedigree.relationship("s1","s3"),0);
        assertEquals(.5,pedigree.relationship("s3","s4"),0);
        assertEquals(0,pedigree.relationship("s1","s5"),0);
        for (String rows : List.of("a\t0\t0\na\t0\t0\n", "a\tb\t0\nb\ta\t0\n", "a\tmissing\t0\n", "a\t\t0\n")) {
            Path file = root.resolve("bad.tsv"); Files.writeString(file,"ID\tSIRE\tDAM\n"+rows);
            assertThrows(IllegalArgumentException.class,()->QMixedModel.readPedigree(file));
        }
        var rows = Files.readAllLines(ROOT.resolve("pedigree.tsv")); var reordered = new ArrayList<>(rows.subList(1,rows.size()));
        Collections.reverse(reordered); reordered.add(0,rows.get(0));
        Path file = root.resolve("reordered.tsv"); Files.write(file,reordered);
        assertEquals(.5,QMixedModel.readPedigree(file).relationship("s3","s4"),0);
    }
    @Test void nullScoreIsScaleInvariantAndDegeneracyIsReported() throws Exception {
        Data d = data(); var model = model("combined","residual"); var state = model.fitNull(d.y());
        double[] g = state.center(d.g()); var expected = state.score(g,state.project(g));
        double[] transformed = d.g().clone(); for (int i = 0; i < transformed.length; i++) transformed[i] = 10 + 3*transformed[i];
        g = state.center(transformed); var actual = state.score(g,state.project(g));
        assertEquals(expected.statistic(),actual.statistic(),1e-12);
        assertEquals(expected.beta()/3,actual.beta(),1e-12);
        double[] constant = new double[d.g().length]; Arrays.fill(constant,2);
        g = state.center(constant); assertEquals("monomorphic-or-collinear",state.score(g,state.project(g)).status());
        assertEquals("monomorphic-or-collinear",model.refit(d.y(),constant,state).status());
    }
    @Test void fullPipelinePreservesAlignmentAndCheckpointBytes(@TempDir Path root) throws Exception {
        Data d = data();
        StringBuilder genotypes = new StringBuilder("ID"), traits = new StringBuilder("ID");
        for (String id : d.ids()) genotypes.append(',').append(id);
        for (int i = d.ids().length-1; i >= 0; i--) traits.append(',').append(d.ids()[i]);
        genotypes.append("\ng"); for (double v : d.g()) genotypes.append(',').append(v);
        genotypes.append("\nconstant"); for (String id : d.ids()) genotypes.append(",2"); genotypes.append('\n');
        traits.append("\ny"); for (int i = d.y().length-1; i >= 0; i--) traits.append(',').append(d.y()[i]); traits.append('\n');
        Path gp = root.resolve("g.csv"), yp = root.resolve("y.csv"), output = root.resolve("out.csv");
        Files.writeString(gp,genotypes); Files.writeString(yp,traits);
        Map<String,String> values = new LinkedHashMap<>();
        values.put("genotype_file",gp.toString()); values.put("expression_file",yp.toString()); values.put("output_file",output.toString());
        values.put("covariate_file",ROOT.resolve("data.csv").toAbsolutePath().toString()); values.put("covariate_fixed","x");
        values.put("covariate_random","group"); values.put("pedigree_file",ROOT.resolve("pedigree.tsv").toAbsolutePath().toString());
        values.put("genotype_id_column","id"); values.put("expression_id_column","id"); values.put("keep_checkpoints","true");
        values.put("mixed_model","score"); values.put("genotype_block_rows","1");
        var config = new QeQTLAnalysisConfig(values);
        String old = System.getProperty(CpuBackend.BLAS_PROPERTY);
        try {
            System.setProperty(CpuBackend.BLAS_PROPERTY,"java");
            run(config,gp,yp);
            byte[] first = Files.readAllBytes(output);
            var lines = Files.readAllLines(output); assertEquals(3,lines.size());
            assertTrue(lines.get(1).startsWith("\"g\",\"y\",48,1.0,"));
            String[] fields = lines.get(1).split(",");
            assertEquals(-.288385915578612,Double.parseDouble(fields[7]),3e-6);
            assertTrue(lines.get(2).contains("monomorphic-or-collinear"));
            config.set("resume","true"); run(config,gp,yp); assertArrayEquals(first,Files.readAllBytes(output));
            Files.writeString(gp,genotypes.toString().replace("\ng,", "\ng_changed,"));
            assertThrows(java.io.IOException.class,()->run(config,gp,yp));
            assertArrayEquals(first,Files.readAllBytes(output));
        } finally {
            if (old == null) System.clearProperty(CpuBackend.BLAS_PROPERTY); else System.setProperty(CpuBackend.BLAS_PROPERTY,old);
            QeQTLAnalysis.config = null; QeQTLAnalysis.mContexts = null;
        }
    }
    static void run(QeQTLAnalysisConfig config, Path g, Path y) throws Exception {
        QeQTLAnalysis.config = config; QeQTLAnalysis.profiler = new QeQTLProfiler(false);
        QeQTLAnalysis.gpuPrecision = GpuPrecision.FP64;
        try (var context = new CpuBackend().discoverGpuDevices().get(0).openContext()) {
            QeQTLAnalysis.mContexts = new GpuContext[]{context};
            QeQTLAnalysis.runMatrixAnalysis(new QeQTLAnalysis(),g.toString(),y.toString(),config.getCovariateFilename(),
                new String[]{"x"},null,"none",Double.NaN,0,true,16,1,"csv");
        }
    }
    @Test void failedFitAndInvalidModeCannotProduceResults() throws Exception {
        Data d = data(); var m = model("combined","residual");
        assertThrows(IllegalArgumentException.class,()->m.fitNull(new double[d.y().length]));
        var low = new QMixedModel(d.x(),List.of(org.jlinalg.reml.VarianceComponent.randomIntercept("group",
            Arrays.asList(d.ids()).stream().map(id->Integer.parseInt(id.substring(1))/6).toList()),
            org.jlinalg.reml.VarianceComponent.identity("residual",d.ids().length)),1,"residual");
        assertThrows(IllegalArgumentException.class,()->low.fitNull(d.y()));
        assertThrows(IllegalArgumentException.class,()->QeQTLCommandLine.parse(new String[]{"--pedigree","p.tsv","--precision","fp32"}));
    }
    @Test void responseUnitsDoNotChangeScoreInferenceAndDfOptionsReachTheFitter() throws Exception {
        Data d = data(); var model = model("combined","residual"); var state = model.fitNull(d.y());
        double[] g = state.center(d.g()); var expected = state.score(g,state.project(g));
        for (double scale : new double[]{1e-6,1e6}) {
            double[] y = d.y().clone(); for (int i = 0; i < y.length; i++) y[i] *= scale;
            var scaled = model.fitNull(y); var actual = scaled.score(g,scaled.project(g));
            assertEquals(expected.log10P(),actual.log10P(),3e-7);
            assertEquals(expected.beta(),actual.beta()/scale,3e-7);
        }
        for (String df : List.of("satterthwaite","kenward-roger")) {
            var inference = model("combined",df).refit(d.y(),d.g(),state);
            assertEquals(.138182504098203,inference.beta(),3e-6);
            assertTrue(inference.df()>0 && Double.isFinite(inference.df()));
            assertNotEquals(45,inference.df());
            assertTrue(Double.isFinite(inference.log10P()));
        }
    }
}
