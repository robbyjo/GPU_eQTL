/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.io.TempDir;
import gov.nih.eqtl.*;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;
import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;

class QRelatedSetTest {
    static final Path ROOT = Path.of("test/resources/related-set-reference");
    record Data(String[] ids,double[][] x,double[] y,double[][] g) { }
    static Data data() throws Exception {
        var lines = Files.readAllLines(ROOT.resolve("data.csv")); int n = lines.size()-1;
        String[] ids = new String[n]; double[][] x = new double[n][2], g = new double[3][n]; double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            String[] v = lines.get(i+1).split(","); ids[i] = v[0]; x[i][0] = 1; x[i][1] = Double.parseDouble(v[1]); y[i] = Double.parseDouble(v[2]);
            for (int j = 0; j < 3; j++) g[j][i] = v[j+3].equals("NA") ? Double.NaN : Double.parseDouble(v[j+3]);
        }
        return new Data(ids,x,y,g);
    }
    static QMixedModel model(Data d,boolean pedigree) throws Exception {
        var config = new QeQTLAnalysisConfig(Map.of(pedigree ? "pedigree_file" : "relationship_matrix",
            pedigree ? "test/resources/mixed-reference/pedigree.tsv" : ROOT.resolve("grm.tsv").toString()));
        return QMixedModel.create(config,d.x(),null,null,d.ids());
    }
    static double[] response(Data d,String name) {
        double[] y = d.y().clone();
        if (name.equals("same-mask-new-trait")) for (int i = 0; i < y.length; i++) y[i] *= y[i];
        if (!name.equals("complete")) for (int i : new int[]{2,6,10,19,28}) y[i] = Double.NaN;
        return y;
    }
    static double[][] weighted(Data d,int[] mask) {
        double[] weights = {1,.7,1.3}; double[][] g = new double[3][mask.length];
        for (int j = 0; j < 3; j++) {
            double sum = 0; int count = 0;
            for (int i : mask) if (Double.isFinite(d.g()[j][i])) { sum += d.g()[j][i]; count++; }
            for (int i = 0; i < mask.length; i++) g[j][i] = weights[j]*(Double.isNaN(d.g()[j][mask[i]]) ? sum/count : d.g()[j][mask[i]]);
        }
        return g;
    }
    @Test void pedigreeGrmAndExactMasksMatchIndependentR() throws Exception {
        Data d = data();
        assertTrue(Arrays.stream(d.y()).anyMatch(x -> x == 0)); assertTrue(Arrays.stream(d.y()).anyMatch(x -> x == 1));
        double maxRelative = 0;
        try (var scalar = new QSetMatrixOperations(null)) {
            for (String line : Files.readAllLines(ROOT.resolve("scores.csv")).subList(1,7)) {
                String[] v = line.split(","); boolean related = v[1].equals("TRUE");
                for (boolean pedigree : related ? new boolean[]{false,true} : new boolean[]{false}) {
                    var nullModel = new QSetTraitModel(response(d,v[0]),d.x(),related ? model(d,pedigree) : null);
                    var scores = nullModel.scores(weighted(d,nullModel.samples()),scalar,false);
                    assertEquals(Integer.parseInt(v[2]),nullModel.sampleCount()); assertEquals(Integer.parseInt(v[3]),nullModel.nullDf());
                    double[][] c = scores.covariance();
                    double[] actual = {scores.u()[0],scores.u()[1],scores.u()[2],c[0][0],c[0][1],c[0][2],c[1][1],c[1][2],c[2][2]};
                    for (int i = 0; i < actual.length; i++) {
                        double expected = Double.parseDouble(v[i+4]), error = Math.abs(actual[i]-expected)/Math.max(1,Math.abs(expected));
                        maxRelative = Math.max(maxRelative,error);
                        assertEquals(expected,actual[i],3e-6*Math.max(1,Math.abs(expected)),line+" field="+i);
                    }
                    var burden = nullModel.scores(weighted(d,nullModel.samples()),scalar,true);
                    assertEquals(Double.parseDouble(v[13]),burden.u()[0]*burden.u()[0]/burden.covariance()[0][0],3e-6);
                    assertEquals(Double.parseDouble(v[14]),Arrays.stream(scores.u()).map(x -> x*x).sum(),3e-6*Math.max(1,Double.parseDouble(v[14])));
                }
            }
        }
        System.out.println("RELATED_SET_R_REFERENCE max scaled error="+maxRelative);
    }
    @Test void grmAlignmentAndValidationAreStrict(@TempDir Path directory) throws Exception {
        Path matrix = directory.resolve("grm.tsv");
        Files.writeString(matrix,"ID\tb\ta\tc\nc\t0\t0\t1\na\t0.5\t1\t0\nb\t1\t0.5\t0\n");
        assertArrayEquals(new double[]{1,.5,.5,1},QMixedModel.readRelationship(matrix,new String[]{"a","b"}),0);
        assertThrows(IllegalArgumentException.class,()->QMixedModel.readRelationship(matrix,new String[]{"a","missing"}));
        for (String text : List.of("ID\ta\ta\na\t1\t0\na\t0\t1\n","ID\ta\tb\na\t1\t.5\nb\t.1\t1\n",
            "ID\ta\tb\na\t1\t2\nb\t2\t1\n","ID\ta\tb\na\t1\t0\n")) {
            Files.writeString(matrix,text);
            assertThrows(IllegalArgumentException.class,()->QMixedModel.readRelationship(matrix,new String[]{"a","b"}));
        }
    }
    @Test void cudaSetProductsMatchScalar() throws Exception { hardware(new CudaGpuBackend()); }
    @Test void openClSetProductsMatchScalar() throws Exception { hardware(new JoclGpuBackend()); }
    @Test void nativeCpuSetProductsMatchScalar() throws Exception { hardware(new CpuBackend()); }
    private void hardware(GpuBackend backend) throws Exception {
        List<GpuDevice> devices;
        try { devices = new GpuRuntime(backend).getGpuDevices(true,true); }
        catch (RuntimeException | LinkageError e) { Assumptions.abort("Runtime unavailable: "+e); return; }
        Assumptions.assumeFalse(devices.isEmpty(),"No eligible "+backend.getName()+" device");
        Data d = data();
        try (var context = devices.get(0).openContext(); var accelerated = new QSetMatrixOperations(new GpuContext[]{context});
            var scalar = new QSetMatrixOperations(null)) {
            for (boolean related : new boolean[]{false,true}) {
                var model = new QSetTraitModel(response(d,"pattern"),d.x(),related ? model(d,false) : null);
                for (int size : new int[]{1,17,129}) {
                    double[][] rows = new double[size][model.sampleCount()];
                    for (int i = 0; i < size; i++) for (int j = 0; j < rows[i].length; j++) rows[i][j] = Math.sin((i+1)*(j+1)*.173)+.13*(j%3);
                    for (boolean burden : new boolean[]{false,true}) {
                        var expected = model.scores(rows,scalar,burden); var actual = model.scores(rows,accelerated,burden);
                        assertArrayEquals(expected.u(),actual.u(),2e-10);
                        for (int i = 0; i < expected.u().length; i++) assertArrayEquals(expected.covariance()[i],actual.covariance()[i],2e-9);
                        if (size == 1 || burden) assertEquals(expected.u()[0]*expected.u()[0]/expected.covariance()[0][0],
                            actual.u()[0]*actual.u()[0]/actual.covariance()[0][0],2e-10);
                    }
                }
                var expected = model.scores(weighted(d,model.samples()),scalar,false);
                var actual = model.scores(weighted(d,model.samples()),accelerated,false);
                assertEquals(new QSetScoreDistribution(expected.covariance()).skat(expected.u()).p(),
                    new QSetScoreDistribution(actual.covariance()).skat(actual.u()).p(),2e-11);
                double[][] two = Arrays.copyOf(weighted(d,model.samples()),2);
                var scalarTwo = model.scores(two,scalar,false); var gpuTwo = model.scores(two,accelerated,false);
                var referenceO = new QSetScoreDistribution(scalarTwo.covariance()).omnibus(scalarTwo.u(),
                    new double[]{0,.5,1},false,99,31,1e-6,1048576);
                var acceleratedO = new QSetScoreDistribution(gpuTwo.covariance()).omnibus(gpuTwo.u(),
                    new double[]{0,.5,1},false,99,31,1e-6,1048576);
                assertEquals(referenceO.p(),acceleratedO.p(),referenceO.p()*2e-6);
            }
            assertTrue(accelerated.executedProducts().keySet().stream().anyMatch(k -> k.endsWith("fitted-P-times-genotypes")));
            System.out.println("RELATED_SET_GPU_PASS "+backend.getName()+" "+devices.get(0).getName()+" "+devices.get(0).getDriverVersion());
        }
    }
}
