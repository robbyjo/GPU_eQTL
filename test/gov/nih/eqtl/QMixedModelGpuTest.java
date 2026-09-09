/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.*;
import java.util.List;
import gov.nih.gpu.*;
import gov.nih.gpu.cpu.CpuBackend;
import gov.nih.gpu.cuda.CudaGpuBackend;
import gov.nih.gpu.opencl.JoclGpuBackend;

class QMixedModelGpuTest {
    @Test void cudaProjectionMatchesScalar() throws Exception { check(new CudaGpuBackend()); }
    @Test void openClProjectionMatchesScalar() throws Exception { check(new JoclGpuBackend()); }
    @Test void nativeCpuProjectionMatchesScalar() throws Exception { check(new CpuBackend()); }
    private void check(GpuBackend backend) throws Exception {
        List<GpuDevice> devices;
        try { devices = new GpuRuntime(backend).getGpuDevices(true,true); }
        catch (RuntimeException | LinkageError e) { Assumptions.abort("Runtime unavailable: "+e); return; }
        Assumptions.assumeFalse(devices.isEmpty(),"No eligible "+backend.getName()+" device");
        var data = QMixedModelTest.data(); var state = QMixedModelTest.model("combined","residual").fitNull(data.y());
        String source = "#define BLOCK_SIZE 16\n#define DATATYPE double\n#define N_MIN_1 1\n"
            + "#if defined(cl_khr_fp64)\n#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
            + "#elif defined(cl_amd_fp64)\n#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n#endif\n"+QeQTLAnalysis.eqtlReal;
        try (var pool = new GpuContextPool(new GpuContext[]{devices.get(0).openContext()})) {
            pool.getAllContexts().get(0).compileKernel(source,"eqtlReal",GpuPrecision.FP64);
            for (int count : new int[]{1,17,33}) {
                double[][] variants = new double[count][];
                for (int g = 0; g < count; g++) {
                    double[] row = data.g().clone();
                    for (int i = 0; i < row.length; i++) row[i] += .05*Math.sin((i+1)*(g+1));
                    variants[g] = state.center(row);
                }
                double[][] actual = QMixedModelRunner.project(state,variants,
                    QMixedModelRunner.packProjection(state,data.y().length),pool);
                for (int g = 0; g < count; g++) {
                    double[] expected = state.project(variants[g]);
                    assertArrayEquals(expected,actual[g],3e-12);
                    assertEquals(state.score(variants[g],expected).log10P(),state.score(variants[g],actual[g]).log10P(),1e-11);
                }
            }
            System.out.println("MIXED_PROJECTION_PASS "+backend.getName()+" "+devices.get(0).getName());
        }
    }
}
