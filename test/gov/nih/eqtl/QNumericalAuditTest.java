/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.eqtl.io.QMatrixRowSource.Block;
import gov.nih.eqtl.settest.QSkatReference;

class QNumericalAuditTest {
    @Test
    void mixedModelsCannotSilentlyRunAsFixedEffects() throws Exception {
        assertThrows(IllegalArgumentException.class, () -> QeQTLCommandLine.parse(new String[]{"--random-covariates","Family"}));
        assertEquals("score", QeQTLCommandLine.parse(new String[]{"--pedigree","pedigree.tsv"}).config().getMixedModel());
        assertThrows(IllegalArgumentException.class, () -> QeQTLCommandLine.parse(new String[]{"--pedigree","pedigree.tsv","--mixed-model","none"}));
        assertThrows(IllegalArgumentException.class, () -> new QeQTLAnalysisConfig(
            java.util.Map.of("covariate_random","Family")).requireSupportedAnalysisModel());
    }
    @Test
    void whitenedAssociationMatchesIndependentWeightedLeastSquares(@TempDir Path root) throws Exception {
        // Diagonal GLS is independently solved as weighted centered least squares.
        double[] x = {0, 1, 2, 0, 2, 1};
        double[] y = {1, 3, 2, 4, 7, 2};
        double[] v = {1, 4, 2, 8, 3, 6};
        StringBuilder table = new StringBuilder("id,cohort\n");
        StringBuilder covariance = new StringBuilder("ID");
        String[] ids = new String[x.length];
        boolean[] selected = new boolean[x.length];
        for (int i = 0; i < x.length; i++) {
            ids[i] = "s" + i; selected[i] = true;
            table.append(ids[i]).append(",A\n"); covariance.append('\t').append(ids[i]);
        }
        covariance.append('\n');
        for (int i = 0; i < x.length; i++) {
            covariance.append(ids[i]);
            for (int j = 0; j < x.length; j++) covariance.append('\t').append(i == j ? v[i] : 0);
            covariance.append('\n');
        }
        Files.writeString(root.resolve("cov.csv"), table);
        Files.writeString(root.resolve("v.tsv"), covariance);
        Files.writeString(root.resolve("model.tsv"), "COHORT\tFIXED_COVARIATES\tCOVARIANCE_FILE\nA\t\tv.tsv\n");
        QCohortModel model = QCohortModel.create(QCohortModel.Definitions.load(root.resolve("model.tsv")),
            QCovariateTable.load(root.resolve("cov.csv"), new char[]{','}, "#"), "cohort", selected, ids, null, null);
        double[][] q = new gov.nih.jama.QRDecomposition(model.design()).getQ().getArray();
        var px = QeQTLPreprocessor.prepare(new Block(0, new String[]{"x"}, new double[][]{x.clone()}), q, "x", model.wrap(null));
        var py = QeQTLPreprocessor.prepare(new Block(0, new String[]{"y"}, new double[][]{y.clone()}), q, "y", model.wrap(null));
        double sw = 0, sx = 0, sy = 0;
        for (int i = 0; i < x.length; i++) { sw += 1/v[i]; sx += x[i]/v[i]; sy += y[i]/v[i]; }
        double xx = 0, xy = 0, yy = 0;
        for (int i = 0; i < x.length; i++) {
            double dx = x[i] - sx/sw, dy = y[i] - sy/sw;
            xx += dx*dx/v[i]; xy += dx*dy/v[i]; yy += dy*dy/v[i];
        }
        double r = QeQTLPreprocessor.correlation(px.values()[0], py.values()[0]);
        var actual = QeQTLStatistics.calculate(r, py.standardDeviations()[0], px.standardDeviations()[0], 4, 0);
        assertEquals(xy/xx, actual.effect(), 2e-13);
        assertEquals(xy*xy/(xx*yy), actual.rSquared(), 2e-13);
        assertEquals((xy/xx)/Math.sqrt((yy-xy*xy/xx)/4/xx), actual.tStatistic(), 2e-13);
    }

    @Test
    void kernelDistributionIsInvariantToCommonScaleAndRejectsInvalidSpectra() {
        double reference = QSkatReference.mixturePValue(8, new double[]{2, 2}).pValue();
        for (double scale : new double[]{1e-100, 1e-20, 1e20, 1e100})
            assertEquals(reference, QSkatReference.mixturePValue(8*scale, new double[]{2*scale, 2*scale}).pValue(), 2e-14);
        assertThrows(IllegalArgumentException.class, () -> QSkatReference.mixturePValue(1, new double[]{1, Double.NaN}));
        assertThrows(IllegalArgumentException.class, () -> QSkatReference.mixturePValue(1, new double[]{1, -1}));
    }

    @Test
    void zeroTraitCannotProduceAnOmnibusPValue() {
        assertThrows(IllegalArgumentException.class, () -> QSkatReference.calculateOmnibus(
            new double[][]{{1,-1,1,-1}}, new double[]{1}, new double[4], 3, new double[]{0,1}, 9, 1));
    }
}
