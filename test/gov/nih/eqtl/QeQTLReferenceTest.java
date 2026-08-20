/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.eqtl.io.QDelimitedMatrixSource;
import gov.nih.eqtl.io.QSampleAlignment;
import gov.nih.jama.QRDecomposition;

class QeQTLReferenceTest {
    @Test
    void deterministicFixturePreservesIdsPairsDegreesOfFreedomAndStatistics() throws Exception {
        QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/genotype.csv"), cCommonDelimiter, "#");
        QDelimitedMatrixSource expression = new QDelimitedMatrixSource(
            Path.of("test/resources/eqtl-reference/expression.csv"), cCommonDelimiter, "#");
        QCovariateTable covariates = QCovariateTable.load(
            Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");
        QSampleAlignment alignment = covariates.align(genotype.metadata().sampleIds(),
            expression.metadata().sampleIds(), "genotype_id", "expression_id");
        QCovariateTable.ModelMatrix model = covariates.buildModelMatrix(new String[] { "Age", "Batch" }, null);
        QRDecomposition qr = new QRDecomposition(model.values());
        int regressionDegreesOfFreedom = qr.getRank();
        int errorDegreesOfFreedom = alignment.sampleCount() - regressionDegreesOfFreedom - 1;
        assertEquals(3, regressionDegreesOfFreedom);
        assertEquals(4, errorDegreesOfFreedom);

        QDelimitedMatrixSource.Block genotypeBlock;
        QDelimitedMatrixSource.Block expressionBlock;
        try (QDelimitedMatrixSource.BlockReader reader = genotype.open(alignment.genotypeColumnOrder())) {
            genotypeBlock = reader.readBlock(16);
        }
        try (QDelimitedMatrixSource.BlockReader reader = expression.open(alignment.expressionColumnOrder())) {
            expressionBlock = reader.readBlock(16);
        }
        double[][] covariateQ = qr.getQ().getArray();
        PreparedBlock snps = QeQTLPreprocessor.prepare(genotypeBlock, covariateQ, "Genotype");
        PreparedBlock traits = QeQTLPreprocessor.prepare(expressionBlock, covariateQ, "Expression");

        List<String> identifiers = new ArrayList<>();
        List<Double> rSquared = new ArrayList<>();
        List<Double> effects = new ArrayList<>();
        List<Double> tStatistics = new ArrayList<>();
        List<Double> log10P = new ArrayList<>();
        for (int snp = 0; snp < snps.values().length; snp++) {
            for (int trait = 0; trait < traits.values().length; trait++) {
                double correlation = QeQTLPreprocessor.correlation(snps.values()[snp], traits.values()[trait]);
                QeQTLStatistics.Result result = QeQTLStatistics.calculate(correlation,
                    traits.standardDeviations()[trait], snps.standardDeviations()[snp], errorDegreesOfFreedom, 0);
                identifiers.add(snps.rowIds()[snp] + "/" + traits.rowIds()[trait]);
                rSquared.add(result.rSquared());
                effects.add(result.effect());
                tStatistics.add(result.tStatistic());
                log10P.add(result.log10P());
            }
        }

        assertEquals(6, identifiers.size());
        assertArrayEquals(new String[] { "rs1/iso1", "rs1/iso2", "rs2/iso1", "rs2/iso2", "rs3/iso1", "rs3/iso2" },
            identifiers.toArray(String[]::new));
        assertArrayEquals(new double[] {
            0.33429286679254133, 0.08441598953807962, 0.18924585599861957,
            0.012958263281588318, 0.2735057531379482, 0.13492926768062868 }, toArray(rSquared), 1e-12);
        assertArrayEquals(new double[] {
            0.2594575510480451, -0.07796736407339734, -0.24749432530757393,
            0.03872783242296443, 0.38742926434923264, -0.16272714491582085 }, toArray(effects), 1e-12);
        assertArrayEquals(new double[] {
            1.4172678817212558, -0.6072859782020563, -0.9662702228081002,
            0.229158323805691, 1.2271485511712632, -0.7898729984065591 }, toArray(tStatistics), 1e-12);
        assertArrayEquals(new double[] {
            -0.6394636593437162, -0.23924073411641683, -0.410466712182518,
            -0.08092923510992273, -0.5420393432793345, -0.3244133736673645 }, toArray(log10P), 1e-12);
    }

    private static double[] toArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int i = 0; i < values.size(); i++)
            result[i] = values.get(i);
        return result;
    }
}
