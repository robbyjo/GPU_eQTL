/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.nio.file.Path;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import gov.nih.eqtl.QeQTLPreprocessor.PreparedBlock;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.eqtl.io.QDelimitedMatrixSource;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.eqtl.io.QMissingnessScan;
import gov.nih.eqtl.io.QPolicyMatrixSource;
import gov.nih.eqtl.io.QSampleAlignment;
import gov.nih.eqtl.io.QSampleAlignmentPolicy;
import gov.nih.jama.QRDecomposition;

class QeQTLReferenceTest {
	@Test
	void prefixNormalizationAndCovariateSubsetPreserveReferenceStatistics(
		@TempDir Path directory) throws Exception {
		Path genotypeFile = directory.resolve("genotype-superset.csv");
		Files.writeString(genotypeFile,
			",G3,EXTRA1,G1,G8,G2,G7,G4,G6,G5,EXTRA2\n"
			+ "rs1,1.7,0.5,0.1,1.1,0.4,0.2,1.4,0.8,1.9,1.5\n"
			+ "rs2,0.2,1.0,1.8,0.7,1.4,1.1,0.4,1.6,0.9,0.3\n"
			+ "rs3,1.0,0.2,0.0,1.7,0.3,1.4,0.7,1.2,1.9,1.8\n");
		Path traitFile = directory.resolve("traits-prefixed.csv");
		Files.writeString(traitFile,
			",XT2,XT8,XT1,XT5,XT4,XT7,XT3,XT6\n"
			+ "iso1,2.2,4.9,1.0,3.6,2.9,4.3,2.4,3.8\n"
			+ "iso2,4.0,1.5,4.8,2.8,3.4,1.9,3.7,2.3\n");
		QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
			genotypeFile, cCommonDelimiter, "#");
		QDelimitedMatrixSource traits = new QDelimitedMatrixSource(
			traitFile, cCommonDelimiter, "#");
		QCovariateTable covariates = QCovariateTable.load(
			Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");
		QSampleAlignment alignment = covariates.align(genotype.metadata().sampleIds(),
			traits.metadata().sampleIds(), "genotype_id", "expression_id",
			QSampleAlignmentPolicy.COVARIATE_SUBSET, null, "X");
		QRDecomposition qr = new QRDecomposition(covariates.buildModelMatrix(
			new String[] {"Age", "Batch"}, null).values());
		QMatrixRowSource.Block genotypeBlock;
		QMatrixRowSource.Block traitBlock;
		try (QMatrixRowSource.BlockReader reader = genotype.open(alignment.genotypeColumnOrder())) {
			genotypeBlock = reader.readBlock(16);
		}
		try (QMatrixRowSource.BlockReader reader = traits.open(alignment.expressionColumnOrder())) {
			traitBlock = reader.readBlock(16);
		}
		PreparedBlock snps = QeQTLPreprocessor.prepare(
			genotypeBlock, qr.getQ().getArray(), "Genotype");
		PreparedBlock preparedTraits = QeQTLPreprocessor.prepare(
			traitBlock, qr.getQ().getArray(), "Expression");
		List<Double> rSquared = new ArrayList<>();
		for (double[] snp : snps.values())
			for (double[] trait : preparedTraits.values()) {
				double correlation = QeQTLPreprocessor.correlation(snp, trait);
				rSquared.add(correlation * correlation);
			}
		assertArrayEquals(new String[] {"rs1", "rs2", "rs3"}, snps.rowIds());
		assertArrayEquals(new String[] {"iso1", "iso2"}, preparedTraits.rowIds());
		assertEquals(2, alignment.genotypeExtraSampleCount());
		assertEquals(8, alignment.expressionIdsStripped());
		assertArrayEquals(new double[] {
			0.33429286679254133, 0.08441598953807962, 0.18924585599861957,
			0.012958263281588318, 0.2735057531379482, 0.13492926768062868},
			toArray(rSquared), 1e-12);
	}

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

        QMatrixRowSource.Block genotypeBlock;
        QMatrixRowSource.Block expressionBlock;
        try (QMatrixRowSource.BlockReader reader = genotype.open(alignment.genotypeColumnOrder())) {
            genotypeBlock = reader.readBlock(16);
        }
        try (QMatrixRowSource.BlockReader reader = expression.open(alignment.expressionColumnOrder())) {
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

	@Test
	void exactTraitMissingnessPatternsRecomputeSampleSetRankAndStatistics() throws Exception {
		QDelimitedMatrixSource genotype = new QDelimitedMatrixSource(
			Path.of("test/resources/eqtl-reference/genotype.csv"), cCommonDelimiter, "#");
		QDelimitedMatrixSource expression = new QDelimitedMatrixSource(
			Path.of("test/resources/missing-reference/expression.csv"), cCommonDelimiter, "#");
		QCovariateTable covariates = QCovariateTable.load(
			Path.of("test/resources/eqtl-reference/covariates.csv"), cCommonDelimiter, "#");
		QSampleAlignment alignment = covariates.align(genotype.metadata().sampleIds(),
			expression.metadata().sampleIds(), "genotype_id", "expression_id");
		QCovariateTable.ModelMatrix model = covariates.buildModelMatrix(
			new String[] {"Age", "Batch"}, null);
		QMissingnessScan scan = QMissingnessScan.scan("trait", expression,
			alignment.expressionColumnOrder());
		assertEquals(2, scan.patterns().size());

		List<String> identifiers = new ArrayList<>();
		List<Double> rSquared = new ArrayList<>();
		List<Double> effects = new ArrayList<>();
		List<Double> tStatistics = new ArrayList<>();
		List<Integer> sampleCounts = new ArrayList<>();
		List<Integer> residualDf = new ArrayList<>();
		for (QMissingnessScan.Pattern pattern : scan.patterns()) {
			int[] observed = complement(pattern.missingSamples(), alignment.sampleCount());
			double[][] selectedModel = selectRows(model.values(), observed);
			QRDecomposition qr = new QRDecomposition(selectedModel);
			assertEquals(3, qr.getRank());
			int dfe = observed.length - qr.getRank() - 1;
			BitSet selectedTraits = new BitSet((int) scan.rowCount());
			for (long row : pattern.rowIndices())
				selectedTraits.set((int) row);
			QMatrixRowSource selectedExpression = new QPolicyMatrixSource(expression, scan,
				QMissingValuePolicy.ERROR, selectedTraits);
			QMatrixRowSource.Block genotypeBlock;
			QMatrixRowSource.Block expressionBlock;
			try (QMatrixRowSource.BlockReader reader = genotype.open(
				select(alignment.genotypeColumnOrder(), observed))) {
				genotypeBlock = reader.readBlock(16);
			}
			try (QMatrixRowSource.BlockReader reader = selectedExpression.open(
				select(alignment.expressionColumnOrder(), observed))) {
				expressionBlock = reader.readBlock(16);
			}
			PreparedBlock snps = QeQTLPreprocessor.prepare(genotypeBlock,
				qr.getQ().getArray(), "Genotype");
			PreparedBlock traits = QeQTLPreprocessor.prepare(expressionBlock,
				qr.getQ().getArray(), "Expression");
			for (int snp = 0; snp < snps.values().length; snp++) {
				for (int trait = 0; trait < traits.values().length; trait++) {
					double correlation = QeQTLPreprocessor.correlation(
						snps.values()[snp], traits.values()[trait]);
					QeQTLStatistics.Result result = QeQTLStatistics.calculate(correlation,
						traits.standardDeviations()[trait], snps.standardDeviations()[snp], dfe, 0);
					identifiers.add(snps.rowIds()[snp] + "/" + traits.rowIds()[trait]);
					rSquared.add(result.rSquared());
					effects.add(result.effect());
					tStatistics.add(result.tStatistic());
					sampleCounts.add(observed.length);
					residualDf.add(dfe);
				}
			}
		}

		assertArrayEquals(new String[] {"rs1/iso1", "rs2/iso1", "rs3/iso1",
			"rs1/iso2", "rs2/iso2", "rs3/iso2"}, identifiers.toArray(String[]::new));
		assertArrayEquals(new Integer[] {8, 8, 8, 7, 7, 7}, sampleCounts.toArray(Integer[]::new));
		assertArrayEquals(new Integer[] {4, 4, 4, 3, 3, 3}, residualDf.toArray(Integer[]::new));
		assertArrayEquals(new double[] {0.33429286679254133, 0.18924585599861957,
			0.2735057531379482, 0.22047136051758234, 0.03308959550638096,
			0.35070209723441076}, toArray(rSquared), 1e-12);
		assertArrayEquals(new double[] {0.2594575510480451, -0.24749432530757393,
			0.38742926434923264, -0.11012349969527942, 0.05316436390087009,
			-0.23150143426852185}, toArray(effects), 1e-12);
		assertArrayEquals(new double[] {1.4172678817212558, -0.9662702228081002,
			1.2271485511712632, -0.921129477588173, 0.32041529561311494,
			-1.272939634032138}, toArray(tStatistics), 1e-12);
	}

	private static int[] complement(BitSet excluded, int size) {
		int[] result = new int[size - excluded.cardinality()];
		int output = 0;
		for (int i = excluded.nextClearBit(0); i < size; i = excluded.nextClearBit(i + 1))
			result[output++] = i;
		return result;
	}

	private static int[] select(int[] values, int[] indices) {
		int[] result = new int[indices.length];
		for (int i = 0; i < indices.length; i++)
			result[i] = values[indices[i]];
		return result;
	}

	private static double[][] selectRows(double[][] values, int[] indices) {
		double[][] result = new double[indices.length][];
		for (int i = 0; i < indices.length; i++)
			result[i] = values[indices[i]].clone();
		return result;
	}

    private static double[] toArray(List<Double> values) {
        double[] result = new double[values.size()];
        for (int i = 0; i < values.size(); i++)
            result[i] = values.get(i);
        return result;
    }
}
