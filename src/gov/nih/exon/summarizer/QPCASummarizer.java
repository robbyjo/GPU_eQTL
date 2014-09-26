package gov.nih.exon.summarizer;

import gov.nih.jama.Matrix;
import gov.nih.jama.SingularValueDecomposition;

import static qstats.QStatsUtils.calcCorrelationMatrixTransposed;

/**
 * Summarizer based on Correlation PCA (only first PC)
 * @author Roby Joehanes
 *
 */
public class QPCASummarizer implements ISummarizer {
	@Override
	public double[] summarize(double[][] allData) {
		int dim = allData.length;
		if (dim == 0)
			return null;
		int
			p = allData[0].length,
			sectionStartIdx[] = new int[] {0}, //{0, 40, 80},
			numSections = sectionStartIdx.length;
		double result[] = new double[p];

		for (int sectionNo = 0; sectionNo < numSections; sectionNo++) {
			int
				curStart = sectionStartIdx[sectionNo],
				curP = (sectionNo == numSections - 1 ? p : sectionStartIdx[sectionNo + 1]) - curStart;
			double data[][];
			if (numSections > 1) {
				data = new double[dim][curP];
				for (int i = 0; i < dim; i++)
					System.arraycopy(allData[i], curStart, data[i], 0, curP);
			} else
				data = allData;
			double corr[][] = calcCorrelationMatrixTransposed(data);
			data = null;
			Matrix matrix = new Matrix(corr);
			SingularValueDecomposition svd = new SingularValueDecomposition(matrix);
			corr = null;

			double
				uMatrix[][] = svd.getU().getArray(),
				d[] = svd.getSingularValues(),
				lambda = d[0];
			for (int i = 0; i < curP; i++) // pick only first PC
				result[curStart + i] = uMatrix[i][0] * lambda;
			svd.purge();
			svd = null; matrix = null;
			/*
			double sum = 0;
			curP = d.length;
			for (int i = 0; i < curP; i++)
				sum += d[i];
			System.out.print(lambda / sum + "\t"); // Print proportion for diagnostics
			//*/
		}
		//System.out.println();
		return result;
	}

	@Override
	public String getPrefix() {
		return "pca";
	}
}
