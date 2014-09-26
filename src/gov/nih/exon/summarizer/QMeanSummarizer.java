package gov.nih.exon.summarizer;

import static qmath.QMathConstants.kUndefinedValue;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QMeanSummarizer implements ISummarizer {
	@Override
	public double[] summarize(double[][] data) {
		int dim = data.length;
		if (dim == 0)
			return null;
		int
			n = data[0].length,
			numData;
		double
			result[] = new double[n],
			sum;
		for (int i = 0; i < n; i++) {
			sum = 0; numData = 0;
			for (int j = 0; j < dim; j++) {
				double val = data[j][i];
				if (val != kUndefinedValue) {
					sum += val; numData++;
				}
			}
			result[i] = sum / numData;
		}
		return result;
	}

	@Override
	public String getPrefix() {
		return "mean";
	}
}
