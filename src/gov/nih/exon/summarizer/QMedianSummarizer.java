package gov.nih.exon.summarizer;

import static java.util.Arrays.sort;
import static qmath.QMathConstants.kUndefinedValue;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QMedianSummarizer implements ISummarizer {
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
			temp[] = new double[dim];
		for (int i = 0; i < n; i++) {
			numData = 0;
			for (int j = 0; j < dim; j++) {
				double val = data[j][i];
				if (val != kUndefinedValue)
					temp[numData++] = val;
			}
			switch (numData) {
				case 0: result[i] = kUndefinedValue; break;
				case 1: result[i] = temp[0]; break;
				case 2: result[i] = (temp[0] + temp[1]) / 2; break;
				default:
					sort(temp, 0, numData);
					int med = numData/2;
					//if (med + 1 >= numData)
					//	System.out.println("Gotcha!");
					result[i] = numData % 1 == 1 ? temp[med] : (temp[med] + temp[med + 1])/2;
			}
		}
		return result;
	}

	@Override
	public String getPrefix() {
		return "median";
	}
}
