package gov.nih.exon.coulter;

import static gov.nih.utils.QDataUtils.shortenArray;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterRSqResultFilter implements ICoulterResultFilter
{
	private double mRSq;

	public QCoulterRSqResultFilter(double rsq)
	{	mRSq = rsq; }

	@Override
	public double[][] filter(double[][] exonData, double[][] result) {
		int
			numRows = result.length,
			numCols = result[0].length,
			rsqCol = numCols - 1,
			numQualifiedRows = 0;
		double[][] selected = new double[numRows][];

		for (int rowNo = 0; rowNo < numRows; rowNo++)
			if (result[rowNo][rsqCol] >= mRSq)
				selected[numQualifiedRows++] = exonData[rowNo];
		return shortenArray(selected, numQualifiedRows);
	}
}
