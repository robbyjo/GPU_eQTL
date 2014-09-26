package gov.nih.exon.coulter;

import gov.nih.parallel.ISimpleParallelJob;

import java.util.Set;
import java.util.TreeSet;

import jdistlib.rng.MersenneTwister;

import static gov.nih.exon.coulter.QCoulterDataAnalysis.analyze;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.test;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterCrossValidationJob implements ISimpleParallelJob {
	protected double
		mFilteredExons[][],
		mCoulterData[][],
		mXVResults[][],
		mRSqResult[];
	protected int mNumTest;
	private int[] mTrainingIdx, mTestIdx;
	private Set<Integer> mTestSet = new TreeSet<Integer>();
	private MersenneTwister mRandom = new MersenneTwister();

	public QCoulterCrossValidationJob(double[][] filteredExons, double[][] coulterData, double[][] results, int numTestedInd) {
		mFilteredExons = filteredExons;
		mCoulterData = coulterData;
		mXVResults = results;
		mNumTest = numTestedInd;
		mTestIdx = new int[mNumTest];
		mTrainingIdx = new int[coulterData.length - mNumTest];
		mRSqResult = new double[filteredExons.length];
	}

	@Override
	public boolean execute(int testNo) {
		int numSubjects = mCoulterData.length;
		mTestSet.clear();
		while (mTestSet.size() < mNumTest)
			mTestSet.add(mRandom.nextInt(numSubjects));
		for (int subjectNo = 0, idx = 0, idx2 = 0; subjectNo < numSubjects; subjectNo++)
			if (!mTestSet.contains(subjectNo))
				mTrainingIdx[idx++] = subjectNo;
			else
				mTestIdx[idx2++] = subjectNo;
		double[][] trainingResult = analyze(mFilteredExons, mCoulterData, mRSqResult, mTrainingIdx, false);
		double[] testResult = test(mFilteredExons, mCoulterData, trainingResult, mTestIdx);
		mXVResults[testNo] = testResult;
		return true;
	}

	@Override
	public QCoulterCrossValidationJob clone() {
		return new QCoulterCrossValidationJob(mFilteredExons, mCoulterData, mXVResults, mNumTest);
	}
}
