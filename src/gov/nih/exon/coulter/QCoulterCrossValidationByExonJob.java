package gov.nih.exon.coulter;

import gov.nih.parallel.ISimpleParallelJob;

import java.util.Set;
import java.util.TreeSet;

import jdistlib.rng.MersenneTwister;

import static gov.nih.utils.QDataUtils.shortenArray;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.analyze;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.test;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterCrossValidationByExonJob implements ISimpleParallelJob {
	protected double
		mExonData[][],
		mCoulterData[][],
		mXVResults[][][],
		mRSq[],
		mRSqResult[];
	protected int
		mNumTest,
		mTrainingIdx[],
		mTestIdx[],
		mAvgExonCount[],
		mExonCount[][]; // numRSq x numExons
	private Set<Integer> mTestSet = new TreeSet<Integer>();
	private MersenneTwister mRandom = new MersenneTwister();

	public QCoulterCrossValidationByExonJob(double[][] exonData, double[][] coulterData, double[] rsq,
		double[][][] xvResults, int[][] exonCount, int[] avgExonCount, int numTestedInd) {
		mExonData = exonData;
		mCoulterData = coulterData;
		mXVResults = xvResults;
		mRSq = rsq;
		mExonCount = exonCount;
		mAvgExonCount = avgExonCount;
		mNumTest = numTestedInd;
		mTestIdx = new int[mNumTest];
		mTrainingIdx = new int[coulterData.length - mNumTest];
		mRSqResult = new double[exonData.length];
	}

	@Override
	public boolean execute(int testNo) {
		int
			numSubjects = mCoulterData.length,
			numExons = mExonData.length,
			numRSq = mRSq.length;
		mTestSet.clear();
		while (mTestSet.size() < mNumTest)
			mTestSet.add(mRandom.nextInt(numSubjects));
		for (int subjectNo = 0, idx = 0, idx2 = 0; subjectNo < numSubjects; subjectNo++)
			if (!mTestSet.contains(subjectNo))
				mTrainingIdx[idx++] = subjectNo;
			else
				mTestIdx[idx2++] = subjectNo;
		double[][] trainingResult = analyze(mExonData, mCoulterData, mRSqResult, mTrainingIdx, false);
		for (int rsqNo = 0; rsqNo < numRSq; rsqNo++) {
			double
				curRSq = mRSq[rsqNo],
				filteredResult[][] = new double[numExons][],
				filteredExons[][] = new double[numExons][];
			int
				exonIdx = 0,
				curExonCount[] = mExonCount[rsqNo];
			for (int exonNo = 0; exonNo < numExons; exonNo++) {
				double[]
					curExon = mExonData[exonNo],
					curResult = trainingResult[exonNo];
				if (mRSqResult[exonNo] >= curRSq) {
					filteredExons[exonIdx] = curExon;
					filteredResult[exonIdx] = curResult;
					exonIdx++;
					synchronized(curExonCount) {
						curExonCount[exonNo]++;
					}
				}
			}
			synchronized(mAvgExonCount) {
				mAvgExonCount[rsqNo] += exonIdx;
			}
			filteredExons = shortenArray(filteredExons, exonIdx);
			filteredResult = shortenArray(filteredResult, exonIdx);
			double[] testResult = test(filteredExons, mCoulterData, filteredResult, mTestIdx);
			mXVResults[rsqNo][testNo] = testResult;
		}
		testNo++;
		if (testNo % 100 == 0)
			System.out.println("Iteration " + testNo);
		return true;
	}

	@Override
	public QCoulterCrossValidationByExonJob clone() {
		QCoulterCrossValidationByExonJob job = new QCoulterCrossValidationByExonJob(mExonData, mCoulterData,
			mRSq, mXVResults, mExonCount, mAvgExonCount, mNumTest);
		return job;
	}
}
