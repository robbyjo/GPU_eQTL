package gov.nih.exon.coulter;

import gov.nih.parallel.ISimpleParallelJob;
import jama.QRDecomposition;
import qstats.glm.QGLMResult;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterAnalysisJob implements ISimpleParallelJob
{
	protected double[][] mCoulter, mExon, mResult;
	protected double[] mRSq;
	protected QRDecomposition mQR;
	protected int[] mSubjectIdx;

	// Cache
	protected double[][] mY, mTempX, mBetas;
	protected QGLMResult mGLM;

	public QCoulterAnalysisJob(double[][] exon, double[][] coulter, double[][] result, double[] rsq, int[] subjectIdx) {
		init(new QRDecomposition(coulter), exon, coulter, result, rsq, subjectIdx);
	}

	private QCoulterAnalysisJob(QRDecomposition qr, double[][] exon, double[][] coulter, double[][] result, double[] rsq, int[] subjectIdx) {
		init(qr, exon, coulter, result, rsq, subjectIdx);
	}

	private void init(QRDecomposition qr, double[][] exon, double[][] coulter, double[][] result, double[] rsq, int[] subjectIdx) {
		mExon = exon;
		mCoulter = coulter;
		mResult = result;
		mSubjectIdx = subjectIdx;
		mQR = qr;
		mRSq = rsq;

		// Initialize caches
		int numSubjects = subjectIdx.length;
		mTempX = new double[1][numSubjects];
		mY = new double[1][numSubjects];
		mBetas = new double[1][];
		mGLM = new QGLMResult(mCoulter, mY, qr, false);
	}

	@Override
	public boolean execute(int iterNo) {
		// Shrinkage
		//mResult[iterNo] = calcShrinkageBetas(mExon[iterNo], mCoulter, mEta, mThreshold, mLODThreshold, mCritValue);

		// Old regression routine
		//mResult[iterNo] = qstats.QStatsUtils.multipleRegressionFindBetaCorrectX(mCoulter, mExon[iterNo]).mBeta;

		// New regression routine
		int numSubjects = mSubjectIdx.length;
		double[]
			oldY = mExon[iterNo],
			newY = mY[0];
		for (int i = 0; i < numSubjects; i++)
			newY[i] = oldY[mSubjectIdx[i]];

		int idx = mGLM.mNumVars; // * 2;
		double[] betas = mBetas[0] = new double[idx];

		mGLM.update(mY, mTempX, mBetas);
		mGLM.mY = mY;
		mGLM.mBetas = mBetas;
		mGLM.calcMoreStats();
		//mGLM.calcFBetas();
		//mGLM.calcPOfFBetas();
		//System.arraycopy(mGLM.mPFBetas, 0, betas, mGLM.mNumVars, mGLM.mNumVars);
		mResult[iterNo] = betas;
		mRSq[iterNo] = mGLM.mRSq;
		return true;
	}

	@Override
	public QCoulterAnalysisJob clone() {
		return new QCoulterAnalysisJob(mQR, mExon, mCoulter, mResult, mRSq, mSubjectIdx);
	}
}
