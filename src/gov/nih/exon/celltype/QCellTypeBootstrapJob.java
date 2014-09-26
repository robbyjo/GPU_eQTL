package gov.nih.exon.celltype;


import gov.nih.exon.reader.QExonData;
import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.IJobOwner;
import gov.nih.parallel.QSynchronizedCounter;

import static java.lang.Math.abs;
import static java.util.Arrays.sort;

import static qstats.QStatsUtils.findQuantile;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCellTypeBootstrapJob extends Thread implements IGenericParallelJob {
	protected QSynchronizedCounter mLock, mCounter;
	protected boolean mIsCanceled = false;
	protected double[][] mResult;
	protected double[] mAlphas;
	protected int[][] mPermutation;
	protected QExonData mData;
	static final int kNumInds = 40, kNumInds2 = 2*kNumInds;

	// Many of the codes are copy and paste from my eQTL Hotspot code
	public QCellTypeBootstrapJob(QExonData data, int[][] permutation, double[][] result, double[] alphas,
		QSynchronizedCounter lock, QSynchronizedCounter counter) {
		mData = data;
		mResult = result;
		mPermutation = permutation;
		mAlphas = alphas;
		mLock = lock;
		mCounter = counter;
	}

	@Override
	public void cancel()
	{	mIsCanceled = true; }

	@Override
	public void execute()
	{
		int exonNo = mCounter.next();
		if (exonNo < 0)
			return;

		double[][] data = mData.getData();
		int
			numIter = mPermutation.length,
			numCols = data[0].length, // Should be 120
			shuffledIdx[];
		double
			sumCL,
			sumPAX,
			sumPBMC,
			pax_cl[] = new double[numIter],
			pbmc_cl[] = new double[numIter],
			pbmc_pax[] = new double[numIter],
			cellTypes[][] = new double[][] { pax_cl, pbmc_cl, pbmc_pax};
		int
			numCellTypes = cellTypes.length,
			numAlphas = mAlphas.length;
		bailOut:
		do {
			double[]
				curExon = data[exonNo],
				curResult = mResult[exonNo];
			for (int iterNo = 0; iterNo < numIter; iterNo++) {
				sumCL = sumPAX = sumPBMC = 0;
				shuffledIdx = mPermutation[iterNo];
				for (int i = 0; i < kNumInds; i++)
					sumCL += curExon[shuffledIdx[i]];
				for (int i = kNumInds; i < kNumInds2; i++)
					sumPAX += curExon[shuffledIdx[i]];
				for (int i = kNumInds2; i < numCols; i++)
					sumPBMC += curExon[shuffledIdx[i]];
				sumCL /= kNumInds;
				sumPAX /= kNumInds;
				sumPBMC /= kNumInds;
				pax_cl[iterNo] = abs(sumPAX - sumCL);
				pbmc_cl[iterNo] = abs(sumPBMC - sumCL);
				pbmc_pax[iterNo] = abs(sumPBMC - sumPAX);
				if (mIsCanceled)
					break bailOut;
			}
			sort(pax_cl);
			sort(pbmc_cl);
			sort(pbmc_pax);
			for (int cellTypeNo = 0, idx = 0; cellTypeNo < numCellTypes; cellTypeNo++) {
				double[] curType = cellTypes[cellTypeNo];
				for (int alphaNo = 0; alphaNo < numAlphas; alphaNo++, idx++) {
					curResult[idx] = findQuantile(curType, mAlphas[alphaNo]);
				}
			}
			exonNo = mCounter.next();
		} while (exonNo >= 0);

		//System.out.println("Thread #" + mID + " done");
		synchronized (mLock) {
			mLock.next();
			mLock.notifyAll();
		}
	}

	@Override
	public void run()
	{	execute(); }

	/**
	 * Set the owner of this job
	 * @param owner
	 */
	@Override
	public void setOwner(IJobOwner owner)
	{	}

	/**
	 * Get the owner of this job
	 * @return
	 */
	@Override
	public IJobOwner getOwner()
	{	return null; }
}