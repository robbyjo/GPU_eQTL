package gov.nih.exon.summarizer;

import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.IJobOwner;
import gov.nih.parallel.QSynchronizedCounter;

import java.util.Map;


/**
 * 
 * @author Roby Joehanes
 *
 */
public class QSummarizerJob extends Thread implements IGenericParallelJob {
	protected QSynchronizedCounter mLock, mCounter;
	protected boolean mIsCanceled = false;
	protected double[][] mResult;
	protected Map<String, double[][]> mGeneData;
	protected String[] mGeneIDs;
	protected ISummarizer mSummarizer;

	// Many of the codes are copy and paste from my eQTL Hotspot code
	public QSummarizerJob(Map<String, double[][]> geneData, String[] geneIDs, ISummarizer summarizer, double[][] result,
		QSynchronizedCounter lock, QSynchronizedCounter counter) {
		mGeneData = geneData;
		mGeneIDs = geneIDs;
		mSummarizer = summarizer;
		mResult = result;
		mLock = lock;
		mCounter = counter;
	}

	@Override
	public void cancel()
	{	mIsCanceled = true; }

	@Override
	public void execute()
	{
		int geneNo = mCounter.next();
		if (geneNo < 0)
			return;
		
		do {
			String geneID = mGeneIDs[geneNo];
			double[][] data = mGeneData.get(geneID);
			double[] summary;
			if (data.length == 1) {
				int length = data[0].length;
				summary = new double[length];
				System.arraycopy(data[0], 0, summary, 0, length);
			} else
				summary = mSummarizer.summarize(data);
			mResult[geneNo] = summary;
			if (mIsCanceled)
				break;
			geneNo = mCounter.next();
		} while (geneNo >= 0);

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