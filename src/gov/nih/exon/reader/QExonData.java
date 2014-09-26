package gov.nih.exon.reader;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QExonData
{
	protected int[] mExonIDs;
	protected String[] mGeneIDs;
	protected double[][] mData;

	public QExonData(int[] exonIDs, String[] geneIDs, double[][] data)
	{
		mExonIDs = exonIDs;
		mGeneIDs = geneIDs;
		mData = data;
	}

	public int[] getExonIDs()
	{	return mExonIDs; }

	public String[] getGeneIDs()
	{	return mGeneIDs; }

	public double[][] getData()
	{	return mData; }

	public void purge()
	{
		mExonIDs = null;
		mGeneIDs = null;
		mData = null;
	}
}
