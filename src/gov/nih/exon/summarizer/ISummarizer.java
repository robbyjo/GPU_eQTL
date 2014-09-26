package gov.nih.exon.summarizer;

/**
 * Summarizing multidimensional data into 1
 * @author Roby Joehanes
 *
 */
public interface ISummarizer {
	public double[] summarize(double[][] data);
	public String getPrefix();
}
