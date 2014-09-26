package gov.nih.exon.coulter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public interface ICoulterResultFilter {
	public double[][] filter(double[][] exonData, double[][] result);
}
