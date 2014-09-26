package gov.nih.exon.snp;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QTrimodalGaussianDistribution {
	public double p0, mu0, sigma0Sq, p1, mu1, sigma1Sq, mu2, sigma2Sq; 

	public QTrimodalGaussianDistribution() {
		p0 = 0.25; mu0 = -1; sigma0Sq = 2; p1 = 0.5; mu1 = 0; sigma1Sq = 2; mu2 = 1; sigma2Sq = 2;
	}

	public QTrimodalGaussianDistribution(double _p0, double _mu0, double _sigmaSq0, double _p1, double _mu1, double _sigmaSq1, double _mu2, double _sigmaSq2) {
		p0 = _p0; mu0 = _mu0; sigma0Sq = _sigmaSq0; p1 = _p1; mu1 = _mu1; sigma1Sq = _sigmaSq1; mu2 = _mu2; sigma2Sq = _sigmaSq2;
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append(String.format("p0 = %f, μ0 = %f, σ0^2 = %f\n", p0, mu0, sigma0Sq));
		buf.append(String.format("p1 = %f, μ1 = %f, σ1^2 = %f\n", p1, mu1, sigma1Sq));
		buf.append(String.format("p2 = %f, μ2 = %f, σ2^2 = %f\n", 1-p0-p1, mu2, sigma2Sq));
		return buf.toString();
	}
}
