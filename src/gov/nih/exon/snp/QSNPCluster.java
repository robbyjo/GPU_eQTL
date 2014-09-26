package gov.nih.exon.snp;

import gov.nih.exon.annotation.QAnnotation;
import qstats.glm.QGLM;
import qstats.glm.QGLMResult;
import gov.nih.parallel.ISimpleParallelJob;
import gov.nih.table.QTableData;
import gov.nih.utils.QStringUtils;

import static java.lang.Math.*;
import static qstats.QStatsUtils.calcNormalPDF;
import static gov.nih.parallel.QParallelUtils.executeSimpleParallelJob;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableData;
import static gov.nih.utils.QFileUtils.writeText;
import static gov.nih.utils.QStringUtils.sTab;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;

/**
 * 
 * @author Roby Joehanes
 */
public class QSNPCluster implements ISimpleParallelJob {
	private static final double defaultEta = 1e-8;

	protected double residuals[][], results[][], eta;
	protected QTrimodalGaussianDistribution param;

	public QSNPCluster(double[][] _residuals, double[][] _results, QTrimodalGaussianDistribution _param) {
		this (_residuals, _results, _param, defaultEta);
	}

	public QSNPCluster(double[][] _residuals, double[][] _results, QTrimodalGaussianDistribution _param, double _eta) {
		residuals = _residuals; results = _results; param = _param; eta = _eta;
	}

	private static final QTrimodalGaussianDistribution estimateParam(double[] curResidual) {
		QTrimodalGaussianDistribution par = new QTrimodalGaussianDistribution();
		double min = curResidual[0], max = min, mean = min, var = min * min;
		int n = curResidual.length;
		for (int i = 1; i < n; i++) {
			double val = curResidual[i];
			mean += val;
			mean += val * val;
			if (val < min) min = val;
			if (val > max) max = val;
		}
		var = (var - mean * mean / n) / (n - 1);
		mean /= n;
		double sd = abs(mean-min)/2;
		par.mu0 = min; par.sigma0Sq = sd * sd;
		par.mu1 = mean; par.sigma1Sq = var;
		sd = abs(max-mean)/2;
		par.mu2 = max; par.sigma2Sq = sd * sd;

		min = mean - sd; max = mean + sd;
		int p0 = 1, p1 = 1;
		for (int i = 0; i < n; i++) {
			double val = curResidual[i];
			if (val < min) p0++;
			else if (val < max) p1++;
		}
		double dn = (n+2);
		par.p0 = p0 / dn; par.p1 = p1 / dn;
		return new QTrimodalGaussianDistribution();
	}

	@Override
	public boolean execute(int iterNo) {
		double[] curResidual = residuals[iterNo];
		QTrimodalGaussianDistribution res = estimateTrimodalGaussDistribution(curResidual, estimateParam(curResidual), eta);
		double[] curResult = results[iterNo];
		curResult[0] = Double.isNaN(res.p0) ? 0 : res.p0;
		curResult[1] = Double.isNaN(res.mu0) ? 0 : res.mu0;
		curResult[2] = Double.isNaN(res.sigma0Sq) ? 0 : res.sigma0Sq;
		curResult[3] = Double.isNaN(res.p1) ? 0 : res.p1;
		curResult[4] = Double.isNaN(res.mu1) ? 0 : res.mu1;
		curResult[5] = Double.isNaN(res.sigma1Sq) ? 0 : res.sigma1Sq;
		curResult[6] = Double.isNaN(res.mu2) ? 0 : res.mu2;
		curResult[7] = Double.isNaN(res.sigma2Sq) ? 0 : res.sigma2Sq;
		testMembership(curResidual, curResult);
		if ((iterNo % 10000) == 0) System.out.println(iterNo);
		return true;
	}

	@Override
	public QSNPCluster clone() {
		return new QSNPCluster(residuals, results, param, eta);
	}

	public static final QTrimodalGaussianDistribution estimateTrimodalGaussDistribution(double[] vec, QTrimodalGaussianDistribution param, double eta) {
		int n = vec.length;
		double
			oldLogLik, newLogLik = 0, sigma0Sq = param.sigma0Sq, sigma1Sq = param.sigma1Sq, sigma2Sq = param.sigma2Sq,
			p0 = param.p0, p1 = param.p1, p2 = 1-p0-p1, mu0 = param.mu0, mu1 = param.mu1, mu2 = param.mu2, sum0, sum1, sum2, sumvec0, sumvec1, sumvec2,
			member0[] = new double[n], member1[] = new double[n], member2[] = new double[n];

		do {
			// Compute membership and likelihood
			oldLogLik = newLogLik;
			newLogLik = sum0 = sum1 = sum2 = sumvec0 = sumvec1 = sumvec2 = 0;
			for (int i = 0; i < n; i++) {
				double v = vec[i],
					v0 = sigma0Sq > eta ? p0 * calcNormalPDF(v, mu0, sigma0Sq) : 0,
					v1 = sigma1Sq > eta ? p1 * calcNormalPDF(v, mu1, sigma1Sq) : 0,
					v2 = sigma2Sq > eta ? p2 * calcNormalPDF(v, mu2, sigma2Sq) : 0,
					sum = v0 + v1 + v2;
				sum0 += member0[i] = v0 = v0 / sum;
				sum1 += member1[i] = v1 = v1 / sum;
				sum2 += member2[i] = v2 = v2 / sum;
				sumvec0 += v0 * v; sumvec1 += v1 * v; sumvec2 += v2 * v;
				newLogLik += log(sum);
			}
			if (abs(oldLogLik - newLogLik) <= eta) break;
			// Update parameters
			p0 = sum0 / n; p1 = sum1 / n; p2 = 1-p0-p1;
			mu0 = sumvec0 / sum0; mu1 = sumvec1 / sum1; mu2 = sumvec2 / sum2;
			sumvec0 = sumvec1 = sumvec2 = 0;
			for (int i = 0; i < n; i++) {
				double diff, v = vec[i];
				diff = v - mu0; sumvec0 += member0[i] * diff * diff;
				diff = v - mu1; sumvec1 += member1[i] * diff * diff;
				diff = v - mu2; sumvec2 += member2[i] * diff * diff;
			}
			sigma0Sq = sumvec0 / sum0; sigma1Sq = sumvec1 / sum1; sigma2Sq = sumvec2 / sum2; //iterNo++;
			//if (iterNo % 1000 == 0)
			//	System.out.println(String.format("p0 = %f, μ0 = %f, σ0^2 = %f, p1 = %f, μ1 = %f, σ1^2 = %f, p2 = %f, μ2 = %f, σ2^2 = %f", p0, mu0, sigma0Sq, p1, mu1, sigma1Sq, 1-p0-p1, mu2, sigma2Sq));
		} while (true);
		//System.out.println(String.format("p0 = %f, μ0 = %f, σ0^2 = %f, p1 = %f, μ1 = %f, σ1^2 = %f, p2 = %f, μ2 = %f, σ2^2 = %f", p0, mu0, sigma0Sq, p1, mu1, sigma1Sq, 1-p0-p1, mu2, sigma2Sq));
		//System.out.println(iterNo);
		return new QTrimodalGaussianDistribution(p0, mu0, sigma0Sq, p1, mu1, sigma1Sq, mu2, sigma2Sq);
	}

	public static final double[][] clusterSNPs(double[][] residuals, QAnnotation[] annots) {
		QTrimodalGaussianDistribution param = new QTrimodalGaussianDistribution();
		int nrow = residuals.length;
		double[][] results = new double[nrow][12];
		executeSimpleParallelJob(new QSNPCluster(residuals, results, param), nrow);
		return results;
	}

	static final double[][] calcMembership(double[] vec, double p0, double mu0, double sigma0Sq, double p1, double mu1, double sigma1Sq, double p2, double mu2, double sigma2Sq) {
		int n = vec.length;
		double member0[] = new double[n], member1[] = new double[n], member2[] = new double[n];
		for (int i = 0; i < n; i++) {
			double v = vec[i],
				v0 = sigma0Sq > 0 ? p0 * calcNormalPDF(v, mu0, sigma0Sq) : 0,
				v1 = sigma1Sq > 0 ? p1 * calcNormalPDF(v, mu1, sigma1Sq) : 0,
				v2 = sigma2Sq > 0 ? p2 * calcNormalPDF(v, mu2, sigma2Sq) : 0,
				sum = v0 + v1 + v2;
			member0[i] = v0 / sum;
			member1[i] = v1 / sum;
			member2[i] = v2 / sum;
		}
		return new double[][] { member0, member1, member2 };
	}

	public static final void testMembership(double[] residuals, double[] d) {
		double p0 = d[0], mu0 = d[1], sigma0Sq = d[2], p1 = d[3], mu1 = d[4], sigma1Sq = d[5], p2, mu2 = d[6], sigma2Sq = d[7];
		boolean
			has0 = !Double.isNaN(p0) & !Double.isNaN(mu0) & !Double.isNaN(sigma0Sq) & p0 > 0 & sigma0Sq > 0,
			has1 = !Double.isNaN(p1) & !Double.isNaN(mu1) & !Double.isNaN(sigma1Sq) & p1 > 0 & sigma1Sq > 0,
			has2 = !Double.isNaN(mu2) & !Double.isNaN(sigma2Sq) & sigma2Sq > 0;
		int numClusters = 0;

		if (has0) numClusters++; else p0 = mu0 = sigma0Sq = 0;
		if (has1) numClusters++; else p1 = mu1 = sigma1Sq = 0;
		if (has2) { numClusters++; p2 = 1-p0-p1; } else p2 = mu2 = sigma2Sq = 0;
		d[8] = numClusters;
		if (numClusters <= 1) { d[9] = d[10] = d[11] = Double.NaN; return; }
		double[][] members = calcMembership(residuals, p0, mu0, sigma0Sq, p1, mu1, sigma1Sq, p2, mu2, sigma2Sq);
		QGLMResult result = null;
		if (numClusters == 3) {
			double[] member0 = null, member1 = null;
			if (p0 <= p1 && p1 <= p2) {
				member0 = members[1];
				member1 = members[0];
			} else if (p0 <= p2 && p2 <= p1) {
				member0 = members[2];
				member1 = members[0];
			} else if (p1 <= p0 && p0 <= p2) {
				member0 = members[0];
				member1 = members[1];
			} else if (p1 <= p2 && p2 <= p0) {
				member0 = members[2];
				member1 = members[1];
			} else if (p2 <= p0 && p0 <= p1) {
				member0 = members[0];
				member1 = members[2];
			} else { // p2 <= p1 && p1 <= p0
				member0 = members[1];
				member1 = members[2];
			}
			result = QGLM.multipleRegressionTransposed(new double[][] {member0, member1}, residuals, true);
		} else {
			double[] member0 = has0? (has1? (p0 <= p1? members[0] : members[1]) : (p0 <= p2? members[0] : members[2])) : (p1 <= p2? members[1] : members[2]);
			result = QGLM.multipleRegressionTransposed(new double[][] {member0}, residuals, true);
		}
		result.calcPValues();
		result.calcPOfFBetas();
		d[9] = result.mP_of_F;
		d[10] = result.mPFBetas[1];
		if (numClusters == 3)
			d[11] = result.mPFBetas[2];
	}

	public static final void clusterSNPs(String dataFile, String annotFile, String outputFile) {
		try {
			System.out.println("Data file = " + dataFile);
			System.out.println("Annotation file = " + annotFile);
			System.out.println("Output file = " + outputFile);
			long t1 = System.currentTimeMillis();
			QTableData tbl = readDelimitedFileAsTableData(dataFile, cCommaDelimiter, "#", true, true);
			long t2 = System.currentTimeMillis();
			double[][] data = tbl.extractData();
			System.out.println(data.length + " x " + data[0].length);
			System.out.println(String.format("Time = %d ms", (t2 - t1)));

			t1 = System.currentTimeMillis();
			QAnnotation[] annots = QAnnotation.load(annotFile);
			t2 = System.currentTimeMillis();
			System.out.println(String.format("Annotation of %d exons loaded.", annots.length));
			System.out.println(String.format("Time = %d ms", (t2 - t1)));

			t1 = System.currentTimeMillis();
			//estimateTrimodalGaussDistribution(data[0], new QTrimodalGaussianDistribution(), defaultEta);
			double[][] result = clusterSNPs(data, annots);
			t2 = System.currentTimeMillis();
			System.out.println(String.format("Time = %d ms", (t2 - t1)));

			//String[] resultColNames = new String[] { "p0", "mu0", "sigma0Sq", "p1", "mu1", "sigma1Sq", "mu2", "sigma2Sq"};
			String[] resultColNames = new String[] { "p0", "mu0", "sigma0Sq", "p1", "mu1", "sigma1Sq", "mu2", "sigma2Sq", "numClusters", "pF", "pClust1", "pClust2"};
			String str = QStringUtils.toDelimitedString(result, tbl.getAllRowNames(), resultColNames, sTab);
			writeText(str, outputFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		clusterSNPs(args[0], args[1], args[2]);
		//testMembership(args[0], args[1], args[2]);
		System.exit(0);
	}
}
