package gov.nih.exon.net;

import qstats.glm.QGLM;
import qstats.glm.QGLMResult;
import qstats.glm.QGLMResultBasic;
import gov.nih.parallel.ISimpleParallelJob;
import gov.nih.table.QTableData;
import gov.nih.utils.QFileUtils;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;

import static qstats.QStatsUtils.calcPOfF;
import static gov.nih.parallel.QParallelUtils.*;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableData;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QGeneRegNetPilot implements ISimpleParallelJob {
	enum EWhichComparison {
		PAX_CL,
		PBMC_CL,
		PBMC_PAX
	}
	private double[][]
		pvals, data;
	private double[][] result;
	private EWhichComparison whichOne;

	private double designMat1[][], designMat2[][], y[], yy[][] = new double[1][];
	private int idx1_start, idx1_stop, idx2_start, idx2_stop;

	public QGeneRegNetPilot(double[][] dt, double[][] rs, double[][] ps, EWhichComparison which) {
		data = dt;
		result = rs;
		pvals = ps;
		whichOne = which;
		designMat1 = new double[70][4];
		designMat2 = new double[70][3];
		yy[0] = y = new double[70];
		for (int i = 0; i < 70; i++)
			designMat1[i][0] = designMat2[i][0] = 1;
		for (int i = 0; i < 35; i++)
			designMat1[i][2] = designMat2[i][2] = 1;
		switch(which) {
			case PAX_CL: idx1_start = 35; idx1_stop = 70; idx2_start = 0; idx2_stop = 35; break;
			case PBMC_CL: idx1_start = 70; idx1_stop = 105; idx2_start = 0; idx2_stop = 35; break;
			case PBMC_PAX: idx1_start = 70; idx1_stop = 105; idx2_start = 35; idx2_stop = 70; break;
		}
	}

	@Override
	public boolean execute(int iterNo) {
		int numGenes = data.length;
		result[iterNo][iterNo] = 1;
		pvals[iterNo][iterNo] = 1;
		double[] x = data[iterNo];
		QGLMResult glm1 = null, glm2 = null;
		for (int i = idx1_start; i < idx1_stop; i++) {
			int j = i - idx1_start;
			designMat1[j][1] = designMat1[j][3] = designMat2[j][1] = x[i];
		}
		for (int i = idx2_start; i < idx2_stop; i++) {
			int j = 35 + i - idx2_start;
			designMat1[j][1] = designMat2[j][1] = x[i];
		}
		for (int geneNo = 0; geneNo < numGenes; geneNo++) {
			if (geneNo == iterNo) continue;
			x = data[geneNo];
			for (int i = idx1_start; i < idx1_stop; i++)
				y[i - idx1_start] = x[i];
			for (int i = idx2_start; i < idx2_stop; i++)
				y[35 + i - idx2_start] = x[i];
			if (glm1 == null || glm2 == null) {
				glm1 = QGLM.multipleRegression(designMat1, y);
				glm2 = QGLM.multipleRegression(designMat2, y);
			} else {
				glm1 = glm1.recalculate(yy);
				glm2 = glm2.recalculate(yy);
			}
			glm1.calcMoreStats();
			glm2.calcMoreStats();
			QGLMResultBasic glm_diff = glm1.calcDiffStatistic(glm2);
			result[geneNo][iterNo] = (float) glm_diff.mRSq;
			pvals[geneNo][iterNo] = 1 - calcPOfF(glm_diff.mF, glm_diff.mDFReg, glm_diff.mDFErr); 
		}
		return true;
	}

	@Override
	public QGeneRegNetPilot clone() {
		return new QGeneRegNetPilot(data, result, pvals, whichOne);
	}

	static final int[] countGraduation(double[][] result, double[] graduation) {
		int len = graduation.length, n = result.length;
		int[] counts = new int[graduation.length];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) continue;
				double val = result[i][j];
				for (int k = 0; k < len; k++)
					if (val >= graduation[k]) counts[k]++;
			}
		}
		return counts;
	}

	public static void main(String[] args) {
		//annotateNet("/Users/joehanesr/Desktop/cc-resid-revised1-gene.txt",
		//		"/Users/joehanesr/Desktop/cc-resid-revised1-gene-diff.txt",
		//		"/Users/joehanesr/Desktop/annotation-gene.txt");
		try {
			String
				fileName = args[0],
				ext, baseName;
			int pos = fileName.lastIndexOf('.');
			if (pos >= 0) {
				ext = fileName.substring(pos + 1);
				baseName = fileName.substring(0, pos);
			} else {
				baseName = fileName;
				ext = "";
			}

			long mem1, mem2, time1, time2; // Benchmark purposes

			mem1 = QSystemUtils.usedMemoryAfterGC();
			QTableData dataTable = readDelimitedFileAsTableData(fileName, cCommaDelimiter, "#", true, true); // num genes x num samples
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));

			double[][] data = dataTable.extractData();
			int
				numSubjects = data[0].length,
				numGenes = data.length;
			System.out.println("Data loaded = " + numGenes + " x " + numSubjects);
			double[][]
				result = new double[numGenes][numGenes],
				pvals = new double[numGenes][numGenes];
			double[] rsq_graduation = new double[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
			int[] rsq_counts;
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));

			System.out.println("PAX-CL");
			time1 = System.currentTimeMillis();
			executeSimpleParallelJob(new QGeneRegNetPilot(data, result, pvals, EWhichComparison.PAX_CL), numGenes);
			time2 = System.currentTimeMillis();
			System.out.println("Regression time = " + (time2 - time1));
			time1 = System.currentTimeMillis();
			QFileUtils.writeText(result, baseName + "-rsq-pax-cl." + ext);
			QFileUtils.writeText(pvals, baseName + "-pvals-pax-cl." + ext);
			time2 = System.currentTimeMillis();
			System.out.println("Saving time = " + (time2 - time1));
			rsq_counts = countGraduation(result, rsq_graduation);
			System.out.println(QStringUtils.toString(rsq_graduation));
			System.out.println(QStringUtils.toString(rsq_counts));

			System.out.println("PBMC-CL");
			time1 = System.currentTimeMillis();
			executeSimpleParallelJob(new QGeneRegNetPilot(data, result, pvals, EWhichComparison.PBMC_CL), numGenes);
			time2 = System.currentTimeMillis();
			System.out.println("Regression time = " + (time2 - time1));
			time1 = System.currentTimeMillis();
			QFileUtils.writeText(result, baseName + "-rsq-pbmc-cl." + ext);
			QFileUtils.writeText(pvals, baseName + "-pvals-pbmc-cl." + ext);
			time2 = System.currentTimeMillis();
			System.out.println("Saving time = " + (time2 - time1));
			rsq_counts = countGraduation(result, rsq_graduation);
			System.out.println(QStringUtils.toString(rsq_graduation));
			System.out.println(QStringUtils.toString(rsq_counts));

			System.out.println("PBMC-PAX");
			time1 = System.currentTimeMillis();
			executeSimpleParallelJob(new QGeneRegNetPilot(data, result, pvals, EWhichComparison.PBMC_PAX), numGenes);
			time2 = System.currentTimeMillis();
			System.out.println("Regression time = " + (time2 - time1));
			time1 = System.currentTimeMillis();
			QFileUtils.writeText(result, baseName + "-rsq-pbmc-pax." + ext);
			QFileUtils.writeText(pvals, baseName + "-pvals-pbmc-pax." + ext);
			time2 = System.currentTimeMillis();
			System.out.println("Saving time = " + (time2 - time1));
			rsq_counts = countGraduation(result, rsq_graduation);
			System.out.println(QStringUtils.toString(rsq_graduation));
			System.out.println(QStringUtils.toString(rsq_counts));

			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
