package gov.nih.exon.coulter;

import jama.QRDecomposition;
import qstats.glm.QGLMResult;
import gov.nih.table.QTableDataExtra;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;

import static java.lang.Math.sqrt;
import static qmath.QMatrixUtils.seq;
import static qstats.QStatsUtils.calcVariance;
import static gov.nih.parallel.QParallelUtils.executeJobSerially;
import static gov.nih.parallel.QParallelUtils.executeSimpleParallelJob;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableDataExtra;

import static gov.nih.exon.coulter.QCoulterDataUtils.saveResult;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kBAP;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kEOP;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kLYP;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kMOP;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kNEP;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kNumResultColumns;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kPLT;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kRBC;
import static gov.nih.exon.coulter.QCoulterDataAnalysis.kWBC;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;
import static gov.nih.utils.QStringUtils.cTabDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterSim {

	static final double[][] analyze(double[][] exonData, double[][] coulterData, double[] rsq, int[] subjectIdx, boolean parallel)
	{
		int
			numExons = exonData.length,
			numSubjects = subjectIdx.length;
		double[][]
			X = new double[numSubjects][kNumResultColumns],
			result = new double[numExons][];

		for (int i = 0; i < numSubjects; i++) {
			double
				curData[] = coulterData[subjectIdx[i]],
				curX[] = X[i],
				nep = curData[kNEP] / 100.0,
				lyp = curData[kLYP] / 100.0,
				mop = curData[kMOP] / 100.0,
				eop = curData[kEOP] / 100.0,
				bap = curData[kBAP] / 100.0;
			curX[0] = 1;
			curX[1] = nep; // Ignored
			curX[1] = lyp;
			curX[2] = mop;
			curX[3] = eop;
			curX[4] = bap;
			curX[5] = curData[kWBC];
			curX[6] = curData[kRBC];
			curX[7] = curData[kPLT];
			/*
			curX[8] = curData[kHGB];
			curX[9] = curData[kHCT];
			curX[10] = curData[kMCV];
			curX[11] = curData[kMCH];
			curX[12] = curData[kMCHC];
			curX[13] = curData[kRDW];
			curX[14] = curData[kMPV];
			//*/
		}
		QCoulterAnalysisJob job = new QCoulterAnalysisJob(exonData, X, result, rsq, subjectIdx);
		if (parallel)
			executeSimpleParallelJob(job, numExons);
		else
			executeJobSerially(job, numExons);
		return result;
	}

	static double[][] simulateRiskFactorsInExprData(double[][] exonData, double[][] coulterData, double[] rsq,
		int[] selectedExons, int[] affectedSubjectIdx, double stdevRatio) {
		int
			numExons = exonData.length,
			numSelectedExons = selectedExons.length,
			numSubjects = exonData[0].length,
			numAffectedSubjects = affectedSubjectIdx.length;
		double newExons[][] = new double[numExons][];
		boolean[] affectedSubjects = new boolean[numSubjects];

		for (int subjectNo = 0; subjectNo < numAffectedSubjects; subjectNo++)
			affectedSubjects[affectedSubjectIdx[subjectNo]] = true;

		for (int exonNo = 0; exonNo < numSelectedExons; exonNo++) {
			int exonIdx = selectedExons[exonNo];
			double
				curExon[] = exonData[exonIdx],
				newCurExon[] = new double[numSubjects],
				curStdev = stdevRatio * sqrt(calcVariance(curExon));
			newExons[exonIdx] = newCurExon;
			for (int subjectNo = 0; subjectNo < numSubjects; subjectNo++)
				if (affectedSubjects[subjectNo])
					newCurExon[subjectNo] = curExon[subjectNo] + curStdev;
		}
		for (int exonNo = 0; exonNo < numExons; exonNo++)
			if (newExons[exonNo] == null)
				newExons[exonNo] = exonData[exonNo];
		return newExons;
	}

	static final double[][] predictCoulterData(double[][] exonData, double[][] results) {
		int
			numSelectedExons = exonData.length,
			numTestSubjects = exonData[0].length,
			subjectIdx[] = new int[numSelectedExons];
		double
			selectedExons[][] = new double[numTestSubjects][numSelectedExons],
			prediction[][] = new double[numTestSubjects][],
			rsq[] = new double[numTestSubjects];

		for (int exonNo = 0; exonNo < numSelectedExons; exonNo++) {
			subjectIdx[exonNo] = exonNo;
			for (int indNo = 0; indNo < numTestSubjects; indNo++)
				selectedExons[indNo][exonNo] = exonData[exonNo][indNo];
		}

		// The test set is pretty small. Let's dispense with the parallel algorithm.
		QCoulterAnalysisJob job = new QCoulterAnalysisJob(selectedExons, results, prediction, rsq, subjectIdx);
		executeJobSerially(job, numTestSubjects);
		//executeSimpleParallelJob(job, numTestSubjects);

		for (int indNo = 0; indNo < numTestSubjects; indNo++) {
			double pred[] = prediction[indNo];
			pred[0] = 1 - (pred[1] + pred[2] + pred[3] + pred[4]); // Recover NEU proportion
		}
		return prediction;
	}

	static final double[] predictRiskFactorFromExonData(double[][] exonData, int[] affectedSubjectIdx) {
		int
			numExons = exonData.length,
			numSubjects = exonData[0].length,
			numAffectedSubjects = affectedSubjectIdx.length;
		double[] exonPVals = new double[numExons];
		double[][]
			affectedSubjectArray = new double[numSubjects][2],
			tempY = new double[1][],
			tempX = new double[1][numSubjects],
			beta = new double[1][];

		for (int subjectNo = 0; subjectNo < numAffectedSubjects; subjectNo++)
			affectedSubjectArray[affectedSubjectIdx[subjectNo]][1] = 1;
		for (int subjectNo = 0; subjectNo < numSubjects; subjectNo++)
			affectedSubjectArray[subjectNo][0] = 1;

		QRDecomposition qr = new QRDecomposition(affectedSubjectArray);
		QGLMResult glm = new QGLMResult(affectedSubjectArray, tempY, qr, false);
		for (int exonNo = 0; exonNo < numExons; exonNo++) {
			tempY[0] = exonData[exonNo];
			beta[0] = new double[2];
			glm.update(tempY, tempX, beta);
			glm.mY = tempY;
			glm.mBetas = beta;
			glm.calcMoreStats();
			glm.calcFBetas();
			glm.calcPOfFBetas();
			exonPVals[exonNo] = glm.mPFBetas[1];
		}
		return exonPVals; // Later on we can slice it or dice it with FDR
	}

	static final double[] predictRiskFactorFromUpdatedCoulterData(double[][] updatedCoulter, int[] affectedSubjectIdx) {
		int
			numComponents = updatedCoulter[0].length,
			numSubjects = updatedCoulter.length,
			numAffectedSubjects = affectedSubjectIdx.length;
		double[] pvals = new double[numComponents];
		double[][]
			affectedSubjectArray = new double[numSubjects][2],
			tempY = new double[1][numSubjects],
			tempX = new double[1][numSubjects],
			beta = new double[1][];
	
		for (int subjectNo = 0; subjectNo < numAffectedSubjects; subjectNo++)
			affectedSubjectArray[affectedSubjectIdx[subjectNo]][1] = 1;
		for (int subjectNo = 0; subjectNo < numSubjects; subjectNo++)
			affectedSubjectArray[subjectNo][0] = 1;
	
		QRDecomposition qr = new QRDecomposition(affectedSubjectArray);
		QGLMResult glm = new QGLMResult(affectedSubjectArray, tempY, qr, false);
		double[] y = tempY[0];
		for (int compNo = 0; compNo < numComponents; compNo++) {
			for (int subjectNo = 0; subjectNo < numSubjects; subjectNo++)
				y[subjectNo] = updatedCoulter[subjectNo][compNo];
			beta[0] = new double[2];
			glm.update(tempY, tempX, beta);
			glm.mY = tempY;
			glm.mBetas = beta;
			glm.calcMoreStats();
			glm.calcFBetas();
			glm.calcPOfFBetas();
			pvals[compNo] = glm.mPFBetas[1];
		}
		return pvals;
	}

	public static void main(String[] args) {
		try {
			long time1, time2, mem1, mem2;
			String
				fileName = args[0],
				coulterFilename = args[1],
				ext, baseName;
			int pos = fileName.lastIndexOf('.');
			if (pos >= 0) {
				ext = fileName.substring(pos + 1);
				baseName = fileName.substring(0, pos);
			} else {
				baseName = fileName;
				ext = "";
			}
			mem1 = QSystemUtils.usedMemoryAfterGC();
			time1 = System.currentTimeMillis();
			QTableDataExtra
				data = readDelimitedFileAsTableDataExtra(fileName, cCommaDelimiter, "#", 1, 1, true),
				coulterTable = readDelimitedFileAsTableDataExtra(coulterFilename, cTabDelimiter, "#", 2, 1, true);
			time2 = System.currentTimeMillis();
			mem2 = QSystemUtils.usedMemoryAfterGC();
			double[][]
				coulterData = coulterTable.getData(),
				exonData = data.getData();
			//double[][] trimmedCoulterData = rearrangeCoulterData(coulterData, testedColumns);
			assert(coulterData != null);

			System.out.println("Memory used = " + (mem2 - mem1));
			System.out.println("Load time = " + (time2 - time1));

			/*
			int
				numSubjects = coulterData.length,
				//testedColumns[] = new int[] { kNEP, kLYP, kMOP, kEOP, kBAP, kWBC, kRBC, kPLT },
				subjectIdx[] = new int[numSubjects];
			double[] rsqResult = new double[exonData.length];
			for (int i = 0; i < numSubjects; i++)
				subjectIdx[i] = i;
			time1 = System.currentTimeMillis();
			double[][] result = analyze(exonData, coulterData, rsqResult, subjectIdx, true);
			time2 = System.currentTimeMillis();
			System.out.println("Analysis time = " + (time2 - time1));
			String[] exonIDs = data.getAllRowNames()[0];
			saveResult(baseName + "-coef." + ext, result, exonIDs);
			System.exit(0);
			//*/

			double[] rsq = seq(0.4, 0.86, 0.025);
			//double[][] xvResult = simpleXV(1000, exonData, coulterData, rsq);
			int[][] exonCount = new int[rsq.length][exonData.length];
			int[] avgExonCount = new int[rsq.length];
			double[][] xvResult = null;// exonBasedXV(1000, exonData, coulterData, exonCount, avgExonCount, rsq);
			saveResult(baseName + "-xv33-15." + ext, xvResult, null);
			saveResult(baseName + "-exoncount33-15." + ext, exonCount, null);
			System.out.println(QStringUtils.toString(rsq));
			System.out.println(QStringUtils.toString(avgExonCount));
			//String[] exonIDs = data.getAllRowNames()[0];
			//saveResult(baseName + "-result." + ext, result, exonIDs);
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
