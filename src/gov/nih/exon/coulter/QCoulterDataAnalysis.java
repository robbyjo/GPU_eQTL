package gov.nih.exon.coulter;

import gov.nih.table.QTableDataExtra;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;

import static qmath.QMatrixUtils.meanColumns;
import static qmath.QMatrixUtils.seq;
import static gov.nih.parallel.QParallelUtils.executeJobSerially;
import static gov.nih.parallel.QParallelUtils.executeSimpleParallelJob;
import static gov.nih.utils.QFileUtils.*;

import static gov.nih.exon.coulter.QCoulterDataUtils.saveResult;
import static gov.nih.utils.QStringUtils.cTabDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterDataAnalysis
{
	public static final int
		kWBC = 0,
		kRBC = 1,
		kHGB = 2,
		kHCT = 3,
		kMCV = 4,
		kMCH = 5,
		kMCHC = 6,
		kRDW = 7,
		kPLT = 8,
		kMPV = 9,
		kNEP = 10,
		kLYP = 11,
		kMOP = 12,
		kEOP = 13,
		kBAP = 14,
		kNEC = 15,
		kLYC = 16,
		kMOC = 17,
		kEOC = 18,
		kBAC = 19,
		kNumResultColumns = 8;

	public static final double
		kDefaultPenalty = 5,
		kDetectionThreshold = 1e-6,
		kDefaultLODThreshold = 2.5,
		kDefaultCriticalValue = 0.0001,
		kTestProportion = 1.0 / 3;

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

	static final double[] test(double[][] exonData, double[][] coulterData, double[][] results, int[] testIdx) {
		// results is of numExons x 13 dimension
		int
			numSelectedExons = exonData.length,
			numTestSubjects = testIdx.length,
			subjectIdx[] = new int[numSelectedExons];
		double
			selectedExons[][] = new double[numTestSubjects][numSelectedExons],
			prediction[][] = new double[numTestSubjects][],
			rsq[] = new double[numTestSubjects];

		for (int exonNo = 0; exonNo < numSelectedExons; exonNo++) {
			subjectIdx[exonNo] = exonNo;
			for (int indNo = 0; indNo < numTestSubjects; indNo++) {
				int indIdx = testIdx[indNo];
				selectedExons[indNo][exonNo] = exonData[exonNo][indIdx];
			}
		}

		// The test set is pretty small. Let's dispense with the parallel algorithm.
		QCoulterAnalysisJob job = new QCoulterAnalysisJob(selectedExons, results, prediction, rsq, subjectIdx);
		executeJobSerially(job, numTestSubjects);
		//executeSimpleParallelJob(job, numTestSubjects);

		int[] testedCol = new int[] { 0, 0, 0, 0, 0, kWBC, kRBC, kPLT}; //, kHGB, kHCT, kMCV, kMCH, kMCHC, kRDW, kMPV };
		double[] summaryStats = new double[kNumResultColumns * 3];
		double[] meanCols = meanColumns(coulterData);
		for (int colNo = 0; colNo < 5; colNo++)
			meanCols[kNEP + colNo] = meanCols[kNEP + colNo] / 100;
		for (int indNo = 0; indNo < numTestSubjects; indNo++) {
			double
				obsData[] = coulterData[testIdx[indNo]],
				pred[] = prediction[indNo],
				curObs,
				curPred,
				diff;
			pred[0] = 1 - (pred[1] + pred[2] + pred[3] + pred[4]);
			for (int colNo = 0; colNo < 5; colNo++) {
				curObs = obsData[kNEP + colNo] / 100.0;
				curPred = pred[colNo];
				if (curPred < 0)
					curPred = 0;
				diff = (curPred - curObs);
				summaryStats[colNo] += diff * diff;
				diff = (meanCols[kNEP + colNo] - curObs);
				summaryStats[colNo + kNumResultColumns] += curObs; //diff*diff;
			}
			//*
			for (int colNo = 5; colNo < kNumResultColumns; colNo++)
			{
				curObs = obsData[testedCol[colNo]];
				curPred = pred[colNo];
				if (curPred < 0)
					curPred = 0;
				diff = (curPred - curObs);
				summaryStats[colNo] += diff * diff;
				diff = (meanCols[testedCol[colNo]] - curObs);
				summaryStats[colNo + kNumResultColumns] += curObs; //diff*diff;
			}
			//*/
		}

		/*
		for (int colNo = 0; colNo < kNumResultColumns; colNo++) {
			double meanObs = summaryStats[colNo + kNumResultColumns];// / numTestSubjects;
			summaryStats[colNo + kNumResultColumns] = summaryStats[colNo] / meanObs;
			summaryStats[colNo] = Math.sqrt(summaryStats[colNo] / (numTestSubjects));
			summaryStats[colNo + 2 * kNumResultColumns] = meanObs;
		}
		//*/
		for (int colNo = 0; colNo < kNumResultColumns; colNo++) {
			double meanObs = summaryStats[colNo + kNumResultColumns] / numTestSubjects;
			summaryStats[colNo] = Math.sqrt(summaryStats[colNo] / (numTestSubjects));
			summaryStats[colNo + kNumResultColumns] = summaryStats[colNo] / meanObs;
			summaryStats[colNo + 2 * kNumResultColumns] = meanObs;
		}
		return summaryStats;
	}

	static final double[] doCrossValidation(int numCVs, double[][] filteredExons, double[][] coulterData, double testProportion) {
		int numTests = (int) Math.round(coulterData.length * testProportion);
		double[][] xvResults = new double[numCVs][];
		QCoulterCrossValidationJob job = new QCoulterCrossValidationJob(filteredExons, coulterData, xvResults, numTests);
		if (filteredExons.length > 100)
			executeSimpleParallelJob(job, numCVs);
		else
			executeJobSerially(job, numCVs);
		return meanColumns(xvResults);
	}

	static final double[][] rearrangeCoulterData(double[][] coulterData, int[] inclCols) {
		int
			numRows = coulterData.length,
			numInclCols = inclCols.length;
		double[][] X = new double[numRows][numInclCols];

		for (int i = 0; i < numRows; i++) {
			double
				curData[] = coulterData[i],
				curX[] = X[i];
			for (int j = 0; j < numInclCols; j++) {
				int colNo = inclCols[j];
				// Normalize cell type properties
				curX[j] = colNo >= kNEP && colNo <= kBAP ? curData[colNo] / 100.0 : curData[colNo];
			}
		}
		return X;
	}

	static final double[][] simpleXV(int numIter, double[][] exonData, double[][] coulterData, double[] rsqResult, double[] rsq) {
		int
			numRSq = rsq.length,
			numSubjects = coulterData.length,
			//testedColumns[] = new int[] { kNEP, kLYP, kMOP, kEOP, kBAP, kWBC, kRBC, kPLT },
			subjectIdx[] = new int[numSubjects];
		for (int i = 0; i < numSubjects; i++)
			subjectIdx[i] = i;

		long time1, time2;
		time1 = System.currentTimeMillis();
		double[][] result = analyze(exonData, coulterData, rsqResult, subjectIdx, true);
		time2 = System.currentTimeMillis();
		System.out.println("Analysis time = " + (time2 - time1));
		time1 = System.currentTimeMillis();
		double[][] xvResult = new double[numRSq][kNumResultColumns * 4 + 2];
		for (int rsqNo = 0; rsqNo < numRSq; rsqNo++) {
			double
				curRSq = rsq[rsqNo],
				curXV[] = xvResult[rsqNo];
			curXV[0] = curRSq;
			ICoulterResultFilter filter = new QCoulterRSqResultFilter(curRSq);
			double[][] filteredExons = filter.filter(exonData, result);
			curXV[1] = filteredExons.length;

			double[][] filteredResults = filter.filter(result, result);
			double[] testResult = test(filteredExons, coulterData, filteredResults, subjectIdx);
			System.arraycopy(testResult, 0, curXV, 2, 2 * kNumResultColumns);
			//System.out.println(QStringUtils.toString(testResult));
			long t1 = System.currentTimeMillis();
			double[] summary = doCrossValidation(numIter, filteredExons, coulterData, kTestProportion);
			long t2 = System.currentTimeMillis();
			System.out.println("Time for XV-ing R^2 = " + curRSq + ": " + (t2 - t1));
			System.arraycopy(summary, 0, curXV, 2 * kNumResultColumns + 2, 2 * kNumResultColumns);
		}
		time2 = System.currentTimeMillis();
		System.out.println("1000x cross validation time = " + (time2 - time1));
		return xvResult;
	}

	static final double[][] exonBasedXV(int numIter, double[][] exonData, double[][] coulterData,
		int[][] exonCount, int[] avgExonCount, double[] rsq) {
		int
			numRSq = rsq.length,
			numTestedInd = (int) Math.floor(coulterData.length * kTestProportion);
		double[][][] xvResults = new double[numRSq][numIter][];
		double[][] result = new double[numRSq][];
		QCoulterCrossValidationByExonJob job = new QCoulterCrossValidationByExonJob(exonData, coulterData, rsq, xvResults,
			exonCount, avgExonCount, numTestedInd);
		long time1, time2;
		time1 = System.currentTimeMillis();
		executeSimpleParallelJob(job, numIter);
		for (int rsqNo = 0; rsqNo < numRSq; rsqNo++)
			result[rsqNo] = meanColumns(xvResults[rsqNo]);
		time2 = System.currentTimeMillis();
		System.out.println("1000x cross validation time (by exons) = " + (time2 - time1));
		xvResults = null;
		return result;
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
				data = readDelimitedFileAsTableDataExtra(fileName, cTabDelimiter, "#", 1, 1, true),
				coulterTable = readDelimitedFileAsTableDataExtra(coulterFilename, cTabDelimiter, "#", 2, 1, true);
			time2 = System.currentTimeMillis();
			mem2 = QSystemUtils.usedMemoryAfterGC();
			double[][]
				coulterData = coulterTable.getData(),
				exonData = data.getData();
			//double[][] trimmedCoulterData = rearrangeCoulterData(coulterData, testedColumns);

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
			double[][] xvResult = exonBasedXV(1000, exonData, coulterData, exonCount, avgExonCount, rsq);
			saveResult(baseName + "-xv33-8." + ext, xvResult, null);
			saveResult(baseName + "-exoncount33-8." + ext, exonCount, null);
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
