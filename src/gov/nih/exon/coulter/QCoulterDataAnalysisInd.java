package gov.nih.exon.coulter;

import gov.nih.jama.QRDecomposition;

import java.util.Set;
import java.util.TreeSet;

import qtables.QTableDataExtra;
import gov.nih.parallel.ISimpleParallelJob;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;
import jdistlib.rng.MersenneTwister;

import static qmath.QMatrixUtils.meanColumns;
import static qmath.QMatrixUtils.seq;
import static gov.nih.parallel.QParallelUtils.*;
import static gov.nih.utils.QDataUtils.shortenArray;
import static qutils.QFileUtils.readDelimitedFileAsTableDataExtra;

import static gov.nih.exon.coulter.QCoulterDataAnalysis.*;
import static gov.nih.exon.coulter.QCoulterDataUtils.saveResult;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;
import static gov.nih.utils.QStringUtils.cTabDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterDataAnalysisInd implements ISimpleParallelJob {
	protected double
		mExonData[][],
		mCoulterData[][], // numSubjects x numColumns
		mXVResults[][][], // numRSq x numXVIter x numResult
		mRSq[];
	protected int
		mNumTest,
		mTrainingIdx[],
		mTestIdx[],
		mAvgExonCount[][],
		mExonCount[][][]; // numColumns x numRSq x numExons
	private Set<Integer> mTestSet = new TreeSet<Integer>();
	private MersenneTwister mRandom = new MersenneTwister();

	public QCoulterDataAnalysisInd(double[][] exonData, double[][] coulterData, double[] rsq,
		double[][][] xvResults, int[][][] exonCount, int[][] avgExonCount, int numTestedInd) {
		mExonData = exonData;
		mCoulterData = coulterData;
		mXVResults = xvResults;
		mRSq = rsq;
		mExonCount = exonCount;
		mAvgExonCount = avgExonCount;
		mNumTest = numTestedInd;
		mTestIdx = new int[mNumTest];
		mTrainingIdx = new int[coulterData.length - mNumTest];
	}

	private double[][] fillCache(int[] selectedIdx) {
		int
			numColumns = mCoulterData[0].length,
			numSubjects = selectedIdx.length;
		double
			sumXCache[] = new double[numColumns],
			sumXSqCache[] = new double[numColumns],
			cache[][] = new double[][] { sumXCache, sumXSqCache };
		for (int colNo = 0; colNo < numColumns; colNo++) {
			double sumX = 0, sumXSq = 0;
			for (int subNo = 0; subNo < numSubjects; subNo++) {
				double val = mCoulterData[selectedIdx[subNo]][colNo];
				sumX += val;
				sumXSq += val * val;
			}
			sumXCache[colNo] = sumX;
			sumXSqCache[colNo] = sumXSq;
		}
		return cache;
	}

	private double[][] analyze(int[] selectedIdx) {
		double
			cache[][] = fillCache(selectedIdx),
			sumXCache[] = cache[0],
			sumXSqCache[] = cache[1];
		int
			numColumns = mCoulterData[0].length,
			numSubjects = selectedIdx.length,
			numExons = mExonData.length;
		double[][] result = new double[numExons][numColumns * 3];

		for (int exonNo = 0; exonNo < numExons; exonNo++) {
			double sumY = 0, sumYSq = 0, curExon[] = mExonData[exonNo], curResult[] = result[exonNo];
			for (int subNo = 0; subNo < numSubjects; subNo++) {
				double val = curExon[selectedIdx[subNo]];
				sumY += val;
				sumYSq += val * val;
			}
			double yBar = sumY / numSubjects;

			for (int colNo = 0; colNo < numColumns; colNo++) {
				double sumXY = 0;
				for (int subNo = 0; subNo < numSubjects; subNo++) {
					int idx = selectedIdx[subNo];
					double
						valX = mCoulterData[idx][colNo],
						valY = curExon[idx];
					sumXY += valX * valY;
				}
				double
					sumX = sumXCache[colNo],
					xBar = sumX / numSubjects,
					sXY = sumXY - sumX * yBar,
					sXX = sumXSqCache[colNo] - sumX * xBar,
					beta1 = sXY / sXX,
					beta0 = yBar - beta1 * xBar,
					rsq = (beta1 * sXY) / (sumYSq - sumY * yBar);
				curResult[colNo] = beta0;
				curResult[colNo + numColumns] = beta1;
				curResult[colNo + 2 * numColumns] = rsq;
			}
		}
		return result;
	}

	static final void test(double[][] exonData, double[][] coulterData, double[][] results, int[] testIdx, double[] testResult, int colNo) {
		// results is of numExons x 13 dimension
		int
			numSelectedExons = exonData.length,
			numTestSubjects = testIdx.length,
			numColumns = coulterData[0].length,
			nextCol = colNo + numColumns;
		double[][]
			selectedExons = new double[numTestSubjects][numSelectedExons],
			Xresult = new double[numSelectedExons][2],
			prediction = new double[numTestSubjects][2];

		for (int exonNo = 0; exonNo < numSelectedExons; exonNo++) {
			Xresult[exonNo][0] = results[exonNo][colNo];
			Xresult[exonNo][1] = results[exonNo][nextCol];
			for (int indNo = 0; indNo < numTestSubjects; indNo++) {
				int indIdx = testIdx[indNo];
				selectedExons[indNo][exonNo] = exonData[exonNo][indIdx];
			}
		}
		// Analyze
		QRDecomposition qr = new QRDecomposition(Xresult);
		double[][]
			tempY = new double[1][],
			tempX = new double[1][numSelectedExons],
			tempBeta = new double[1][];
		for (int indNo = 0; indNo < numTestSubjects; indNo++) {
			tempY[0] = selectedExons[indNo];
			tempBeta[0] = prediction[indNo];
			qr.solveTranspose(tempY, tempX, tempBeta);
		}

		double sumSq = 0, sumObs = 0;
		for (int indNo = 0; indNo < numTestSubjects; indNo++) {
			double
				obs = coulterData[indNo][colNo],
				pred = prediction[indNo][1];
			if (pred < 0)
				pred = 0;
			double diff = pred - obs;
			sumSq += diff * diff;
			sumObs += obs;
		}
		sumSq = Math.sqrt(sumSq / (numTestSubjects - 1));
		sumObs = sumObs / numTestSubjects;
		testResult[colNo] += sumSq / sumObs;
		testResult[nextCol] += sumObs;
	}

	@Override
	public boolean execute(int iterNo) {
		int
			numSubjects = mCoulterData.length,
			numColumns = mCoulterData[0].length,
			numExons = mExonData.length,
			numRSq = mRSq.length;
		mTestSet.clear();
		while (mTestSet.size() < mNumTest)
			mTestSet.add(mRandom.nextInt(numSubjects));
		for (int subjectNo = 0, idx = 0, idx2 = 0; subjectNo < numSubjects; subjectNo++)
			if (!mTestSet.contains(subjectNo))
				mTrainingIdx[idx++] = subjectNo;
			else
				mTestIdx[idx2++] = subjectNo;
		double[][] trainingResult = analyze(mTrainingIdx);
		double[][] testResult = new double[numRSq][numColumns * 2];
		for (int rsqNo = 0; rsqNo < numRSq; rsqNo++)
			mXVResults[rsqNo][iterNo] = testResult[rsqNo];

		for (int colNo = 0; colNo < numColumns; colNo++) {
			int rsqIdx = colNo + 2 * numColumns;
			int[][] curColExonCount = mExonCount[colNo];
			for (int rsqNo = 0; rsqNo < numRSq; rsqNo++) {
				double
					curRSq = mRSq[rsqNo],
					filteredResult[][] = new double[numExons][],
					filteredExons[][] = new double[numExons][];
				int
					exonIdx = 0,
					curExonCount[] = curColExonCount[rsqNo];
				for (int exonNo = 0; exonNo < numExons; exonNo++) {
					double[]
						curExon = mExonData[exonNo],
						curResult = trainingResult[exonNo];
					if (curResult[rsqIdx] >= curRSq) {
						filteredExons[exonIdx] = curExon;
						filteredResult[exonIdx] = curResult;
						exonIdx++;
						synchronized(curExonCount) {
							curExonCount[exonNo]++;
						}
					}
				}
				if (exonIdx > 0) {
					synchronized(mAvgExonCount) {
						mAvgExonCount[colNo][rsqNo] += exonIdx;
					}
					filteredExons = shortenArray(filteredExons, exonIdx);
					filteredResult = shortenArray(filteredResult, exonIdx);
					test(filteredExons, mCoulterData, filteredResult, mTestIdx, testResult[rsqNo], colNo);
				}
				//else
				//	throw new RuntimeException();
			}
		}
		iterNo++;
		if (iterNo % 100 == 0)
			System.out.println("Iteration " + iterNo);
		return true;
	}

	@Override
	public QCoulterDataAnalysisInd clone() {
		QCoulterDataAnalysisInd job = new QCoulterDataAnalysisInd(mExonData, mCoulterData,
			mRSq, mXVResults, mExonCount, mAvgExonCount, mNumTest);
		return job;
	}

	static final double[][] exonBasedXV(int numIter, double[][] exonData, double[][] coulterData,
		int[][][] exonCount, int[][] avgExonCount, double[] rsq, double kTestProportion) {
		int
			numRSq = rsq.length,
			numTestedInd = (int) Math.floor(coulterData.length * kTestProportion);
		double[][][] xvResults = new double[numRSq][numIter][];
		double[][] result = new double[numRSq][];
		QCoulterDataAnalysisInd job = new QCoulterDataAnalysisInd(exonData, coulterData, rsq, xvResults,
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

			System.out.println("Memory used = " + (mem2 - mem1));
			System.out.println("Load time = " + (time2 - time1));
			coulterData = rearrangeCoulterData(coulterData, new int[] { kNEP, kLYP, kMOP, kEOP, kBAP, kWBC, kRBC, kPLT});

			int numCoulterCols = coulterData[0].length;
			double[] rsq = seq(0.2, 0.6, 0.025);
			//double[][] xvResult = simpleXV(1000, exonData, coulterData, rsq);
			int[][][] exonCount = new int[numCoulterCols][rsq.length][exonData.length];
			int[][] avgExonCount = new int[numCoulterCols][rsq.length];
			double[][] xvResult = exonBasedXV(1000, exonData, coulterData, exonCount, avgExonCount, rsq, 1.0/3);
			saveResult(baseName + "-xvind-33-8." + ext, xvResult, null);
			for (int i = 0; i < numCoulterCols; i++)
				saveResult(baseName + "-exoncountind-"+i+"-33-8." + ext, exonCount[i], null);
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
