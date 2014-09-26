package gov.nih.exon.celltype;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import jdistlib.rng.MersenneTwister;

import gov.nih.exon.reader.QExonData;
import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.QSynchronizedCounter;
import gov.nih.parallel.QThreadPool;
import gov.nih.utils.QSystemUtils;

import static gov.nih.exon.reader.QExonDataReader.readExonDataFileUngrouped;
import static gov.nih.utils.QFileUtils.writeText;
import static gov.nih.utils.QStringUtils.sLn;
import static gov.nih.utils.QStringUtils.sTab;
import static gov.nih.utils.QStringUtils.cTabDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCellTypeBootstrap
{
	private static final int[][] createPermutation(int numIterations, int numCols, int numReps) {
		MersenneTwister random = new MersenneTwister();
		int[][] result = new int[numIterations][numCols * numReps];
		for (int iterNo = 0; iterNo < numIterations; iterNo++) {
			int[] curPermutation = result[iterNo];
			for (int colNo = 0; colNo < numCols; colNo++) {
				int item = random.nextInt(numCols);
				for (int repNo = 0; repNo < numReps; repNo++) {
					curPermutation[repNo*numCols + colNo] = item + repNo*numCols;
				}
			}
		}
		return result;
	}

	public static final double[][] runBootstrap(QExonData exonData, int numIterations, double[] alphas)
	{
		// The parallel pattern is copied from my QeQTLAnalysisHotspot class
		double[][] data = exonData.getData();
		int
			numExons = data.length,
			numCols = data[0].length,
			numAlphas = alphas.length,
			numCellTypes = 3,
			numSubjects = numCols / numCellTypes,
			numThreads = QSystemUtils.kNumCPUCores;
		int[][] permutation = createPermutation(numIterations, numSubjects, numCellTypes);
		double[] newAlphas = new double[numAlphas];
		for (int alphaNo = 0; alphaNo < numAlphas; alphaNo++) {
			double alpha = alphas[alphaNo];
			newAlphas[alphaNo] = 1 - alpha;
		}
		double[][] result = new double[numExons][numAlphas * numCellTypes];
		List<IGenericParallelJob> jobList = new ArrayList<IGenericParallelJob>();
		QSynchronizedCounter
			counter = new QSynchronizedCounter(0, numIterations),
			lock = new QSynchronizedCounter(0, numThreads);
		synchronized(lock) {
			for (int i = 0; i < numThreads; i++) {
				QCellTypeBootstrapJob job = new QCellTypeBootstrapJob(exonData, permutation, result, newAlphas, lock, counter);
				jobList.add(job);
			}
			QThreadPool.mDefaultThreadPool.addAllJobs(jobList);
			while (lock.hasNext()) {
				try	{ lock.wait(); }
				catch (Exception e) {} // Interrupted exception is ignored
				//System.out.println(lock);
			}
		}
		return result;
	}

	public static void writeResult(QExonData exonData, double[][] result, double[] alphas, String saveFile, String delimiter) throws IOException {
		String[] cellTypes = new String[] { "PAX_CL_", "PBMC_CL_", "PBMC_PAX_"};
		int
			numCellTypes = cellTypes.length,
			numAlphas = alphas.length;
		double[] newAlphas = new double[numAlphas];
		for (int alphaNo = 0; alphaNo < numAlphas; alphaNo++) {
			double alpha = alphas[alphaNo];
			newAlphas[alphaNo] = 1 - alpha;
		}

		StringBuilder buf = new StringBuilder();
		buf.append("ExonID");
		buf.append(delimiter);
		buf.append("GeneID");
		for (int cellTypeNo = 0; cellTypeNo < numCellTypes; cellTypeNo++) {
			String cellTypeStr = cellTypes[cellTypeNo];
			for (int alphaNo = 0; alphaNo < numAlphas; alphaNo++) {
				buf.append(delimiter);
				buf.append(cellTypeStr + newAlphas[alphaNo]);
			}
		}
		buf.append(sLn);

		String[] geneIDs = exonData.getGeneIDs();
		int[] exonIDs = exonData.getExonIDs();
		int
			numExons = exonIDs.length,
			numCols = result[0].length;
		for (int exonNo = 0; exonNo < numExons; exonNo++) {
			double[] curResult = result[exonNo];
			buf.append(exonIDs[exonNo]); buf.append(delimiter);
			buf.append(geneIDs[exonNo]);
			for (int colNo = 0; colNo < numCols; colNo++) {
				buf.append(delimiter);
				buf.append(curResult[colNo]);
			}
			buf.append(sLn);
		}
		writeText(buf.toString(), saveFile);
	}

	public static void main(String[] args) {
		try {
			double[] alphas = new double[] { 0.05, 0.01, 0.001, 0.0001 };
			int numIterations = 100000;
			long time1, time2, mem1, mem2;
			String fileName = args[0], ext, baseName;
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
			QExonData data = readExonDataFileUngrouped(fileName, cTabDelimiter, true);
			time2 = System.currentTimeMillis();
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));
			System.out.println("Load time = " + (time2 - time1));
			time1 = System.currentTimeMillis();
			double[][] result = runBootstrap(data, numIterations, alphas);
			writeResult(data, result, alphas, baseName + "_meancutoff_" + ext, sTab);
			time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
