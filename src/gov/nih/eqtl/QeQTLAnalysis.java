/*
 * Roby Joehanes
 * 
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package gov.nih.eqtl;

import java.io.File;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.BitSet;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import net.sourceforge.jdistlib.T;
import gov.nih.eqtl.datastructure.QGeneticSNPData;
import gov.nih.eqtl.datastructure.QSNPData;
import gov.nih.eqtl.datastructure.QSNPDataInt;
import gov.nih.eqtl.io.QPlinkLoader;
import gov.nih.eqtl.io.QCovariateTable;
import gov.nih.eqtl.io.QBinaryMatrixCache;
import gov.nih.eqtl.io.QDelimitedMatrixSource;
import gov.nih.eqtl.io.QLocalPatternImputedSource;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.eqtl.io.QMissingnessReport;
import gov.nih.eqtl.io.QMissingnessScan;
import gov.nih.eqtl.io.QPolicyMatrixSource;
import gov.nih.eqtl.io.QSampleAlignment;
import gov.nih.eqtl.io.QVariantMatrixSource;
import gov.nih.eqtl.datastructure.QGeneExpressionData;
import gov.nih.eqtl.datastructure.QSNPDataReal;
import gov.nih.jama.QRDecomposition;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuPrecision;
import gov.nih.gpu.GpuRuntime;
import gov.nih.gpu.GpuTuning;
import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.IJobOwner;
import gov.nih.parallel.QSynchronizedCounter;
import gov.nih.table.QTableData;
import gov.nih.utils.EPlatform;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;
import gov.nih.utils.matrix.EMultiplicationMode;
import gov.nih.utils.matrix.QMatrixUtils;
import static java.lang.Math.*;
import static gov.nih.utils.QStatsUtils.calcStdDevAndStandardize;
import static gov.nih.utils.QStringUtils.cCommonDelimiter;
import static gov.nih.utils.QStringUtils.sLn;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QeQTLAnalysis implements IJobOwner
{
	static ExecutorService threadPool;
	static boolean DEBUG = false;
	static boolean simplifyResult = false, rsqOnly = false;
	static final double
		kMB = 1024 * 1024,
		kGB = 1024 * kMB;

	static final int
		kExitCodeNormal = 0,
		kExitCodeErrorInvalidParam = -1,
		kExitCodeErrorPlatformNot64Bit = -2,
		kExitCodeErrorInitOpenCLFailure = -3,
		kExitCodeError64bitGPUNotFound = -4,
		kExitCodeErrorCantLoadFile = -5,
		kExitCodeErrorNotEnoughMem = -6,
		kExitCodeErrorWrongCovarSpec = -7,
		kExitCodeErrorCovarNotFound = -8,
		kExitCodeErrorNumIndGenoCovarNotMatch = -9,
		kExitCodeErrorNumIndGenoExprNotMatch = -10,
		kExitCodeErrorCovarMissingValues = -11,
		kExitCodeErrorGenoMissingValues = -12,
		kExitCodeErrorAnalysisFailure = -13;

	// Global variables
	static GpuRuntime GPU_RUNTIME;
	static GpuContext[] mContexts = null;
	static GpuPrecision gpuPrecision = GpuPrecision.FP64;
	static QResidualizationMode residualizationMode = QResidualizationMode.AUTO;
	static QeQTLAnalysisConfig config = null;
	static QeQTLProfiler profiler = new QeQTLProfiler(false);

	public static final String eqtlCat =
		"#define CPY 4" + sLn +
		"#define AS(i, j) As[j + i * (BLOCK_SIZE*CPY+1)]" + sLn +
		"#define BS(i, j) Bs[j + i * (BLOCK_SIZE*CPY+1)]" + sLn +
		"__kernel void eqtlCat(__global DATATYPE* C, __global DATATYPE* A, __global DATATYPE* B, __local DATATYPE* As, __local DATATYPE* Bs, const int wA, const int wB) {" +
		"int bx = get_group_id(0), by = get_group_id(1), tx = get_local_id(0), ty = get_local_id(1)," +
		"aBegin = wA * BLOCK_SIZE * by, aEnd = aBegin + wA, bStep = BLOCK_SIZE * wB, aOffset = wA * ty + tx, bOffset = wB * ty + tx, ssq = (tx % COL_PER_EQTL) == (COL_PER_EQTL-1);" +
		"DATATYPE cval = 0, temp;" + sLn +
		"#pragma unroll" + sLn +
		"for (int a = aBegin, b = BLOCK_SIZE * bx; a < aEnd;)  {" + sLn +
		"	#pragma unroll" + sLn +
		"	for (int k = 0; k < CPY; ++k) {" + sLn +
		"		AS(ty,tx+k*BLOCK_SIZE) = A[a + aOffset];" + sLn +
		"		BS(tx,ty+k*BLOCK_SIZE) = B[b + bOffset]; a += BLOCK_SIZE; b += bStep;" + sLn +
		"	}" + sLn +
		"	barrier(CLK_LOCAL_MEM_FENCE);" + sLn +
		"	#pragma unroll" + sLn +
		"	for (int k = 0; k < BLOCK_SIZE*CPY; ++k) { temp = AS(ty,k); cval += (ssq ? temp : 1) * temp * BS(tx,k); } " + sLn +
		"	barrier(CLK_LOCAL_MEM_FENCE); }" + sLn +
		"C[BLOCK_SIZE * (wB * by + bx) + bOffset] = cval; }";

	public static final String eqtlReal =
		"#define CPY 4" + sLn +
		"#define AS(i, j) As[j + i * (BLOCK_SIZE*CPY+1)]" + sLn +
		"#define BS(i, j) Bs[j + i * (BLOCK_SIZE*CPY+1)]" + sLn +
		"__kernel void eqtlReal(__global DATATYPE* C, __global DATATYPE* A, __global DATATYPE* B, __local DATATYPE* As, __local DATATYPE* Bs, const int wA, const int wB) {" +
		"int bx = get_group_id(0), by = get_group_id(1), tx = get_local_id(0), ty = get_local_id(1)," +
		"aBegin = wA * BLOCK_SIZE * by, aEnd = aBegin + wA, bStep = BLOCK_SIZE * wB, aOffset = wA * ty + tx, bOffset = wB * ty + tx;" +
		"DATATYPE cval = 0;" + sLn +
		"#pragma unroll" + sLn +
		"for (int a = aBegin, b = BLOCK_SIZE * bx; a < aEnd;)  {" + sLn +
		"	#pragma unroll" + sLn +
		"	for (int k = 0; k < CPY; ++k) {" + sLn +
		"		AS(ty,tx+k*BLOCK_SIZE) = A[a + aOffset];" + sLn +
		"		BS(tx,ty+k*BLOCK_SIZE) = B[b + bOffset]; a += BLOCK_SIZE; b += bStep;" + sLn +
		"	}" + sLn +
		"	barrier(CLK_LOCAL_MEM_FENCE);" + sLn +
		"	#pragma unroll" + sLn +
		"	for (int k = 0; k < BLOCK_SIZE*CPY; ++k) cval += AS(ty,k) * BS(tx,k);" + sLn +
		"	barrier(CLK_LOCAL_MEM_FENCE); }" + sLn +
		"C[BLOCK_SIZE * (wB * by + bx) + bOffset] = cval / N_MIN_1; }";

	static final int roundUpNearestMultiple(int number, int multiple) {
		int rem = number % multiple;
		if (rem > 0)
			number = ((number / multiple) * multiple) + multiple;
		return number;
	}

	static final int roundDownNearestMultiple(int number, int multiple) {
		int rem = number % multiple;
		return number - rem;
	}

	public void eSNPAnalysis(QGeneticSNPData popn, QGeneExpressionData expDataTbl, double[][] covarQ, double rsq0,
		int dfo, int dfe, boolean isAdditive)
	{
		eSNPAnalysisInternal(popn, expDataTbl, covarQ, null, null, rsq0, dfo, dfe, isAdditive);
	}

	public void eSNPAnalysisPrepared(QGeneticSNPData popn, QGeneExpressionData expDataTbl,
		double[] snpDataSD, double[] expDataSD, double rsq0, int dfo, int dfe, boolean isAdditive)
	{
		if (snpDataSD == null || expDataSD == null)
			throw new IllegalArgumentException("Prepared standard deviations must not be null");
		eSNPAnalysisInternal(popn, expDataTbl, null, expDataSD, snpDataSD,
			rsq0, dfo, dfe, isAdditive);
	}

	private void eSNPAnalysisInternal(QGeneticSNPData popn, QGeneExpressionData expDataTbl,
		double[][] covarQ, double[] preparedExpDataSD, double[] preparedSnpDataSD, double rsq0,
		int dfo, int dfe, boolean isAdditive)
	{
		//long mem1 = QSystemUtils.usedMemoryAfterGC();
		double[][] expData = expDataTbl.getData();
		int
			numSNPs = popn.getNumSNPs(),
			numInds = popn.getNumIndividuals(),
			numETraits = expData.length;
		if (numInds != expData[0].length)
			throw new RuntimeException();
		boolean isCategoricalSNP = popn.getSNPs().get(0) instanceof QSNPDataInt;
		if (isCategoricalSNP)
			throw new RuntimeException("This branch is broken. Sorry.");

		double[][] snpData = new double[numSNPs][];
		List<QSNPData> snps = popn.getSNPs();
		for (int i = 0; i < numSNPs; i++)
			snpData[i] = snps.get(i).getSNPValues();

		long time1, time2;
		// Clean the expression data to a residual
		if (preparedExpDataSD == null && covarQ != null) {
			System.out.println("Taking residuals from expression data.");
			time1 = System.currentTimeMillis();
			double[][] resid = QMatrixUtils.parallelMatrixMultiplication(expData, covarQ, null, 1, numETraits, numInds, EMultiplicationMode.XMinusXYYt);
			System.arraycopy(resid, 0, expData, 0, resid.length);
			time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			resid = null;
			QSystemUtils.runGCAggressively();
			System.out.println("Taking residuals from SNP data.");
			time1 = System.currentTimeMillis();
			resid = QMatrixUtils.parallelMatrixMultiplication(snpData, covarQ, null, 1, numSNPs, numInds, EMultiplicationMode.XMinusXYYt);
			for (int i = 0; i < numSNPs; i++)
				System.arraycopy(resid[i], 0, snpData[i], 0, numInds);
			time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			resid = null;
			QSystemUtils.runGCAggressively();
		}

		// Compute SD
		time1 = System.currentTimeMillis();
		double[] expDataSD = preparedExpDataSD, snpDataSD = preparedSnpDataSD;
		if (expDataSD == null) {
			expDataSD = new double[numETraits];
			snpDataSD = new double[numSNPs];
			for (int i = 0; i < numETraits; i++)
				expDataSD[i] = calcStdDevAndStandardize(expData[i]);
			for (int i = 0; i < numSNPs; i++)
				snpDataSD[i] = calcStdDevAndStandardize(snpData[i]);
		} else if (expDataSD.length != numETraits || snpDataSD.length != numSNPs) {
			throw new IllegalArgumentException("Prepared standard-deviation counts do not match the matrices");
		}
		time2 = System.currentTimeMillis();
		System.out.println("Standardizing time = " + (time2 - time1));

		int colPereQTL = isCategoricalSNP ? (isAdditive ? 3 : 4) : 1;

		time1 = System.currentTimeMillis();
		final int localBlockSize = 16;
		String
			program,
			kernelName,
			hdr = kernelHeader(localBlockSize, gpuPrecision);

		if (isCategoricalSNP) {
			program = eqtlCat;
			kernelName = "eqtlCat";
			hdr += "#define COL_PER_EQTL " + colPereQTL + sLn;
		} else {
			program = eqtlReal;
			kernelName = "eqtlReal";
			hdr += "#define N_MIN_1 " + (numInds - 1) + sLn;
		}

		GpuContextPool gpuContextsPool = new GpuContextPool(mContexts);
		List<GpuContext> gpuContexts = gpuContextsPool.getAllContexts();
		int numDevices = gpuContexts.size();
		try {
			for (int i = 0; i < numDevices; i++)
			{
				GpuContext ctx = gpuContexts.get(i);
				ctx.compileKernel(hdr + program, kernelName, gpuPrecision);
			}
		} catch (RuntimeException e) {
			gpuContextsPool.close();
			throw e;
		}
		time2 = System.currentTimeMillis();
		System.out.println("Kernel compile time = " + (time2 - time1));

		int globalBlockSize = config.getBlockSize();
		int numThreads = min(config.getNumThreads(), Runtime.getRuntime().availableProcessors() + 1);
		System.out.println("Num threads = " + numThreads);
		threadPool = Executors.newFixedThreadPool(numThreads);

		//time1 = System.currentTimeMillis();
		int
			numETraitsPerBlock = globalBlockSize,
			numESNPsPerBlock = globalBlockSize / colPereQTL,
			nrow = roundUpNearestMultiple(numInds, localBlockSize);
		//int[] snpCounts = new int[numSNPs];
		//int[][] eTraitIndices = new int[numSNPs][numStoredETraitsPerESNP];
		//double[][][] snpResults = new double[numSNPs][numStoredETraitsPerESNP][];
		//QeQTLSNPAnalysisResults results = new QeQTLSNPAnalysisResults(snpCounts, eTraitIndices, snpResults);

		//time2 = System.currentTimeMillis();
		//System.out.println("Setup time = " + (time2 - time1));

		Writer fw = null;
		try {
			fw = new PrintWriter(new FileOutputStream(config.getOutputFilename()), true);
			if (rsqOnly) {
				fw.write("Rs_ID,ProbesetID,RSq,Dir");
			} else {
				fw.write("Rs_ID,ProbesetID,RSq,Fx,T,log10P");
			}
			fw.write(sLn);
		} catch (Exception e) {
			threadPool.shutdownNow();
			gpuContextsPool.close();
			throw new RuntimeException("Cannot open the analysis output file", e);
		}

		try {
			int numIters = (int) ceil(numSNPs * 1.0 / numESNPsPerBlock);
			QSynchronizedCounter counter = new QSynchronizedCounter(0, numIters);
			List<Future<?>> jobs = new ArrayList<Future<?>>(numThreads);
			if (isCategoricalSNP) {
				//for (int i = 0; i < numThreads; i++)
				//jobs.add(threadPool.submit(new QeQTLSNPJobCat(popn, expData, covarQ, gpuContextsPool, results,
				//	numETraitsPerBlock, numESNPsPerBlock, localBlockSize, rsq0, isAdditive, counter)));
			} else {
				for (int i = 0; i < numThreads; i++)
					jobs.add(threadPool.submit(new QeQTLSNPJobReal(popn, expDataTbl, expDataSD, snpDataSD, gpuContextsPool,
						numETraitsPerBlock, numESNPsPerBlock, localBlockSize, dfo, dfe, rsq0,
						isAdditive, gpuPrecision, fw, counter)));
			}
			threadPool.shutdown();
			for (Future<?> job : jobs) {
				job.get();
			}
			if (!threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS)) {
				throw new RuntimeException("Timed out while waiting for eQTL workers");
			}
			fw.flush();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			threadPool.shutdownNow();
			throw new RuntimeException("Interrupted while waiting for eQTL workers", e);
		} catch (Exception e) {
			threadPool.shutdownNow();
			throw new RuntimeException("eQTL worker failed", e);
		} finally {
			gpuContextsPool.close();
			if (fw != null)
				try {
					fw.close();
				} catch (IOException e) {
				}
		}
	}

	public void eSNPAnalysisStreamed(QBinaryMatrixCache genotypeCache,
		QBinaryMatrixCache expressionCache, double rsq0, int dfo, int dfe,
		int genotypeRowsPerBlock, int expressionRowsPerBlock)
	{
		eSNPAnalysisStreamed(genotypeCache, expressionCache, rsq0, dfo, dfe,
			genotypeRowsPerBlock, expressionRowsPerBlock,
			Path.of(config.getOutputFilename()),
			config.getCheckpointDirectory() == null
				? Path.of(config.getOutputFilename() + ".checkpoint")
				: Path.of(config.getCheckpointDirectory()),
			config.getResume(), config.getKeepCheckpoints(), true, false);
	}

	private void eSNPAnalysisStreamed(QBinaryMatrixCache genotypeCache,
		QBinaryMatrixCache expressionCache, double rsq0, int dfo, int dfe,
		int genotypeRowsPerBlock, int expressionRowsPerBlock,
		Path outputPath, Path checkpointDirectory, boolean resume,
		boolean keepCheckpoints, boolean closeGpuContexts, boolean includeSampleStatistics)
	{
		final int localBlockSize = 16;
		int numInds = genotypeCache.sampleCount();
		if (expressionCache.sampleCount() != numInds)
			throw new IllegalArgumentException("Prepared cache sample counts differ");
		int totalBlocks = (int) ceil(genotypeCache.rowCount() * 1.0 / genotypeRowsPerBlock);
		String analysisSignature = QBinaryMatrixCache.analysisSignature(genotypeCache, expressionCache,
			genotypeRowsPerBlock, expressionRowsPerBlock, dfo, dfe, rsq0,
			simplifyResult, rsqOnly, gpuPrecision, includeSampleStatistics);
		QAnalysisCheckpoint checkpoint;
		try {
			checkpoint = QAnalysisCheckpoint.open(checkpointDirectory, analysisSignature, totalBlocks,
				resume, keepCheckpoints);
		} catch (IOException e) {
			throw new RuntimeException("Cannot initialize analysis checkpoint", e);
		}
		int completedAtStart = checkpoint.completedCount();
		if (completedAtStart > 0)
			System.out.println("Completed checkpoint blocks found: " + completedAtStart + " / " + totalBlocks);
		String outputHeader = rsqOnly ? "Rs_ID,ProbesetID,RSq,Dir" : "Rs_ID,ProbesetID,RSq,Fx,T,log10P";
		if (includeSampleStatistics)
			outputHeader += ",N,DF";
		if (completedAtStart == totalBlocks) {
			try {
				long assemblyStarted = profiler.start();
				checkpoint.assemble(outputPath, outputHeader);
				profiler.record(QeQTLProfiler.Phase.OUTPUT_ASSEMBLY, assemblyStarted, totalBlocks, 0);
				return;
			} catch (IOException e) {
				throw new RuntimeException("Cannot assemble completed checkpoint", e);
			}
		}

		String header = kernelHeader(localBlockSize, gpuPrecision)
			+ "#define N_MIN_1 " + (numInds - 1) + sLn;

		GpuContextPool contextPool = closeGpuContexts
			? new GpuContextPool(mContexts) : GpuContextPool.borrowed(mContexts);
		try {
			long compileStarted = profiler.start();
			for (GpuContext context : contextPool.getAllContexts()) {
				context.setProfilingEnabled(profiler.isEnabled());
				context.compileKernel(header + eqtlReal, "eqtlReal", gpuPrecision);
			}
			profiler.record(QeQTLProfiler.Phase.KERNEL_COMPILE, compileStarted,
				contextPool.getAllContexts().size(), 0);
		} catch (RuntimeException e) {
			contextPool.close();
			throw e;
		}

		int requestedThreads = config.getNumThreads();
		int maxThreads = Runtime.getRuntime().availableProcessors() + 1;
		int numThreads = min(requestedThreads, maxThreads);
		System.out.println("Num threads = " + numThreads);
		threadPool = Executors.newFixedThreadPool(numThreads);
		try {
			Deque<Future<?>> pending = new ArrayDeque<Future<?>>();
			for (int blockNumber = 0; blockNumber < totalBlocks; blockNumber++) {
				if (checkpoint.isComplete(blockNumber))
					continue;
				long rowOffset = (long) blockNumber * genotypeRowsPerBlock;
				long cacheReadStarted = profiler.start();
				QeQTLPreprocessor.PreparedBlock prepared = genotypeCache.readBlock(rowOffset, genotypeRowsPerBlock);
				profiler.record(QeQTLProfiler.Phase.GENOTYPE_CACHE_READ, cacheReadStarted,
					prepared.values().length, (long) prepared.values().length * numInds * Double.BYTES);
				pending.addLast(threadPool.submit(new QeQTLStreamedJobReal(prepared, expressionCache,
					checkpoint, blockNumber, contextPool, genotypeRowsPerBlock,
					expressionRowsPerBlock, localBlockSize, dfo, dfe, rsq0, gpuPrecision, profiler,
					includeSampleStatistics)));
				if (pending.size() >= numThreads)
					pending.removeFirst().get();
			}
			while (!pending.isEmpty())
				pending.removeFirst().get();
			threadPool.shutdown();
			if (!threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS))
				throw new RuntimeException("Timed out while waiting for streamed eQTL workers");
			long assemblyStarted = profiler.start();
			checkpoint.assemble(outputPath, outputHeader);
			profiler.record(QeQTLProfiler.Phase.OUTPUT_ASSEMBLY, assemblyStarted, totalBlocks, 0);
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			threadPool.shutdownNow();
			throw new RuntimeException("Interrupted while waiting for streamed eQTL workers", e);
		} catch (Exception e) {
			threadPool.shutdownNow();
			throw new RuntimeException("Streamed eQTL worker failed", e);
		} finally {
			contextPool.close();
		}
	}

	@Override
	public void reportCompletion(IGenericParallelJob job)
	{ }

	@Override
	public boolean isCanceled()
	{	return false; }

	static final void dumpToString(QeQTLSNPAnalysisResults results, QGeneticSNPData popn, QGeneExpressionData tbl, String filename)
	{
		int[] snpCounts = results.mSNPCounts;
		int[][] eTraitIndices = results.mETraitIndices;
		double[][][] snpResults = results.mSNPResults;
		String[] geneID = tbl.getGeneIDs();
		List<QSNPData> snpList =  popn.getSNPs();
		int numSNPs = snpCounts.length;
		Writer fw = null;
		try {
			fw = new PrintWriter(new FileOutputStream(filename), true);
			fw.write("Rs_ID,ProbesetID,RSq,Fx,T,log10P");
			fw.write(sLn);
			for (int snpNo = 0; snpNo < numSNPs; snpNo++)
			{
				int numETraitsInCurSNP = snpCounts[snpNo];
				int[] curETraitIndices = eTraitIndices[snpNo];
				double[][] curSNPResults = snpResults[snpNo];
				if (numETraitsInCurSNP == 0 || curETraitIndices == null)
					continue;
				String snpID = snpList.get(snpNo).getID();
				for (int eTraitNo = 0; eTraitNo < numETraitsInCurSNP; eTraitNo++)
				{
					String probesetID = geneID[curETraitIndices[eTraitNo]];
					double[] curResults = curSNPResults[eTraitNo];
					fw.write(snpID + "," + probesetID + "," + curResults[0]);
					for (int i = 1; i < curResults.length; i++)
						fw.write("," + curResults[i]);
					fw.write(sLn);
				}
			}
			fw.flush();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (fw != null)
				try {
					fw.close();
				} catch (IOException e) {
				}
		}
	}

	private static final void initGPUs()
	{
		GpuRuntime gpuRuntime = getGpuRuntime();
		System.out.println("GPU backend: " + gpuRuntime.getBackend().getName());
		System.out.println("GPU runtime: " + gpuRuntime.getBackend().getRuntimeDescription());
		List<GpuDevice> devices;
		try {
			devices = gpuRuntime.getGpuDevices(true, gpuPrecision.requiresDoublePrecision());
		} catch (Throwable e) {
			System.err.println("Cannot initialize the GPU backend: " + e.getMessage());
			System.exit(kExitCodeErrorInitOpenCLFailure);
			return;
		}
		int numDevices = devices.size();
		if (numDevices == 0) {
			System.err.println("Cannot find an available GPU supporting " + gpuPrecision.optionName() + ".");
			List<GpuDevice> allGpuDevices = gpuRuntime.getGpuDevices(true, false);
			if (gpuPrecision == GpuPrecision.FP64 && !allGpuDevices.isEmpty()) {
				System.err.println("The following GPUs were detected but do not report double-precision support:");
				for (GpuDevice device : allGpuDevices) {
					System.err.println("  " + device.getName() + " (" + device.getComputeApiVersion() + ")");
				}
			} else {
				System.err.println("No GPU was reported by the system OpenCL ICD drivers.");
			}
			System.exit(kExitCodeError64bitGPUNotFound);
			return;
		}
		System.out.println("Found " + numDevices + " suitable GPU" + (numDevices > 1 ? "s": "") + " in this machine:");
		boolean hasHostUnifiedMemory = false;
		for (GpuDevice device : devices) {
			System.out.println(device.getName() + " [" + device.getVendor() + ", " + device.getComputeApiVersion()
				+ ", VRAM=" + String.format("%.2f GiB", device.getGlobalMemoryBytes() / kGB) + "]");
			hasHostUnifiedMemory = hasHostUnifiedMemory | device.hasUnifiedMemory();
		}
		if (hasHostUnifiedMemory) {
			System.err.println("NOTE: At least one GPU reports unified host memory. Performance depends on the device and driver; this is not necessarily an error.");
		}
		mContexts = new GpuContext[numDevices];
		try {
			for (int i = 0; i < numDevices; i++) {
				mContexts[i] = devices.get(i).openContext();
			}
		} catch (RuntimeException e) {
			for (GpuContext context : mContexts) {
				if (context != null) context.close();
			}
			throw e;
		}
		String devName = null;
		boolean allDevicesTheSame = true;
		for (GpuContext context : mContexts) {
			GpuDevice dev = context.getDevice();
			if (devName == null) {
				devName = dev.getName();
			} else if (!devName.equalsIgnoreCase(dev.getName())) {
				allDevicesTheSame = false;
			}
		}
		if (!allDevicesTheSame) {
			System.err.println("WARNING: The detected GPUs are not identical. Performance may be constrained by the weakest GPU.");
		} else {
			System.out.println("All the detected GPUs appear to be identical.");
		}
	}

	private static final void printUsage()
	{
		System.err.print(QeQTLCommandLine.usage());
	}
	private static final void dumpGPUInfo()
	{
		System.out.println(getGpuRuntime().describeGpuDevices());
		System.exit(kExitCodeNormal);
	}

	private static GpuRuntime getGpuRuntime()
	{
		if (GPU_RUNTIME == null)
			GPU_RUNTIME = GpuRuntime.createDefault();
		return GPU_RUNTIME;
	}

	static final int getDefaultBlockSize(int numInds, int numSNPs, int numETraits)
	{
		if (mContexts == null || mContexts.length == 0)
			throw new IllegalStateException("GPU contexts are required for automatic block sizing");
		List<GpuDevice> devices = new ArrayList<GpuDevice>(mContexts.length);
		for (GpuContext context : mContexts)
			devices.add(context.getDevice());
		GpuTuning.BlockRecommendation recommendation = GpuTuning.recommendBlockSize(
			devices, numInds, numSNPs, numETraits, gpuPrecision);
		System.out.println("Automatic block-size limiter: " + recommendation.limitingDevice()
			+ "; estimated device buffers="
			+ String.format("%.2f GiB", recommendation.estimatedDeviceBytes() / kGB));
		return recommendation.blockSize();
	}

	private static int configureThreadCount(long requiredIterations, int numDevices, boolean streamed)
	{
		return configureThreadCount(requiredIterations, numDevices, streamed, Long.MAX_VALUE, 1);
	}

	private static int configureThreadCount(long requiredIterations, int numDevices, boolean streamed,
		long availableHeapBytes, long estimatedWorkerBytes)
	{
		int threads = config.getNumThreads();
		if (threads == 0) {
			threads = GpuTuning.recommendThreadCount(QSystemUtils.kNumCPUCores,
				numDevices, requiredIterations, streamed, availableHeapBytes, estimatedWorkerBytes);
			config.setNumThreads(threads);
			System.out.println("The thread count was not specified; using " + threads
				+ " (" + numDevices + " GPU context" + (numDevices == 1 ? "" : "s") + ")");
		}
		if (threads > QSystemUtils.kNumCPUCores)
			System.err.println("WARNING: num_threads exceeds the available CPU cores ("
				+ QSystemUtils.kNumCPUCores + ") and may reduce performance.");
		if (threads < Math.min((long) numDevices, requiredIterations))
			System.err.println("WARNING: num_threads is lower than the number of usable GPUs; some GPUs may remain idle.");
		if (streamed && estimatedWorkerBytes > 1
			&& estimatedWorkerBytes * (double) threads > availableHeapBytes * 0.75)
			System.err.println("WARNING: explicit num_threads may require more than 75% of the available JVM heap; "
				+ "reduce --threads or streamed block rows if the run exhausts memory.");
		return threads;
	}

	private static String kernelHeader(int localBlockSize, GpuPrecision precision)
	{
		String header = "#define BLOCK_SIZE " + localBlockSize + sLn
			+ "#define DATATYPE " + (precision == GpuPrecision.FP32 ? "float" : "double") + sLn;
		if (precision == GpuPrecision.FP64)
			header += "#if defined(cl_khr_fp64)" + sLn
				+ "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" + sLn
				+ "#elif defined(cl_amd_fp64)" + sLn
				+ "#pragma OPENCL EXTENSION cl_amd_fp64 : enable" + sLn + "#endif" + sLn;
		return header;
	}

	private static void validateBlockCapacity(int blockSize)
	{
		if (blockSize < 16 || blockSize % 16 != 0)
			throw new IllegalArgumentException("block_size must be a positive multiple of 16");
		if ((long) blockSize * blockSize > Integer.MAX_VALUE)
			throw new IllegalArgumentException("block_size is too large for a Java result array: " + blockSize);
	}

	static final void checkCovars(Set<String> allColNames, String[] covars)
	{
		try {
			for (String str : covars) {
				if (str.contains("*")) {
					for (String sstr : str.split("\\*"))
						if (!allColNames.contains(sstr))
							throw new RuntimeException(sstr);
				} else
					if (!allColNames.contains(str))
						throw new RuntimeException(str);
			}
		} catch (Exception e) {
			System.err.println("Error: Covariate '" + e.getMessage() + "' is not found in the covariate file");
			System.exit(kExitCodeErrorCovarNotFound);
		}
	}

	private static void runMatrixAnalysis(QeQTLAnalysis plugin, String genotypeFilename,
		String expressionFilename, String covariateFilename, String[] fixedCovariates,
		String[] factorCovariates, String thresholdType, double threshold, int dfOffset,
		boolean isAdditive, int configuredBlockSize, int numDevices, String genotypeFormat) throws Exception
	{
		long metadataStarted = profiler.start();
		System.out.println("Scanning matrix metadata and identifiers...");
		QDataType predictorType = config.getPredictorDataType();
		QDataType traitType = config.getTraitDataType();
		QMissingValuePolicy predictorMissing = config.getPredictorMissingPolicy();
		QMissingValuePolicy traitMissing = config.getTraitMissingPolicy();
		String covariateMissing = config.getCovariateMissingPolicy();
		if (!genotypeFormat.equals("csv") && !predictorType.isGenotype())
			throw new IllegalArgumentException("VCF/BCF input requires --predictor-type genotype");
		if (predictorMissing == QMissingValuePolicy.LOCAL_PATTERN && !predictorType.isGenotype())
			throw new IllegalArgumentException("local-pattern imputation is only valid for --predictor-type genotype");
		if (predictorMissing == QMissingValuePolicy.PATTERN)
			throw new IllegalArgumentException("Pattern-wise predictor deletion is not yet supported; use mean, local-pattern, error, zero, or exclude-row");
		if (!(covariateMissing.equals("error") || covariateMissing.equals("complete-samples")))
			throw new IllegalArgumentException("covariate_missing must be error or complete-samples");
		QMatrixRowSource genotypeSource;
		if (genotypeFormat.equals("csv")) {
			genotypeSource = new QDelimitedMatrixSource(Path.of(genotypeFilename), cCommonDelimiter, "#");
		} else {
			if (!isAdditive)
				throw new IllegalArgumentException("VCF/BCF input currently supports the additive model only");
			Path qcOutput = config.getVariantQcOutputFilename() == null
				? Path.of(config.getOutputFilename() == null
					? genotypeFilename + ".variants.tsv" : config.getOutputFilename() + ".variants.tsv")
				: Path.of(config.getVariantQcOutputFilename());
			if (config.getOutputFilename() != null && qcOutput.toAbsolutePath().normalize().equals(
				Path.of(config.getOutputFilename()).toAbsolutePath().normalize()))
				throw new IllegalArgumentException("Variant QC output must differ from association output");
			QVariantMatrixSource.Options variantOptions = new QVariantMatrixSource.Options(
				QVariantMatrixSource.Format.parse(genotypeFormat, Path.of(genotypeFilename)),
				QVariantMatrixSource.GenotypeField.parse(config.getGenotypeField()),
				predictorMissing == QMissingValuePolicy.EXCLUDE_ROW
					? QVariantMatrixSource.MissingPolicy.EXCLUDE_VARIANT
					: QVariantMatrixSource.MissingPolicy.PRESERVE,
				QVariantMatrixSource.MultiallelicPolicy.parse(config.getMultiallelicPolicy()),
				config.getMinimumMaf(), config.getMinimumMac(), qcOutput);
			QVariantMatrixSource variantSource = new QVariantMatrixSource(
				Path.of(genotypeFilename), variantOptions);
			genotypeSource = variantSource;
			QVariantMatrixSource.Summary variantSummary = variantSource.summary();
			System.out.println("Variant genotype field = " + variantSource.selectedField());
			System.out.println("Variant QC: input=" + variantSummary.inputRecords()
				+ ", included=" + variantSummary.includedVariants()
				+ ", monomorphic=" + variantSummary.monomorphicVariants()
				+ ", singletons=" + variantSummary.singletons()
				+ ", doubletons=" + variantSummary.doubletons());
			System.out.println("Variant annotation/QC output: " + qcOutput.toAbsolutePath().normalize());
		}
		QMatrixRowSource expressionSource = new QDelimitedMatrixSource(
			Path.of(expressionFilename), cCommonDelimiter, "#");
		if (genotypeSource.metadata().rowCount() > Integer.MAX_VALUE
			|| expressionSource.metadata().rowCount() > Integer.MAX_VALUE)
			throw new IllegalArgumentException("This release supports at most " + Integer.MAX_VALUE + " rows per matrix");

		System.out.println("Declared matrix types: predictor=" + predictorType.optionName()
			+ ", trait=" + traitType.optionName());
		System.out.println("Missing-value policies: predictor=" + predictorMissing.optionName()
			+ ", trait=" + traitMissing.optionName() + ", covariates=" + covariateMissing);

		QSampleAlignment fullAlignment;
		QSampleAlignment alignment;
		QCovariateTable covariates = null;
		QCovariateTable.Missingness covariateMissingness = null;
		boolean[] retainedCovariateRows = null;
		String[] canonicalSampleIds;
		double[][] covariateQ = null;
		double[][] covariateModel = null;
		int covariateRank = 1;
		if (covariateFilename != null) {
			System.out.println("Loading and validating covariate data...");
			covariates = QCovariateTable.load(Path.of(covariateFilename), cCommonDelimiter, "#");
			fullAlignment = covariates.align(genotypeSource.metadata().sampleIds(),
				expressionSource.metadata().sampleIds(), config.getGenotypeIdColumn(), config.getExpressionIdColumn());
			canonicalSampleIds = reorder(genotypeSource.metadata().sampleIds(), fullAlignment.genotypeColumnOrder());
			covariateMissingness = covariates.inspectMissingness(fixedCovariates);
			retainedCovariateRows = covariateMissingness.completeRows();
			if (covariateMissing.equals("error")) {
				java.util.Arrays.fill(retainedCovariateRows, true);
				alignment = fullAlignment;
			} else {
				alignment = fullAlignment.select(retainedCovariateRows);
				if (alignment.sampleCount() < fullAlignment.sampleCount())
					System.err.println("WARNING: Removed " + (fullAlignment.sampleCount() - alignment.sampleCount())
						+ " sample(s) missing a selected covariate.");
			}
		} else {
			fullAlignment = QCovariateTable.alignDirectly(genotypeSource.metadata().sampleIds(),
				expressionSource.metadata().sampleIds());
			alignment = fullAlignment;
			canonicalSampleIds = reorder(genotypeSource.metadata().sampleIds(), alignment.genotypeColumnOrder());
		}

		System.out.println("Inspecting matrix missingness...");
		QMissingnessScan predictorScan = QMissingnessScan.scan("predictor", genotypeSource,
			alignment.genotypeColumnOrder());
		QMissingnessScan traitScan = QMissingnessScan.scan("trait", expressionSource,
			alignment.expressionColumnOrder());
		Path missingnessOutput = missingnessOutput(genotypeFilename);
		if (config.getOutputFilename() != null && missingnessOutput.equals(
			Path.of(config.getOutputFilename()).toAbsolutePath().normalize()))
			throw new IllegalArgumentException("Missingness QC output must differ from association output");
		QMissingnessReport.write(missingnessOutput, predictorScan, predictorType, predictorMissing,
			traitScan, traitType, traitMissing, covariateMissingness, covariateMissing, canonicalSampleIds);
		canonicalSampleIds = reorder(genotypeSource.metadata().sampleIds(), alignment.genotypeColumnOrder());
		System.out.println("Missingness QC output: " + missingnessOutput);
		System.out.println("Missing values: predictor=" + predictorScan.totalMissingValues()
			+ ", trait=" + traitScan.totalMissingValues()
			+ ", selected covariates=" + (covariateMissingness == null ? 0 : covariateMissingness.totalMissing()));
		if (config.getInspectMissingness()) {
			System.out.println("Missingness inspection completed; no association analysis was run.");
			return;
		}
		if (covariateMissingness != null && covariateMissing.equals("error")
			&& covariateMissingness.totalMissing() > 0)
			throw new IllegalArgumentException("Selected covariates contain missing values; see " + missingnessOutput);
		if (predictorMissing == QMissingValuePolicy.ERROR && predictorScan.hasMissingValues())
			throw new IllegalArgumentException("Predictor matrix contains missing values; see " + missingnessOutput);
		if (traitMissing == QMissingValuePolicy.ERROR && traitScan.hasMissingValues())
			throw new IllegalArgumentException("Trait matrix contains missing values; see " + missingnessOutput);

		if (predictorMissing == QMissingValuePolicy.LOCAL_PATTERN) {
			System.err.println("WARNING: local-pattern is a nearest flanking-genotype proxy; it is not "
				+ "phasing or reference-panel imputation. Flanks per side=" + config.getPredictorFlankCount());
			genotypeSource = new QLocalPatternImputedSource(genotypeSource, predictorScan,
				config.getPredictorFlankCount());
		} else {
			genotypeSource = new QPolicyMatrixSource(genotypeSource, predictorScan, predictorMissing);
		}
		if (!(traitMissing == QMissingValuePolicy.PATTERN && traitScan.hasMissingValues()))
			expressionSource = new QPolicyMatrixSource(expressionSource, traitScan, traitMissing);

		if (covariates != null) {
			QCovariateTable.ModelMatrix model = covariates.buildModelMatrix(
				fixedCovariates, factorCovariates, retainedCovariateRows);
			if (model.automaticFactors().length > 0)
				System.out.println("Automatically categorical covariates: " + String.join(", ", model.automaticFactors()));
			System.out.println("Covariate model columns: " + String.join(", ", model.columnNames()));
			QRDecomposition decomposition = new QRDecomposition(model.values());
			if (!decomposition.isFullRank())
				throw new IllegalArgumentException("Covariate model is rank deficient: rank "
					+ decomposition.getRank() + " for " + decomposition.getNumColumns() + " columns");
			covariateRank = decomposition.getRank();
			covariateQ = decomposition.getQ().getArray();
			covariateModel = model.values();
		}

		long genotypeRowCount = genotypeSource.metadata().rowCount();
		long expressionRowCount = traitMissing == QMissingValuePolicy.PATTERN && traitScan.hasMissingValues()
			? traitScan.rowCount() : expressionSource.metadata().rowCount();
		int numIndividuals = alignment.sampleCount();
		int numSnps = (int) genotypeRowCount;
		int numTraits = (int) expressionRowCount;
		System.out.println("Sample alignment: genotype via " + alignment.genotypeIdColumn()
			+ ", expression via " + alignment.expressionIdColumn());
		System.out.println("Samples reordered: genotype=" + alignment.genotypeReorderedCount()
			+ ", expression=" + alignment.expressionReorderedCount());
		System.out.println("Number of SNPs = " + numSnps);
		System.out.println("Number of e-traits = " + numTraits);
		System.out.println("Number of aligned individuals = " + numIndividuals);

		int regressionDegreesOfFreedom = covariateRank + (isAdditive ? 1 : 2) - 1;
		int errorDegreesOfFreedom = numIndividuals - regressionDegreesOfFreedom - 1;
		if (errorDegreesOfFreedom - dfOffset <= 0)
			throw new IllegalArgumentException("Non-positive residual degrees of freedom after df_offset");
		System.out.println("Degrees of Freedom for Regression = " + regressionDegreesOfFreedom);
		System.out.println("Degrees of Freedom for Errors = " + errorDegreesOfFreedom);
		if (dfOffset != 0)
			System.out.println("Degrees of Freedom for Offset = " + dfOffset);
		profiler.record(QeQTLProfiler.Phase.METADATA_AND_ALIGNMENT, metadataStarted,
			numIndividuals, 0);
		if (config.getValidateOnly()) {
			if (configuredBlockSize <= 0 || config.getNumThreads() == 0)
				System.out.println("Automatic block/thread tuning is deferred until a real GPU analysis run.");
			System.out.println("Validation completed successfully; --validate-only requested, so no analysis was run.");
			return;
		}

		int globalBlockSize = configuredBlockSize;
		if (globalBlockSize <= 0) {
			globalBlockSize = getDefaultBlockSize(numIndividuals, numSnps, numTraits);
			if (globalBlockSize <= 0)
				throw new IllegalArgumentException("Cannot infer a positive block size; specify --block-size");
			config.setBlockSize(globalBlockSize);
			System.out.println("The block size was not specified; using " + globalBlockSize);
		} else {
			System.out.println("Block size = " + globalBlockSize);
		}
		validateBlockCapacity(globalBlockSize);

		long numRequiredIterations = ((numSnps + (long) globalBlockSize - 1) / globalBlockSize)
			* ((numTraits + (long) globalBlockSize - 1) / globalBlockSize);
		System.out.println("Given the block size, " + numRequiredIterations + " iteration(s) are needed");
		boolean dynamicTraitPatterns = traitMissing == QMissingValuePolicy.PATTERN
			&& traitScan.hasMissingValues();
		boolean stream = dynamicTraitPatterns || config.getGenotypeBlockRows() > 0
			|| config.getExpressionBlockRows() > 0;
		int genotypeRows = 0;
		int expressionRows = 0;
		long schedulableJobs = numRequiredIterations;
		if (stream) {
			genotypeRows = config.getGenotypeBlockRows() > 0 ? config.getGenotypeBlockRows() : globalBlockSize;
			expressionRows = config.getExpressionBlockRows() > 0 ? config.getExpressionBlockRows() : globalBlockSize;
			genotypeRows = roundUpNearestMultiple(genotypeRows, 16);
			expressionRows = roundUpNearestMultiple(expressionRows, 16);
			if (genotypeRows > globalBlockSize || expressionRows > globalBlockSize)
				throw new IllegalArgumentException("Streamed row blocks must not exceed block_size ("
					+ globalBlockSize + ")");
			schedulableJobs = (numSnps + (long) genotypeRows - 1) / genotypeRows;
		}
		if (stream) {
			Runtime runtime = Runtime.getRuntime();
			long availableHeap = runtime.maxMemory() - (runtime.totalMemory() - runtime.freeMemory());
			long estimatedWorkerBytes = GpuTuning.estimateStreamedWorkerBytes(numIndividuals,
				genotypeRows, expressionRows, gpuPrecision);
			System.out.println("Estimated streamed memory per worker = "
				+ String.format("%.2f MiB", estimatedWorkerBytes / kMB));
			configureThreadCount(schedulableJobs, numDevices, true,
				availableHeap, estimatedWorkerBytes);
		} else {
			configureThreadCount(schedulableJobs, numDevices, false);
		}

		double rSquaredThreshold = thresholdType.equals("none") ? 0.0 : threshold;
		if (thresholdType.equals("pval")) {
			double tValue = T.quantile(threshold / 2.0, errorDegreesOfFreedom - dfOffset, false, false);
			tValue *= tValue;
			rSquaredThreshold = tValue / (errorDegreesOfFreedom + tValue);
			System.out.println("P-value of " + threshold + " corresponds to R^2 of " + rSquaredThreshold);
		}

		long start = System.currentTimeMillis();
		long analysisStarted = profiler.start();
		if (dynamicTraitPatterns) {
			if (!isAdditive)
				throw new IllegalArgumentException("Pattern-wise trait deletion currently supports the additive model only");
			if (config.getResume() || config.getKeepCheckpoints())
				throw new IllegalArgumentException("Pattern-wise trait deletion does not yet support --resume or --keep-checkpoints");
			System.err.println("WARNING: Exact pattern-wise trait deletion will run one GPU pass per distinct "
				+ "missingness pattern and can be very slow. Patterns=" + traitScan.patterns().size());
			runTraitPatternAnalysis(plugin, genotypeSource, expressionSource, traitScan,
				alignment, covariateModel, thresholdType, threshold, dfOffset, genotypeRows,
				expressionRows, covariateFilename != null);
			profiler.record(QeQTLProfiler.Phase.ANALYSIS_WALL, analysisStarted,
				(long) numSnps * numTraits, 0);
			System.out.println("Total analysis time (in seconds) = "
				+ (System.currentTimeMillis() - start) / 1000.0);
			return;
		}
		QGpuResidualizer gpuResidualizer = null;
		if (covariateQ != null && residualizationMode != QResidualizationMode.CPU) {
			gpuResidualizer = new QGpuResidualizer(mContexts, covariateQ, gpuPrecision, profiler);
			System.out.println("Fixed-effect residualization = GPU (" + gpuPrecision.optionName()
				+ ", " + mContexts.length + " context" + (mContexts.length == 1 ? "" : "s") + ")");
		} else if (covariateQ != null) {
			System.out.println("Fixed-effect residualization = CPU FP64");
		}
		try {
		if (stream) {
			System.out.println("Bounded-RAM matrix mode: genotype blocks=" + genotypeRows
				+ " rows, expression blocks=" + expressionRows + " rows");
			Path outputPath = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
			Path cacheDirectory = config.getCacheDirectory() == null
				? outputPath.getParent().resolve(".gpu-eqtl-cache")
				: Path.of(config.getCacheDirectory());
			long signatureStarted = profiler.start();
			String preprocessingTag = gpuResidualizer == null ? null : gpuResidualizer.cacheSignatureTag();
			String genotypeSignature = QBinaryMatrixCache.signature("Genotype", genotypeSource,
				alignment.genotypeColumnOrder(), covariateQ, preprocessingTag);
			String expressionSignature = QBinaryMatrixCache.signature("Expression", expressionSource,
				alignment.expressionColumnOrder(), covariateQ, preprocessingTag);
			profiler.record(QeQTLProfiler.Phase.CACHE_SIGNATURES, signatureStarted, 2, 0);
			long genotypeCacheStarted = profiler.start();
			try (QBinaryMatrixCache genotypeCache = QBinaryMatrixCache.openOrBuild(cacheDirectory,
				"Genotype", genotypeSignature, genotypeSource, alignment.genotypeColumnOrder(), covariateQ,
				genotypeRows, config.getRebuildCache(), gpuResidualizer)) {
				profiler.record(QeQTLProfiler.Phase.GENOTYPE_CACHE_OPEN_OR_BUILD, genotypeCacheStarted,
					genotypeCache.rowCount(), java.nio.file.Files.size(genotypeCache.path()));
				long expressionCacheStarted = profiler.start();
				try (QBinaryMatrixCache expressionCache = QBinaryMatrixCache.openOrBuild(cacheDirectory,
				"Expression", expressionSignature, expressionSource, alignment.expressionColumnOrder(), covariateQ,
				expressionRows, config.getRebuildCache(), gpuResidualizer)) {
					profiler.record(QeQTLProfiler.Phase.EXPRESSION_CACHE_OPEN_OR_BUILD, expressionCacheStarted,
						expressionCache.rowCount(), java.nio.file.Files.size(expressionCache.path()));
					if (gpuResidualizer != null) {
						gpuResidualizer.close();
						gpuResidualizer = null;
					}
					plugin.eSNPAnalysisStreamed(genotypeCache, expressionCache, rSquaredThreshold, dfOffset,
						errorDegreesOfFreedom, genotypeRows, expressionRows);
				}
			}
		} else {
			System.out.println("Loading aligned matrices into RAM (set --genotype-block-rows or --expression-block-rows to stream).");
			QMatrixRowSource.Block genotypeBlock;
			QMatrixRowSource.Block expressionBlock;
			try (QMatrixRowSource.BlockReader reader = genotypeSource.open(alignment.genotypeColumnOrder())) {
				genotypeBlock = reader.readBlock(numSnps);
			}
			try (QMatrixRowSource.BlockReader reader = expressionSource.open(alignment.expressionColumnOrder())) {
				expressionBlock = reader.readBlock(numTraits);
			}
			QeQTLPreprocessor.PreparedBlock preparedGenotypes = QeQTLPreprocessor.prepare(
				genotypeBlock, covariateQ, "Genotype", gpuResidualizer);
			QeQTLPreprocessor.PreparedBlock preparedExpressions = QeQTLPreprocessor.prepare(
				expressionBlock, covariateQ, "Expression", gpuResidualizer);
			if (gpuResidualizer != null) {
				gpuResidualizer.close();
				gpuResidualizer = null;
			}
			QGeneticSNPData genotypeData = toGenotypeData(new QMatrixRowSource.Block(
				preparedGenotypes.rowOffset(), preparedGenotypes.rowIds(), preparedGenotypes.values()));
			QGeneExpressionData expressionData = new QGeneExpressionData(preparedExpressions.values(),
				preparedExpressions.rowIds(), canonicalSampleIds);
			plugin.eSNPAnalysisPrepared(genotypeData, expressionData,
				preparedGenotypes.standardDeviations(), preparedExpressions.standardDeviations(),
				rSquaredThreshold, dfOffset, errorDegreesOfFreedom, isAdditive);
		}
		} finally {
			if (gpuResidualizer != null)
				gpuResidualizer.close();
		}
		profiler.record(QeQTLProfiler.Phase.ANALYSIS_WALL, analysisStarted,
			(long) numSnps * numTraits, 0);
		System.out.println("Total analysis time (in seconds) = " + (System.currentTimeMillis() - start) / 1000.0);
	}

	private static void runTraitPatternAnalysis(QeQTLAnalysis plugin,
		QMatrixRowSource predictorSource, QMatrixRowSource rawTraitSource,
		QMissingnessScan traitScan,
		QSampleAlignment alignment, double[][] covariateModel,
		String thresholdType, double threshold, int dfOffset,
		int predictorRowsPerBlock, int traitRowsPerBlock, boolean hasCovariates) throws Exception
	{
		Path output = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
		Path parent = output.getParent();
		if (parent == null)
			parent = Path.of(".").toAbsolutePath().normalize();
		Files.createDirectories(parent);
		Path work = Files.createTempDirectory(parent, ".gpu-eqtl-pattern-").toAbsolutePath().normalize();
		Path partialOutput = Files.createTempFile(parent, output.getFileName().toString(), ".partial");
		boolean complete = false;
		try (BufferedWriter combined = Files.newBufferedWriter(partialOutput, StandardCharsets.UTF_8)) {
			combined.write(rsqOnly ? "Rs_ID,ProbesetID,RSq,Dir,N,DF"
				: "Rs_ID,ProbesetID,RSq,Fx,T,log10P,N,DF");
			combined.newLine();
			int patternNumber = 0;
			for (QMissingnessScan.Pattern pattern : traitScan.patterns()) {
				patternNumber++;
				BitSet missing = pattern.missingSamples();
				int[] observed = complement(missing, alignment.sampleCount());
				int[] predictorColumns = select(alignment.genotypeColumnOrder(), observed);
				int[] traitColumns = select(alignment.expressionColumnOrder(), observed);
				double[][] patternQ = null;
				int covariateRank = 1;
				if (hasCovariates) {
					double[][] selectedModel = selectRows(covariateModel, observed);
					QRDecomposition decomposition = new QRDecomposition(selectedModel);
					if (!decomposition.isFullRank())
						throw new IllegalArgumentException("Covariate model is rank deficient for trait missingness pattern "
							+ pattern.id() + " (" + pattern.rowIndices().length + " trait row(s), N="
							+ observed.length + ")");
					covariateRank = decomposition.getRank();
					patternQ = decomposition.getQ().getArray();
				}
				int regressionDegreesOfFreedom = covariateRank;
				int errorDegreesOfFreedom = observed.length - regressionDegreesOfFreedom - 1;
				if (errorDegreesOfFreedom - dfOffset <= 0)
					throw new IllegalArgumentException("Non-positive residual degrees of freedom for trait missingness pattern "
						+ pattern.id() + " (N=" + observed.length + ", covariate rank=" + covariateRank + ")");
				double patternThreshold = rSquaredThreshold(thresholdType, threshold,
					errorDegreesOfFreedom, dfOffset);
				BitSet selectedTraits = new BitSet((int) traitScan.rowCount());
				for (long row : pattern.rowIndices())
					selectedTraits.set(Math.toIntExact(row));
				QMatrixRowSource traitSource = new QPolicyMatrixSource(rawTraitSource,
					traitScan, QMissingValuePolicy.ERROR, selectedTraits);
				String suffix = "Pattern" + String.format("%06d", pattern.id());
				System.out.println("Trait pattern " + patternNumber + "/" + traitScan.patterns().size()
					+ ": traits=" + pattern.rowIndices().length + ", N=" + observed.length
					+ ", DF=" + (errorDegreesOfFreedom - dfOffset));

				QGpuResidualizer residualizer = null;
				if (patternQ != null && residualizationMode != QResidualizationMode.CPU)
					residualizer = new QGpuResidualizer(mContexts, patternQ, gpuPrecision, profiler);
				Path predictorCachePath = null;
				Path traitCachePath = null;
				try {
					String preprocessingTag = residualizer == null ? null : residualizer.cacheSignatureTag();
					String predictorKind = "Predictor" + suffix;
					String traitKind = "Trait" + suffix;
					String predictorSignature = QBinaryMatrixCache.signature(predictorKind,
						predictorSource, predictorColumns, patternQ, preprocessingTag);
					String traitSignature = QBinaryMatrixCache.signature(traitKind,
						traitSource, traitColumns, patternQ, preprocessingTag);
					try (QBinaryMatrixCache predictorCache = QBinaryMatrixCache.openOrBuild(work,
						predictorKind, predictorSignature, predictorSource, predictorColumns, patternQ,
						predictorRowsPerBlock, true, residualizer);
						 QBinaryMatrixCache traitCache = QBinaryMatrixCache.openOrBuild(work,
						traitKind, traitSignature, traitSource, traitColumns, patternQ,
						traitRowsPerBlock, true, residualizer)) {
						predictorCachePath = predictorCache.path();
						traitCachePath = traitCache.path();
						if (residualizer != null) {
							residualizer.close();
							residualizer = null;
						}
						Path groupOutput = work.resolve(suffix + ".csv");
						Path checkpoint = work.resolve(suffix + ".checkpoint");
						plugin.eSNPAnalysisStreamed(predictorCache, traitCache, patternThreshold,
							dfOffset, errorDegreesOfFreedom, predictorRowsPerBlock, traitRowsPerBlock,
							groupOutput, checkpoint, false, false, false, true);
						appendWithoutHeader(groupOutput, combined);
						Files.deleteIfExists(groupOutput);
					}
				} finally {
					if (residualizer != null)
						residualizer.close();
					if (predictorCachePath != null)
						Files.deleteIfExists(predictorCachePath);
					if (traitCachePath != null)
						Files.deleteIfExists(traitCachePath);
				}
			}
			combined.flush();
			complete = true;
		} finally {
			for (GpuContext context : mContexts)
				if (context != null)
					context.close();
			cleanupPatternWorkDirectory(work, parent);
			if (!complete)
				Files.deleteIfExists(partialOutput);
		}
		if (complete)
			moveAtomically(partialOutput, output);
	}

	private static double rSquaredThreshold(String type, double threshold,
		int errorDegreesOfFreedom, int dfOffset) {
		if (type.equals("none"))
			return 0;
		if (!type.equals("pval"))
			return threshold;
		double tValue = T.quantile(threshold / 2.0,
			errorDegreesOfFreedom - dfOffset, false, false);
		tValue *= tValue;
		return tValue / (errorDegreesOfFreedom + tValue);
	}

	private static int[] complement(BitSet excluded, int size) {
		int[] selected = new int[size - excluded.cardinality()];
		int output = 0;
		for (int index = excluded.nextClearBit(0); index < size; index = excluded.nextClearBit(index + 1))
			selected[output++] = index;
		return selected;
	}

	private static int[] select(int[] values, int[] indices) {
		int[] selected = new int[indices.length];
		for (int i = 0; i < indices.length; i++)
			selected[i] = values[indices[i]];
		return selected;
	}

	private static double[][] selectRows(double[][] values, int[] indices) {
		if (values == null)
			throw new IllegalArgumentException("Covariate model is required");
		double[][] selected = new double[indices.length][];
		for (int i = 0; i < indices.length; i++)
			selected[i] = values[indices[i]].clone();
		return selected;
	}

	private static void appendWithoutHeader(Path input, BufferedWriter output) throws IOException {
		try (BufferedReader reader = Files.newBufferedReader(input, StandardCharsets.UTF_8)) {
			if (reader.readLine() == null)
				throw new IOException("Pattern result has no header: " + input);
			String line;
			while ((line = reader.readLine()) != null) {
				output.write(line);
				output.newLine();
			}
		}
	}

	private static void cleanupPatternWorkDirectory(Path work, Path expectedParent) throws IOException {
		Path normalized = work.toAbsolutePath().normalize();
		if (!normalized.startsWith(expectedParent.toAbsolutePath().normalize())
			|| !normalized.getFileName().toString().startsWith(".gpu-eqtl-pattern-"))
			throw new IOException("Refusing to clean unexpected pattern work directory " + normalized);
		try (java.util.stream.Stream<Path> paths = Files.walk(normalized)) {
			for (Path path : paths.sorted(java.util.Comparator.reverseOrder()).toList())
				Files.deleteIfExists(path);
		}
	}

	private static void moveAtomically(Path source, Path target) throws IOException {
		try {
			Files.move(source, target, StandardCopyOption.ATOMIC_MOVE,
				StandardCopyOption.REPLACE_EXISTING);
		} catch (AtomicMoveNotSupportedException e) {
			Files.move(source, target, StandardCopyOption.REPLACE_EXISTING);
		}
	}

	private static QGeneticSNPData toGenotypeData(QMatrixRowSource.Block block)
	{
		QGeneticSNPData result = new QGeneticSNPData();
		for (int row = 0; row < block.rowCount(); row++) {
			QSNPDataReal snp = new QSNPDataReal();
			snp.setID(block.rowIds()[row]);
			snp.setSNPValues(block.values()[row]);
			result.addSNP(snp);
		}
		return result;
	}

	private static void validateFinite(QMatrixRowSource.Block block, String matrixName)
	{
		for (int row = 0; row < block.rowCount(); row++)
			for (double value : block.values()[row])
				if (value == gov.nih.utils.QDataUtils.kUndefinedValue || !Double.isFinite(value))
					throw new IllegalArgumentException(matrixName + " row '" + block.rowIds()[row]
						+ "' contains a missing or non-finite value");
	}

	private static Path missingnessOutput(String predictorFilename)
	{
		String configured = config.getMissingnessQcOutputFilename();
		if (configured != null)
			return Path.of(configured).toAbsolutePath().normalize();
		String base = config.getOutputFilename() == null ? predictorFilename : config.getOutputFilename();
		return Path.of(base + ".missingness.tsv").toAbsolutePath().normalize();
	}

	private static String[] reorder(String[] values, int[] order)
	{
		String[] result = new String[order.length];
		for (int i = 0; i < order.length; i++)
			result[i] = values[order[i]];
		return result;
	}

	private static String resolveGenotypeFormat(String configured, String filename)
	{
		String format = configured == null ? "auto" : configured.toLowerCase(java.util.Locale.ROOT);
		if (!format.equals("auto")) {
			if (format.equals("vcf.gz") || format.equals("vcfgz"))
				return "vcf";
			return format;
		}
		String name = filename == null ? "" : filename.toLowerCase(java.util.Locale.ROOT);
		if (name.endsWith(".bcf"))
			return "bcf";
		if (name.endsWith(".vcf") || name.endsWith(".vcf.gz"))
			return "vcf";
		return "csv";
	}

	public static void main(String[] args)
	{
		long timea = System.currentTimeMillis();
		System.out.println("GPU eQTL analysis software version 2.0. By: Roby Joehanes");
		if (args == null || args.length < 1) {
			printUsage();
			System.exit(kExitCodeErrorInvalidParam);
		}

		QeQTLCommandLine.Result commandLine;
		try {
			commandLine = QeQTLCommandLine.parse(args);
		} catch (Exception e) {
			System.err.println("ERROR: " + e.getMessage());
			printUsage();
			System.exit(kExitCodeErrorInvalidParam);
			return;
		}
		if (commandLine.help()) {
			System.out.print(QeQTLCommandLine.usage());
			return;
		}
		if (commandLine.gpuBackend() != null)
			System.setProperty("eqtl.gpu.backend", commandLine.gpuBackend());
		config = commandLine.config();
		profiler = new QeQTLProfiler(config.getProfile());
		try {
			gpuPrecision = config.getGpuPrecision();
			residualizationMode = config.getResidualizationMode();
		} catch (IllegalArgumentException e) {
			System.err.println("ERROR: " + e.getMessage());
			System.exit(kExitCodeErrorInvalidParam);
			return;
		}
		if (commandLine.printGpuInfo())
			dumpGPUInfo();
		DEBUG = commandLine.debug();
		System.out.println("GPU matrix-product precision = " + gpuPrecision.optionName());
		System.out.println("Fixed-effect residualization mode = " + residualizationMode.optionName());
		if (gpuPrecision == GpuPrecision.FP32)
			System.err.println("NOTE: FP32 is enabled explicitly. GPU matrix products and GPU covariate "
				+ "projection are approximate; CPU projection and statistical calculations remain FP64.");
		if (config.get("lambda") != null)
			System.err.println("NOTE: The obsolete lambda setting is ignored; automatic sizing uses GPU memory limits.");

		if (!EPlatform.is64Bit()) {
			System.err.println("ERROR: Your operating system / platform seems to be 32-bit, not 64-bit.");
			System.err.println("This program will not run properly in 32-bit platform due to severe limitation to the amount data that can be processed."
				+ " You may have a 64-bit machine, but you will need a 64-bit operating system to run in it. You also need 64-bit Java."
				+ " Please contact your system administrator.");
			System.exit(kExitCodeErrorPlatformNot64Bit);
		}

		int numDevices;
		if (config.getValidateOnly() || config.getInspectMissingness()) {
			numDevices = 1;
			System.out.println("Validation/inspection mode: GPU initialization is skipped.");
		} else {
			System.out.println("Initializing GPUs...");
			initGPUs();
			numDevices = mContexts.length;
		}
		String
			famFilename = config.getFamilyFilename(),
			genoFilename = config.getGenotypeFilename(),
			exprFilename = config.getExpressionFilename(),
			covarFilename = config.getCovariateFilename(),
			pedigreeFilename = config.getPedigreeFilename(),
			outputFilename = config.getOutputFilename(),
			genoFileformat = config.getGenotypeFileFormat(),
			exprFileformat = config.getExpressionFileFormat(),
			genoModel = config.getGenotypeModel(),
			thresholdType = config.getThresholdType(),
			covarFixed[] = config.getFixedCovariates(),
			covarRandom[] = config.getRandomCovariates(),
			covarFactor[] = config.getFactorCovariates();
		genoFileformat = resolveGenotypeFormat(genoFileformat, genoFilename);
		simplifyResult = config.getSimplifyOutput();
		rsqOnly = config.getRSqOutput();
		double t0 = config.getThresholdValue();
		boolean isAdditive = true;
		int
			globalBlockSize = config.getBlockSize(),
			dfOffset = config.getDFOffset();

		if (famFilename != null)
			System.out.println("Family file: " + famFilename);
		System.out.println("Genotype file: " + genoFilename);
		System.out.println("Expression file: " + exprFilename);
		if (covarFilename != null) {
			System.out.println("Covariate file: " + covarFilename);
			if (covarFixed == null && covarRandom == null) {
				System.err.println("ERROR: Fixed or random covariates are not mentioned, but the covariate file is specified.");
				System.exit(kExitCodeErrorWrongCovarSpec);
			} else if (covarRandom != null) {
				System.err.println("WARNING: Random covariates are currently ignored.");
			}
			if (covarFixed == null) {
				System.err.println("ERROR: At present, you will need to specify at least one fixed covariate with a specified covariate file.");
				System.exit(kExitCodeErrorWrongCovarSpec);
			}
			System.out.println("Fixed covariates: " + QStringUtils.toString(covarFixed));
		} else {
			if (covarFixed != null || covarRandom != null) {
				System.err.println("ERROR: Fixed or random covariates are mentioned without specifying the covariate file.");
				System.exit(kExitCodeErrorWrongCovarSpec);
			}
		}
		if (pedigreeFilename != null)
			System.out.println("Pedigree file: " + pedigreeFilename);
		System.out.println("Output file: " + outputFilename);
		System.out.println("Threshold type: " + thresholdType);
		System.out.println("Threshold: " + t0);
		if (!(thresholdType.equals("none") || thresholdType.equals("pval") || thresholdType.equals("rsq"))) {
			System.err.println("ERROR: Threshold type has to be either 'none' or 'rsq' or 'pval'.");
			System.exit(kExitCodeErrorInvalidParam);
		}
		if (!(genoModel.equals("additive") || genoModel.equals("full"))) {
			System.err.println("ERROR: Genotype model has to be either 'additive' or 'full'.");
			System.exit(kExitCodeErrorInvalidParam);
		}
		isAdditive = genoModel.equals("additive");
		System.out.println("Linear model = " + (isAdditive? "Additive" : "Full"));

		if (genoFilename == null || !new File(genoFilename).exists()) {
			System.err.println(String.format("ERROR: Genotype file %s does not exist!", genoFilename));
			System.exit(kExitCodeErrorCantLoadFile);
		}
		if (exprFilename == null || !new File(exprFilename).exists()) {
			System.err.println(String.format("ERROR: Expression file %s does not exist!", exprFilename));
			System.exit(kExitCodeErrorCantLoadFile);
		}
		if (covarFilename != null && !new File(covarFilename).exists()) {
			System.err.println(String.format("ERROR: Covariate file %s does not exist!", covarFilename));
			System.exit(kExitCodeErrorCantLoadFile);
		}
		if (famFilename != null && !new File(famFilename).exists()) {
			System.err.println(String.format("ERROR: Family file %s does not exist!", famFilename));
			System.exit(kExitCodeErrorCantLoadFile);
		}
		if (pedigreeFilename != null && !new File(pedigreeFilename).exists()) {
			System.err.println(String.format("ERROR: Pedigree file %s does not exist!", pedigreeFilename));
			System.exit(kExitCodeErrorCantLoadFile);
		}

		long availMem = QSystemUtils.getMemoryStatistics(true).mAvailable;
		System.out.println("Available memory (MB) = " + String.format("%5.2f", (availMem / kMB)));

		QeQTLAnalysis plugin = new QeQTLAnalysis();
		boolean modernGenotypeFormat = genoFileformat.equalsIgnoreCase("csv")
			|| genoFileformat.equalsIgnoreCase("vcf") || genoFileformat.equalsIgnoreCase("bcf");
		if (modernGenotypeFormat && exprFileformat.equalsIgnoreCase("csv")
			&& (config.getGenotypeFileHeader() || !genoFileformat.equalsIgnoreCase("csv"))) {
			try {
				runMatrixAnalysis(plugin, genoFilename, exprFilename, covarFilename, covarFixed,
					covarFactor, thresholdType, t0, dfOffset, isAdditive, globalBlockSize, numDevices,
					genoFileformat.toLowerCase(java.util.Locale.ROOT));
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(kExitCodeErrorAnalysisFailure);
				return;
			}
			long timeb = System.currentTimeMillis();
			System.out.println("Overall time (in seconds) = " + (timeb - timea) / 1000.0);
			finishProfiling();
			System.exit(kExitCodeNormal);
			return;
		}
		if (config.getValidateOnly()) {
			System.err.println("ERROR: --validate-only requires VCF/BCF or headered CSV genotype input and CSV expression input.");
			System.exit(kExitCodeErrorInvalidParam);
			return;
		}
		double[][] covarMatrix =  null;
		int covarRank = 1; // At least we have the intercept
		BitSet missing = null;
		try
		{
			if (covarFilename != null && covarFixed != null) {
				System.out.println("Loading covariate data...");
				QTableData covarTable = QTableData.load(covarFilename, cCommonDelimiter, "#", true, true, false);
				System.out.println("Covariate table is loaded, " + covarTable.getNumberOfRows() + " x " + covarTable.getNumberOfColumns());
				missing = new BitSet(covarTable.getNumberOfRows());
				Set<String> allColNames = new HashSet<String>();
				for (String str : covarTable.getAllColumnNames())
					allColNames.add(str);
				checkCovars(allColNames, covarFixed);
				if (covarRandom != null)
					checkCovars(allColNames, covarRandom);
				if (covarFactor != null)
					checkCovars(allColNames, covarFactor);
				covarMatrix = covarTable.buildModelMatrix(covarFixed, covarFactor, true, missing);
				if (DEBUG) {
					System.out.println(QStringUtils.toString(covarMatrix));
					System.out.flush();
				}
				if (!missing.isEmpty()) {
					System.err.println("ERROR: Missing values are found in the covariate data! Subset your data first to exclude those with missing values.");
					System.exit(kExitCodeErrorCovarMissingValues);
				}
				QRDecomposition qr = new QRDecomposition(covarMatrix);
				covarMatrix = qr.getQ().getArray();
				covarRank = qr.getRank();
			}
			// At this point covarMatrix is the Q matrix of the QR decomposition
			System.out.println("Loading genotype data...");
			long time1 = System.currentTimeMillis();
			long gb1 = QSystemUtils.usedMemoryAfterGC();
			//System.out.println(new File(".").getAbsolutePath());
			QGeneticSNPData snpData = genoFileformat.equalsIgnoreCase("tped") ?
				new QPlinkLoader().load(famFilename, genoFilename) :
				QGeneticSNPData.load(genoFilename, cCommonDelimiter, "#", true, config.getGenotypeFileHeader());
			long time2 = System.currentTimeMillis();
			System.out.println("Genetic data load time (in seconds) = " + (time2 - time1) / 1000.0);
			assert (snpData != null);
			int
				numInds = snpData.getNumIndividuals(),
				numSNPs = snpData.getNumSNPs();
			time1 = System.currentTimeMillis();
			System.out.println("Number of SNPs = " + numSNPs);
			System.out.println("Number of individuals in SNP data = " + numInds);
			if (covarMatrix != null && covarMatrix.length != numInds) {
				System.err.println("ERROR: The number of individuals in the SNP data does not match with that of covariate data! " + numInds + " (SNP) vs. " + covarMatrix.length + " (Covar)");
				System.exit(kExitCodeErrorNumIndGenoCovarNotMatch);
			}
			if (covarMatrix != null) {
				for (QSNPData snp: snpData.getSNPs())
					if (snp.getNumMissing() > 0) {
						System.err.println("ERROR: Adjustment with covariate data cannot have missing SNPs. If the SNP data are peppered with missing values, the GPU speed-up will no longer work correctly. You should consider imputation.");
						System.exit(kExitCodeErrorGenoMissingValues);
					}
			}
			System.out.println("Loading expression data...");
			QGeneExpressionData exprData = QGeneExpressionData.load(exprFilename, cCommonDelimiter, "#", true, true);
			time2 = System.currentTimeMillis();
			int numETraits = exprData.getNumberOfRows();
			long gb2 = QSystemUtils.usedMemoryAfterGC();

			// A lot of sanity checks
			System.out.println("Expression data load time (in seconds) = " + (time2 - time1) / 1000.0);
			System.out.println("Memory use for loading (in MB) = " + String.format("%5.2f", (gb2 - gb1) / kMB));
			System.out.println("Total memory use after loading (in MB) = " + String.format("%5.2f", gb2 / kMB));
			System.out.println("Number of e-traits = " + numETraits);
			System.out.println("Number of individuals in expression data = " + exprData.getNumberOfColumns());
			if (exprData.getNumberOfColumns() != numInds) {
				System.err.println("ERROR: The number of individuals in the SNP data does not match with that of expression data! " + numInds + " (SNP) vs. " + exprData.getNumberOfColumns() + " (Expr)");
				System.exit(kExitCodeErrorNumIndGenoExprNotMatch);
			}

			if (globalBlockSize <= 0) {
				globalBlockSize = getDefaultBlockSize(numInds, numSNPs, numETraits);
				System.out.println("The block_size parameter is not specified or is invalid, trying " + globalBlockSize);
				config.setBlockSize(globalBlockSize);
			} else {
				System.out.println("Block size = " + globalBlockSize);
			}
			validateBlockCapacity(globalBlockSize);

			long numReqdIter = ((numSNPs + (long) globalBlockSize - 1) / globalBlockSize)
				* ((numETraits + (long) globalBlockSize - 1) / globalBlockSize);
			System.out.println("Given the block size, " + numReqdIter + " iteration(s) are needed");
			configureThreadCount(numReqdIter, numDevices, false);

			int
				dfr = covarRank + (isAdditive ? 1 : 2) - 1,
				dfe = numInds - dfr - 1;
			System.out.println("Degrees of Freedom for Regression = " + dfr);
			System.out.println("Degrees of Freedom for Errors = " + dfe);
			if (dfOffset != 0) {
				System.out.println("Degrees of Freedom for Offset = " + dfOffset);
			}
			long numPairs = ((long) numSNPs) * numETraits;
			if (thresholdType.equals("none") && numPairs * 8 < gb2) {
				System.out.println("WARNING: You chose to store all of the results. You should be aware that there are " + numPairs + " pairs to be analyzed.");
				System.out.println("You will need at least " + String.format("%5.2f", (numPairs * 8) / kMB) + " MB of RAM,");
				System.out.println("But you appear to have only " + String.format("%5.2f", (gb2 - gb1) / kMB) + " MB of free RAM left.");
				System.out.println("At some point, this program will exhaust all memory and will throw an out of memory error. Be warned.");
			}

			/*
			 * R^2 = T^2 / (DF + T^2)
			 * T = sqrt(DF / (1/R^2 - 1))
			 */
			if (thresholdType.equals("pval")) {
				double tval = T.quantile(t0/2.0, dfe - dfOffset, false, false);
				tval *= tval;
				double rsq = tval / (dfe + tval);
				System.out.println("P-value of " + t0 + " corresponds to R^2 of " + rsq);
				t0 = rsq;
			}

			// Analysis begins
			time1 = System.currentTimeMillis();
			plugin.eSNPAnalysis(snpData, exprData, covarMatrix, t0, dfOffset, dfe, isAdditive);
			time2 = System.currentTimeMillis();
			System.out.println("Total analysis time (in seconds) = " + (time2 - time1) / 1000.0);
			//time1 = System.currentTimeMillis();
			//dumpToString(results, snpData, exprData, outputFilename);
			//time2 = System.currentTimeMillis();
			//System.out.println("Total write time (in seconds) = " + (time2 - time1) / 1000.0);
			//System.out.println(data.mPopn);
			//System.out.println(data);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(kExitCodeErrorAnalysisFailure);
			return;
		}
		long timeb = System.currentTimeMillis();
		System.out.println("Overall time (in seconds) = " + (timeb - timea) / 1000.0);
		System.exit(kExitCodeNormal);
	}

	private static void finishProfiling()
	{
		if (!profiler.isEnabled())
			return;
		profiler.printSummary(System.out);
		String output = config.getProfileOutputFilename();
		if (output != null) {
			try {
				profiler.writeCsv(Path.of(output));
				System.out.println("Profiling CSV: " + Path.of(output).toAbsolutePath().normalize());
			} catch (IOException e) {
				throw new RuntimeException("Cannot write profiling CSV", e);
			}
		}
	}
}
