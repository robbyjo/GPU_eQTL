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
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Deque;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.List;
import java.util.Locale;
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
import gov.nih.eqtl.io.QGenomicRegions;
import gov.nih.eqtl.io.QInMemoryPreparedMatrix;
import gov.nih.eqtl.io.QLocalPatternImputedSource;
import gov.nih.eqtl.io.QMatrixRowSource;
import gov.nih.eqtl.io.QMissingnessReport;
import gov.nih.eqtl.io.QMissingnessScan;
import gov.nih.eqtl.io.QPolicyMatrixSource;
import gov.nih.eqtl.io.QPatternVariantSource;
import gov.nih.eqtl.io.QRawMatrixCache;
import gov.nih.eqtl.io.QPreparedMatrix;
import gov.nih.eqtl.io.QSampleAlignment;
import gov.nih.eqtl.io.QSampleAlignmentPolicy;
import gov.nih.eqtl.io.QVariantMatrixSource;
import gov.nih.eqtl.settest.QSetTestMethod;
import gov.nih.eqtl.settest.QSetTestPolicy;
import gov.nih.eqtl.settest.QSetTestPolicy.FailurePolicy;
import gov.nih.eqtl.settest.QSetTestRunner;
import gov.nih.eqtl.settest.QVariantSetTable;
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
	private static final String TEST_FAIL_AFTER_TRAIT_PATTERN_PROPERTY =
		"eqtl.test.trait.pattern.fail.after";
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
			QAnalysisProgress progress = new QAnalysisProgress("Association",
				(long) numSNPs * numETraits, 0);
			try {
			List<Future<?>> jobs = new ArrayList<Future<?>>(numThreads);
			if (isCategoricalSNP) {
				//for (int i = 0; i < numThreads; i++)
				//jobs.add(threadPool.submit(new QeQTLSNPJobCat(popn, expData, covarQ, gpuContextsPool, results,
				//	numETraitsPerBlock, numESNPsPerBlock, localBlockSize, rsq0, isAdditive, counter)));
			} else {
				for (int i = 0; i < numThreads; i++)
					jobs.add(threadPool.submit(new QeQTLSNPJobReal(popn, expDataTbl, expDataSD, snpDataSD, gpuContextsPool,
						numETraitsPerBlock, numESNPsPerBlock, localBlockSize, dfo, dfe, rsq0,
						isAdditive, gpuPrecision, fw, counter, progress)));
			}
			threadPool.shutdown();
			for (Future<?> job : jobs) {
				job.get();
			}
			if (!threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS)) {
				throw new RuntimeException("Timed out while waiting for eQTL workers");
			}
			progress.complete();
			fw.flush();
			} finally {
				progress.close();
			}
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
		QPreparedMatrix expressionCache, double rsq0, int dfo, int dfe,
		int genotypeRowsPerBlock, int expressionRowsPerBlock)
	{
		eSNPAnalysisStreamed(genotypeCache, expressionCache, rsq0, dfo, dfe,
			genotypeRowsPerBlock, expressionRowsPerBlock,
			Path.of(config.getOutputFilename()),
			config.getCheckpointDirectory() == null
				? Path.of(config.getOutputFilename() + ".checkpoint")
				: Path.of(config.getCheckpointDirectory()),
			config.getResume(), config.getKeepCheckpoints(), true, false, "Association");
	}

	private void eSNPAnalysisStreamed(QBinaryMatrixCache genotypeCache,
		QPreparedMatrix expressionCache, double rsq0, int dfo, int dfe,
		int genotypeRowsPerBlock, int expressionRowsPerBlock,
		Path outputPath, Path checkpointDirectory, boolean resume,
		boolean keepCheckpoints, boolean closeGpuContexts, boolean includeSampleStatistics,
		String progressLabel)
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
		long totalComparisons = genotypeCache.rowCount() * expressionCache.rowCount();
		long completedComparisons = 0;
		for (int blockNumber = 0; blockNumber < totalBlocks; blockNumber++) {
			if (checkpoint.isComplete(blockNumber)) {
				long rows = Math.min(genotypeRowsPerBlock,
					genotypeCache.rowCount() - (long) blockNumber * genotypeRowsPerBlock);
				completedComparisons += rows * expressionCache.rowCount();
			}
		}
		QAnalysisProgress progress = new QAnalysisProgress(progressLabel,
			totalComparisons, completedComparisons);
		if (completedAtStart > 0)
			System.out.println("Completed checkpoint blocks found: " + completedAtStart + " / " + totalBlocks);
		String outputHeader = rsqOnly ? "Rs_ID,ProbesetID,RSq,Dir" : "Rs_ID,ProbesetID,RSq,Fx,T,log10P";
		if (includeSampleStatistics)
			outputHeader += ",N,DF";
		if (completedAtStart == totalBlocks) {
			try {
				progress.complete();
				long assemblyStarted = profiler.start();
				checkpoint.assemble(outputPath, outputHeader);
				profiler.record(QeQTLProfiler.Phase.OUTPUT_ASSEMBLY, assemblyStarted, totalBlocks, 0);
				return;
			} catch (IOException e) {
				throw new RuntimeException("Cannot assemble completed checkpoint", e);
			} finally {
				progress.close();
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
			progress.close();
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
					progress, includeSampleStatistics)));
				if (pending.size() >= numThreads)
					pending.removeFirst().get();
			}
			while (!pending.isEmpty())
				pending.removeFirst().get();
			threadPool.shutdown();
			if (!threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS))
				throw new RuntimeException("Timed out while waiting for streamed eQTL workers");
			progress.complete();
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
			progress.close();
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
		System.out.println("Compute backend: " + gpuRuntime.getBackend().getName());
		System.out.println("Compute runtime: " + gpuRuntime.getBackend().getRuntimeDescription());
		List<GpuDevice> devices;
		try {
			devices = gpuRuntime.getGpuDevices(true, gpuPrecision.requiresDoublePrecision());
		} catch (Throwable e) {
			System.err.println("Cannot initialize the compute backend: " + e.getMessage());
			System.exit(kExitCodeErrorInitOpenCLFailure);
			return;
		}
		int numDevices = devices.size();
		if (numDevices == 0) {
			System.err.println("Cannot find an available compute device supporting " + gpuPrecision.optionName() + ".");
			List<GpuDevice> allGpuDevices = gpuRuntime.getGpuDevices(true, false);
			if (gpuPrecision == GpuPrecision.FP64 && !allGpuDevices.isEmpty()) {
				System.err.println("The following devices were detected but do not report double-precision support:");
				for (GpuDevice device : allGpuDevices) {
					System.err.println("  " + device.getName() + " (" + device.getComputeApiVersion() + ")");
				}
			} else {
				System.err.println("No usable GPU or CPU compute backend was reported.");
			}
			System.exit(kExitCodeError64bitGPUNotFound);
			return;
		}
		boolean cpuOnly = devices.stream().allMatch(device -> "cpu".equalsIgnoreCase(device.getBackendName()));
		if (cpuOnly && "auto".equalsIgnoreCase(gpuRuntime.getBackend().getName()))
			System.err.println("WARNING: No suitable GPU was selected; automatically falling back to CPU. "
				+ "The results remain compatible, but a large analysis may be substantially slower.");
		else if (cpuOnly)
			System.out.println("CPU compute was selected explicitly; results use the same association pipeline.");
		System.out.println("Found " + numDevices + " suitable compute device" + (numDevices > 1 ? "s": "") + ":");
		boolean hasHostUnifiedMemory = false;
		for (GpuDevice device : devices) {
			String memoryLabel = "cpu".equalsIgnoreCase(device.getBackendName()) ? "JVM heap limit=" : "VRAM=";
			System.out.println(device.getName() + " [" + device.getVendor() + ", " + device.getComputeApiVersion()
				+ ", " + memoryLabel + String.format("%.2f GiB", device.getGlobalMemoryBytes() / kGB) + "]");
			hasHostUnifiedMemory |= device.hasUnifiedMemory()
				&& !"cpu".equalsIgnoreCase(device.getBackendName());
		}
		if (hasHostUnifiedMemory) {
			System.err.println("NOTE: At least one GPU reports unified host memory. Performance depends on the device and driver; this is not necessarily an error.");
		}
		mContexts = new GpuContext[numDevices];
		try {
			for (int i = 0; i < numDevices; i++) {
				mContexts[i] = devices.get(i).openContext();
				if ("cpu".equalsIgnoreCase(devices.get(i).getBackendName()))
					System.out.println("CPU matrix engine: " + devices.get(i).getComputeApiVersion());
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
			System.err.println("WARNING: The detected compute devices are not identical. Performance may be constrained by the weakest device.");
		} else {
			System.out.println("All the selected compute devices appear to be identical.");
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
			throw new IllegalStateException("Compute contexts are required for automatic block sizing");
		List<GpuDevice> devices = new ArrayList<GpuDevice>(mContexts.length);
		for (GpuContext context : mContexts)
			devices.add(context.getDevice());
		boolean cpu = devices.stream().allMatch(device -> "cpu".equalsIgnoreCase(device.getBackendName()));
		GpuTuning.BlockRecommendation recommendation;
		if (cpu) {
			Runtime runtime = Runtime.getRuntime();
			long availableHeap = runtime.maxMemory() - (runtime.totalMemory() - runtime.freeMemory());
			recommendation = GpuTuning.recommendCpuBlockSize(numInds, numSNPs, numETraits,
				gpuPrecision, availableHeap);
		} else {
			recommendation = GpuTuning.recommendBlockSize(
				devices, numInds, numSNPs, numETraits, gpuPrecision);
		}
		System.out.println("Automatic block-size limiter: " + recommendation.limitingDevice()
			+ "; estimated " + (cpu ? "host working buffers=" : "device buffers=")
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
		boolean cpu = mContexts != null && mContexts.length > 0
			&& Arrays.stream(mContexts).allMatch(context ->
				"cpu".equalsIgnoreCase(context.getDevice().getBackendName()));
		if (threads == 0) {
			threads = cpu
				? GpuTuning.recommendCpuWorkerCount(requiredIterations, streamed,
					availableHeapBytes, estimatedWorkerBytes)
				: GpuTuning.recommendThreadCount(QSystemUtils.kNumCPUCores,
					numDevices, requiredIterations, streamed, availableHeapBytes, estimatedWorkerBytes);
			config.setNumThreads(threads);
			System.out.println("The thread count was not specified; using " + threads
				+ " (" + numDevices + " compute context" + (numDevices == 1 ? "" : "s") + ")");
		}
		if (threads > QSystemUtils.kNumCPUCores)
			System.err.println("WARNING: num_threads exceeds the available CPU cores ("
				+ QSystemUtils.kNumCPUCores + ") and may reduce performance.");
		if (!cpu && threads < Math.min((long) numDevices, requiredIterations))
			System.err.println("WARNING: num_threads is lower than the number of usable GPUs; some GPUs may remain idle.");
		if (cpu && threads > 2)
			System.err.println("WARNING: The CPU backend has one exclusive BLAS context; more than two pipeline "
				+ "workers usually add memory pressure without increasing matrix-product throughput.");
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

	static void runMatrixAnalysis(QeQTLAnalysis plugin, String genotypeFilename,
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
		QSampleAlignmentPolicy sampleAlignmentPolicy = config.getSampleAlignmentPolicy();
		QTraitCacheMode traitCacheMode = config.getTraitCacheMode();
		boolean preprocessOnly = config.getPreprocessOnly();
		String cohortModelFilename = config.getCohortModelFilename();
		String cohortColumn = config.getCohortColumn();
		if ((cohortModelFilename == null) != (cohortColumn == null))
			throw new IllegalArgumentException("--cohort-model and --cohort-column must be specified together");
		if (cohortModelFilename != null && covariateFilename == null)
			throw new IllegalArgumentException("Cohort-aware analysis requires --covariates");
		if (cohortModelFilename != null && config.getAnalysisMethod().isSetTest())
			throw new IllegalArgumentException("Cohort-aware covariance is not yet supported for burden/SKAT/SKAT-O analyses");
		QCohortModel.Definitions cohortDefinitions = cohortModelFilename == null ? null
			: QCohortModel.Definitions.load(Path.of(cohortModelFilename));
		if (preprocessOnly && (config.getValidateOnly() || config.getInspectMissingness()))
			throw new IllegalArgumentException("--preprocess-only cannot be combined with --validate-only or --inspect-missingness");
		if (!preprocessOnly && !config.getValidateOnly() && !config.getInspectMissingness()
			&& config.getOutputFilename() == null)
			throw new IllegalArgumentException("--output is required for association analysis");
		if (!genotypeFormat.equals("csv") && !predictorType.isGenotype())
			throw new IllegalArgumentException("VCF/BCF input requires --predictor-type genotype");
		if (predictorMissing == QMissingValuePolicy.LOCAL_PATTERN && !predictorType.isGenotype())
			throw new IllegalArgumentException("local-pattern imputation is only valid for --predictor-type genotype");
		if (predictorMissing == QMissingValuePolicy.PATTERN)
			throw new IllegalArgumentException("Pattern-wise predictor deletion is not yet supported; use mean, local-pattern, error, zero, or exclude-row");
		if (!(covariateMissing.equals("error") || covariateMissing.equals("complete-samples")))
			throw new IllegalArgumentException("covariate_missing must be error or complete-samples");
		QVariantMatrixSource.FrequencyScope frequencyScope =
			QVariantMatrixSource.FrequencyScope.parse(config.getFrequencyScope());
		boolean regionsRequested = (config.getRegions() != null && !config.getRegions().isBlank())
			|| config.getRegionsFilename() != null;
		boolean windowOptionsPresent = config.get("window_size") != null
			|| config.get("window_stride") != null;
		int windowSize = config.getWindowSize();
		int windowStride = config.getWindowStride();
		if (windowOptionsPresent) {
			if (windowSize < 1)
				throw new IllegalArgumentException("--window-size must be positive when sliding windows are requested");
			if (windowStride < 1 || windowStride > windowSize)
				throw new IllegalArgumentException("--window-stride must be positive and no larger than --window-size");
			if (!config.getAnalysisMethod().isSetTest())
				throw new IllegalArgumentException("--window-size/--window-stride require --analysis burden, skat, or skat-o");
			if (config.getVariantSetsFilename() != null || regionsRequested)
				throw new IllegalArgumentException("Automatic sliding windows cannot be combined with --variant-sets, --region, or --regions-file");
		}
		if (genotypeFormat.equals("csv") && regionsRequested)
			throw new IllegalArgumentException("--region/--regions-file requires VCF or BCF genotype input");
		if (genotypeFormat.equals("csv") && frequencyScope != QVariantMatrixSource.FrequencyScope.ALIGNED)
			throw new IllegalArgumentException("--frequency-scope pattern currently requires VCF or BCF genotype input");
		if (genotypeFormat.equals("csv") && preprocessOnly)
			throw new IllegalArgumentException("--preprocess-only requires VCF or BCF genotype input");
		QMatrixRowSource genotypeSource;
		QVariantMatrixSource variantSource = null;
		Path variantQcOutput = null;
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
				config.getMinimumMaf(), config.getMinimumMac(), qcOutput,
				config.getVariantQcThreads(), config.getVariantQcCheckpointDirectory() == null
					? null : Path.of(config.getVariantQcCheckpointDirectory()),
				config.getVariantIndexFilename() == null ? null : Path.of(config.getVariantIndexFilename()),
				config.getRegions(),
				config.getRegionsFilename() == null ? null : Path.of(config.getRegionsFilename()),
				QGenomicRegions.Coordinates.parse(config.getRegionCoordinates()), frequencyScope);
			variantQcOutput = qcOutput.toAbsolutePath().normalize();
			double minimumMaf = config.getMinimumMaf();
			double minimumMac = config.getMinimumMac();
			System.out.println("Variant frequency filters: MAF "
				+ (minimumMaf == 0 ? "disabled" : ">= " + minimumMaf)
				+ ", MAC " + (minimumMac == 0 ? "disabled" : ">= " + minimumMac
					+ (config.get("min_mac") == null ? " (default)" : ""))
				+ ", scope=" + frequencyScope.name().toLowerCase(Locale.ROOT));
			variantSource = QVariantMatrixSource.openForAlignment(
				Path.of(genotypeFilename), variantOptions);
			genotypeSource = variantSource;
		}
		QMatrixRowSource expressionSource = new QDelimitedMatrixSource(
			Path.of(expressionFilename), cCommonDelimiter, "#");
		String[] sourceGenotypeSampleIds = genotypeSource.metadata().sampleIds();
		String[] sourceExpressionSampleIds = expressionSource.metadata().sampleIds();

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
		QCohortModel cohortModel = null;
		double[][] covariateQ = null;
		double[][] covariateModel = null;
		int covariateRank = 1;
		if (covariateFilename != null) {
			System.out.println("Loading and validating covariate data...");
			covariates = QCovariateTable.load(Path.of(covariateFilename), cCommonDelimiter, "#");
			fullAlignment = covariates.align(sourceGenotypeSampleIds,
				sourceExpressionSampleIds, config.getGenotypeIdColumn(),
				config.getExpressionIdColumn(), sampleAlignmentPolicy,
				config.getGenotypeIdStripPrefix(), config.getExpressionIdStripPrefix());
			canonicalSampleIds = reorder(sourceGenotypeSampleIds, fullAlignment.genotypeColumnOrder());
			String[] missingnessTerms = cohortDefinitions == null ? fixedCovariates
				: cohortDefinitions.missingnessTerms(cohortColumn, fixedCovariates);
			covariateMissingness = covariates.inspectMissingness(missingnessTerms);
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
			if (sampleAlignmentPolicy != QSampleAlignmentPolicy.STRICT)
				throw new IllegalArgumentException("covariate-subset alignment requires --covariates and explicit canonical rows");
			fullAlignment = QCovariateTable.alignDirectly(sourceGenotypeSampleIds,
				sourceExpressionSampleIds, config.getGenotypeIdStripPrefix(),
				config.getExpressionIdStripPrefix());
			alignment = fullAlignment;
			canonicalSampleIds = reorder(sourceGenotypeSampleIds, alignment.genotypeColumnOrder());
		}
		if (variantSource != null) {
			String qcThreadSelection = config.getVariantQcThreads() == 0 ? "automatic" : "explicit";
			System.out.println("Scanning aligned-sample variant QC for " + alignment.sampleCount()
				+ " sample(s) in " + Path.of(genotypeFilename).toAbsolutePath().normalize()
				+ "; progress will be reported periodically...");
			System.out.println("Variant QC workers = " + variantSource.qcThreadCount()
				+ " (" + qcThreadSelection + "; ordered output)");
			System.out.flush();
			long variantQcStarted = profiler.start();
			variantSource.selectAnalysisSamples(alignment.genotypeColumnOrder());
			QVariantMatrixSource.Summary variantSummary = variantSource.summary();
			profiler.record(QeQTLProfiler.Phase.VARIANT_QC, variantQcStarted,
				variantSummary.inputRecords(), 0);
			System.out.println("Variant genotype field = " + variantSource.selectedField());
			System.out.println("Variant QC samples = " + variantSource.analysisSampleCount()
				+ " aligned sample(s)");
			System.out.println("Variant QC: input=" + variantSummary.inputRecords()
				+ ", included=" + variantSummary.includedVariants()
				+ ", monomorphic=" + variantSummary.monomorphicVariants()
				+ ", singletons=" + variantSummary.singletons()
				+ ", doubletons=" + variantSummary.doubletons());
			System.out.println("Variant annotation/QC output: " + variantQcOutput);
			if (variantSource.variantIndex() != null)
				System.out.println("Variant index: " + variantSource.variantIndex()
					+ (variantSource.regions() == null ? " (available for indexed resume)" : " (indexed region access)"));
			if (variantSource.regions() != null) {
				System.out.println("Genomic region sets = " + variantSummary.regionSets()
					+ "; merged indexed queries = " + variantSource.regions().queryIntervals().size());
				if (!variantSummary.emptyRegionSets().isEmpty())
					System.err.println("WARNING: Region sets with no source variant: "
						+ String.join(", ", variantSummary.emptyRegionSets()));
			}
			System.out.println("Variant QC checkpoint: " + variantSource.qcCheckpointDirectory()
				+ (variantSource.resumedQcRecords() > 0
					? " (reused " + variantSource.resumedQcRecords() + " records)" : ""));
		}
		if (genotypeSource.metadata().rowCount() > Integer.MAX_VALUE
			|| expressionSource.metadata().rowCount() > Integer.MAX_VALUE)
			throw new IllegalArgumentException("This release supports at most " + Integer.MAX_VALUE + " rows per matrix");
		int[] predictorColumnOrder = alignment.genotypeColumnOrder();
		QRawMatrixCache preprocessedPredictors = null;
		if (variantSource != null) {
			Path cacheRoot = matrixCacheDirectory(genotypeFilename);
			String rawSignature = QRawMatrixCache.signature(variantSource,
				alignment.genotypeColumnOrder());
			if (preprocessOnly) {
				int rowsPerBlock = config.getGenotypeBlockRows() > 0
					? config.getGenotypeBlockRows() : 2048;
				preprocessedPredictors = QRawMatrixCache.openOrBuild(
					cacheRoot.resolve("aligned-raw"), rawSignature, variantSource,
					alignment.genotypeColumnOrder(), rowsPerBlock, config.getRebuildCache());
			} else {
				preprocessedPredictors = QRawMatrixCache.openIfPresent(
					cacheRoot.resolve("aligned-raw"), rawSignature, variantSource);
			}
			if (preprocessedPredictors != null) {
				genotypeSource = preprocessedPredictors;
				predictorColumnOrder = identity(alignment.sampleCount());
			}
		}
		if (windowOptionsPresent && preprocessedPredictors == null) {
			Path cacheRoot = matrixCacheDirectory(genotypeFilename);
			String rawSignature = QRawMatrixCache.signature(genotypeSource,
				predictorColumnOrder);
			int rowsPerBlock = config.getGenotypeBlockRows() > 0
				? config.getGenotypeBlockRows() : 2048;
			preprocessedPredictors = QRawMatrixCache.openOrBuild(
				cacheRoot.resolve("aligned-raw"), rawSignature, genotypeSource,
				predictorColumnOrder, rowsPerBlock, config.getRebuildCache());
			genotypeSource = preprocessedPredictors;
			predictorColumnOrder = identity(alignment.sampleCount());
		}
		QMatrixRowSource unfilledPredictorSource = genotypeSource;

		System.out.println("Inspecting matrix missingness...");
		Path missingnessCacheDirectory = matrixCacheDirectory(genotypeFilename).resolve("missingness");
		QMissingnessScan predictorScan = QMissingnessScan.scanOrLoad("predictor", genotypeSource,
			predictorColumnOrder, missingnessCacheDirectory, config.getRebuildCache());
		QMissingnessScan traitScan = QMissingnessScan.scanOrLoad("trait", expressionSource,
			alignment.expressionColumnOrder(), missingnessCacheDirectory, config.getRebuildCache());
		Path missingnessOutput = missingnessOutput(genotypeFilename);
		if (config.getOutputFilename() != null && missingnessOutput.equals(
			Path.of(config.getOutputFilename()).toAbsolutePath().normalize()))
			throw new IllegalArgumentException("Missingness QC output must differ from association output");
		QMissingnessReport.write(missingnessOutput, predictorScan, predictorType, predictorMissing,
			traitScan, traitType, traitMissing, covariateMissingness, covariateMissing,
			canonicalSampleIds, fullAlignment);
		canonicalSampleIds = reorder(sourceGenotypeSampleIds, alignment.genotypeColumnOrder());
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

		boolean dynamicTraitPatterns = traitMissing == QMissingValuePolicy.PATTERN
			&& (traitScan.hasMissingValues()
				|| frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN);
		if (frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN
			&& traitMissing != QMissingValuePolicy.PATTERN)
			throw new IllegalArgumentException("--frequency-scope pattern requires --trait-missing pattern");
		if (!preprocessOnly && predictorMissing == QMissingValuePolicy.LOCAL_PATTERN) {
			System.err.println("WARNING: local-pattern is a nearest flanking-genotype proxy; it is not "
				+ "phasing or reference-panel imputation. Flanks per side=" + config.getPredictorFlankCount());
			genotypeSource = new QLocalPatternImputedSource(genotypeSource, predictorScan,
				config.getPredictorFlankCount());
		} else if (!preprocessOnly && !(dynamicTraitPatterns && variantSource != null)) {
			genotypeSource = new QPolicyMatrixSource(genotypeSource, predictorScan, predictorMissing);
		}
		if (!preprocessOnly && !(traitMissing == QMissingValuePolicy.PATTERN && traitScan.hasMissingValues()))
			expressionSource = new QPolicyMatrixSource(expressionSource, traitScan, traitMissing);

		if (covariates != null) {
			double[][] modelValues;
			String[] modelColumnNames;
			if (cohortDefinitions == null) {
				QCovariateTable.ModelMatrix model = covariates.buildModelMatrix(
					fixedCovariates, factorCovariates, retainedCovariateRows);
				if (model.automaticFactors().length > 0)
					System.out.println("Automatically categorical covariates: " + String.join(", ", model.automaticFactors()));
				modelValues = model.values();
				modelColumnNames = model.columnNames();
			} else {
				cohortModel = QCohortModel.create(cohortDefinitions, covariates, cohortColumn,
					retainedCovariateRows, canonicalSampleIds, fixedCovariates, factorCovariates);
				modelValues = cohortModel.design();
				modelColumnNames = cohortModel.columnNames();
				Path audit = config.getOutputFilename() == null
					? Path.of(cohortModelFilename + ".audit.tsv")
					: Path.of(config.getOutputFilename() + ".cohorts.tsv");
				cohortModel.writeAudit(audit);
				System.out.println("Cohort model audit: " + audit.toAbsolutePath().normalize());
			}
			System.out.println("Covariate model columns: " + String.join(", ", modelColumnNames));
			QRDecomposition decomposition = new QRDecomposition(modelValues);
			if (!decomposition.isFullRank())
				throw new IllegalArgumentException("Covariate model is rank deficient: rank "
					+ decomposition.getRank() + " for " + decomposition.getNumColumns() + " columns");
			covariateRank = decomposition.getRank();
			covariateQ = decomposition.getQ().getArray();
			covariateModel = modelValues;
		}
		QSetTestMethod setTestMethod = config.getAnalysisMethod();
		if (setTestMethod.isSetTest()) {
			QVariantSetTable setDefinitions;
			Path explicitSetPath = config.getVariantSetsFilename() == null
				? null : Path.of(config.getVariantSetsFilename());
			if (windowOptionsPresent) {
				int scanRows = config.getGenotypeBlockRows() > 0
					? config.getGenotypeBlockRows() : 1024;
				setDefinitions = QVariantSetTable.fromSlidingWindows(unfilledPredictorSource,
					predictorColumnOrder, windowSize, windowStride, scanRows);
				System.out.println("Set definitions: " + setDefinitions.sets().size()
					+ " nonempty automatic sliding windows; size=" + windowSize
					+ " bp, stride=" + windowStride + " bp, one-based grid anchored at 1");
			} else if (explicitSetPath != null) {
				setDefinitions = QVariantSetTable.load(explicitSetPath);
			} else if (variantSource != null && variantSource.regions() != null) {
				setDefinitions = QVariantSetTable.fromVariantQc(variantQcOutput,
					variantSource.regions().setIds());
				System.out.println("Set definitions: exact ALT-effect memberships adapted from aligned VCF/BCF region QC");
			} else {
				throw new IllegalArgumentException("--variant-sets or --window-size is required for --analysis "
					+ setTestMethod.optionName() + " unless indexed VCF/BCF regions define the sets");
			}
			if (traitMissing == QMissingValuePolicy.PATTERN && traitScan.hasMissingValues())
				throw new IllegalArgumentException("Set-test production execution currently requires complete, mean-filled, zero-filled, or excluded traits; exact trait patterns are not yet supported");
			if (gpuPrecision != GpuPrecision.FP64)
				throw new IllegalArgumentException("Set-test production execution currently requires --precision fp64");
			double[][] setDesign = covariateModel;
			if (setDesign == null) {
				setDesign = new double[alignment.sampleCount()][1];
				for (double[] row : setDesign) row[0] = 1.0;
			}
			QSetTestPolicy setPolicy = new QSetTestPolicy(config.getSetMinimumMaf(),
				config.getSetMaximumMaf(), config.getSetMinimumMac(), config.getSetMaximumMac(),
				predictorMissing, FailurePolicy.parse(config.getSetAbsentVariantPolicy(),
					FailurePolicy.ERROR), FailurePolicy.parse(config.getSetDegeneratePolicy(),
					FailurePolicy.ERROR));
			Path setOutput = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
			Path setAudit = config.getSetAuditFilename() == null
				? Path.of(config.getOutputFilename() + ".sets.tsv").toAbsolutePath().normalize()
				: Path.of(config.getSetAuditFilename()).toAbsolutePath().normalize();
			Path setCheckpoint = config.getCheckpointDirectory() == null
				? Path.of(config.getOutputFilename() + ".checkpoint").toAbsolutePath().normalize()
				: Path.of(config.getCheckpointDirectory()).toAbsolutePath().normalize();
			int setTraitRows = config.getExpressionBlockRows() > 0
				? config.getExpressionBlockRows() : 1024;
			QSetTestRunner.run(unfilledPredictorSource, predictorColumnOrder, expressionSource,
				alignment.expressionColumnOrder(), setDesign, setDefinitions,
				new QSetTestRunner.Options(setTestMethod,
					explicitSetPath, setOutput, setAudit, setCheckpoint,
					setTraitRows, setPolicy, config.getResume(), config.getKeepCheckpoints(),
					thresholdType, threshold, config.getSetBlockSize(), config.getSkatORhoGrid(),
					config.getSkatOSimulations(), config.getSkatOSeed()));
			return;
		}

		long genotypeRowCount = genotypeSource.metadata().rowCount();
		long expressionRowCount = traitMissing == QMissingValuePolicy.PATTERN && traitScan.hasMissingValues()
			? traitScan.rowCount() : expressionSource.metadata().rowCount();
		int numIndividuals = alignment.sampleCount();
		int numSnps = (int) genotypeRowCount;
		int numTraits = (int) expressionRowCount;
		System.out.println("Sample alignment: genotype via " + alignment.genotypeIdColumn()
			+ ", expression via " + alignment.expressionIdColumn()
			+ ", policy=" + alignment.policy().optionName());
		System.out.println("Samples reordered: genotype=" + alignment.genotypeReorderedCount()
			+ ", expression=" + alignment.expressionReorderedCount());
		if (fullAlignment.genotypeIdsStripped() > 0 || fullAlignment.expressionIdsStripped() > 0)
			System.out.println("Sample-ID prefixes stripped: genotype=" + fullAlignment.genotypeIdsStripped()
				+ " ('" + fullAlignment.genotypeIdStripPrefix() + "'), expression="
				+ fullAlignment.expressionIdsStripped() + " ('"
				+ fullAlignment.expressionIdStripPrefix() + "')");
		if (fullAlignment.genotypeExtraSampleCount() > 0 || fullAlignment.expressionExtraSampleCount() > 0)
			System.out.println("Matrix-only samples excluded by covariate-subset alignment: genotype="
				+ fullAlignment.genotypeExtraSampleCount() + ", expression="
				+ fullAlignment.expressionExtraSampleCount());
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
		if (preprocessOnly) {
			System.out.println("Preprocessed aligned variant cache: " + preprocessedPredictors.path());
			System.out.println("Preprocessing completed successfully; association analysis was not run.");
			System.out.println("Reuse this cache by keeping the same --cache-dir, sample alignment, "
				+ "regions, genotype field, variant-inclusion policy, and MAF/MAC settings.");
			return;
		}
		if (config.getValidateOnly()) {
			if (configuredBlockSize <= 0 || config.getNumThreads() == 0)
				System.out.println("Automatic block/thread tuning is deferred until a real analysis run.");
			System.out.println("Validation completed successfully; --validate-only requested, so no analysis was run.");
			return;
		}

		int globalBlockSize = configuredBlockSize;
		if (globalBlockSize <= 0) {
			globalBlockSize = getDefaultBlockSize(numIndividuals, numSnps, numTraits);
			if (globalBlockSize <= 0)
				throw new IllegalArgumentException("Cannot infer a positive block size; specify --block-size");
			int explicitStreamCapacity = Math.max(config.getGenotypeBlockRows(),
				config.getExpressionBlockRows());
			if (explicitStreamCapacity > globalBlockSize) {
				globalBlockSize = roundUpNearestMultiple(explicitStreamCapacity, 16);
				System.out.println("Automatic block capacity was raised to " + globalBlockSize
					+ " to honor the explicitly requested streamed row blocks");
			}
			config.setBlockSize(globalBlockSize);
			System.out.println("The block size was not specified; using " + globalBlockSize);
		} else {
			System.out.println("Block size = " + globalBlockSize);
		}
		validateBlockCapacity(globalBlockSize);

		long numRequiredIterations = ((numSnps + (long) globalBlockSize - 1) / globalBlockSize)
			* ((numTraits + (long) globalBlockSize - 1) / globalBlockSize);
		System.out.println("Given the block size, " + numRequiredIterations + " iteration(s) are needed");
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
			if (cohortModel != null && cohortModel.hasWhitening())
				throw new IllegalArgumentException("Cohort covariance whitening currently requires complete or filled traits; it cannot be combined with exact trait-pattern deletion");
			if (!isAdditive)
				throw new IllegalArgumentException("Pattern-wise trait deletion currently supports the additive model only");
			int maximumPatterns = config.getMaximumTraitPatterns();
			QTraitPatternScheduler requestedScheduler = config.getTraitPatternScheduler();
			QTraitPatternScheduler patternScheduler = requestedScheduler == QTraitPatternScheduler.AUTO
				? ((maximumPatterns == 0 || traitScan.patterns().size() <= maximumPatterns)
					? QTraitPatternScheduler.PATTERN : QTraitPatternScheduler.GENOTYPE)
				: requestedScheduler;
			System.out.println("Trait-pattern scheduler = " + patternScheduler.optionName()
				+ (requestedScheduler == QTraitPatternScheduler.AUTO ? " (automatic)" : " (explicit)"));
			if (patternScheduler == QTraitPatternScheduler.PATTERN) {
				System.err.println("WARNING: Pattern-outer exact deletion runs one predictor preparation/compute "
					+ "pass per distinct missingness pattern. Completed groups are restartable. Patterns="
					+ traitScan.patterns().size());
				validateTraitPatternWorkload(traitScan, genotypeRowCount, covariateModel,
					dfOffset, maximumPatterns);
				if (config.getUnestimableTraitPolicy() == QUnestimableTraitPolicy.SKIP)
					throw new IllegalArgumentException("--unestimable-trait-patterns skip currently requires "
						+ "--trait-pattern-scheduler genotype");
			} else {
				if (gpuPrecision != GpuPrecision.FP64)
					throw new IllegalArgumentException("Genotype-outer exact deletion currently requires --precision fp64");
				System.out.println("Genotype-outer exact deletion prepares each predictor block once and "
					+ "uses QR-equivalent pattern sufficient statistics; result order is genotype block, "
					+ "trait block, variant, then original estimable trait row.");
				if (predictorScan.hasMissingValues())
					System.out.println("Predictor missing values are handled inside each exact trait mask using "
						+ predictorMissing.optionName() + " policy.");
				if (frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN)
					System.out.println("Pattern-specific MAF/MAC filters are applied independently inside "
						+ "each exact trait mask.");
			}
			if (variantSource != null) {
				Path cacheRoot = matrixCacheDirectory(genotypeFilename);
				QRawMatrixCache rawPredictors = preprocessedPredictors;
				boolean closeRawPredictors = false;
				if (rawPredictors == null) {
					String rawSignature = QRawMatrixCache.signature(variantSource,
						alignment.genotypeColumnOrder());
					rawPredictors = QRawMatrixCache.openOrBuild(
						cacheRoot.resolve("aligned-raw"), rawSignature, variantSource,
						alignment.genotypeColumnOrder(), genotypeRows, config.getRebuildCache());
					closeRawPredictors = true;
				}
				try {
					if (patternScheduler == QTraitPatternScheduler.GENOTYPE)
						runGenotypeOuterTraitPatternAnalysis(rawPredictors, identity(alignment.sampleCount()),
							expressionSource, traitScan, alignment, covariateModel, thresholdType,
							threshold, dfOffset, genotypeRows, expressionRows, predictorMissing,
							frequencyScope);
					else
						runTraitPatternAnalysis(plugin, rawPredictors, expressionSource, traitScan,
							alignment, covariateModel, thresholdType, threshold, dfOffset, genotypeRows,
							expressionRows, covariateFilename != null, true, true,
							predictorMissing, frequencyScope, config.getMinimumMaf(), config.getMinimumMac());
				} finally {
					if (closeRawPredictors)
						rawPredictors.close();
				}
			} else {
				if (patternScheduler == QTraitPatternScheduler.GENOTYPE) {
					QMatrixRowSource genotypeOuterSource = predictorMissing == QMissingValuePolicy.MEAN
						|| predictorMissing == QMissingValuePolicy.ZERO
						? unfilledPredictorSource : genotypeSource;
				runGenotypeOuterTraitPatternAnalysis(genotypeOuterSource, predictorColumnOrder,
						expressionSource, traitScan, alignment, covariateModel, thresholdType,
						threshold, dfOffset, genotypeRows, expressionRows, predictorMissing,
						frequencyScope);
				} else
					runTraitPatternAnalysis(plugin, genotypeSource, expressionSource, traitScan,
						alignment, covariateModel, thresholdType, threshold, dfOffset, genotypeRows,
						expressionRows, covariateFilename != null, false, false,
						predictorMissing, frequencyScope, config.getMinimumMaf(), config.getMinimumMac());
			}
			profiler.record(QeQTLProfiler.Phase.ANALYSIS_WALL, analysisStarted,
				(long) numSnps * numTraits, 0);
			System.out.println("Total analysis time (in seconds) = "
				+ (System.currentTimeMillis() - start) / 1000.0);
			return;
		}
		QGpuResidualizer gpuResidualizer = null;
		if (covariateQ != null && residualizationMode != QResidualizationMode.CPU) {
			gpuResidualizer = new QGpuResidualizer(mContexts, covariateQ, gpuPrecision, profiler);
			System.out.println("Fixed-effect residualization = selected compute backend (" + gpuPrecision.optionName()
				+ ", " + mContexts.length + " context" + (mContexts.length == 1 ? "" : "s") + ")");
		} else if (covariateQ != null) {
			System.out.println("Fixed-effect residualization = CPU FP64");
		}
		QeQTLPreprocessor.Residualizer preprocessingResidualizer = cohortModel == null
			? gpuResidualizer : cohortModel.wrap(gpuResidualizer);
		try {
		if (stream) {
			System.out.println("Bounded-RAM matrix mode: genotype blocks=" + genotypeRows
				+ " rows, expression blocks=" + expressionRows + " rows");
			Path outputPath = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
			Path cacheDirectory = config.getCacheDirectory() == null
				? outputPath.getParent().resolve(".gpu-eqtl-cache")
				: Path.of(config.getCacheDirectory());
			long signatureStarted = profiler.start();
			String preprocessingTag = preprocessingResidualizer == null ? null
				: preprocessingResidualizer.cacheSignatureTag();
			String genotypeSignature = QBinaryMatrixCache.signature("Genotype", genotypeSource,
				predictorColumnOrder, covariateQ, preprocessingTag);
			String expressionSignature = QBinaryMatrixCache.signature("Expression", expressionSource,
				alignment.expressionColumnOrder(), covariateQ, preprocessingTag);
			profiler.record(QeQTLProfiler.Phase.CACHE_SIGNATURES, signatureStarted, 2, 0);
			long genotypeCacheStarted = profiler.start();
			try (QBinaryMatrixCache genotypeCache = QBinaryMatrixCache.openOrBuild(cacheDirectory,
				"Genotype", genotypeSignature, genotypeSource, predictorColumnOrder, covariateQ,
				genotypeRows, config.getRebuildCache(), preprocessingResidualizer)) {
				profiler.record(QeQTLProfiler.Phase.GENOTYPE_CACHE_OPEN_OR_BUILD, genotypeCacheStarted,
					genotypeCache.rowCount(), java.nio.file.Files.size(genotypeCache.path()));
				long expressionCacheStarted = profiler.start();
				try (QBinaryMatrixCache expressionCache = QBinaryMatrixCache.openOrBuild(cacheDirectory,
				"Expression", expressionSignature, expressionSource, alignment.expressionColumnOrder(), covariateQ,
				expressionRows, config.getRebuildCache(), preprocessingResidualizer)) {
					profiler.record(QeQTLProfiler.Phase.EXPRESSION_CACHE_OPEN_OR_BUILD, expressionCacheStarted,
						expressionCache.rowCount(), java.nio.file.Files.size(expressionCache.path()));
					if (gpuResidualizer != null) {
						gpuResidualizer.close();
						gpuResidualizer = null;
					}
					QPreparedMatrix preparedTraits = selectPreparedTraitMatrix(expressionCache,
						expressionRows, numIndividuals, genotypeRows, traitCacheMode);
					try {
						plugin.eSNPAnalysisStreamed(genotypeCache, preparedTraits, rSquaredThreshold, dfOffset,
							errorDegreesOfFreedom, genotypeRows, expressionRows);
					} finally {
						closeMemoryPreparedMatrix(preparedTraits);
					}
				}
			}
		} else {
			System.out.println("Loading aligned matrices into RAM (set --genotype-block-rows or --expression-block-rows to stream).");
			QMatrixRowSource.Block genotypeBlock;
			QMatrixRowSource.Block expressionBlock;
			try (QMatrixRowSource.BlockReader reader = genotypeSource.open(predictorColumnOrder)) {
				genotypeBlock = reader.readBlock(numSnps);
			}
			try (QMatrixRowSource.BlockReader reader = expressionSource.open(alignment.expressionColumnOrder())) {
				expressionBlock = reader.readBlock(numTraits);
			}
			QeQTLPreprocessor.PreparedBlock preparedGenotypes = QeQTLPreprocessor.prepare(
				genotypeBlock, covariateQ, "Genotype", preprocessingResidualizer);
			QeQTLPreprocessor.PreparedBlock preparedExpressions = QeQTLPreprocessor.prepare(
				expressionBlock, covariateQ, "Expression", preprocessingResidualizer);
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

	private static void runGenotypeOuterTraitPatternAnalysis(QMatrixRowSource predictorSource,
		int[] predictorColumns, QMatrixRowSource rawTraitSource, QMissingnessScan traitScan,
		QSampleAlignment alignment, double[][] covariateModel, String thresholdType,
		double threshold, int dfOffset, int predictorRowsPerBlock,
		int traitRowsPerBlock, QMissingValuePolicy predictorMissing,
		QVariantMatrixSource.FrequencyScope frequencyScope) throws Exception {
		Path output = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
		Path parent = output.getParent() == null
			? Path.of(".").toAbsolutePath().normalize() : output.getParent();
		Path cacheRoot = config.getCacheDirectory() == null
			? parent.resolve(".gpu-eqtl-cache")
			: Path.of(config.getCacheDirectory()).toAbsolutePath().normalize();
		QTraitPatternModelSet models = QTraitPatternModelSet.create(traitScan,
			covariateModel, dfOffset, thresholdType, threshold,
			config.getUnestimableTraitPolicy());
		Path audit = Path.of(output.toString() + ".trait-patterns.tsv");
		models.writeAudit(audit, dfOffset, QTraitPatternScheduler.GENOTYPE.optionName());
		System.out.println("Trait-pattern audit: " + audit);
		if (models.excludedTraitRows() > 0
			&& config.getUnestimableTraitPolicy() == QUnestimableTraitPolicy.ERROR)
			throw new IllegalArgumentException("Unestimable trait missingness patterns are present; first: "
				+ models.firstExclusion() + ". The complete audit was written to " + audit
				+ ". Use --unestimable-trait-patterns skip to exclude and audit them.");
		if (models.excludedTraitRows() > 0)
			System.err.println("WARNING: Excluding " + models.excludedTraitRows()
				+ " unestimable trait row(s); see " + audit);
		if (models.estimableTraitRows() == 0)
			throw new IllegalArgumentException("No estimable trait rows remain after pattern/rank checks");

		try (QPatternPreparedTraitSource preparedSource = new QPatternPreparedTraitSource(
			rawTraitSource, alignment.expressionColumnOrder(), models);
			 QBinaryMatrixCache traitDiskCache = QBinaryMatrixCache.openOrBuildPrepared(
				cacheRoot.resolve("trait-pattern-global"), "TraitPatternsGlobal",
				preparedSource.signature(), rawTraitSource.metadata().path(), preparedSource,
				traitRowsPerBlock, config.getRebuildCache())) {
			QPreparedMatrix preparedTraits = selectPreparedTraitMatrix(traitDiskCache,
				traitRowsPerBlock, alignment.sampleCount(), predictorRowsPerBlock,
				config.getTraitCacheMode());
			try {
				 runGenotypeOuterBlocks(predictorSource, predictorColumns, preparedTraits,
					models, output, predictorRowsPerBlock, traitRowsPerBlock, dfOffset,
					predictorMissing, config.getMinimumMaf(), config.getMinimumMac(),
					frequencyScope);
			} finally {
				closeMemoryPreparedMatrix(preparedTraits);
			}
		}
	}

	private static void runGenotypeOuterBlocks(QMatrixRowSource predictorSource,
		int[] predictorColumns, QPreparedMatrix traits, QTraitPatternModelSet models,
		Path output, int predictorRowsPerBlock, int traitRowsPerBlock, int dfOffset,
		QMissingValuePolicy predictorMissing, double minimumMaf, double minimumMac,
		QVariantMatrixSource.FrequencyScope frequencyScope)
		throws Exception {
		int totalBlocks = (int) ((predictorSource.metadata().rowCount()
			+ predictorRowsPerBlock - 1) / predictorRowsPerBlock);
		Path checkpointDirectory = config.getCheckpointDirectory() == null
			? Path.of(output.toString() + ".checkpoint")
			: Path.of(config.getCheckpointDirectory()).toAbsolutePath().normalize();
		String signature = genotypeOuterCheckpointSignature(predictorSource,
			predictorColumns, traits, models, predictorRowsPerBlock, traitRowsPerBlock, dfOffset,
			predictorMissing, frequencyScope, minimumMaf, minimumMac);
		QAnalysisCheckpoint checkpoint = QAnalysisCheckpoint.open(checkpointDirectory,
			signature, totalBlocks, config.getResume(), config.getKeepCheckpoints());
		Path patternQcDirectory = Path.of(checkpointDirectory.toString() + ".pattern-qc");
		String patternQcSignature = signature + ":pattern-qc-v2";
		QGenotypeOuterPatternQcCheckpoint patternQc =
			QGenotypeOuterPatternQcCheckpoint.open(patternQcDirectory, patternQcSignature,
				totalBlocks, models.models().length, config.getResume(),
				config.getKeepCheckpoints());
		long completedComparisons = 0;
		for (int block = 0; block < totalBlocks; block++)
			if (checkpoint.isComplete(block) && patternQc.isComplete(block)) {
				completedComparisons = Math.addExact(completedComparisons,
					patternQc.activeComparisons(block, models));
			}
		QAnalysisProgress progress = QAnalysisProgress.dynamic(
			"Association genotype-outer trait patterns", completedComparisons);
		GpuContextPool contextPool = new GpuContextPool(mContexts);
		try {
			String header = kernelHeader(16, GpuPrecision.FP64) + "#define N_MIN_1 1" + sLn;
			long compileStarted = profiler.start();
			for (GpuContext context : contextPool.getAllContexts()) {
				context.setProfilingEnabled(profiler.isEnabled());
				context.compileKernel(header + eqtlReal, "eqtlReal", GpuPrecision.FP64);
			}
			profiler.record(QeQTLProfiler.Phase.KERNEL_COMPILE, compileStarted,
				contextPool.getAllContexts().size(), 0);
			int numThreads = min(config.getNumThreads(),
				min(Runtime.getRuntime().availableProcessors() + 1,
					contextPool.getAllContexts().size()));
			System.out.println("Num threads = " + numThreads);
			threadPool = Executors.newFixedThreadPool(numThreads);
			Deque<Future<?>> pending = new ArrayDeque<>();
			try (QMatrixRowSource.BlockReader reader = predictorSource.open(predictorColumns)) {
				for (int block = 0; block < totalBlocks; block++) {
					QMatrixRowSource.Block raw = reader.readBlock(predictorRowsPerBlock);
					if (raw == null || raw.rowOffset() != (long) block * predictorRowsPerBlock)
						throw new IOException("Predictor source ended or reordered during genotype-outer analysis");
					if (checkpoint.isComplete(block) && patternQc.isComplete(block)) continue;
					QGenotypeOuterPatternJob.QMatrixBlock submitted =
						new QGenotypeOuterPatternJob.QMatrixBlock(raw.rowOffset(), raw.rowIds(), raw.values());
					pending.addLast(threadPool.submit(new QGenotypeOuterPatternJob(submitted,
						traits, models, checkpoint, block, patternQc, contextPool,
						predictorRowsPerBlock, traitRowsPerBlock, predictorMissing, dfOffset,
						profiler, progress, minimumMaf, minimumMac,
						frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN)));
					if (pending.size() >= numThreads) pending.removeFirst().get();
				}
				if (reader.readBlock(1) != null)
					throw new IOException("Predictor source grew during genotype-outer analysis");
			}
			while (!pending.isEmpty()) pending.removeFirst().get();
			threadPool.shutdown();
			if (!threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS))
				throw new RuntimeException("Timed out waiting for genotype-outer workers");
			progress.complete();
			maybeFailBeforeGenotypeOuterAssembly();
			String outputHeader = rsqOnly ? "Rs_ID,ProbesetID,RSq,Dir,N,DF"
				: "Rs_ID,ProbesetID,RSq,Fx,T,log10P,N,DF";
			patternQc.assemble(Path.of(output.toString() + ".pattern-variant-qc.tsv"),
				models, minimumMaf, minimumMac,
				frequencyScope.name().toLowerCase(Locale.ROOT),
				frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN);
			checkpoint.assemble(output, outputHeader);
			patternQc.finishSuccess();
		} catch (Exception e) {
			if (threadPool != null) {
				threadPool.shutdownNow();
				try {
					if (!threadPool.awaitTermination(1, TimeUnit.MINUTES))
						e.addSuppressed(new IOException(
							"Genotype-outer workers did not stop within one minute"));
				} catch (InterruptedException interrupted) {
					Thread.currentThread().interrupt();
					e.addSuppressed(interrupted);
				}
			}
			throw e;
		} finally {
			progress.close();
			contextPool.close();
		}
	}

	private static void maybeFailBeforeGenotypeOuterAssembly() throws IOException {
		if (Boolean.parseBoolean(System.getProperty(
			"eqtl.test.genotype.outer.fail.before.assembly", "false")))
			throw new IOException("Injected test failure before genotype-outer assembly");
	}

	private static String genotypeOuterCheckpointSignature(QMatrixRowSource predictorSource,
		int[] predictorColumns, QPreparedMatrix traits, QTraitPatternModelSet models,
		int predictorRows, int traitRows, int dfOffset, QMissingValuePolicy predictorMissing,
		QVariantMatrixSource.FrequencyScope frequencyScope, double minimumMaf,
		double minimumMac)
		throws IOException {
		MessageDigest digest;
		try { digest = MessageDigest.getInstance("SHA-256"); }
		catch (NoSuchAlgorithmException e) { throw new IllegalStateException(e); }
		patternSignatureUpdate(digest, "gpu-eqtl-genotype-outer-pattern-v2");
		patternSignatureUpdate(digest, QBinaryMatrixCache.signature("PatternGenotypeRoot",
			predictorSource, predictorColumns, null));
		patternSignatureUpdate(digest, traits.signature());
		patternSignatureUpdate(digest, models.signature());
		patternSignatureUpdate(digest, predictorRows);
		patternSignatureUpdate(digest, traitRows);
		patternSignatureUpdate(digest, dfOffset);
		patternSignatureUpdate(digest, predictorMissing.optionName());
		patternSignatureUpdate(digest, frequencyScope.name());
		patternSignatureUpdate(digest, Double.doubleToLongBits(minimumMaf));
		patternSignatureUpdate(digest, Double.doubleToLongBits(minimumMac));
		patternSignatureUpdate(digest, simplifyResult ? 1 : 0);
		patternSignatureUpdate(digest, rsqOnly ? 1 : 0);
		patternSignatureUpdate(digest, GpuPrecision.FP64.optionName());
		return HexFormat.of().formatHex(digest.digest());
	}

	private static void validateTraitPatternWorkload(QMissingnessScan traitScan,
		long predictorRows, double[][] covariateModel, int dfOffset, int maximumPatterns) {
		int patternCount = traitScan.patterns().size();
		long observedAcrossPatterns = 0;
		int minimumObserved = traitScan.sampleCount();
		int maximumObserved = 0;
		int unestimablePatterns = 0;
		long unestimableTraitRows = 0;
		int designColumns = covariateModel == null ? 1 : covariateModel[0].length;
		for (QMissingnessScan.Pattern pattern : traitScan.patterns()) {
			int observed = traitScan.sampleCount() - pattern.missingSamples().cardinality();
			observedAcrossPatterns = Math.addExact(observedAcrossPatterns, observed);
			minimumObserved = Math.min(minimumObserved, observed);
			maximumObserved = Math.max(maximumObserved, observed);
			if (observed - designColumns - 1 - dfOffset <= 0) {
				unestimablePatterns++;
				unestimableTraitRows = Math.addExact(unestimableTraitRows,
					pattern.rowIndices().length);
			}
		}
		double fullCohortPasses = observedAcrossPatterns / (double) traitScan.sampleCount();
		double preparedBytes = predictorRows * (double) observedAcrossPatterns * Double.BYTES;
		System.out.println("Trait-pattern preflight: patterns=" + patternCount
			+ ", observed N range=" + minimumObserved + "-" + maximumObserved
			+ ", estimated predictor-preparation work="
			+ String.format(Locale.ROOT, "%.1f full-cohort pass(es)", fullCohortPasses)
			+ ", upper-bound prepared numeric data=" + formatBinaryBytes(preparedBytes));
		if (unestimablePatterns > 0)
			System.err.println("WARNING: " + unestimablePatterns + " trait pattern(s) containing "
				+ unestimableTraitRows + " trait row(s) have too few complete samples for "
				+ designColumns + " covariate-model column(s), one tested predictor, and df_offset="
				+ dfOffset + ".");
		boolean tooManyPatterns = maximumPatterns > 0 && patternCount > maximumPatterns;
		if (tooManyPatterns || unestimablePatterns > 0) {
			StringBuilder reason = new StringBuilder("Trait-pattern preflight failed: ");
			if (tooManyPatterns)
				reason.append(patternCount).append(" patterns exceed --max-trait-patterns ")
					.append(maximumPatterns).append(". ");
			if (unestimablePatterns > 0)
				reason.append(unestimablePatterns).append(" patterns are not estimable with the selected covariates. ");
			reason.append("Review the missingness QC report. Use trait mean/zero imputation, prefilter "
				+ "unestimable traits, or explicitly raise --max-trait-patterns (0 disables only the "
				+ "pattern-count limit) after reviewing the reported work estimate. Restartability "
				+ "preserves completed groups but does not reduce this workload.");
			throw new IllegalArgumentException(reason.toString());
		}
	}

	private static String formatBinaryBytes(double bytes) {
		if (!Double.isFinite(bytes))
			return "more than addressable storage";
		String[] units = {"B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB"};
		int unit = 0;
		while (bytes >= 1024 && unit < units.length - 1) {
			bytes /= 1024;
			unit++;
		}
		return String.format(Locale.ROOT, "%.2f %s", bytes, units[unit]);
	}

	private static void runTraitPatternAnalysis(QeQTLAnalysis plugin,
		QMatrixRowSource predictorSource, QMatrixRowSource rawTraitSource,
		QMissingnessScan traitScan,
		QSampleAlignment alignment, double[][] covariateModel,
		String thresholdType, double threshold, int dfOffset,
		int predictorRowsPerBlock, int traitRowsPerBlock, boolean hasCovariates,
		boolean variantPredictor, boolean predictorColumnsAreAligned,
		QMissingValuePolicy predictorMissing,
		QVariantMatrixSource.FrequencyScope frequencyScope,
		double minimumMaf, double minimumMac) throws Exception
	{
		Path output = Path.of(config.getOutputFilename()).toAbsolutePath().normalize();
		Path parent = output.getParent();
		if (parent == null)
			parent = Path.of(".").toAbsolutePath().normalize();
		Files.createDirectories(parent);
		Path cacheRoot = config.getCacheDirectory() == null
			? parent.resolve(".gpu-eqtl-cache").toAbsolutePath().normalize()
			: Path.of(config.getCacheDirectory()).toAbsolutePath().normalize();
		Path preparedCacheRoot = cacheRoot.resolve("trait-patterns");
		Path statisticsCacheRoot = cacheRoot.resolve("pattern-variant-statistics");
		Files.createDirectories(preparedCacheRoot);
		Path patternQcOutput = variantPredictor
			? Path.of(output.toString() + ".pattern-variant-qc.tsv") : null;
		Path checkpointDirectory = config.getCheckpointDirectory() == null
			? Path.of(output.toString() + ".checkpoint")
			: Path.of(config.getCheckpointDirectory()).toAbsolutePath().normalize();
		try {
			String checkpointSignature = traitPatternCheckpointSignature(predictorSource,
				rawTraitSource, traitScan, alignment, covariateModel, thresholdType, threshold,
				dfOffset, predictorRowsPerBlock, traitRowsPerBlock, hasCovariates,
				variantPredictor, predictorColumnsAreAligned, predictorMissing,
				frequencyScope, minimumMaf, minimumMac);
			QTraitPatternCheckpoint checkpoint = QTraitPatternCheckpoint.open(checkpointDirectory,
				checkpointSignature, traitScan.patterns().size(), config.getResume(),
				config.getKeepCheckpoints());
			int completedAtStart = checkpoint.completedResultCount();
			if (completedAtStart > 0)
				System.out.println("Completed trait-pattern groups found: " + completedAtStart
					+ " / " + traitScan.patterns().size());
			String outputHeader = rsqOnly ? "Rs_ID,ProbesetID,RSq,Dir,N,DF"
				: "Rs_ID,ProbesetID,RSq,Fx,T,log10P,N,DF";
			String patternQcHeader = "pattern_id\ttrait_rows\tobserved_samples\tinput_variants"
				+ "\tincluded_variants\tmonomorphic_variants\tbelow_min_maf\tbelow_min_mac"
				+ "\tno_call_variants\tmissing_genotypes\tfrequency_scope\tstatistics_cache\treused";
			for (int patternIndex = 0; patternIndex < traitScan.patterns().size(); patternIndex++) {
				QMissingnessScan.Pattern pattern = traitScan.patterns().get(patternIndex);
				int patternNumber = patternIndex + 1;
				if (checkpoint.isResultComplete(patternIndex)
					&& (!variantPredictor || checkpoint.isQcComplete(patternIndex))) {
					System.out.println("Reusing completed trait pattern " + patternNumber + "/"
						+ traitScan.patterns().size() + " (pattern_id=" + pattern.id() + ")");
					continue;
				}
				BitSet missing = pattern.missingSamples();
				int[] observed = complement(missing, alignment.sampleCount());
				int[] predictorColumns = predictorColumnsAreAligned
					? observed.clone() : select(alignment.genotypeColumnOrder(), observed);
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
				QMatrixRowSource effectivePredictorSource = predictorSource;
				QPatternVariantSource patternVariantSource = null;
				if (variantPredictor) {
					QMissingValuePolicy patternPolicy = predictorMissing == QMissingValuePolicy.MEAN
						|| predictorMissing == QMissingValuePolicy.ZERO ? predictorMissing
						: QMissingValuePolicy.ERROR;
					patternVariantSource = QPatternVariantSource.openOrBuild(statisticsCacheRoot,
						suffix, predictorSource, predictorColumns, patternPolicy,
						frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN,
						minimumMaf, minimumMac, config.getRebuildCache());
					effectivePredictorSource = patternVariantSource;
					QPatternVariantSource.Summary stats = patternVariantSource.summary();
					String patternQcLine = pattern.id() + "\t" + pattern.rowIndices().length + "\t"
						+ observed.length + "\t" + stats.inputVariants() + "\t"
						+ stats.includedVariants() + "\t" + stats.monomorphicVariants() + "\t"
						+ stats.belowMinimumMaf() + "\t" + stats.belowMinimumMac() + "\t"
						+ stats.noCallVariants() + "\t" + stats.missingGenotypes() + "\t"
						+ frequencyScope.name().toLowerCase(Locale.ROOT) + "\t"
						+ stats.cachePath() + "\t" + stats.reused();
					if (!checkpoint.isQcComplete(patternIndex))
						checkpoint.commitQcLine(patternIndex, patternQcLine);
					System.out.println("Pattern-specific variants: input=" + stats.inputVariants()
						+ ", included=" + stats.includedVariants()
						+ ", monomorphic=" + stats.monomorphicVariants()
						+ (frequencyScope == QVariantMatrixSource.FrequencyScope.PATTERN
							? ", below MAF=" + stats.belowMinimumMaf()
								+ ", below MAC=" + stats.belowMinimumMac() : ""));
					if (stats.includedVariants() == 0) {
						System.err.println("WARNING: No variable variants remain for trait pattern "
							+ pattern.id() + "; its traits are skipped.");
						if (!checkpoint.isResultComplete(patternIndex))
							checkpoint.commitEmptyResult(patternIndex);
						maybeFailAfterTraitPattern(patternNumber);
						continue;
					}
				}
				if (checkpoint.isResultComplete(patternIndex))
					continue;
				Path groupOutput = checkpoint.groupOutput(patternIndex);
				if (Files.isRegularFile(groupOutput)) {
					System.out.println("Recovering assembled trait-pattern group " + patternNumber
						+ " from " + groupOutput);
					checkpoint.commitResultFromOutput(patternIndex);
					maybeFailAfterTraitPattern(patternNumber);
					continue;
				}

				QGpuResidualizer residualizer = null;
				if (patternQ != null && residualizationMode != QResidualizationMode.CPU)
					residualizer = new QGpuResidualizer(mContexts, patternQ, gpuPrecision, profiler);
				try {
					String preprocessingTag = residualizer == null ? null : residualizer.cacheSignatureTag();
					String predictorKind = "Predictor" + suffix;
					String traitKind = "Trait" + suffix;
					String predictorSignature = QBinaryMatrixCache.signature(predictorKind,
						effectivePredictorSource, predictorColumns, patternQ, preprocessingTag);
					String traitSignature = QBinaryMatrixCache.signature(traitKind,
						traitSource, traitColumns, patternQ, preprocessingTag);
					try (QBinaryMatrixCache predictorCache = QBinaryMatrixCache.openOrBuild(preparedCacheRoot,
						predictorKind, predictorSignature, effectivePredictorSource, predictorColumns, patternQ,
						predictorRowsPerBlock, config.getRebuildCache(), residualizer);
						 QBinaryMatrixCache traitCache = QBinaryMatrixCache.openOrBuild(preparedCacheRoot,
						traitKind, traitSignature, traitSource, traitColumns, patternQ,
						traitRowsPerBlock, config.getRebuildCache(), residualizer)) {
						if (residualizer != null) {
							residualizer.close();
							residualizer = null;
						}
						Path blockCheckpoint = checkpoint.blockCheckpointDirectory(patternIndex);
						QPreparedMatrix preparedTraits = selectPreparedTraitMatrix(traitCache,
							traitRowsPerBlock, observed.length, predictorRowsPerBlock,
							config.getTraitCacheMode());
						try {
							plugin.eSNPAnalysisStreamed(predictorCache, preparedTraits, patternThreshold,
								dfOffset, errorDegreesOfFreedom, predictorRowsPerBlock, traitRowsPerBlock,
								groupOutput, blockCheckpoint, Files.exists(blockCheckpoint),
								config.getKeepCheckpoints(), false, true,
								"Association trait pattern " + pattern.id());
						} finally {
							closeMemoryPreparedMatrix(preparedTraits);
						}
						checkpoint.commitResultFromOutput(patternIndex);
						maybeFailAfterTraitPattern(patternNumber);
					}
				} finally {
					if (residualizer != null)
						residualizer.close();
				}
			}
			checkpoint.assembleResults(output, outputHeader);
			if (patternQcOutput != null) {
				checkpoint.assembleQc(patternQcOutput, patternQcHeader);
				System.out.println("Pattern-specific variant QC summary: " + patternQcOutput);
			}
			checkpoint.finishSuccess();
		} finally {
			for (GpuContext context : mContexts)
				if (context != null)
					context.close();
		}
	}

	private static String traitPatternCheckpointSignature(QMatrixRowSource predictorSource,
		QMatrixRowSource traitSource, QMissingnessScan traitScan, QSampleAlignment alignment,
		double[][] covariateModel, String thresholdType, double threshold, int dfOffset,
		int predictorRowsPerBlock, int traitRowsPerBlock, boolean hasCovariates,
		boolean variantPredictor, boolean predictorColumnsAreAligned,
		QMissingValuePolicy predictorMissing, QVariantMatrixSource.FrequencyScope frequencyScope,
		double minimumMaf, double minimumMac) throws IOException {
		MessageDigest digest;
		try {
			digest = MessageDigest.getInstance("SHA-256");
		} catch (NoSuchAlgorithmException e) {
			throw new IllegalStateException("SHA-256 is unavailable", e);
		}
		patternSignatureUpdate(digest, "gpu-eqtl-trait-pattern-checkpoint-v1");
		int[] predictorColumns = predictorColumnsAreAligned ? identity(alignment.sampleCount())
			: alignment.genotypeColumnOrder();
		patternSignatureUpdate(digest, QBinaryMatrixCache.signature("TraitPatternPredictorRoot",
			predictorSource, predictorColumns, null));
		patternSignatureUpdate(digest, QBinaryMatrixCache.signature("TraitPatternTraitRoot",
			traitSource, alignment.expressionColumnOrder(), null));
		for (int value : alignment.genotypeColumnOrder()) patternSignatureUpdate(digest, value);
		for (int value : alignment.expressionColumnOrder()) patternSignatureUpdate(digest, value);
		if (covariateModel == null) {
			patternSignatureUpdate(digest, -1);
		} else {
			patternSignatureUpdate(digest, covariateModel.length);
			patternSignatureUpdate(digest, covariateModel[0].length);
			for (double[] row : covariateModel)
				for (double value : row)
					patternSignatureUpdate(digest, Double.doubleToLongBits(value));
		}
		patternSignatureUpdate(digest, thresholdType);
		patternSignatureUpdate(digest, Double.doubleToLongBits(threshold));
		patternSignatureUpdate(digest, dfOffset);
		patternSignatureUpdate(digest, predictorRowsPerBlock);
		patternSignatureUpdate(digest, traitRowsPerBlock);
		patternSignatureUpdate(digest, hasCovariates ? 1 : 0);
		patternSignatureUpdate(digest, variantPredictor ? 1 : 0);
		patternSignatureUpdate(digest, predictorColumnsAreAligned ? 1 : 0);
		patternSignatureUpdate(digest, predictorMissing.optionName());
		patternSignatureUpdate(digest, frequencyScope.name());
		patternSignatureUpdate(digest, Double.doubleToLongBits(minimumMaf));
		patternSignatureUpdate(digest, Double.doubleToLongBits(minimumMac));
		patternSignatureUpdate(digest, simplifyResult ? 1 : 0);
		patternSignatureUpdate(digest, rsqOnly ? 1 : 0);
		patternSignatureUpdate(digest, gpuPrecision.optionName());
		patternSignatureUpdate(digest, residualizationMode.name());
		patternSignatureUpdate(digest, traitScan.patterns().size());
		for (QMissingnessScan.Pattern pattern : traitScan.patterns()) {
			patternSignatureUpdate(digest, pattern.id());
			patternSignatureUpdate(digest, pattern.rowIndices().length);
			for (long row : pattern.rowIndices()) patternSignatureUpdate(digest, row);
			long[] missingWords = pattern.missingSamples().toLongArray();
			patternSignatureUpdate(digest, missingWords.length);
			for (long word : missingWords)
				patternSignatureUpdate(digest, word);
		}
		return HexFormat.of().formatHex(digest.digest());
	}

	private static void patternSignatureUpdate(MessageDigest digest, String value) {
		digest.update(value.getBytes(StandardCharsets.UTF_8));
		digest.update((byte) 0);
	}

	private static void patternSignatureUpdate(MessageDigest digest, long value) {
		patternSignatureUpdate(digest, Long.toString(value));
	}

	private static void maybeFailAfterTraitPattern(int completedPatternNumber) throws IOException {
		String configured = System.getProperty(TEST_FAIL_AFTER_TRAIT_PATTERN_PROPERTY);
		if (configured == null || configured.isBlank())
			return;
		int failAfter;
		try {
			failAfter = Integer.parseInt(configured.trim());
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException(TEST_FAIL_AFTER_TRAIT_PATTERN_PROPERTY
				+ " must be a positive integer", e);
		}
		if (failAfter <= 0)
			throw new IllegalArgumentException(TEST_FAIL_AFTER_TRAIT_PATTERN_PROPERTY
				+ " must be a positive integer");
		if (completedPatternNumber >= failAfter)
			throw new IOException("Injected test failure after trait pattern "
				+ completedPatternNumber);
	}

	static double rSquaredThreshold(String type, double threshold,
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

	private static QPreparedMatrix selectPreparedTraitMatrix(QBinaryMatrixCache diskCache,
		int traitRowsPerBlock, int sampleCount, int predictorRowsPerBlock,
		QTraitCacheMode mode) throws IOException {
		if (mode == QTraitCacheMode.DISK) {
			System.out.println("Prepared trait cache residency = disk (residualization/standardization was performed once)");
			return diskCache;
		}
		long estimate = QInMemoryPreparedMatrix.estimateResidentBytes(
			diskCache.rowCount(), diskCache.sampleCount());
		Runtime runtime = Runtime.getRuntime();
		long used = runtime.totalMemory() - runtime.freeMemory();
		long available = Math.max(0, runtime.maxMemory() - used);
		long workerBytes = GpuTuning.estimateStreamedWorkerBytes(sampleCount,
			predictorRowsPerBlock, traitRowsPerBlock, gpuPrecision);
		long workerReserve = saturatedMultiply(workerBytes, Math.max(1, config.getNumThreads()));
		long safetyReserve = Math.max(256L * 1024 * 1024, runtime.maxMemory() / 10);
		long required = saturatedAdd(estimate, saturatedAdd(workerReserve, safetyReserve));
		boolean safe = required <= available;
		if (!safe && mode == QTraitCacheMode.AUTO) {
			System.out.println("Prepared trait cache residency = disk (memory estimate "
				+ humanBytes(estimate) + "; available JVM heap " + humanBytes(available)
				+ " after reserving workers/headroom)");
			return diskCache;
		}
		if (!safe)
			throw new IllegalArgumentException("--trait-cache memory needs about "
				+ humanBytes(required) + " of currently available JVM heap including worker/headroom reserves, but "
				+ humanBytes(available) + " is available; increase -Xmx, reduce block rows/threads, or use --trait-cache disk");
		System.out.println("Prepared trait cache residency = memory (estimated " + humanBytes(estimate)
			+ "; residualization/standardization was performed once)");
		long started = profiler.start();
		QInMemoryPreparedMatrix memory = QInMemoryPreparedMatrix.load(diskCache, traitRowsPerBlock);
		profiler.record(QeQTLProfiler.Phase.TRAIT_CACHE_MEMORY_LOAD, started,
			diskCache.rowCount(), saturatedMultiply(saturatedMultiply(
				diskCache.sampleCount(), diskCache.rowCount()), Double.BYTES));
		return memory;
	}

	private static void closeMemoryPreparedMatrix(QPreparedMatrix prepared) {
		if (prepared instanceof QInMemoryPreparedMatrix memory)
			memory.close();
	}

	private static long saturatedMultiply(long left, long right) {
		try {
			return Math.multiplyExact(left, right);
		} catch (ArithmeticException e) {
			return Long.MAX_VALUE;
		}
	}

	private static long saturatedAdd(long left, long right) {
		try {
			return Math.addExact(left, right);
		} catch (ArithmeticException e) {
			return Long.MAX_VALUE;
		}
	}

	private static String humanBytes(long bytes) {
		if (bytes == Long.MAX_VALUE)
			return "more than the addressable heap";
		double gib = bytes / (1024.0 * 1024 * 1024);
		return String.format(java.util.Locale.ROOT, "%.2f GiB", gib);
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

	private static Path matrixCacheDirectory(String predictorFilename)
	{
		if (config.getCacheDirectory() != null)
			return Path.of(config.getCacheDirectory()).toAbsolutePath().normalize();
		String anchor = config.getOutputFilename() == null
			? predictorFilename : config.getOutputFilename();
		Path parent = Path.of(anchor).toAbsolutePath().normalize().getParent();
		if (parent == null)
			parent = Path.of(".").toAbsolutePath().normalize();
		return parent.resolve(".gpu-eqtl-cache").toAbsolutePath().normalize();
	}

	private static int[] identity(int count)
	{
		int[] result = new int[count];
		for (int i = 0; i < count; i++)
			result[i] = i;
		return result;
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
		if (config.getInspectCache() || config.getPruneCache() || config.getApplyCachePrune()) {
			try {
				QCacheManager.run(new QCacheManager.Options(
					config.getCacheDirectory() == null ? null : Path.of(config.getCacheDirectory()),
					config.getCacheReportFilename() == null ? null : Path.of(config.getCacheReportFilename()),
					config.getPruneCache(), config.getApplyCachePrune(),
					config.getCachePruneOlderThanDays()));
				return;
			} catch (Exception error) {
				System.err.println("ERROR: " + error.getMessage());
				System.exit(kExitCodeErrorInvalidParam);
				return;
			}
		}
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
		System.out.println("Matrix-product precision = " + gpuPrecision.optionName());
		System.out.println("Fixed-effect residualization mode = " + residualizationMode.optionName());
		if (gpuPrecision == GpuPrecision.FP32)
			System.err.println("NOTE: FP32 is enabled explicitly. Backend matrix products and accelerated covariate "
				+ "projection are approximate; explicit CPU projection and statistical calculations remain FP64.");
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
		if (config.getValidateOnly() || config.getInspectMissingness() || config.getPreprocessOnly()) {
			numDevices = 1;
			System.out.println("Validation/preprocessing mode: compute-backend initialization is skipped.");
		} else {
			System.out.println("Initializing compute backend...");
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
			boolean cohortAware = config.getCohortModelFilename() != null;
			if (covarFixed == null && covarRandom == null && !cohortAware) {
				System.err.println("ERROR: Fixed or random covariates are not mentioned, but the covariate file is specified.");
				System.exit(kExitCodeErrorWrongCovarSpec);
			} else if (covarRandom != null) {
				System.err.println("WARNING: Random covariates are currently ignored.");
			}
			if (covarFixed == null && !cohortAware) {
				System.err.println("ERROR: At present, you will need to specify at least one fixed covariate with a specified covariate file.");
				System.exit(kExitCodeErrorWrongCovarSpec);
			}
			System.out.println("Fixed covariates: " + (covarFixed == null ? "<per cohort>" : QStringUtils.toString(covarFixed)));
		} else {
			if (covarFixed != null || covarRandom != null) {
				System.err.println("ERROR: Fixed or random covariates are mentioned without specifying the covariate file.");
				System.exit(kExitCodeErrorWrongCovarSpec);
			}
		}
		if (pedigreeFilename != null)
			System.out.println("Pedigree file: " + pedigreeFilename);
		System.out.println("Output file: " + (outputFilename == null ? "<not requested>" : outputFilename));
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
		if (config.getPreprocessOnly()) {
			System.err.println("ERROR: --preprocess-only requires VCF or BCF genotype input and CSV expression input.");
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
