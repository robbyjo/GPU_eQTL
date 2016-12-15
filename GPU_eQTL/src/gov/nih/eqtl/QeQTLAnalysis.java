/*
 * Roby Joehanes
 * 
 * Copyright 2009 Roby Joehanes
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
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.bridj.BridJ;

import com.nativelibs4java.opencl.CLDevice;

import jdistlib.T;
import gov.nih.eqtl.datastructure.QGeneticSNPData;
import gov.nih.eqtl.datastructure.QSNPData;
import gov.nih.eqtl.datastructure.QSNPDataInt;
import gov.nih.eqtl.io.QPlinkLoader;
import gov.nih.eqtl.datastructure.QGeneExpressionData;
import gov.nih.jama.QRDecomposition;
import gov.nih.opencl.QOCLContext;
import gov.nih.opencl.QOCLContextPool;
import gov.nih.opencl.QOpenCL;
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
		kGB = 1024 * kMB,
		kLambda = 0.75;

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
		kExitCodeErrorGenoMissingValues = -12;

	// Global variables
	static QOCLContext[] mContexts = null;
	static long mMinAllocMemSize = 0;
	static QeQTLAnalysisConfig config = null;

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
		int dfo, int dfe, double lambda, boolean isAdditive)
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
		if (covarQ != null) {
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
		double[] expDataSD = new double[numETraits], snpDataSD = new double[numSNPs];
		for (int i = 0; i < numETraits; i++)
			expDataSD[i] = calcStdDevAndStandardize(expData[i]);
		for (int i = 0; i < numSNPs; i++)
			snpDataSD[i] = calcStdDevAndStandardize(snpData[i]);
		time2 = System.currentTimeMillis();
		System.out.println("Standardizing time = " + (time2 - time1));

		int colPereQTL = isCategoricalSNP ? (isAdditive ? 3 : 4) : 1;

		time1 = System.currentTimeMillis();
		final int
			localBlockSize = 16,
			sizeof_data = 8;
		String
			program,
			kernelName,
			hdr = "#define BLOCK_SIZE " + localBlockSize + sLn + "#define DATATYPE double" + sLn;

		if (isCategoricalSNP) {
			program = eqtlCat;
			kernelName = "eqtlCat";
			hdr += "#define COL_PER_EQTL " + colPereQTL + sLn;
		} else {
			program = eqtlReal;
			kernelName = "eqtlReal";
			hdr += "#define N_MIN_1 " + (numInds - 1) + sLn;
		}

		QOCLContextPool gpuContextsPool = new QOCLContextPool(mContexts);
		List<QOCLContext> gpuContexts = gpuContextsPool.getAllContexts();
		int numDevices = gpuContexts.size();
		for (int i = 0; i < numDevices; i++)
		{
			try {
				QOCLContext ctx = gpuContexts.get(i);
				String
					dblExt = ctx.getDoubleExtensionHeader(),
					temp = hdr;
				if (dblExt != null) {
					temp += "#pragma OPENCL EXTENSION " + dblExt + //$NON-NLS-1$
						": enable" + sLn; //$NON-NLS-1$
				} else {
					// Should never happen
					throw new RuntimeException("ERROR: GPU not capable of double precision!");
				}
				ctx.setKernelFromProgram(temp + program, kernelName);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		time2 = System.currentTimeMillis();
		System.out.println("Kernel compile time = " + (time2 - time1));

		int globalBlockSize = config.getBlockSize();
		int numThreads = config.getNumThreads();
		System.out.println("Num threads = " + numThreads);
		int maxNumThreads = Runtime.getRuntime().availableProcessors() + 1;
		threadPool = Executors.newFixedThreadPool(numThreads > maxNumThreads ? maxNumThreads : numThreads);

		//time1 = System.currentTimeMillis();
		int
			numETraitsPerBlock = globalBlockSize,
			numESNPsPerBlock = globalBlockSize / colPereQTL,
			nrow = roundUpNearestMultiple(numInds, localBlockSize);
		System.out.println("Estimating buffer size");
		long
			reqdCacheSize = 4 * sizeof_data * numThreads * numETraitsPerBlock * nrow,
			availMem = QSystemUtils.getMemoryStatistics(false).mAvailable,
			mem = availMem - reqdCacheSize;

		if (mem <= 0) {
			System.err.println("WARNING: Not enough memory for storage. This program will fail at some point. Try reducing the values for num_threads and block_size.");
			mem = availMem;
		}

		//System.out.println("Allocating memory");
		//int numStoredETraitsPerESNP = 32; // (int) ceil((lambda * min(mem, availMem) * 0.5) / (numSNPs * (sizeof_data + 4)));
		//System.out.println("Number stored e-traits per eSNP = " + numStoredETraitsPerESNP);

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
			e.printStackTrace();
		} finally {
		}

		try {
			int numIters = (int) ceil(numSNPs * 1.0 / numESNPsPerBlock);
			QSynchronizedCounter counter = new QSynchronizedCounter(0, numIters);
			if (isCategoricalSNP) {
				//for (int i = 0; i < numThreads; i++)
				//threadPool.execute(new QeQTLSNPJobCat(popn, expData, covarQ, gpuContextsPool, results,
				//	numETraitsPerBlock, numESNPsPerBlock, localBlockSize, rsq0, isAdditive, counter);
			} else {
				for (int i = 0; i < numThreads; i++)
					threadPool.execute(new QeQTLSNPJobReal(popn, expDataTbl, expDataSD, snpDataSD, gpuContextsPool,
						numETraitsPerBlock, numESNPsPerBlock, localBlockSize, dfo, dfe, rsq0, isAdditive, fw, counter));
			}
			threadPool.shutdown();
			try {
				threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
			} catch (InterruptedException exx) {}
			fw.flush();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			gpuContextsPool.dispose();
			if (fw != null)
				try {
					fw.close();
				} catch (IOException e) {
				}
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
		if (EPlatform.isLinux())
		{
			String[] paths = config.getLibraryPaths();
			if (paths == null) {
				BridJ.addLibraryPath("/usr/lib64");
				BridJ.addLibraryPath("/usr/lib64/nvidia");
			} else {
				for (String path: paths)
					BridJ.addLibraryPath(path);
			}
			int numNVidiaDevices = QOpenCL.getNumNVidiaDevicesInLinux();
			if (numNVidiaDevices > 0) {
				if (QOpenCL.checkNVidiaDevicePermissionInLinux())
					System.out.println("Device permission appears to be okay.");
				else {
					System.out.println("WARNING: The device permission appears to be not okay. Most likely error will result soon.");
					System.out.println("All NVidia devices (/dev/nvidia*) must be globally read-writable.");
				}
			}
		}
		System.out.println("Detected OpenCL library: " + QOpenCL.getDetectedOpenCLLib());
		try {
			QOpenCL.init();
		} catch (Throwable e) {
			System.err.println("Cannot initialize OpenCL library. Please contact your system administrator.");
			System.exit(kExitCodeErrorInitOpenCLFailure);
		}
		CLDevice[] devices = QOpenCL.getAllGPUDevices(true, true);
		int numDevices = devices == null ? 0 : devices.length;
		if (numDevices == 0) {
			System.err.println("Cannot find any GPUs capable of 64-bit floating point arithmetic.");
			devices = QOpenCL.getAllGPUDevices(true, false);
			numDevices = devices == null ? 0 : devices.length;
			if (numDevices > 0) {
				System.err.println(numDevices + " GPU" + (numDevices > 1 ? "s": "") + " installed in this machine "
					+ (numDevices > 1 ? "are": "is") + " only capable of 32-bit floating point arithmetic:");
				for (CLDevice device: devices)
					System.err.println(device.getName());
				System.err.println("For NVidia cards, you will need cards with Compute Capability version 2.0 or higher:");
				System.err.println("https://developer.nvidia.com/cuda-gpus");
				System.err.println("For ATI/AMD cards, the availability of double-precision capabilities may VARY by operating system." +
				" For example Radeon HD 4850 is capable of double-precision in Windows, but not in Mac.");
				System.err.println("Contact your system administrators.");
			}
			System.exit(kExitCodeError64bitGPUNotFound);
		}
		System.out.println("Found " + numDevices + " suitable GPU" + (numDevices > 1 ? "s": "") + " in this machine:");
		boolean hasHostUnifiedMemory = false;
		for (CLDevice device: devices) {
			System.out.println(device.getName());
			hasHostUnifiedMemory = hasHostUnifiedMemory | device.isHostUnifiedMemory();
		}
		if (hasHostUnifiedMemory) {
			System.err.println("WARNING: Some of the GPUs have a unified memory access with host."+
				" This is often the case with the low-end GPUs. Such memory architecture will significantly reduce the performance of the computation" +
				" due to extensive memory bus transfer bottleneck.");
		}
		mContexts = QOpenCL.getAllGPUContexts(true, true);
		assert(numDevices == mContexts.length); // Sanity check. If not, something is seriously wrong
		// Detect whether all the detected GPUs are equivalent
		// Note: Only the device name is checked.
		String devName = null;
		boolean allDevicesTheSame = true;
		for (QOCLContext context: mContexts) {
			CLDevice dev = context.getDevice();
			if (devName == null) {
				devName = dev.getName();
			} else {
				if (!devName.equalsIgnoreCase(dev.getName())) {
					allDevicesTheSame = false;
					break;
				}
			}
			if (mMinAllocMemSize == 0)
				mMinAllocMemSize = dev.getMaxMemAllocSize();
			else
				mMinAllocMemSize = min(mMinAllocMemSize, dev.getMaxMemAllocSize());
		}
		if (!allDevicesTheSame) {
			System.err.println("WARNING: These GPUs may not be identical in their features and capabilities."+
				" Most likely, the performance will be constrained by the weakest GPU.");
		} else {
			System.out.println("All the detected GPUs appear to be identical.");
		}
		if (mMinAllocMemSize <= 0 && config.getBlockSize() <= 0) {
			System.err.println("WARNING: The GPU(s) do not report their RAM sizes correctly." +
				" This program will attempt to use a conservative setting that may or may not work." +
				" You may want to manually specify block_size entry in your input.");
			System.err.println("A good rule of thumb for specifying block_size is b = -p + 0.5*sqrt(p*p + c/2), " +
				"where p is the number of individuals rounded up to a multiple of 16 " +
				"and c is the amount of memory that can be allocated, and round b down to a multiple of 512.");
			System.err.println("At the moment, the detected c is " + mMinAllocMemSize);
		}
	}

	private static final void printUsage()
	{
		System.err.println("ERROR: No parameters specified!");
		System.err.println("Usage: (program_name) --printgpuinfo will print GPU information installed in this machine.");
		System.err.println("(program_name) (inifile) will run the eQTL analysis with parameters specified in the inifile.");
		System.exit(kExitCodeErrorInvalidParam);
	}

	private static final void dumpGPUInfo()
	{
		System.out.println(QOpenCL.dumpAllPlatformInformation(false, true));
		System.exit(kExitCodeNormal);
	}

	static final int getDefaultBlockSize(int numInds, int numSNPs, int numETraits)
	{
		int
			n = roundUpNearestMultiple(numInds, QOpenCL.kDefaultAlignment),
			n_snps = roundUpNearestMultiple(numSNPs, QOpenCL.kDefaultAlignment),
			n_etraits = roundUpNearestMultiple(numETraits, QOpenCL.kDefaultAlignment);
		n = roundDownNearestMultiple((int) (sqrt(n * n + mMinAllocMemSize / 2.0) - n), 512);
		return min(n, min(n_snps, n_etraits));
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

	public static void main(String[] args)
	{
		long timea = System.currentTimeMillis();
		System.out.println("GPU eQTL analysis software version 1.0. By: Roby Joehanes");
		if (args == null || args.length < 1)
			printUsage();

		List<String> params = new ArrayList<String>();
		for (int i = 0; i < args.length; i++)
			params.add(args[i]);
		if (params.contains("--printgpuinfo"))
			dumpGPUInfo();

		if (params.contains("--debug"))
			DEBUG = true;

		if (!EPlatform.is64Bit()) {
			System.err.println("ERROR: Your operating system / platform seems to be 32-bit, not 64-bit.");
			System.err.println("This program will not run properly in 32-bit platform due to severe limitation to the amount data that can be processed."
				+ " You may have a 64-bit machine, but you will need a 64-bit operating system to run in it. You also need 64-bit Java."
				+ " Please contact your system administrator.");
			System.exit(kExitCodeErrorPlatformNot64Bit);
		}

		try {
			config = QeQTLAnalysisConfig.loadConfig(args[0]);
		} catch (Exception e) {
			System.err.println("ERROR: Failed to load the input configuration file!");
			System.exit(kExitCodeErrorCantLoadFile);
		}
		System.out.println("Initializing GPUs...");
		initGPUs();
		int numDevices = mContexts.length;
		String
			famFilename = config.getFamilyFilename(),
			genoFilename = config.getGenotypeFilename(),
			exprFilename = config.getExpressionFilename(),
			covarFilename = config.getCovariateFilename(),
			pedigreeFilename = config.getPedigreeFilename(),
			outputFilename = config.getOutputFilename(),
			genoFileformat = config.getGenotypeFileFormat(),
			genoModel = config.getGenotypeModel(),
			thresholdType = config.getThresholdType(),
			covarFixed[] = config.getFixedCovariates(),
			covarRandom[] = config.getRandomCovariates(),
			covarFactor[] = config.getFactorCovariates();
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

			int numReqdIter = (int) (ceil(numSNPs * 1.0 / globalBlockSize) * ceil(numETraits * 1.0 / globalBlockSize));
			System.out.println("Given the block size, " + numReqdIter + " iteration(s) are needed");

			int
				recMaxNumThreads = min(QSystemUtils.kNumCPUCores, min(numReqdIter, numDevices * 2)),
				numThreads = config.getNumThreads();
			if (numThreads == 0) {
				numThreads = recMaxNumThreads;
				System.out.println("The num_threads parameter is not specified, trying " + numThreads);
				config.setNumThreads(recMaxNumThreads);
			}
			if (numThreads > QSystemUtils.kNumCPUCores) {
				System.err.println("WARNING: Specifying num_threads value that is greater than that is available will degrade performance."
				+ "The number of cores available in the system is " + QSystemUtils.kNumCPUCores);
			}

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
			plugin.eSNPAnalysis(snpData, exprData, covarMatrix, t0, dfOffset, dfe, kLambda, isAdditive);
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
		}
		long timeb = System.currentTimeMillis();
		System.out.println("Overall time (in seconds) = " + (timeb - timea) / 1000.0);
		System.exit(kExitCodeNormal);
	}
}
