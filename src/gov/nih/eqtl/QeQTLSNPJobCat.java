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

import java.util.Arrays;

import gov.nih.eqtl.datastructure.QGeneticSNPData;
import gov.nih.gpu.GpuContext;
import gov.nih.gpu.GpuContextPool;
import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.IJobOwner;
import gov.nih.parallel.QSynchronizedCounter;
import gov.nih.utils.QDataUtils;
import gov.nih.utils.matrix.QMatrixUtils;
//import gov.nih.utils.QStatsUtils;
import static java.lang.Math.min;
import static gov.nih.gpu.GpuRuntime.DEFAULT_ALIGNMENT;

public class QeQTLSNPJobCat implements IGenericParallelJob, Runnable {
	protected QGeneticSNPData popn;
	protected GpuContextPool contextPool;
	protected double expData[][], covarQ[][], flatETraits[], flatSNPs[], RSq0;
	protected int numETraitsPerBlock, numETraits, numSNPsPerBlock, numSNPs, numInds, nrow, blockSize, snpCol;
	protected int[][] genotypeCount;
	protected int snpCounts[], eTraitIndices[][], colPereQTL, sumYOffset, sumYSqOffset;
	protected double snpResults[][][];
	protected boolean isAdditive;
	protected static final int fillIns[][] = { {0, 1, 0, -1}, {0, 0, 1, 0}, {0, 1, 1, 1} }; // Fill in matrix
	static final int sizeof_data = 8;
	protected QeQTLSNPAnalysisResults results;
	protected QSynchronizedCounter counter;

	/**
	 * Round up the number to the nearest multiple
	 * @param number
	 * @param multiple
	 * @return
	 */
	static final int roundUpNearestMultiple(int number, int multiple) {
		int rem = number % multiple;
		if (rem > 0)
			number = ((number / multiple) * multiple) + multiple;
		return number;
	}

	public QeQTLSNPJobCat(QGeneticSNPData popn, double[][] expData, double[][] covarQ, GpuContextPool contextPool,
		QeQTLSNPAnalysisResults results, int numETraitsPerBlock, int numSNPsPerBlock, int blockSize, double RSq0, boolean isAdditive,
		QSynchronizedCounter ct)
	{
		this.popn = popn;
		this.expData = expData;
		this.contextPool = contextPool;
		this.numETraitsPerBlock = numETraitsPerBlock;
		this.numSNPsPerBlock = numSNPsPerBlock;
		this.isAdditive = isAdditive;
		if (isAdditive) {
			colPereQTL = 3;
			sumYOffset = 1;
			sumYSqOffset = 2;
		} else {
			colPereQTL = 4;
			sumYOffset = 2;
			sumYSqOffset = 3;
		}
		numSNPs = popn.getNumSNPs();
		numInds = popn.getNumIndividuals();
		numETraits = expData.length;
		snpCol = numSNPsPerBlock * colPereQTL;
		nrow = roundUpNearestMultiple(numInds, DEFAULT_ALIGNMENT);
		this.blockSize = blockSize;
		flatETraits = new double[numETraitsPerBlock * nrow];
		flatSNPs = new double[snpCol * nrow];
		genotypeCount = new int[numSNPsPerBlock][4];
		this.RSq0 = RSq0;
		this.results = results;
		this.snpCounts = results.mSNPCounts;
		this.eTraitIndices = results.mETraitIndices;
		this.snpResults = results.mSNPResults;
	}

	@Override
	public void execute() {
		int iterNo = counter.next();
		if (iterNo >= 0)
		{
			do {
				execute(iterNo);
				iterNo = counter.next();
			} while (iterNo >= 0);
		}
	}

	public boolean execute(int eSNPBlockNo)
	{
		int
			curSNPOffset = eSNPBlockNo * numSNPsPerBlock,
			curNumSNPs = min(numSNPsPerBlock, numSNPs - curSNPOffset);
			//curSNPCol = 4 * curNumSNPs,
			//curSNPColAligned = roundUpNearestMultiple(curSNPCol, kDefaultAlignment);
		long nElems = numETraitsPerBlock * snpCol;
		if (curNumSNPs < numSNPsPerBlock)
			Arrays.fill(flatSNPs, 0);

		double[] xyResult;
		for (int k = 0; k < curNumSNPs; k++)
		{
			int[]
				snps = popn.getSNPs().get(curSNPOffset+k).getSNPCodes(),
				count = genotypeCount[k];
			count[0] = count[1] = count[2] = count[3] = 0;
			int multiple_k = colPereQTL * k;
			if (isAdditive) {
				for (int l = 0; l < numInds; l++)
				{
					int snpCode = snps[l] + 1, offset = l * snpCol + multiple_k;
					flatSNPs[offset] = fillIns[0][snpCode];
					flatSNPs[offset + 1] = flatSNPs[offset + 2] = fillIns[2][snpCode];
					// Include count genotype here
					count[snpCode]++;
				}
				count[2] = 0; // Set heterozygote count as zero for additive models
			} else {
				for (int l = 0; l < numInds; l++)
				{
					int snpCode = snps[l] + 1, offset = l * snpCol + multiple_k;
					flatSNPs[offset] = fillIns[0][snpCode];
					flatSNPs[offset + 1] = fillIns[1][snpCode];
					flatSNPs[offset + 2] = flatSNPs[offset + 3] = fillIns[2][snpCode];
					// Include count genotype here
					count[snpCode]++;
				}
			}
		}
		//double[] snpSDs = null;
		if (covarQ != null) {
			QMatrixUtils.calculateXMinusYYtXInPlace(flatSNPs, covarQ);
			//snpSDs = QStatsUtils.calcStdDevAndStandardizeByColumn(flatSNPs, numInds, curNumSNPs);
		}

		for (int curETraitOffset = 0; curETraitOffset < numETraits; curETraitOffset += numETraitsPerBlock)
		{
			int
				curNumETraits = min(numETraitsPerBlock, numETraits - curETraitOffset);
				//curNumETraitsAligned = roundUpNearestMultiple(curNumETraits, kDefaultAlignment);
			if (curNumETraits < numETraitsPerBlock)
				Arrays.fill(flatETraits, 0);
			//System.out.println("curSNPCol = " + curSNPCol + ", curNumETraits = " + curNumETraits + ", numInds = " + numInds + ", nrow = " + nrow + ", blockSize = " + blockSize);
			for (int k = 0; k < curNumETraits; k++)
				System.arraycopy(expData[curETraitOffset+k], 0, flatETraits, nrow * k, numInds);
			//eTraitOffset[deviceNo] = i;
			//snpOffset[deviceNo] = j;
			GpuContext gpuContext = contextPool.reserveContext();
			try {
				long localMemoryBytes = (long) (blockSize + 1) * (4L * blockSize) * sizeof_data;
				int activeSnpColumns = roundUpNearestMultiple(curNumSNPs * colPereQTL, blockSize);
				int activeETraits = roundUpNearestMultiple(curNumETraits, blockSize);
				xyResult = gpuContext.executeDoubleKernel(
					flatETraits,
					flatSNPs,
					(int) nElems,
					localMemoryBytes,
					nrow,
					snpCol,
					new long[] { activeSnpColumns, activeETraits },
					new long[] { blockSize, blockSize });
			} finally {
				contextPool.releaseContext(gpuContext);
			}

			for (int snpNo = 0; snpNo < curNumSNPs; snpNo++)
			{
				int
					curSNPNo = snpNo + curSNPOffset,
					curCount[] = genotypeCount[snpNo],
					numData = numInds - curCount[0];
				for (int eTraitNo = 0; eTraitNo < curNumETraits; eTraitNo++)
				{
					int
						curOffset = eTraitNo * numETraitsPerBlock + snpNo * colPereQTL,
						curETraitNo = eTraitNo + curETraitOffset;
					// Since we're dealing with 0s and 1s, sumX1 = sumX1Sq, and sumX2 = sumX2Sq
					// Also, due to orthogonality, sumX1X2 = 0
					double
						sumX1Y = xyResult[curOffset],
						sumX2Y = xyResult[curOffset + 1],
						sumY = xyResult[curOffset + sumYOffset],
						sumYSq = xyResult[curOffset + sumYSqOffset],
						sumX1 = curCount[1] - curCount[3],
						sumX1Sq = curCount[1] + curCount[3],
						sumX2Sq = curCount[2],
						s_YY = sumYSq - sumY * sumY / numData, // = SST
						ssRegr,
						//F,
						//MSE,
						RSq, sumX2, s_X1Y, s_X1X1, s_X2Y, s_X2X2, s_X1X2, D, beta1 = 0, beta2 = 0;
					if (sumX1Sq > 0 && sumX2Sq > 0)
					{
						sumX2 = sumX2Sq;
						s_X1Y = sumX1Y - sumX1 * sumY / numData;
						s_X1X1 = sumX1Sq - sumX1 * sumX1 / numData;
						s_X2Y = sumX2Y - sumX2 * sumY / numData;
						s_X2X2 = sumX2Sq - sumX2 * sumX2 / numData;
						s_X1X2 = 0 - sumX1 * sumX2 / numData;
						D = (s_X1X1 * s_X2X2) - (s_X1X2 * s_X1X2);
						if (D == 0)
						{
							beta1 = s_X1Y / s_X1X1;
							beta2 = s_X2Y / s_X2X2;
							ssRegr = (s_X1X1 == 0) ? beta2 * s_X2Y : beta1 * s_X1Y;
						} else {
							beta1 = (((s_X1Y * s_X2X2) - (s_X2Y * s_X1X2)) / D);
							beta2 = (((s_X2Y * s_X1X1) - (s_X1Y * s_X1X2)) / D);
							//beta0 = (sumY - beta1 * sumX1Sq - beta2 * sumX2Sq) / numData,
							ssRegr = beta1 * s_X1Y + beta2 * s_X2Y; // Sum of squares of regression (SSR)
						}
						//MSE = (s_YY - ssRegr) / (numData - 3); // Mean Square Error (MSE) = (SST-SSR) / DF_E
						//if (MSE < 0) // Gotta be very close to zero
						//	MSE = -MSE;
						//F = ssRegr / (MSE * 2); // F = MSR / MSE, where MSR = SSR / DF_R
					} else {
						if (sumX1Sq > 0)
						{
							s_X1Y = sumX1Y - sumX1 * sumY / numData;
							s_X1X1 = sumX1Sq - sumX1 * sumX1 / numData;
							beta1 = s_X1Y / s_X1X1;
							ssRegr = (s_X1X1 == 0) ? 0 : beta1 * s_X1Y;
						} else if (sumX2Sq > 0)
						{
							sumX2 = sumX2Sq;
							s_X2Y = sumX2Y - sumX2 * sumY / numData;
							s_X2X2 = sumX1Sq - sumX2 * sumX2 / numData;
							beta2 = s_X2Y / s_X2X2;
							ssRegr = (s_X2X2 == 0) ? 0 : beta2 * s_X2Y;
						} else
							ssRegr = 0;
						//MSE = (s_YY - ssRegr) / (numData - 2); // Mean Square Error (MSE) = (SST-SSR) / DF_E
						//if (MSE < 0) // Gotta be very close to zero
						//	MSE = -MSE;
						//F = ssRegr / MSE; // F = MSR / MSE, where MSR = SSR / DF_R
					}
					RSq = ssRegr / s_YY;
					//if (RSq > 1)
					//	System.out.println("Here!");
					//if (Double.isInfinite(RSq))
					//{
						//System.out.println(String.format("sX1Y = %f, sX2Y = %f, sY = %f, sYSq = %f, sX1 = %f, sX1Sq = %f, sX2Sq = %f",
						//	sumX1Y, sumX2Y, sumY, sumYSq, sumX1, sumX1Sq, sumX2Sq));
					//} else
					if (RSq >= RSq0)
					{
						int idx = snpCounts[curSNPNo];
						if (isAdditive) {
							snpResults[curSNPNo][idx] = new double [] {RSq, beta1};
						} else {
							snpResults[curSNPNo][idx] = new double [] {RSq, beta1, beta2};
						}
						eTraitIndices[curSNPNo][idx] = curETraitNo;
						idx++;
						if (idx == eTraitIndices[curSNPNo].length) // Enlarge storage array if exhausted
						{
							int newLen = (idx * 3)/2;
							int[] newIndices = new int[newLen];
							double[][] newResult = new double[newLen][];
							System.arraycopy(snpResults[curSNPNo], 0, newResult, 0, idx);
							System.arraycopy(eTraitIndices[curSNPNo], 0, newIndices, 0, idx);
							snpResults[curSNPNo] = newResult;
							eTraitIndices[curSNPNo] = newIndices;
						}
						snpCounts[curSNPNo] = idx;
					}
				}
			}
		}
		for (int snpNo = 0; snpNo < curNumSNPs; snpNo++)
		{
			int
				curSNPNo = snpNo + curSNPOffset,
				idx = snpCounts[curSNPNo];
			if (idx < eTraitIndices[curSNPNo].length)
			{
				if (idx > 0)
				{
					eTraitIndices[curSNPNo] = QDataUtils.shortenArray(eTraitIndices[curSNPNo], idx);
					snpResults[curSNPNo] = QDataUtils.shortenArray(snpResults[curSNPNo], idx);
				} else {
					eTraitIndices[curSNPNo] = null;
					snpResults[curSNPNo] = null;
				}
			}
		}

		return true;
	}

	@Override
	public void setOwner(IJobOwner owner) {
	}

	@Override
	public IJobOwner getOwner() {
		return null;
	}

	@Override
	public void cancel() {
	}

	@Override
	public void run() {
		execute();
	}
}
