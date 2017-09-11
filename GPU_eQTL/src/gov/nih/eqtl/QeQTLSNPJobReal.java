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

import java.io.IOException;
import java.io.Writer;
import java.nio.DoubleBuffer;
import java.util.Arrays;
import java.util.List;

import jdistlib.T;

import org.bridj.Pointer;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem;
import com.nativelibs4java.opencl.CLQueue;

import gov.nih.eqtl.datastructure.QGeneExpressionData;
import gov.nih.eqtl.datastructure.QGeneticSNPData;
import gov.nih.eqtl.datastructure.QSNPData;
import gov.nih.opencl.QOCLContext;
import gov.nih.opencl.QOCLContextPool;
import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.IJobOwner;
import gov.nih.parallel.QSynchronizedCounter;
import gov.nih.utils.QSystemUtils;

import static java.lang.Math.log;
import static java.lang.Math.log10;
import static java.lang.Math.min;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;
import static gov.nih.opencl.QOpenCL.kDefaultAlignment;
import static gov.nih.utils.QStringUtils.sLn;

public class QeQTLSNPJobReal implements IGenericParallelJob, Runnable {
	static final double
		log10Of2 = log10(2),
		logOf10 = log(10);
	protected QGeneticSNPData popn;
	protected QGeneExpressionData expDataTbl;
	protected QOCLContextPool contextPool;
	protected CLEvent gpuEvents;
	protected double expData[][], expSD[], snpSD[], flatETraits[], flatSNPs[], RSq0;
	protected int numETraitsPerBlock, numETraits, numSNPsPerBlock, numSNPs, numInds, nrow, blockSize, snpCol;
	protected CLBuffer<Double> bufferA = null, bufferB = null, bufferC = null;
	protected int dfo, dfe, netdfe;
	protected boolean isAdditive;
	protected QSynchronizedCounter counter;
	static final int sizeof_data = 8;
	protected Writer fw;

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

	public QeQTLSNPJobReal(QGeneticSNPData popn, QGeneExpressionData expDataTbl, double[] expSD, double[] snpSD, QOCLContextPool contextPool,
		int numETraitsPerBlock, int numSNPsPerBlock, int blockSize,
		int dfo, int dfe, double RSq0, boolean isAdditive, Writer fw, QSynchronizedCounter ct)
	{
		this.popn = popn;
		this.expDataTbl = expDataTbl;
		this.expData = expDataTbl.getData();
		this.contextPool = contextPool;
		this.numETraitsPerBlock = numETraitsPerBlock;
		this.numSNPsPerBlock = numSNPsPerBlock;
		this.isAdditive = isAdditive;
		this.expSD = expSD;
		this.snpSD = snpSD;
		this.dfe = dfe;
		this.dfo = dfo;
		netdfe = dfe - dfo;
		this.fw = fw;
		numSNPs = popn.getNumSNPs();
		numInds = popn.getNumIndividuals();
		numETraits = expData.length;
		snpCol = numSNPsPerBlock;
		nrow = roundUpNearestMultiple(numInds, kDefaultAlignment);
		this.blockSize = blockSize;
		flatETraits = new double[numETraitsPerBlock * nrow];
		flatSNPs = new double[snpCol * nrow];
		this.RSq0 = RSq0;
		counter = ct;
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

		//long time1, time2;
		//time1 = System.currentTimeMillis();
		List<QSNPData> snpList = popn.getSNPs();
		String[] geneID = expDataTbl.getGeneIDs();
		double[] xyResult;
		for (int k = 0; k < curNumSNPs; k++)
		{
			double[] snps = snpList.get(curSNPOffset+k).getSNPValues();
			for (int l = 0; l < numInds; l++)
				flatSNPs[l * snpCol + k] = snps[l] - 1;
		}
		//time2 = System.currentTimeMillis();
		//System.out.println("SNP filling time = " + (time2 - time1));

		for (int curETraitOffset = 0; curETraitOffset < numETraits; curETraitOffset += numETraitsPerBlock)
		{
			System.out.println(curSNPOffset+","+curETraitOffset);
			//time1 = System.currentTimeMillis();
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
			QOCLContext oclContext = contextPool.reserveContext();
			synchronized (oclContext) {
				CLContext context = oclContext.getContext();
				CLKernel kernel = oclContext.getKernel();
				CLQueue queue = oclContext.getQueue();
				bufferA = context.createDoubleBuffer(CLMem.Usage.Input, DoubleBuffer.wrap(flatETraits), true);
				bufferB = context.createDoubleBuffer(CLMem.Usage.Input, DoubleBuffer.wrap(flatSNPs), true);
				bufferC = context.createDoubleBuffer(CLMem.Usage.Output, nElems);
				queue.finish();
				kernel.setArg(0,  bufferC);
				kernel.setArg(1,  bufferA);
				kernel.setArg(2,  bufferB);
				kernel.setLocalArg(3, (blockSize+1) * (4*blockSize) * sizeof_data);
				kernel.setLocalArg(4, (blockSize+1) * (4*blockSize) * sizeof_data);
				kernel.setArg(5, nrow);
				kernel.setArg(6, snpCol);
				gpuEvents = kernel.enqueueNDRange(queue, new int[] { snpCol, numETraitsPerBlock }, new int[] {blockSize, blockSize});
				queue.flush();
				queue.finish();
				Pointer<Double> outBuffer = bufferC.read(queue, gpuEvents);
				queue.finish();
				xyResult = outBuffer.getDoubles();
				bufferA.release();
				bufferB.release();
				bufferC.release();
			}
			contextPool.releaseContext(oclContext);
			//time2 = System.currentTimeMillis();
			//System.out.println("GPU run time = " + (time2 - time1));

			//time1 = System.currentTimeMillis();
			for (int snpNo = 0; snpNo < curNumSNPs; snpNo++)
			{
				int curSNPNo = snpNo + curSNPOffset;
				String snpID = snpList.get(curSNPNo).getID();
				for (int eTraitNo = 0; eTraitNo < curNumETraits; eTraitNo++)
				{
					int
						curOffset = eTraitNo * numETraitsPerBlock + snpNo,
						curETraitNo = eTraitNo + curETraitOffset;
					String probesetID = geneID[curETraitNo];
					double
						rawEffect = xyResult[curOffset],
						RSq = rawEffect * rawEffect;
					if (RSq >= RSq0)
					{
						if (QeQTLAnalysis.rsqOnly) {
							if (QeQTLAnalysis.simplifyResult) {
								RSq = Math.round(RSq * 10000) / 10000.0;
							}
							synchronized(fw) {
								try {
									fw.write(snpID + "," + probesetID + "," + RSq + (rawEffect < 0 ? ",-" : ",+") + sLn);
									fw.flush();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						} else {
							double
								realEffect = rawEffect * expSD[curETraitNo] / snpSD[curSNPNo],
								t = sqrt(RSq*dfe / (1 - RSq)),
								p = log10Of2 + T.cumulative(t, netdfe, false, true) / logOf10;
							if (QeQTLAnalysis.simplifyResult) {
								RSq = Math.round(RSq * 10000) / 10000.0;
								realEffect = Math.round(realEffect * 10000) / 10000.0;
								t = Math.round(t * 10000) / 10000.0;
								p = Math.round(p * 10000) / 10000.0;
							}
							synchronized(fw) {
								try {
									if (QeQTLAnalysis.pvalOnly) {
										fw.write(p + sLn);
									} else {
										fw.write(snpID + "," + probesetID + "," + RSq + "," + realEffect +
											"," + (signum(rawEffect) * t) + "," + p + sLn);
									}
									fw.flush();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						}
					}
				}
			}
			xyResult = null;
			QSystemUtils.runGC();
			//time2 = System.currentTimeMillis();
			//System.out.println("Post-processing time = " + (time2 - time1));
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
