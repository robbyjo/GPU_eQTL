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
package gov.nih.opencl;

import java.nio.FloatBuffer;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem;
import com.nativelibs4java.opencl.CLQueue;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QOCLContext
{
	public static final String
		sDoubleExtAMD = "cl_amd_fp64", //$NON-NLS-1$
		sDoubleExtStandard = "cl_khr_fp64"; //$NON-NLS-1$
	protected CLDevice mDevice;
	protected CLContext mContext;
	protected CLQueue mQueue;
	protected CLKernel mKernel;
	protected long
		mMaxWorkGroupSize,
		mMaxWorkItemSizes[];
	protected int[]
		mLocalBlockSize,
		mGlobalBlockSize;
	protected String mDoubleHeader = null;

	public QOCLContext(CLDevice device)
	{
		mDevice = device;
		mContext = device.getPlatform().createContext(null, device);
		mQueue = mContext.createDefaultQueue();
		mMaxWorkGroupSize = device.getMaxWorkGroupSize();
		mMaxWorkItemSizes = device.getMaxWorkItemSizes();
		if (device.isDoubleSupported())
			mDoubleHeader = device.isDoubleSupportedAMD() ? sDoubleExtAMD : sDoubleExtStandard;
		mLocalBlockSize = detectLocalBlockSize();
		mGlobalBlockSize = detectGlobalBlockSize();
	}

	protected int[] detectLocalBlockSize()
	{
		String dummyKernel = "__kernel void temp(__global float* A) { int bx = get_group_id(0), by = get_group_id(1); A[bx*128+by] = bx * by; }"; //$NON-NLS-1$
		CLKernel kernel = createKernel(dummyKernel, "temp"); //$NON-NLS-1$
		int[]
			localBlockSizes = new int[] { 128, 128 },
			globalBlockSizes = new int[] { 128, 128 };
		CLBuffer<Float> bufferC = mContext.createFloatBuffer(CLMem.Usage.Output, FloatBuffer.wrap(new float[128*128]), false);
		do {
			try {
				kernel.setArg(0,  bufferC);
				kernel.enqueueNDRange(mQueue, globalBlockSizes, localBlockSizes);
				mQueue.flush();
				mQueue.finish();
				break;
			} catch (Throwable e) {
				if (localBlockSizes[0] == 2) {
					localBlockSizes = null;
					break;
				}
				localBlockSizes[0] = localBlockSizes[0] / 2;
				localBlockSizes[1] = localBlockSizes[1] / 2;
			}
		} while (true);
		bufferC.release();
		return localBlockSizes;
	}

	protected int[] detectGlobalBlockSize()
	{
		String dummyKernel = "__kernel void temp(__global float* A) { int bx = get_group_id(0), by = get_group_id(1); A[(bx&0x7F)*(by&0x7F)] = bx * by; }"; //$NON-NLS-1$
		CLKernel kernel = createKernel(dummyKernel, "temp"); //$NON-NLS-1$
		int
			curGuess = 128,
			minGuess = 0,
			maxGuess = 0;
		int[]
			localBlockSizes = mLocalBlockSize != null ? mLocalBlockSize : detectLocalBlockSize(),
			globalBlockSizes = new int[] { curGuess, curGuess };
		do {
			try {
				globalBlockSizes[0] = globalBlockSizes[1] = curGuess;
				CLBuffer<Float> bufferC = mContext.createFloatBuffer(CLMem.Usage.Output, FloatBuffer.wrap(new float[128*128]), false);
				kernel.setArg(0,  bufferC);
				kernel.enqueueNDRange(mQueue, globalBlockSizes, localBlockSizes);
				mQueue.flush();
				mQueue.finish();
				bufferC.release();
				if (curGuess >= 32*1024) // Limit to 32K
					break;
				minGuess = curGuess;
				curGuess = curGuess * 2;
			} catch (Throwable e) {
				maxGuess = curGuess;
				curGuess = (maxGuess - minGuess) / 2;
				// Ensure that curGuess is within 64-byte alignment boundary
				if ((curGuess & 0x3F) > 0)
					curGuess = (curGuess / 64) * 64;
				if (maxGuess - curGuess < 64)
					break;
			}
		} while (true);
		return globalBlockSizes;
	}

	public int[] getLocalBlockSize()
	{	return mLocalBlockSize; }

	public int[] getGlobalBlockSize()
	{	return mGlobalBlockSize; }

	public void setKernel(CLKernel k)
	{	mKernel = k; }

	public CLKernel createKernel(String program, String kernelName)
	{	return mContext.createProgram(program).build().createKernel(kernelName); }

	public CLKernel setKernelFromProgram(String program, String kernelName)
	{
		return mKernel = createKernel(program, kernelName);
	}

	public CLKernel getKernel()
	{	return mKernel; }

	public CLDevice getDevice()
	{	return mDevice; }

	public CLContext getContext()
	{	return mContext; }

	public CLQueue getQueue()
	{	return mQueue; }

	public long getMaxWorkGroupSize()
	{	return mMaxWorkGroupSize; }

	public long[] getMaxWorkItemSizes()
	{	return mMaxWorkItemSizes; }

	public String getDoubleExtensionHeader()
	{	return mDoubleHeader; }

	public void dispose()
	{
		if (mKernel != null)
			mKernel.release();
		mQueue.release();
		mContext.release();
		mDevice.release();
	}
}
