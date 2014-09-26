/*
 * Roby Joehanes
 * 
 * Copyright 2007 Roby Joehanes
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
package gov.nih.utils;

import static java.lang.Math.pow;

/**
 * Memory statistics
 * @author Roby Joehanes
 *
 */
public class QMemoryUsage
{
	private static final double kMegaBytes = 1024 * 1024;
	private static final String sUsedMemory = "UsedMemory"; //$NON-NLS-1$

	public long
		mUsedHeap,
		mUsedNonHeap,
		mUsedBaseLine,
		mAvailable,
		mTotal;

	public QMemoryUsage(long h, long nh, long bl, long t)
	{
		mUsedHeap = h;
		mUsedNonHeap = nh;
		mUsedBaseLine = bl;
		mTotal = t;
		mAvailable = mTotal - mUsedHeap - mUsedNonHeap;
	}

	/**
	 * Rounding up to the desired decimal places
	 * @param val
	 * @param places
	 * @return
	 * @author Roby Joehanes
	 */
	static final double ceil(double val, int places)
	{
		double factor = pow(10, places);
		return Math.ceil(val*factor)/factor;
	}

	public static final String convertToMiB(long b)
	{	return String.valueOf(ceil(b / kMegaBytes, 2)); }

	@Override
	public String toString()
	{
		long total = mUsedHeap + mUsedNonHeap;
		return String.format(sUsedMemory,
			total, convertToMiB(total), mUsedHeap, convertToMiB(mUsedHeap),
			mUsedNonHeap, convertToMiB(mUsedNonHeap), mUsedBaseLine, convertToMiB(mUsedBaseLine),
			mAvailable, convertToMiB(mAvailable), mTotal, convertToMiB(mTotal));
	}
}
