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
package gov.nih.parallel;

/**
 * Synchronized 2D counter. Has a form of an iterator
 * @author Roby Joehanes
 *
 */
public class QSynchronized2DCounter implements ILockable
{
	protected int
		mMaxI,
		mMaxJ,
		mI,
		mJ;

	public QSynchronized2DCounter(int x, int y)
	{
		mMaxI = x; mMaxJ = y;
	}

	/**
	 * Returns <pre>((long) i) << 32 | j</pre>. To extract:
	 * <pre>i = (int) (x >> 32); j = (int) (x & 0xFFFF);</pre>
	 * 
	 * <P>This is 1000x faster than returning new int[] {i, j}.
	 * It's counter-intuitive. :(
	 * @return
	 */
	public synchronized long next()
	{
		if (mI < mMaxI && mJ < mMaxJ)
		{
			long x = (((long) mI) << 32) | mJ;
			mJ++;
			if (mJ == mMaxJ)
			{
				mJ = 0;
				mI++;
			}
			return x;
		}
		return -1l;
	}

	public synchronized boolean hasNext()
	{	return mI < mMaxI && mJ < mMaxJ; }
}
