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
 * Synchronized counter. Has a form of an iterator
 * @author Roby Joehanes
 *
 */
public class QSynchronizedCounter implements ILockable
{
	protected int
		mStart,
		mStop,
		mCurrent;

	public QSynchronizedCounter(int start, int stop)
	{
		mStart = mCurrent = start;
		mStop = stop;
	}

	public synchronized int next()
	{	return mCurrent < mStop ? mCurrent++ : -1; }

	public synchronized boolean hasNext()
	{	return mCurrent < mStop; }

	public int getStart()
	{	return mStart; }

	public int getStop()
	{	return mStop; }

	public void reset()
	{	mCurrent = mStart; }

	@Override
	public String toString()
	{	return mCurrent + "/" + mStop; } //$NON-NLS-1$
}
