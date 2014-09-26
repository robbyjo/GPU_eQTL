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
/**
 * 
 */
package gov.nih.parallel;

import java.util.*;

/**
 * <P>Barrier synchronization pattern.
 * Taken straight out of my CIS 720 class homework
 * with minimal changes.
 * 
 * <P>WARNING: NEVER nest two or more synchronized methods
 * unless you know what you're doing!
 * 
 * @author Roby Joehanes
 */
public class QBarrier
{
	protected int
		mNumThreadEnter = 0,
		mNumThreadExit = 0,
		mTotalThread;
	protected QThreadPool mThreadPool;
	protected Object mLock = new Object();
    protected LinkedList<IGenericParallelJob> mQueue = new LinkedList<IGenericParallelJob>();

	public QBarrier(QThreadPool pool, int max)
	{
		this.mTotalThread = max;
		mThreadPool = pool;
	}

	public synchronized void enter()
	{
		mNumThreadEnter++;
		signal();
	}

	protected synchronized void signal()
	{
		if (!mQueue.isEmpty() && mNumThreadExit < (mNumThreadEnter / mTotalThread) * mTotalThread)
		{
			// Changed the usual template to this since the signal
			// part is always signal all
			while (!mQueue.isEmpty())
			{
				IGenericParallelJob job = mQueue.removeFirst();
				mThreadPool.addJob(job);
				mNumThreadExit++;
			}
		}
	}

	public synchronized boolean leave(IGenericParallelJob job)
	{
		if (!(mNumThreadExit < (mNumThreadEnter / mTotalThread) * mTotalThread))
		{
			mQueue.add(job);
			return false;
		}
		mNumThreadExit++;
		return true;
	}

	public synchronized boolean leaveWithSignal(IGenericParallelJob job)
	{
		if (!(mNumThreadExit < (mNumThreadEnter / mTotalThread) * mTotalThread))
		{
			mQueue.add(job);
			return false;
		}
		mNumThreadExit++;
		signal();
		return true;
	}
}
