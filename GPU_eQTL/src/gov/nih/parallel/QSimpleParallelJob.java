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
 * 
 * @author Roby Joehanes
 *
 */
public class QSimpleParallelJob implements IGenericParallelJob
{
	protected ISimpleParallelJob mJob;
	protected IJobOwner mOwner;
	protected QSynchronizedCounter
		mCounter,
		mLock;
	protected int mID;

	protected QSimpleParallelJob(int id, IJobOwner owner, ISimpleParallelJob job, QSynchronizedCounter counter, QSynchronizedCounter lock)
	{
		setOwner(owner);
		mJob = job;
		mCounter = counter;
		mLock = lock;
		mID = id;
	}

	public void execute()
	{
		//System.out.println("In " + mID);
		int iterNo = mCounter.next();
		if (iterNo >= 0)
		{
			do
			{
				if ((!mJob.execute(iterNo)) || (mOwner != null && mOwner.isCanceled()))
					break;
				iterNo = mCounter.next();
			} while (iterNo >= 0);
		}

		mLock.next();
		synchronized (mLock)
		{
			mLock.notifyAll();
		}
		//System.out.println("Out " + mID);
	}

	public void setOwner(IJobOwner owner)
	{	mOwner = owner; }

	public IJobOwner getOwner()
	{	return mOwner; }

	/**
	 * By default, this is ignored
	 */
	public void cancel()
	{	}
}
