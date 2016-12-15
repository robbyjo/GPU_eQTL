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
public class QParallelJobState
{
	protected QSynchronizedCounter mCounter;
	protected ISimpleParallelJob mJob;
	protected IJobOwner mOwner;
	protected int mNumThreads;

	public QParallelJobState(IJobOwner owner, ISimpleParallelJob job, QSynchronizedCounter ctr, int numThreads)
	{
		mOwner = owner;
		mCounter = ctr;
		mJob = job;
		mNumThreads = numThreads;
	}

	public boolean isResumable()
	{	return mJob != null && mCounter != null && mCounter.hasNext(); }

	public int getNumThreads()
	{	return mNumThreads; }

	public void setNumThreads(int i)
	{
		if (i > 0)
			mNumThreads = i;
	}

	public IJobOwner getOwner()
	{	return mOwner; }

	public ISimpleParallelJob getJob()
	{	return mJob; }

	public QSynchronizedCounter getCounter()
	{	return mCounter; }
}
