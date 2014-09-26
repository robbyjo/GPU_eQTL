/*
 * Roby Joehanes
 * 
 * Copyright 2007 Roby Joehanes.
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

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QParallelUtils
{
	public static final void waitOnLock(ILockable lock)
	{
		do
		{
			synchronized(lock)
			{
				try
				{
					if (!lock.hasNext())
						break;
					lock.wait();
				}
				catch (Exception e) {} // Interrupted exception is ignored
			}
		} while (true);
	}

	/**
	 * Set up a simple parallel job that repeats from <tt>startIter</tt> to <tt>endIter</tt> and run in <tt>numThreads</tt>
	 * threads.
	 * @param owner
	 * @param job
	 * @param startIter
	 * @param endIter
	 * @param numThreads
	 * @return resumable state
	 */
	public static final QParallelJobState executeSimpleParallelJob(IJobOwner owner, ISimpleParallelJob job, int startIter, int endIter, int numThreads)
	{
		List<IGenericParallelJob> jobList = new ArrayList<IGenericParallelJob>();
		QSynchronizedCounter
			counter = new QSynchronizedCounter(startIter, endIter),
			lock = new QSynchronizedCounter(0, numThreads);
		for (int i = 0; i < numThreads; i++) {
			QSimpleParallelJob curJob = new QSimpleParallelJob(i, owner, job.clone(), counter, lock);
			jobList.add(curJob);
		}
		QThreadPool.mDefaultThreadPool.addAllJobs(jobList);
		waitOnLock(lock);
		return counter.hasNext() ? new QParallelJobState(owner, job, counter, numThreads) : null;
	}

	public static final QParallelJobState resumeParallelJob(QParallelJobState state)
	{
		if (state == null || !state.isResumable())
			return null;
		List<IGenericParallelJob> jobList = new ArrayList<IGenericParallelJob>();
		int numThreads = state.getNumThreads();
		IJobOwner owner = state.getOwner();
		ISimpleParallelJob job = state.getJob();
		QSynchronizedCounter
			counter = state.getCounter(),
			lock = new QSynchronizedCounter(0, numThreads);
		for (int i = 0; i < numThreads; i++) {
			QSimpleParallelJob curJob = new QSimpleParallelJob(i, owner, job.clone(), counter, lock);
			jobList.add(curJob);
		}
		QThreadPool.mDefaultThreadPool.addAllJobs(jobList);
		waitOnLock(lock);
		return counter.hasNext() ? new QParallelJobState(owner, job, counter, numThreads) : null;
	}

	/**
	 * Same as <tt>executeSimpleParallelJob(owner, job, numIterations, QThreadPool.kDefaultNumThreads);</tt>
	 * @param owner
	 * @param job
	 * @param numIterations
	 */
	public static final QParallelJobState executeSimpleParallelJob(IJobOwner owner, ISimpleParallelJob job, int numIterations)
	{	return executeSimpleParallelJob(owner, job, 0, numIterations, QThreadPool.kDefaultNumThreads); }

	public static final QParallelJobState executeSimpleParallelJob(IJobOwner owner, ISimpleParallelJob job, int startIter, int endIter)
	{	return executeSimpleParallelJob(owner, job, startIter, endIter, QThreadPool.kDefaultNumThreads); }

	/**
	 * Same as <tt>executeSimpleParallelJob(null, job, numIterations, QThreadPool.kDefaultNumThreads);</tt>
	 * @param job
	 * @param numIterations
	 */
	public static final QParallelJobState executeSimpleParallelJob(ISimpleParallelJob job, int numIterations)
	{	return executeSimpleParallelJob(null, job, 0, numIterations, QThreadPool.kDefaultNumThreads); }

	public static final QParallelJobState executeSimpleParallelJob(ISimpleParallelJob job, int startIter, int endIter)
	{	return executeSimpleParallelJob(null, job, startIter, endIter, QThreadPool.kDefaultNumThreads); }

	/**
	 * Execute the job serially. For debugging
	 * @param job
	 * @param numIterations
	 */
	public static final void executeJobSerially(ISimpleParallelJob job, int numIterations)
	{	executeJobSerially(job, 0, numIterations); }

	public static final void executeJobSerially(ISimpleParallelJob job, int startIter, int endIter)
	{
		for (int iterNo = startIter; iterNo < endIter; iterNo++)
			job.execute(iterNo);
	}
}
