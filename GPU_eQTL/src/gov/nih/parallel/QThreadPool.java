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
/*
 * Created on Aug 12, 2004
 */
package gov.nih.parallel;

import java.util.Collection;
import java.util.LinkedList;

import gov.nih.utils.QSystemUtils;

/**
 * <h3>Why thread pool?</h3>
 * 
 * <P>Spawning threads at will may work in a program if it is simple
 * enough and if the thread terminates quickly that Java and/or
 * operating system, which ultimately handles it, can recycle the
 * threads quickly.
 * 
 * <P>However, if the number of active threads is massive, even if
 * the operating system is efficient, it can slow down the program
 * considerably. On the other hand, we may need to have many
 * tasks run in parallel. A solution to this is to create a
 * limited amount of threads and a class to manage them. This
 * solution is called "thread pooling", which was first conceived
 * in embedded systems.
 * 
 * <h3>What this class do</h3>
 * 
 * <P>This class is a simple implementation of thread pooling.
 * By default, it creates (the number of CPU cores) threads and manage any parallel requests
 * (i.e. parallel jobs) accordingly. All parallel jobs are required
 * to implement IGenericParallelJob. Parallel jobs can be added
 * to the thread pool through addJob or addAllJobs methods.
 * The thread pool will then manage the job in round-robin fashion.
 * 
 * <P>All parallel jobs are required to maintain thread safety,
 * only with regards of deadlock. This is because the current version
 * of Java doesn't allow any deadlock detection and therefore it is
 * virtually impossible to know which thread is deadlocked. If a
 * thread is deadlocked, then it effectively blocks that particular
 * thread. Though it will not block other running threads, it still
 * reduces the number of threads available in the pool.
 * 
 * <P>Threads with race conditions do no harm to the pool, unless if
 * it eventually causes deadlock. Certainly, race conditions should
 * be avoided at all cost.
 * 
 * <P>The jobs must put its parallel task in the execute() method.
 * The thread pool will appoint available threads to invoke
 * execute() method and the assigned thread will wait until the
 * execute() method is finished or go to sleep. If the job is
 * finished, then the job quits the thread pool. However, if the
 * job is just sleeping, then the assigned thread will do some
 * other job and put the sleeping job to the back of the queue.
 * 
 * <P>WARNING: NEVER nest two or more synchronized methods
 * unless you know what you're doing! Often, this will result in
 * deadlock by unwittingly not surrendering the monitor of the
 * first synchronized object.
 * 
 * <h3>Example usage</h3>
 * 
 * <P>Below is an example of a worker thread implementation.
 * Certainly, you can put in the constructor any parameters
 * or values necessary for the jobs to run. However, do NOT
 * invoke execute() method directly.
 * 
 * <pre>public MyParallelJob implements IGenericParallelJob {
 *     public MyParallelJob() {
 *         // Put constructor here.
 *     }
 * 
 *     public boolean cancelJob() {
 *         return false; // By default uncancelable
 *     }
 *
 *     public void execute() {
 *         // Your parallel job goes here.
 *         // Do synchronization routines as necessary.
 *     }
 * }</pre>
 * 
 * <P>First, you can create several new worker objects. It
 * can be more than 8. For example:
 * 
 * <pre>List<MyParallelJob> jobList = new ArrayList<MyParallelJob>();
 * for (int i = 0; i < 1000; i++)
 *     jobList.add(new MyParallelJob());
 * </pre>
 *  
 * <P>Then, add the jobs to the thread pool. You can use the
 * default thread pool at QPoolThread.mDefaultThreadPool or
 * you can create your own thread pool object. If you use the
 * default thread pool, here's the code:
 * 
 * <pre>QThreadPool.addAllJobs(jobList);</pre>
 * 
 * <P>That's it.
 * 
 * @author Roby Joehanes
 */
public class QThreadPool
{
	// Don't make thread pool instance final. This is just in case of
	// flexibility to have performance / versatility tradeoff.
	public static final int
		kDefaultNumThreads = QSystemUtils.kNumCPUCores + 1,
		kMaxNumThreads = 1024;
	public static QThreadPool mDefaultThreadPool = new QThreadPool();

	protected LinkedList<IGenericParallelJob> mJobQueue = new LinkedList<IGenericParallelJob>();
	protected Object mObjectLock = new Object();
	protected boolean mAllowDuplicateJob = false;
	protected QPoolThread[] mThreads;

	/**
	 * Create a thread pool with the default number of threads (i.e. 8).
	 */
	public QThreadPool()
	{	this(kDefaultNumThreads); }

	/**
	 * Create a thread pool with the specified number of threads
	 * @param numThreads
	 */
	public QThreadPool(int numThreads)
	{
		QPoolThread[] threads = mThreads = new QPoolThread[numThreads];
		for (int i = 0; i < numThreads; i++)
		{
			threads[i] = new QPoolThread();
			threads[i].start();
		}
	}

    public synchronized void setAllowDuplicateJob(boolean b)
    {   mAllowDuplicateJob = b; }

    public boolean getAllowDuplicateJob()
    {   return mAllowDuplicateJob; }

    /**
     * Add a single job
     * @param job
     */
	public void addJob(IGenericParallelJob job)
	{
		if (mThreads == null)
			throw new RuntimeException(); // This thread pool is already terminated!
		synchronized(this)
		{
            if (mAllowDuplicateJob || !mJobQueue.contains(job))
                mJobQueue.add(job);
		}
		synchronized(mObjectLock)
		{
			mObjectLock.notify();
		}
   }

    /**
     * Add a list of jobs at once
     * @param jobs
     */
	public void addAllJobs(Collection<? extends IGenericParallelJob> jobs)
	{
		for (IGenericParallelJob job: jobs)
			addJob(job);
   }

	protected IGenericParallelJob getJob()
	{
		IGenericParallelJob job;
		synchronized(mObjectLock)
		{
			while ((job = obtainNextJob()) == null)
			{
				try
				{
					mObjectLock.wait();
				}
				catch (InterruptedException e)
				{
					// Ignore the interruption exception -- the thread has woken up
				}
			}
			return job;
		}
	}

	protected synchronized IGenericParallelJob obtainNextJob()
	{
		if (mJobQueue.size() > 0)
			return mJobQueue.removeFirst();
		return null;
	}

	public int getNumThreads()
	{	return mThreads == null? 0 : mThreads.length; }

	/**
	 * Terminate the thread pool
	 */
	public void terminate()
	{
		int numThreads = mThreads.length;
		for (int i = 0; i < numThreads; i++)
		{
			mThreads[i].terminate();
			mThreads[i] = null;
		}
		mThreads = null;
	}

	/**
	 * Create a barrier pattern for a specified number of parallel jobs.
	 * @see gov.nih.parallel.QBarrier
	 * @param numJobs
	 * @return
	 */
	public QBarrier createBarrier(int numJobs)
	{	return mThreads == null ? null : new QBarrier(this, numJobs); }

	/**
	 * The actual worker thread
	 * @author Roby Joehanes
	 */
	protected class QPoolThread extends Thread
	{
		protected boolean mIsRunning = true;

		@Override
		public void run()
		{
			while (mIsRunning)
			{
				getJob().execute();
			}
		}

		public void terminate()
		{	mIsRunning = false; }
	}
}
