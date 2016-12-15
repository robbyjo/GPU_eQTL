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

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QOCLContextPool
{
	protected List<QOCLContext> mAllContextList = new ArrayList<QOCLContext>();
	protected LinkedList<QOCLContext> mContextList = new LinkedList<QOCLContext>();
	protected Object mObjectLock = new Object();

	public QOCLContextPool(QOCLContext[] contexts)
	{
		for (QOCLContext context: contexts)
		{
			mContextList.add(context);
			mAllContextList.add(context);
		}
	}

	public QOCLContextPool(Collection<QOCLContext> contexts)
	{	mContextList.addAll(contexts); }

	public List<QOCLContext> getAllContexts()
	{	return mAllContextList; }

	public int getNumContexts()
	{	return mAllContextList.size(); }

	public void releaseContext(QOCLContext ctx)
	{
		synchronized(this)
		{
            if (!mContextList.contains(ctx))
            	mContextList.add(ctx);
		}
		synchronized(mObjectLock)
		{
			mObjectLock.notify();
		}
   }

	public QOCLContext reserveContext()
	{
		QOCLContext ctx;
		synchronized(mObjectLock)
		{
			while ((ctx = obtainNextContext()) == null)
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
			return ctx;
		}
	}

	protected synchronized QOCLContext obtainNextContext()
	{
		if (mContextList.size() > 0)
			return mContextList.removeFirst();
		return null;
	}

	public synchronized void dispose()
	{
		for (QOCLContext context: mAllContextList)
			try {
			context.dispose();
			} catch (Exception e) {}
	}
}
