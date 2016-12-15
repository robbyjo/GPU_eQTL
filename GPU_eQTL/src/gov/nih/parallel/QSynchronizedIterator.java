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

import java.util.Iterator;

/**
 * Synchronized wrapper for an iterator
 * @author Roby Joehanes
 */
public class QSynchronizedIterator<T> implements Iterator<T>, ILockable
{
	protected Iterator<T> mIterator;

	public QSynchronizedIterator(Iterator<T> iter)
	{	mIterator = iter; }

	public synchronized boolean hasNext()
	{	return mIterator.hasNext(); }

	public synchronized T next()
	{	return mIterator.hasNext() ? mIterator.next() : null; }

	public synchronized void remove()
	{	mIterator.remove(); }
}
