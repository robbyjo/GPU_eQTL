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
 *
 */
package gov.nih.parallel;

/**
 * <P>Interface for generic parallel job. Any object who wants
 * to make use of parallel programming through our thread pool
 * implementation, QThreadPool, is required to implement this
 * interface.
 * 
 * <P>See QThreadPool comments.
 * 
 * @author Roby Joehanes
 */
public interface IGenericParallelJob
{
	/**
	 * Execute the job
	 */
	public void execute();

	/**
	 * Set the owner of this job
	 * @param owner
	 */
	public void setOwner(IJobOwner owner);

	/**
	 * Get the owner of this job
	 * @return
	 */
	public IJobOwner getOwner();

	/**
	 * Cancel job
	 */
	public void cancel();
}
