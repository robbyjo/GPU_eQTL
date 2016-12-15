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
package gov.nih.eqtl.datastructure;


/**
 * 
 * @author Roby Joehanes
 *
 */
public abstract class QSNPData {
	protected String mID;
	protected int
		mNumIndividuals,
		mChromID,
		mLocation;

	public static final int kMissingCode = -1;

	public QSNPData() {}

	public abstract void setSNPValues(String[] tokens);

	public abstract int[] getSNPCodes();

	/**
	 * Same as getSNPCodes, but returns double values
	 * @return
	 */
	public abstract double[] getSNPValues();

	public abstract int getNumMissing();

	public int getNumIndividuals()
	{	return mNumIndividuals; }

	public double getPctMissing()
	{	return getNumMissing() * 1.0 / mNumIndividuals; }

	public int getChromID()
	{	return mChromID; }

	public void setChromID(int c)
	{	mChromID = c; }

	public int getLocation()
	{	return mLocation; }

	public void setLocation(int i)
	{	mLocation = i; }

	public String getID()
	{	return mID; }

	public void setID(String id)
	{	mID = id; }

	public abstract void purge();
}
