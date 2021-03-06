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
public class QSNPDataReal extends QSNPData {
	protected double[] mSNPData;

	@Override
	public void setSNPValues(String[] tokens)
	{
		mNumIndividuals = tokens.length;
		mSNPData = new double[mNumIndividuals];
		for (int i = 0; i < mNumIndividuals; i++)
			mSNPData[i] = Double.parseDouble(tokens[i]);
	}

	public void setSNPValues(double[] tokens)
	{
		mNumIndividuals = tokens.length;
		mSNPData = tokens;
	}

	@Override
	public int[] getSNPCodes() {
		return null;
	}

	@Override
	public double[] getSNPValues()
	{	return mSNPData; }

	@Override
	public int getNumMissing()
	{	return 0; }

	@Override
	public void purge()
	{	mSNPData = null; }

}
