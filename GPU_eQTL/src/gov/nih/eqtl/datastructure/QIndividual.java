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
 * Created on Feb 27, 2006
 *
 */
package gov.nih.eqtl.datastructure;

/**
 * A class that represents an individual in a population
 * @author Roby Joehanes
 */
public class QIndividual implements Comparable<QIndividual>
{	
	protected QIndividual[] mParents;
	protected QFamily mFamily;
	protected String mID;
	protected ESex mSex;

	protected QIndividual() {}

	public QIndividual(QFamily family, String id, ESex sex)
	{
		setFamily(family);
		setID(id);
		setSex(sex);
	}

	@Override
	public QIndividual clone()
	{
		QIndividual newIndividual = new QIndividual();
		newIndividual.mParents = mParents;
		newIndividual.mID = mID;
		newIndividual.mSex = mSex;
		newIndividual.mFamily = mFamily;
		return newIndividual;
	}

	public int compareTo(QIndividual obj)
	{	return 0; }

	public QIndividual[] getParents() {
		return mParents;
	}

	public void setParents(QIndividual[] parents) {
		mParents = parents;
	}

	public String getID() {
		return mID;
	}

	public void setID(String id) {
		mID = id;
	}

	public ESex getSex() {
		return mSex;
	}

	public void setSex(ESex sex) {
		mSex = sex;
	}

	public QFamily getFamily() {
		return mFamily;
	}

	public void setFamily(QFamily family) {
		mFamily = family;
	}
}
