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
package gov.nih.table;

/**
 * @author Roby Joehanes
 */
public interface IViewableData extends Iterable<double[]>
{
	public IViewableData getBackingStore();
	/**
	 * Filtered column names
	 * @return
	 */
	public String[] getColumnNames();
	/**
	 * All column names
	 * @return
	 */
	public String[] getAllColumnNames();
	public String[] getRowNames();
	public String[] getAllRowNames();

	public String[][] getColumnCategories();

	public String[][] getAllColumnCategories();
	public double[][] extractData();
	public double[][] getUnfilteredData();
	public double[] getColumnWeights();
	public double[] getRowWeights();
	public double[] getAllColumnWeights();
	public double[] getAllRowWeights();
	public int getNumberOfRows();
	public int getNumberOfColumns();
	public IDataIterator<double[]> getDataIterator();
	public IDataIterator<double[]> getAllDataIterator();
	//public IViewableData createView(IDataFilter filter);
	//public IViewableData createView(int columnNo, double value);
	public int count(int columnNo, double value);
	public IViewableData cloneData();
	public String getName();
	public void setName(String s);

	/**
	 * Get rid of everything
	 */
	public void purge();
}
