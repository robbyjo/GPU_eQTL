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
package gov.nih.table;

import java.util.Iterator;

import gov.nih.utils.QDataUtils;

import static qplugin.QPluginConstants.sDefaultDataName;

/**
 * The same as QTableData except it allows multiple row and column names
 * @author Roby Joehanes
 *
 */
public class QTableDataExtra implements IViewableDataExtra
{
	class QTableDataIterator implements IDataIterator<double[]>
	{
		protected int mIdx = 0;

		/* (non-Javadoc)
		 * @see java.util.Iterator#hasNext()
		 */
		public boolean hasNext()
		{	return mData != null && mIdx < mData.length; }

		/* (non-Javadoc)
		 * @see java.util.Iterator#next()
		 */
		public double[] next()
		{	return hasNext() ? mData[mIdx++] : null; }

		public double[] next(int n)
		{
			mIdx += (n - 1);
			if (mIdx >= mData.length)
				return null;
			return hasNext() ? mData[mIdx++] : null;
		}

		public Iterator<double[]> iterator()
		{	return this; }

		/* Unsupported
		 * @see java.util.Iterator#remove()
		 */
		public void remove()
		{	throw new UnsupportedOperationException(); }
	}

	protected static int mCounter = 0;
	protected String mName = sDefaultDataName + (mCounter++); // Title / name for this data
	protected double[][] mData; // row x columns
	protected double[]
		mColumnWeights,
		mRowWeights;
	protected String[][]
		mRowNames, // numRowNames x numRows
		mColNames; // numColNames x numColumns
	protected String[][] mColumnCategories;

	protected String[] mFilenames;

	public QTableDataExtra() {}

	public QTableDataExtra(double[][] data)
	{	this(data, null, null, null, null); }

	public QTableDataExtra(double[][] data, String[][] rowNames, String[][] colNames)
	{	this(data, rowNames, colNames, null, null); }

	public QTableDataExtra(double[][] data, String[][] rowNames, String[][] colNames, double[] rowWeights, double[] colWeights)
	{
		mData = data;
		mRowNames = rowNames;
		mColNames = colNames;
		mColumnCategories = new String[mData[0].length][];
		mRowWeights = rowWeights;
		mColumnWeights = colWeights;
	}

	public void purge()
	{
		mData = null;
		mRowNames = null;
		mColNames = null;
		mRowWeights = null;
		mColumnWeights = null;
		mColumnCategories = null;
		mFilenames = null;
		mName = null;
	}

	public void setName(String n)
	{	mName = n; }

	public String getName()
	{	return mName; }

	public double[][] getData()
	{	return mData; }

	public double[][] getUnfilteredData()
	{	return mData; }

	public double[][] extractData()
	{	return mData; }

	public void setData(double[][] data)
	{	mData = data; }

	public String[][] getRowNames()
	{	return mRowNames; }

	public String[][] getAllRowNames()
	{	return mRowNames; }

	public void setRowNames(String[][] s)
	{	mRowNames = s; }

	public String[][] getColumnNames()
	{	return mColNames; }

	public void setColumnNames(String[][] s)
	{	mColNames = s; }

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getAllColumnNames()
	 */
	public String[][] getAllColumnNames()
	{	return mColNames; }

	public String[] getColumnCategory(int colNo)
	{	return mColumnCategories[colNo]; }

	public void setColumnCategories(String[][] c)
	{	mColumnCategories = c; }

	public void setColumnCategory(int colNo, String[] args)
	{	mColumnCategories[colNo] = args; }

	public String[][] getColumnCategories()
	{	return mColumnCategories; }

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getAllColumnCategories()
	 */
	public String[][] getAllColumnCategories()
	{	return mColumnCategories; }

	public double[] getColumnWeights()
	{	return mColumnWeights; }

	public double[] getAllColumnWeights()
	{	return mColumnWeights; }

	public void setColumnWeights(double[] w)
	{	mColumnWeights = w; }

	public double[] getRowWeights()
	{	return mRowWeights; }

	public double[] getAllRowWeights()
	{	return mRowWeights; }

	public void setRowWeights(double[] w)
	{	mRowWeights = w; }

	public int getNumberOfRows()
	{	return mData.length; }

	public int getNumberOfColumns()
	{	return mData[0].length; }

	public IViewableDataExtra getBackingStore()
	{	return this; }

	/**
	 * Insert a column at the last position
	 * @param columnNames
	 * @param columnCats
	 * @param columnData
	 */
	public void insertColumn(String[] columnNames, String[] columnCats, double[] columnData)
	{
		int
			numRows = columnData.length,
			numCols = mData[0].length;
		if (numRows != mData.length)
			throw new RuntimeException(); // Row mismatch
		String[][] newColNames = null;
		if (mColNames != null)
		{
			assert (mColNames.length == columnNames.length);
			newColNames = new String[mColNames.length][numCols + 1];
			for (int i = 0; i < mColNames.length; i++)
			{
				System.arraycopy(mColNames[i], 0, newColNames[i], 0, numCols);
				newColNames[i][numCols] = columnNames[i];
			}
		}
		String[][] newColCats = new String[numCols + 1][];
		if (mColumnCategories != null)
			System.arraycopy(mColumnCategories, 0, newColCats, 0, numCols);
		newColCats[numCols] = columnCats;

		double[][] newData = new double[numRows][numCols + 1];
		for (int i = 0; i < numRows; i++)
		{
			System.arraycopy(mData[i], 0, newData[i], 0, numCols);
			mData[i][numRows] = columnData[i];
		}

		mData = newData;
		mColNames = newColNames;
		mColumnCategories = newColCats;
	}

	public String[] getFilenames()
	{	return mFilenames; }

	public void setFilenames(String[] filenames)
	{	mFilenames = filenames; }

	public IDataIterator<double[]> getDataIterator()
	{	return new QTableDataIterator(); }

	public IDataIterator<double[]> getAllDataIterator()
	{	return new QTableDataIterator(); }

	/* (non-Javadoc)
	 * @see java.lang.Iterable#iterator()
	 */
	public Iterator<double[]> iterator()
	{	return new QTableDataIterator(); }

	/**
	 * Clone this data. Column / row names and column categories are NOT cloned.
	 */
	public QTableDataExtra cloneData()
	{
		QTableDataExtra newData = new QTableDataExtra(QDataUtils.safeClone(mData), mRowNames, mColNames);
		newData.mColumnCategories = mColumnCategories;
		newData.mColumnWeights = mColumnWeights != null ? mColumnWeights.clone() : null;
		newData.mRowWeights = mRowWeights != null ? mRowWeights.clone() : null;
		return newData;
	}

	public int count(int columnNo, double value)
	{
		int
			count = 0,
			numRows = mData.length;
		for (int i = 0; i < numRows; i++)
			if (mData[i][columnNo] == value)
				count++;
		return count;
	}
}
