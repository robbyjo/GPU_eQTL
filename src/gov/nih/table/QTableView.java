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

import java.util.BitSet;
import java.util.Iterator;

/**
 * Fast table view abstraction. The idea is almost like database view,
 * where only the relevant data indices are stored. The only caveat is
 * that this data structure does NOT monitor consistency nor does it
 * fail-fast when modification occurred.
 * 
 * @author Roby Joehanes
 */
public class QTableView implements IViewableData
{
	class QTableViewIterator implements IDataIterator<double[]>
	{
		protected int
			mIdx = mRowChooser.nextSetBit(0),
			mIterIdx = -1;
		protected IDataIterator<double[]> mParentIterator = mBackingStore.getDataIterator();

		/* (non-Javadoc)
		 * @see java.util.Iterator#hasNext()
		 */
		public boolean hasNext()
		{	return mIdx >= 0; }

		/* (non-Javadoc)
		 * @see java.util.Iterator#next()
		 */
		public double[] next()
		{
			if (!hasNext())
				return null;
			double[] data = mParentIterator.next(mIdx - mIterIdx);
			mIterIdx = mIdx;
			mIdx = mRowChooser.nextSetBit(mIdx + 1);
			return data;
		}

		public double[] next(int n)
		{
			if (!hasNext())
				return null;
			for (int i = 1; i < n && mIdx >= 0; i++)
				mIdx = mRowChooser.nextSetBit(mIdx + 1);
			if (mIdx < 0)
				return null;
			double[] data = mParentIterator.next(mIdx - mIterIdx);
			mIterIdx = mIdx;
			mIdx = mRowChooser.nextSetBit(mIdx + 1);
			return data;
		}

		public Iterator<double[]> iterator()
		{	return this; }

		/* Unsupported
		 * @see java.util.Iterator#remove()
		 */
		public void remove()
		{	throw new UnsupportedOperationException(); }
	}

	protected IViewableData mBackingStore;
	protected BitSet
		mColumnChooser,
		mRowChooser;
	protected String[] mFilenames;

	public QTableView(IViewableData data, BitSet colChooser, BitSet rowChooser)
	{
		mBackingStore = data;
		mColumnChooser = colChooser;
		mRowChooser = rowChooser;
	}

	/*
	 * @see qplugin.IViewableData#getBackingStore()
	 */
	public IViewableData getBackingStore()
	{	return mBackingStore; }

	/*
	 * @see qplugin.IViewableData#getColumnNames()
	 */
	public String[] getColumnNames()
	{
		String
			colNames[] = mBackingStore.getAllColumnNames(),
			newColNames[] = colNames;
		int
			origNumCols = colNames.length,
			curNumCols = getNumberOfColumns();
		if (curNumCols < origNumCols)
		{
			newColNames = new String[curNumCols];
			int idx = 0;
			for (int i = mColumnChooser.nextSetBit(0); i >= 0; i = mColumnChooser.nextSetBit(i+1), idx++)
				newColNames[idx] = colNames[i];
		}
			
		return newColNames;
	}

	public String[] getRowNames()
	{
		String
			rowNames[] = mBackingStore.getAllRowNames(),
			newRowNames[] = rowNames;
		int
			origNumRows = rowNames.length,
			curNumRows = getNumberOfRows();
		if (curNumRows < origNumRows)
		{
			newRowNames = new String[curNumRows];
			int idx = 0;
			for (int i = mRowChooser.nextSetBit(0); i >= 0; i = mRowChooser.nextSetBit(i+1), idx++)
				newRowNames[idx] = rowNames[i];
		}
			
		return newRowNames;
	}

	public double[][] extractData()
	{
		int
			curNumRows = getNumberOfRows(),
			curNumCols = getNumberOfColumns(),
			origNumCols = mBackingStore.getNumberOfColumns(),
			x, y = 0;
		double[][]
			allData = getUnfilteredData(),
			result = new double[curNumRows][];
		if (origNumCols == curNumCols) // Dispense with column filtering
		{
			for (int i = mRowChooser.nextSetBit(0); i >= 0; i = mRowChooser.nextSetBit(i+1), y++)
				result[y] = allData[i];
		}
		else
		{
			for (int i = mRowChooser.nextSetBit(0); i >= 0; i = mRowChooser.nextSetBit(i+1), y++)
			{
				double[]
					di = allData[i],
					ri = result[y] = new double[curNumCols];
				x = 0;
				for (int j = mColumnChooser.nextSetBit(0); j >= 0; j = mColumnChooser.nextSetBit(j+1), x++)
					ri[x] = di[j];
			}
		}
		return result;
	}

	public double[] getColumnWeights()
	{
		double[]
			allWeights = mBackingStore.getAllColumnWeights(),
			result = allWeights;
		int
			origNumCols = allWeights.length,
			curNumCols = getNumberOfColumns();
		if (curNumCols < origNumCols)
		{
			result = new double[curNumCols];
			int idx = 0;
			for (int i = mColumnChooser.nextSetBit(0); i >= 0; i = mColumnChooser.nextSetBit(i+1), idx++)
				result[idx] = allWeights[i];
		}
		return result;
	}

	public double[] getRowWeights()
	{
		double[]
			allWeights = mBackingStore.getAllRowWeights(),
			result = allWeights;
		int
			origNumRows = allWeights.length,
			curNumRows = getNumberOfRows();
		if (curNumRows < origNumRows)
		{
			result = new double[curNumRows];
			int idx = 0;
			for (int i = mRowChooser.nextSetBit(0); i >= 0; i = mRowChooser.nextSetBit(i+1), idx++)
				result[idx] = allWeights[i];
		}
		return result;
	}

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getAllColumnNames()
	 */
	public String[] getAllColumnNames()
	{	return mBackingStore.getAllColumnNames(); }

	public String[] getAllRowNames()
	{	return mBackingStore.getAllRowNames(); }

	public double[][] getUnfilteredData()
	{	return mBackingStore.getUnfilteredData(); }

	public double[] getAllColumnWeights()
	{	return mBackingStore.getAllColumnWeights(); }

	public double[] getAllRowWeights()
	{	return mBackingStore.getAllRowWeights(); }

	public String[][] getColumnCategories()
	{
		String
			colCats[][] = mBackingStore.getAllColumnCategories(),
			newColCats[][] = colCats;
		int
			origNumCols = colCats.length,
			curNumCols = getNumberOfColumns();
		if (curNumCols < origNumCols)
		{
			newColCats = new String[curNumCols][];
			int idx = 0;
			for (int i = mColumnChooser.nextSetBit(0); i >= 0; i = mColumnChooser.nextSetBit(i+1), idx++)
				newColCats[idx] = colCats[i];
		}
		return newColCats;
	}

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getAllColumnCategories()
	 */
	public String[][] getAllColumnCategories()
	{	return mBackingStore.getAllColumnCategories(); }

	/*
	 * @see qplugin.IViewableData#getNumberOfColumns()
	 */
	public int getNumberOfColumns()
	{	return mColumnChooser.cardinality(); }

	/*
	 * @see qplugin.IViewableData#getNumberOfRows()
	 */
	public int getNumberOfRows()
	{	return mRowChooser.cardinality(); }

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getAllDataIterator()
	 */
	public IDataIterator<double[]> getAllDataIterator()
	{	return mBackingStore.getAllDataIterator(); }

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getDataIterator()
	 */
	public IDataIterator<double[]> getDataIterator()
	{	return new QTableViewIterator(); }

	/* (non-Javadoc)
	 * @see java.lang.Iterable#iterator()
	 */
	public Iterator<double[]> iterator()
	{	return getDataIterator(); }

	public String getName()
	{	return mBackingStore.getName(); }

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#createView(qplugin.IDataFilter)
	 */
	public IViewableData createView(IDataFilter filter)
	{
		Object colNames[] = getAllColumnNames();
		String colCats[][] = getAllColumnCategories();
		BitSet
			colChooser = mColumnChooser,
			rowChooser = new BitSet(mRowChooser.cardinality());
		if (!filter.isAllColumnsSelected())
		{
			colChooser = new BitSet(mColumnChooser.cardinality());
			for (int i = mColumnChooser.nextSetBit(0); i >= 0; i = mColumnChooser.nextSetBit(i+1))
			{
				if (filter.isColumnSelected(i, colNames != null ? colNames[i] : null, colCats != null ? colCats[i] : null))
					colChooser.set(i);
			}
		}

		int
			origNumCols = colNames.length,
			curNumCols = getNumberOfColumns(),
			idx = 0;
		if (origNumCols == curNumCols)
		{
			for (double[] data: this)
			{
				if (filter.isRowSelected(idx, data, colNames))
					rowChooser.set(idx);
				idx++;
			}
		}
		else
		{
			for (double[] data: this)
			{
				double newData[] = new double[curNumCols];
				int i = 0;
				for (int j = mColumnChooser.nextSetBit(0); j >= 0; j = mColumnChooser.nextSetBit(j+1), i++)
					newData[i] = data[j];
				if (filter.isRowSelected(idx, newData, colNames))
					rowChooser.set(idx);
				idx++;
			}
		}
		return new QTableView(this, colChooser, rowChooser);
	}

	public IViewableData createView(int columnNo, double value)
	{
		Object colNames[] = getAllColumnNames();
		BitSet
			colChooser = mColumnChooser,
			rowChooser = new BitSet(mRowChooser.cardinality());
	
		int
			origNumCols = colNames.length,
			curNumCols = getNumberOfColumns(),
			idx = 0;
		if (origNumCols == curNumCols)
		{
			for (double[] data: this)
			{
				if (data[columnNo] == value)
					rowChooser.set(idx);
				idx++;
			}
		}
		else
		{
			for (double[] data: this)
			{
				double newData[] = new double[curNumCols];
				int i = 0;
				for (int j = mColumnChooser.nextSetBit(0); j >= 0; j = mColumnChooser.nextSetBit(j+1), i++)
					newData[i] = data[j];
				if (newData[columnNo] == value)
					rowChooser.set(idx);
				idx++;
			}
		}
		return new QTableView(this, colChooser, rowChooser);
	}

	public int count(int columnNo, double value)
	{
		int
			count = 0,
			origNumCols = getAllColumnNames().length,
			curNumCols = getNumberOfColumns();
		if (origNumCols == curNumCols)
		{
			for (double[] data: this)
				if (data[columnNo] == value)
					count++;
		}
		else
		{
			for (double[] data: this)
			{
				double newData[] = new double[curNumCols];
				int i = 0;
				for (int j = mColumnChooser.nextSetBit(0); j >= 0; j = mColumnChooser.nextSetBit(j+1), i++)
					newData[i] = data[j];
				if (data[columnNo] == value)
					count++;
			}
		}
		return count;
	}

	/**
	 * In view, column and row weights are NOT cloned! Row names are also NOT cloned!
	 */
	public QTableData cloneData()
	{
		int
			numRows = getNumberOfRows(),
			numCols = getNumberOfColumns(),
			origNumCols = getAllColumnNames().length,
			rowIdx = 0,
			colIdx;
		double[][] data = new double[numRows][];

		if (origNumCols == numCols)
		{
			for (double[] rowData: this)
				data[rowIdx++] = rowData;
		}
		else
		{
			for (double[] rowData: this)
			{
				double newData[] = data[rowIdx] = new double[numCols];
				colIdx = 0;
				data[rowIdx++] = rowData;
				for (int j = mColumnChooser.nextSetBit(0); j >= 0; j = mColumnChooser.nextSetBit(j+1), colIdx++)
					newData[colIdx] = rowData[j];
			}
		}

		QTableData newTable = new QTableData(data, null, getColumnNames());
		newTable.setColumnCategories(getColumnCategories());
		return newTable;
	}

	public void purge()
	{
		mBackingStore = null;
		mColumnChooser = mRowChooser = null;
		mFilenames = null;
	}

	/**
	 * Set the name of this data
	 */
	public void setName(String name)
	{	mBackingStore.setName(name); }

	/**
	 * Set the file names associated with this data.
	 */
	public void setFilenames(String[] names)
	{	mFilenames = names; }

	/**
	 * Returns the file names associated with this data
	 */
	public String[] getFilenames()
	{	return mFilenames; }
}
