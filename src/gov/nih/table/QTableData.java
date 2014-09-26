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

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import gov.nih.utils.QDataUtils;
import gov.nih.utils.QFileUtils;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.matrix.QMatrixUtils;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableData;

/**
 * Table-based data
 * @author Roby Joehanes
 *
 */
public class QTableData implements IViewableData
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
	protected String mName = "Default" + (mCounter++); // Title / name for this data
	protected double[][] mData; // row x columns
	protected double[]
		mColumnWeights,
		mRowWeights;
	protected String[]
		mRowNames,
		mColNames;
	protected String[][] mColumnCategories;

	protected String[] mFilenames;

	public QTableData() {}

	public QTableData(double[][] data)
	{	this(data, null, null, null, null); }

	public QTableData(double[][] data, String[] rowNames, String[] colNames)
	{	this(data, rowNames, colNames, null, null); }

	public QTableData(double[][] data, String[] rowNames, String[] colNames, double[] rowWeights, double[] colWeights)
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

	public String[] getRowNames()
	{	return mRowNames; }

	public String[] getAllRowNames()
	{	return mRowNames; }

	public void setRowNames(String[] s)
	{	mRowNames = s; }

	public String[] getColumnNames()
	{	return mColNames; }

	public void setColumnNames(String[] s)
	{	mColNames = s; }

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#getAllColumnNames()
	 */
	public String[] getAllColumnNames()
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

	/**
	 * Insert a column at the last position
	 * @param columnName
	 * @param columnCats
	 * @param columnData
	 */
	public void insertColumn(String columnName, String[] columnCats, double[] columnData)
	{
		int
			numRows = columnData.length,
			numCols = mData[0].length;
		if (numRows != mData.length)
			throw new RuntimeException(); // Row mismatch
		String[] newColNames = null;
		if (mColNames != null)
		{
			newColNames = new String[numCols + 1];
			System.arraycopy(mColNames, 0, newColNames, 0, numCols);
			newColNames[numCols] = columnName;
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

	public double[] getColumnData(int columnNo)
	{
		int numRows = mData.length;
		double[] data = new double[numRows];
		for (int i = 0; i < numRows; i++)
			data[i] = mData[i][columnNo];
		return data;
	}

	public String[] getFilenames()
	{	return mFilenames; }

	public void setFilenames(String[] filenames)
	{	mFilenames = filenames; }

	public IViewableData getBackingStore()
	{	return this; }

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
	public QTableData cloneData()
	{
		QTableData newData = new QTableData(QDataUtils.safeClone(mData), mRowNames, mColNames);
		newData.mColumnCategories = mColumnCategories;
		newData.mColumnWeights = mColumnWeights != null ? mColumnWeights.clone() : null;
		newData.mRowWeights = mRowWeights != null ? mRowWeights.clone() : null;
		return newData;
	}

	/* (non-Javadoc)
	 * @see qplugin.IViewableData#createView(qplugin.IDataFilter)
	 */
	public IViewableData createView(IDataFilter filter)
	{
		int
			numCols = mData[0].length,
			numRows = mData.length;
		BitSet
			colChooser = new BitSet(numCols),
			rowChooser = new BitSet(numRows);
		if (!filter.isAllColumnsSelected())
		{
			for (int i = 0; i < numCols; i++)
			{
				if (filter.isColumnSelected(i, mColNames != null ? mColNames[i] : null, mColumnCategories != null ? mColumnCategories[i] : null))
					colChooser.set(i);
			}
		}
		else
			colChooser.set(0, numCols);
	
		for (int i = 0; i < numRows; i++)
		{
			if (filter.isRowSelected(i, mData[i], mColNames))
				rowChooser.set(i);
		}
		return new QTableView(this, colChooser, rowChooser);
	}

	public IViewableData createView(int columnNo, double value)
	{
		int
			numCols = mData[0].length,
			numRows = mData.length;
		BitSet
			colChooser = new BitSet(numCols),
			rowChooser = new BitSet(numRows);
		colChooser.set(0, numCols);
	
		for (int i = 0; i < numRows; i++)
		{
			if (mData[i][columnNo] == value)
				rowChooser.set(i);
		}
		return new QTableView(this, colChooser, rowChooser);
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

	public void swapColumns(int col1, int col2)
	{
		int numRows = mData.length;
		for (int i = 0; i < numRows; i++)
		{
			double
				curRow[] = mData[i],
				temp = curRow[col1];
			curRow[col1] = curRow[col2];
			curRow[col2] = temp;
		}
		if (mColNames != null)
		{
			String temp = mColNames[col1];
			mColNames[col1] = mColNames[col2];
			mColNames[col2] = temp;
		}
		if (mColumnCategories != null)
		{
			String[] temp = mColumnCategories[col1];
			mColumnCategories[col1] = mColumnCategories[col2];
			mColumnCategories[col2] = temp;
		}
		if (mColumnWeights != null)
		{
			double temp = mColumnWeights[col1];
			mColumnWeights[col1] = mColumnWeights[col2];
			mColumnWeights[col2] = temp;
		}
	}

	public void swapRows(int row1, int row2)
	{
		{
			double[] temp = mData[row1];
			mData[row1] = mData[row2];
			mData[row2] = temp;
		}
		if (mRowNames != null)
		{
			String temp = mRowNames[row1];
			mRowNames[row1] = mRowNames[row2];
			mRowNames[row2] = temp;
		}
		if (mRowWeights != null)
		{
			double temp = mRowWeights[row1];
			mRowWeights[row1] = mRowWeights[row2];
			mRowWeights[row2] = temp;
		}
	}

	public static final QTableData load(String filename, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return readDelimitedFileAsTableData(filename, delimiters, commentMarkers, hasRowNames, hasColNames, true); }

	public static final QTableData load(File file, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return readDelimitedFileAsTableData(file, delimiters, commentMarkers, hasRowNames, hasColNames, true); }

	public static final QTableData load(Reader r, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return readDelimitedFileAsTableData(r, delimiters, commentMarkers, hasRowNames, hasColNames, true); }

	/**
	 * The same as QFileUtils.readDelimitedFileAsTableData
	 * @param filename
	 * @param delimiters
	 * @param commentMarkers
	 * @param hasRowNames
	 * @param hasColNames
	 * @param isAllNumeric
	 * @return
	 */
	public static final QTableData load(String filename, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames, boolean isAllNumeric) throws IOException
	{	return readDelimitedFileAsTableData(filename, delimiters, commentMarkers, hasRowNames, hasColNames, isAllNumeric); }

	public static final QTableData load(File file, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames, boolean isAllNumeric) throws IOException
	{	return readDelimitedFileAsTableData(file, delimiters, commentMarkers, hasRowNames, hasColNames, isAllNumeric); }

	public static final QTableData load(Reader r, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames, boolean isAllNumeric) throws IOException
	{	return readDelimitedFileAsTableData(r, delimiters, commentMarkers, hasRowNames, hasColNames, isAllNumeric); }

	public String toDelimitedString(String delimiter, boolean includeRowName, boolean includeColName)
	{	return QStringUtils.toDelimitedString(mData, includeRowName ? mRowNames: null, includeColName? mColNames: null, delimiter); }

	private static final double[][] buildModelMatrix(double[] curData, String[] categories)
	{
		int
			n = curData.length,
			ncats;
		if (categories == null) {
			Set<Double> cats = new HashSet<Double>();
			for (double dbl: curData)
				if (dbl != QFileUtils.kUndefinedValue)
					cats.add(dbl);
			List<Double> catsList = new ArrayList<Double>();
			catsList.addAll(cats);
			double[] newData = new double[n];
			for (int i = 0; i < n; i++)
				newData[i] = catsList.indexOf(curData[i]);
			curData = newData;
			ncats = catsList.size();
			categories = new String[ncats];
			for (int i = 0; i < ncats; i++)
				categories[i] = String.valueOf(catsList.get(i));
		}
		// Always omit the first category (i.e., index 0)
		ncats = categories.length - 1;
		double[][] newData = new double[ncats][n];
		for (int i = 0; i < n; i++) {
			if (curData[i] == QFileUtils.kUndefinedValue)
				continue;
			int v = (int) (curData[i] - 1);
			if (v >= 0)
				newData[v][i] = 1.0;
		}
		return newData;
	}

	public double[][] buildModelMatrix(String[] covars, String[] factorCovars, boolean withIntercept, BitSet missing)
	{
		String[] colNames = getAllColumnNames();
		int
			numCols = colNames.length,
			numRows = getNumberOfRows();
		Set<String> factors = new HashSet<String>();
		Map<String, Integer> idxMap = new HashMap<String, Integer>();
		Map<String, double[][]> idxToDataMap = new HashMap<String, double[][]>();
		if (missing == null)
			missing = new BitSet(numRows);
		double[][] covarMatrix = null;
		if (withIntercept) {
			double[] intercept = new double[numRows];
			Arrays.fill(intercept, 1.0);
			covarMatrix = new double[][] { intercept };
		} else
			covarMatrix = new double[][] {};

		if (factorCovars != null)
			for (String str: factorCovars)
				factors.add(str);
		for (int i = 0; i < numCols; i++)
			idxMap.put(colNames[i], i);

		for (String str : covars) {
			double[][] mtx = null;
			if (str.contains("*")) {
				for (String sstr: str.split("\\*")) {
					double[][] curMtx = idxToDataMap.get(sstr);
					if (curMtx == null) {
						int colIdx = idxMap.get(sstr);
						double[] curData = getColumnData(colIdx);
						for (int i = 0; i < numRows; i++)
							if (curData[i] == QFileUtils.kUndefinedValue)
								missing.set(i);
						String[] colCats = getColumnCategory(colIdx);
						if (factors.contains(sstr) || colCats != null)
							idxToDataMap.put(sstr, buildModelMatrix(curData, colCats));
						else
							idxToDataMap.put(sstr, new double[][] {curData});
					}
					// The first iteration: Simply store in mtx
					if (mtx == null) {
						mtx = curMtx;
						continue;
					}
					// Else, perform cross-product
					double[][] newMtx = new double[mtx.length * curMtx.length][numRows];
					for (int i = 0; i < mtx.length; i++)
						for (int j = 0; j < curMtx.length; j++)
							for (int k = 0; k < numRows; k++)
								newMtx[i*j][k] = mtx[i][k] * curMtx[j][k];
					mtx = newMtx;
				}
			} else {
				mtx = idxToDataMap.get(str);
				if (mtx == null) {
					int colIdx = idxMap.get(str);
					double[] curData = getColumnData(colIdx);
					for (int i = 0; i < numRows; i++)
						if (curData[i] == QFileUtils.kUndefinedValue)
							missing.set(i);
					String[] colCats = getColumnCategory(colIdx);
					if (factors.contains(str) || colCats != null)
						mtx = buildModelMatrix(curData, colCats);
					else
						mtx = new double[][] {curData};
					idxToDataMap.put(str, mtx);
				}
			}
			// Insert the newly formed matrix for that factor into the overall design matrix
			int i = covarMatrix.length;
			double[][] newCovarMatrix = new double[i + mtx.length][];
			System.arraycopy(covarMatrix, 0, newCovarMatrix, 0, i);
			System.arraycopy(mtx, 0, newCovarMatrix, i, mtx.length);
			covarMatrix = newCovarMatrix;
		}
		return QMatrixUtils.transpose(covarMatrix);
	}
}
