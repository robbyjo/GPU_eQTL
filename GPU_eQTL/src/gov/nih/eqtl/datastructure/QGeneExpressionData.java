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
package gov.nih.eqtl.datastructure;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

import com.csvreader.CsvReader;

import gov.nih.table.IDataIterator;
import gov.nih.table.IViewableData;
import gov.nih.utils.QDataUtils;
import gov.nih.utils.QStringUtils;
import static gov.nih.utils.QDataUtils.kUndefinedValue;
import static gov.nih.utils.QFileUtils.sLineFormat;

/**
 * Table-based data
 * @author Roby Joehanes
 *
 */
public class QGeneExpressionData implements IViewableData
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

	protected String[] mGeneIDs;
	protected static int mCounter = 0;
	protected String mName = "Default" + (mCounter++); //$NON-NLS-1$
	protected double[][] mData; // row x columns
	protected double[]
		mColumnWeights,
		mRowWeights;
	protected String[] mColNames;
	protected String[][] mColumnCategories;

	protected String[] mFilenames;

	public QGeneExpressionData() {}

	public QGeneExpressionData(double[][] data)
	{	this(data, null, null, null, null); }

	public QGeneExpressionData(double[][] data, String[] geneIDs, String[] colNames)
	{	this(data, geneIDs, colNames, null, null); }

	public QGeneExpressionData(double[][] data, String[] geneIDs, String[] colNames, double[] rowWeights, double[] colWeights)
	{
		mData = data;
		mGeneIDs = geneIDs;
		mColNames = colNames;
		mColumnCategories = new String[mData[0].length][];
		mRowWeights = rowWeights;
		mColumnWeights = colWeights;
	}

	public void purge()
	{
		mData = null;
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
	{	return null; }

	public String[] getAllRowNames()
	{	return null; }

	public void setRowNames(String[] s)
	{	}

	public String[] getColumnNames()
	{	return mColNames; }

	public void setColumnNames(String[] s)
	{	mColNames = s; }

	public String[] getGeneIDs()
	{	return mGeneIDs; }

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
	public QGeneExpressionData cloneData()
	{
		QGeneExpressionData newData = new QGeneExpressionData(QDataUtils.safeClone(mData), mGeneIDs, mColNames);
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
		if (mGeneIDs != null)
		{
			String temp = mGeneIDs[row1];
			mGeneIDs[row1] = mGeneIDs[row2];
			mGeneIDs[row2] = temp;
		}
		if (mRowWeights != null)
		{
			double temp = mRowWeights[row1];
			mRowWeights[row1] = mRowWeights[row2];
			mRowWeights[row2] = temp;
		}
	}

	public void saveAsBinary(String filename) throws IOException
	{
		DataOutputStream os = new DataOutputStream(new FileOutputStream(filename));
		int
			nrow = mData.length,
			ncol = mData[0].length;
		os.writeInt(nrow);
		os.writeInt(0); // For future expansion
		os.writeInt(ncol);
		os.writeInt(0); // For future expansion
		os.writeInt(mGeneIDs == null ? 0 : 1);
		os.writeInt(mColNames == null ? 0 : 1);
		os.writeInt(0); // For future expansion
		os.writeInt(0); // For future expansion
		if (mGeneIDs != null)
		{
			for (int i = 0; i < nrow; i++) {
				byte[] bytes = mGeneIDs[i].getBytes();
				os.writeInt(bytes.length);
				os.write(bytes);
			}
		}
		if (mColNames != null)
		{
			for (int i = 0; i < ncol; i++)
			{
				byte[] bytes = mColNames[i].getBytes();
				os.writeInt(bytes.length);
				os.write(bytes);
			}
		}
		for (int i = 0; i < nrow; i++)
			for (int j = 0; j < ncol; j++)
				os.writeDouble(mData[i][j]);
		os.flush();
		os.close();
	}

	public static final QGeneExpressionData loadAsBinary(String filename) throws IOException
	{
		DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(filename), 16 * 1024 * 1024));
		String[] geneIDs = null;
		String[] colNames = null;
		double[][] data = null;
		try {
			int nrow = is.readInt();
			is.readInt();
			int ncol = is.readInt();
			is.readInt();
			boolean hasGeneID = is.readInt() != 0;
			boolean hasColumnName = is.readInt() != 0;
			is.readInt();
			is.readInt();
			data = new double[nrow][ncol];

			if (hasGeneID)
			{
				geneIDs = new String[nrow];
				for (int i = 0; i < nrow; i++) {
					int len = is.readInt();
					byte[] bytes = new byte[len];
					is.read(bytes);
					geneIDs[i] = new String(bytes);
				}
			}
			if (hasColumnName)
			{
				colNames = new String[ncol];
				for (int i = 0; i < ncol; i++)
				{
					int len = is.readInt();
					byte[] bytes = new byte[len];
					is.read(bytes);
					colNames[i] = new String(bytes);
				}
			}
			for (int i = 0; i < nrow; i++)
				for (int j = 0; j < ncol; j++)
					data[i][j] = is.readDouble();
		} finally {
			is.close();
		}
		return new QGeneExpressionData(data, geneIDs, colNames);
	}

	/**
	 * The same as QFileUtils.readDelimitedFileAsTableData
	 * @param filename
	 * @param delimiters
	 * @param commentMarkers
	 * @param hasRowNames
	 * @param hasColNames
	 * @return
	 */
	public static final QGeneExpressionData load(String filename, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return load(new FileReader(filename), delimiters, commentMarkers, hasRowNames, hasColNames); }

	public static final QGeneExpressionData loadAsBZip2(String filename, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return load(new InputStreamReader(new BZip2CompressorInputStream(new BufferedInputStream(new FileInputStream(filename)))), delimiters, commentMarkers, hasRowNames, hasColNames); }

	public static final QGeneExpressionData loadAsGZip(String filename, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return load(new InputStreamReader(new GZIPInputStream(new BufferedInputStream(new FileInputStream(filename)))), delimiters, commentMarkers, hasRowNames, hasColNames); }

	public static final QGeneExpressionData load(File file, char[] delimiters, String commentMarkers,
		boolean hasRowNames, boolean hasColNames) throws IOException
	{	return load(new FileReader(file), delimiters, commentMarkers, hasRowNames, hasColNames); }

	public static final QGeneExpressionData load(Reader r, char[] delim, String comments, boolean hasRowNames, boolean hasColNames) throws IOException
	{
		CsvReader reader = null;
		List<double[]> list = new ArrayList<double[]>();
		List<String> rowNames = new ArrayList<String>();
		int lineNo = 0, numTokens = 0;
		String[] tokens = null, colNames = null;
		try
		{
			reader = new CsvReader(r);
			reader.setTrimWhitespace(true);
			reader.setDelimiter(delim);
			if (comments != null) {
				reader.setUseComments(true);
				reader.setComment(comments.charAt(0));
			}
			if (hasColNames) {
				lineNo++;
				if (!reader.readRecord())
					return null;
				colNames = reader.getValues();
				numTokens = colNames.length;
				if (hasRowNames) {
					String[] s = new String[colNames.length - 1];
					System.arraycopy(colNames, 1, s, 0, s.length);
					colNames = s;
				}
			}

			while(reader.readRecord()) {
				tokens = reader.getValues();
				lineNo++;
				int curNumTokens = tokens.length;
				if (numTokens != curNumTokens)
					throw new RuntimeException(String.format(sLineFormat, curNumTokens, numTokens, lineNo));

				double[] curData = null;
				String tok;
				if (hasRowNames)
				{
					curData = new double[numTokens - 1];
					tok = tokens[0].trim();
					rowNames.add(tok);
					for (int i = 1; i < curNumTokens; i++)
					{
						tok = tokens[i].trim();
						curData[i - 1] = "".equals(tok) ? kUndefinedValue : Double.parseDouble(tok);
					}
				}
				else
				{
					curData = new double[numTokens];
					for (int i = 0; i < curNumTokens; i++)
					{
						tok = tokens[i].trim();
						curData[i] = "".equals(tok) ? kUndefinedValue : Double.parseDouble(tok);
					}
				}
				list.add(curData);
			}
		}
		finally
		{
			if (reader != null)
				reader.close();
		}
		int size = list.size();
		if (size == 0)
			return null;
		return new QGeneExpressionData(list.toArray(new double[size][]),
			hasRowNames ? rowNames.toArray(new String[rowNames.size()]) : null,
			hasColNames ? colNames : null);
	}

	public String toDelimitedString(String delimiter, boolean includeRowName, boolean includeColName)
	{	return QStringUtils.toDelimitedString(mData, includeRowName ? mGeneIDs: null, includeColName? mColNames: null, delimiter); }
}
