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
package gov.nih.utils;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;

import com.csvreader.CsvReader;

import gov.nih.table.QTableData;
import static gov.nih.utils.QStringUtils.*;

/**
 * File utilities, salvaged from the old QGene
 * 
 * @author Roby Joehanes
 */
public class QFileUtils
{
	public static final int
		kDefaultXMLIndent = 2,
		kDefaultBufferSize = 4 * 1024 * 1024,
		kUndefinedValue = Integer.MIN_VALUE;

	public static final String
		// General file saving-related error messages / warnings
		sErrNothingToSave = "Error.IO.NothingToSave", //$NON-NLS-1$
		sUnexpectedEOFErr = "Error.IO.Unexpected_EOF", //$NON-NLS-1$
		sFileMissingErr = "Error.IO.File_Missing", //$NON-NLS-1$
		sCantWriteErr = "Error.IO.Cant_Write", //$NON-NLS-1$
		sCantReadErr = "Error.IO.Cant_Read", //$NON-NLS-1$
		sConfirmOverwrite = "Warning.AreYouSureOverwrite", //$NON-NLS-1$

		// Text format
		sTextFormat = "Format.Text", //$NON-NLS-1$
		sTextFormatSuffix = ".txt", //$NON-NLS-1$
		sTextFormatExtension = "*" + sTextFormatSuffix, //$NON-NLS-1$

		sWordDelimiter = ", ", //$NON-NLS-1$
		sDefaultDelim = " \t\n\r\f,", //$NON-NLS-1$
		sLineFormat = "%d vs %d, at line %d", //$NON-NLS-1$
		sQuotes = "\"", // $NON-NLS-1$
		sCommentTokens[] = new String[] { "#", "//" }; //$NON-NLS-1$ //$NON-NLS-2$

	/**
	 * Write an XML Document object to a file name with default indenting (2 spaces)
	 * @param doc The XML document node
	 * @param filename The file name it's going to be saved to
	 */
	public static final void writeToXMLFile(Document doc, String filename) throws IOException
	{	writeToXMLFile(doc, filename,  kDefaultXMLIndent); }

	/**
	 * Write an XML Document object to a writer with default indenting (2 spaces)
	 * @param doc The XML document node
	 * @param w The writer / stream it's going to be saved to
	 * @param ident the number of spaces to ident the XML output
	 */
	public static final void writeToXMLFile(Document doc, Writer w) throws IOException
	{	writeToXMLFile(doc, w, kDefaultXMLIndent); }

	/**
	 * Write an XML Document object to a file name
	 * @param doc The XML document node
	 * @param filename The file name it's going to be saved to
	 * @param ident the number of spaces to ident the XML output
	 */
	public static final void writeToXMLFile(Document doc, String filename, int ident) throws IOException
	{
		// if filename is null, then output to screen instead -- RJ 2004/10/21
		PrintWriter fw = filename != null ? new PrintWriter(filename) :
			new PrintWriter(System.out);
		writeToXMLFile(doc, fw, ident);
		if (filename != null)
			fw.close();
	}

	// Refactored from QDSNRecordParser -- RJ 2004/10/21
	/**
	 * Write an XML Document object to a writer
	 * @param doc The XML document node
	 * @param w The writer / stream it's going to be saved to
	 * @param ident the number of spaces to ident the XML output
	 */
	public static final void writeToXMLFile(Document doc, Writer w, int ident) throws IOException
	{
		DOMSource domSource = new DOMSource(doc.getDocumentElement());
		StringWriter outputstring = new StringWriter();		 
		StreamResult streamResult = new StreamResult(outputstring);
		Transformer xmlserializer = null;
		try
		{	xmlserializer = TransformerFactory.newInstance().newTransformer(); }
		catch (Exception e)
		{} // Should never happen
		xmlserializer.setOutputProperty(OutputKeys.ENCODING,"ISO-8859-1"); //$NON-NLS-1$
		xmlserializer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM,"users.dtd"); //$NON-NLS-1$
		xmlserializer.setOutputProperty(OutputKeys.INDENT,"yes"); //$NON-NLS-1$
		try
		{	xmlserializer.transform(domSource, streamResult); }
		catch (Exception e)
		{	throw new IOException(e.getLocalizedMessage()); }
		QFileUtils.writeText(outputstring.toString(), w); 
	}

	/**
	 * Generic method to write a given text string to a file whose path is specified
	 */
	public static final void writeText(String textOutput, String filename) throws IOException
	{
		writeText(textOutput, filename, false);
	}

	/**
	 * Generic method to write a given text string to a file whose path is specified
	 * @param textOutput what to write
	 * @param filename file name. If null, it will be printed on screen and append flag is ignored
	 * @param append whether to append to the existing file (if any) or not
	 */
	public static final void writeText(String textOutput, String filename, boolean append) throws IOException
	{
		// if filename is null, then output to screen instead -- RJ 2004/10/21
		Writer fw = filename != null ? new PrintWriter(new FileOutputStream(filename, append), true) :
			new PrintWriter(System.out);
		writeText(textOutput, fw);
		if (filename != null)
			fw.close();
	}

	/**
	 * Generic method to write a given text string to a writer stream
	 */
	public static final void writeText(String textOutput, Writer fw) throws IOException
	{
		fw.write(textOutput);
		fw.flush();
	}

	/**
	 * Same as writeText(data, filename, false)
	 * @param data
	 * @param filename
	 * @throws IOException
	 */
	public static final void writeText(double[][] data, String filename) throws IOException
	{	writeText(data, filename, false); }

	/**
	 * Convert data to string, then write it to file name.
	 * @param data
	 * @param filename
	 * @param append Set to true if you want to append to the existing file. If file is not found,
	 * then a new file is created
	 * @throws IOException
	 */
	public static final void writeText(double[][] data, String filename, boolean append) throws IOException
	{
		// if filename is null, then output to screen instead -- RJ 2004/10/21
		Writer fw = filename != null ? new PrintWriter(new FileOutputStream(filename, append), true) :
			new PrintWriter(System.out);
		writeText(data, fw);
		if (filename != null)
			fw.close();
	}

	/**
	 * Convert data to string, then write it.
	 * @param data
	 * @param fw
	 * @throws IOException
	 */
	public static final void writeText(double[][] data, Writer fw) throws IOException
	{
		for (double[] dataRow: data)
		{
			String str = toString(dataRow) + sLn;
			fw.write(str);
		}
	}

	/**
	 * Same as writeText(data, filename, false)
	 * @param data
	 * @param filename
	 * @throws IOException
	 */
	public static final void writeText(float[][] data, String filename) throws IOException
	{	writeText(data, filename, false); }

	/**
	 * Convert data to string, then write it to file name.
	 * @param data
	 * @param filename
	 * @param append Set to true if you want to append to the existing file. If file is not found,
	 * then a new file is created
	 * @throws IOException
	 */
	public static final void writeText(float[][] data, String filename, boolean append) throws IOException
	{
		// if filename is null, then output to screen instead -- RJ 2004/10/21
		Writer fw = filename != null ? new PrintWriter(new FileOutputStream(filename, append), true) :
			new PrintWriter(System.out);
		writeText(data, fw);
		if (filename != null)
			fw.close();
	}

	public static final void writeText(float[][] data, Writer fw) throws IOException
	{
		for (float[] dataRow: data)
		{
			String str = toString(dataRow) + sLn;
			fw.write(str);
		}
	}

	public static final String toString(float[] values)
	{	return toString(values, sWordDelimiter); }

	public static final String toString(float[] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		if (numValues > 0)
		{
			buffer.append(values[0]);
			for (int i = 1; i < numValues; i++)
				buffer.append(delimiter + values[i]);
		}
		return buffer.toString();
	}

	public static final String toString(double[] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert an array of double values into a single,
	 * comma delimited string
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(double[] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder(values.length * 4);
		int numValues = values.length;
		if (numValues > 0)
		{
			buffer.append(values[0]);
			for (int i = 1; i < numValues; i++)
			{
				buffer.append(delimiter);
				buffer.append(values[i]);
			}
		}
		return buffer.toString();
	}

	/**
	 * Read a text file and put it into a string
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static final String readTextFile(String filename) throws IOException
	{
		FileReader reader = new FileReader(filename);
		try {
			return readTextFile(reader);
		} finally {
		if (reader != null)
			reader.close();
		}
	}

	/**
	 * Read a text file and put it into a string
	 * @param f
	 * @return
	 * @throws IOException
	 */
	public static final String readTextFile(File f) throws IOException
	{
		FileReader reader = new FileReader(f);
		try {
			return readTextFile(reader);
		} finally {
		if (reader != null)
			reader.close();
		}
	}

	public static final String readTextFile(InputStream instr) throws IOException
	{	return readTextFile(new InputStreamReader(instr)); }

	/**
	 * Read a text file and put it into a string
	 * @param r
	 * @return
	 * @throws IOException
	 */
	public static final String readTextFile(Reader r) throws IOException
	{
		StringBuilder buf = new StringBuilder(); // new class from Java 1.5
		LineNumberReader reader = null;
		try
		{
			reader = new LineNumberReader(r);
			String s = null;
			while ((s = reader.readLine()) != null)
				buf.append(s + sLn);
		}
		finally
		{
		}
		return buf.toString();
	}

	public static final String[] readTextFileAsArray(String filename, String[] comments) throws IOException
	{	return readTextFileAsArray(new BufferedReader(new FileReader(filename), kDefaultBufferSize), comments); }

	public static final String[] readTextFileAsArray(File file, String[] comments) throws IOException
	{	return readTextFileAsArray(new BufferedReader(new FileReader(file), kDefaultBufferSize), comments); }

	/**
	 * Read a text file as an array of Strings
	 * @param r
	 * @param comments
	 * @return
	 * @throws IOException
	 */
	public static final String[] readTextFileAsArray(Reader r, String[] comments) throws IOException
	{
		List<String> list = readTextFileAsList(r, comments);
		int size = list.size();
		return size == 0 ? null : list.toArray(new String[size]);
	}

	public static final List<String> readTextFileAsList(String filename, String[] comments) throws IOException
	{	return readTextFileAsList(new BufferedReader(new FileReader(filename), kDefaultBufferSize), comments); }

	public static final List<String> readTextFileAsList(File file, String[] comments) throws IOException
	{	return readTextFileAsList(new BufferedReader(new FileReader(file), kDefaultBufferSize), comments); }

	public static final List<String> readTextFileAsList(Reader r, String[] comments) throws IOException
	{
		LineNumberReader reader = null;
		List<String> list = new ArrayList<String>();
		try
		{
			reader = new LineNumberReader(r);
			String s = null;
			do {
				s = reader.readLine();
				if (s == null)
					break;
				if (comments != null)
				{
					String temp = s.trim();
					boolean isComments = false;
					if (temp.length() > 0)
					{
						for (int comNo = 0; comNo < comments.length; comNo++)
							if (temp.startsWith(comments[comNo]))
							{
								isComments = true;
								break;
							}
					}
					else
						isComments = true;
					if (isComments)
						continue;
				}
				list.add(s);
			} while (true);
		}
		finally
		{
			if (reader != null)
				reader.close();
		}
		return list;
	}

	public static final String[][] readDelimitedFile(String fn, char[] delim, String comments) throws IOException
	{	return readDelimitedFile(new BufferedReader(new FileReader(fn), kDefaultBufferSize), delim, comments); }

	public static final String[][] readDelimitedFile(File f, char[] delim, String comments) throws IOException
	{	return readDelimitedFile(new BufferedReader(new FileReader(f), kDefaultBufferSize), delim, comments); }

	/**
	 * Read a delimited file (tab or comma or whatever)
	 * @param r
	 * @param delim The delimiter
	 * @param comments Comment tokens
	 * @return
	 * @throws IOException
	 */
	public static final String[][] readDelimitedFile(Reader r, char[] delim, String comments) throws IOException
	{
		CsvReader reader = null;
		List<String[]> list = new ArrayList<String[]>();
		try
		{
			reader = new CsvReader(r);
			reader.setTrimWhitespace(true);
			reader.setDelimiter(delim);
			if (comments != null) {
				reader.setUseComments(true);
				reader.setComment(comments.charAt(0));
			}
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				int numTokens = tokens.length;
				for (int i = 0; i < numTokens; i++)
					tokens[i] = tokens[i];
				list.add(tokens);
			}
		}
		finally
		{
			if (reader != null)
				reader.close();
		}
		int size = list.size();
		return size == 0 ? null : list.toArray(new String[size][]);
	}

	public static final QTableData readDelimitedFileAsTableData(String filename, char[] delim, String comments,
			boolean hasRowNames, boolean hasColNames) throws IOException
		{	return readDelimitedFileAsTableData(new BufferedReader(new FileReader(filename), kDefaultBufferSize), delim, comments, hasRowNames, hasColNames, true); }

		public static final QTableData readDelimitedFileAsTableData(File f, char[] delim, String comments,
			boolean hasRowNames, boolean hasColNames) throws IOException
		{	return readDelimitedFileAsTableData(new BufferedReader(new FileReader(f), kDefaultBufferSize), delim, comments, hasRowNames, hasColNames, true); }

	public static final QTableData readDelimitedFileAsTableData(String filename, char[] delim, String comments,
		boolean hasRowNames, boolean hasColNames, boolean isAllNumeric) throws IOException
	{	return readDelimitedFileAsTableData(new BufferedReader(new FileReader(filename), kDefaultBufferSize), delim, comments, hasRowNames, hasColNames, isAllNumeric); }

	public static final QTableData readDelimitedFileAsTableData(File f, char[] delim, String comments,
		boolean hasRowNames, boolean hasColNames, boolean isAllNumeric) throws IOException
	{	return readDelimitedFileAsTableData(new BufferedReader(new FileReader(f), kDefaultBufferSize), delim, comments, hasRowNames, hasColNames, isAllNumeric); }

	/**
	 * Read a delimited file (tab or comma or whatever)
	 * @param r
	 * @param delim The delimiter
	 * @param comments Comment tokens
	 * @param hasRowNames set to true if the first column is the row names
	 * @param hasColNames set to true if the first row is the column names
	 * @param isAllNumeric set to true if the data are all numeric (except for the row and column names)
	 * @return
	 * @throws IOException
	 */
	public static final QTableData readDelimitedFileAsTableData(Reader r, char[] delim, String comments,
		boolean hasRowNames, boolean hasColNames, boolean isAllNumeric) throws IOException
	{
		String[][] tbl = readDelimitedFile(r, delim, comments);
		int
			nrow = tbl.length,
			ncol = tbl[0].length;
		if (hasRowNames) ncol--;
		if (hasColNames) nrow--;
		double[][] data = new double[nrow][ncol];
		String[]
			colNames = hasColNames ? new String[ncol] : null,
			rowNames = hasRowNames ? new String[nrow] : null;
		String[][] colcatsArray = null;
		if (hasColNames)
			System.arraycopy(tbl[0], hasRowNames ? 1 : 0, colNames, 0, ncol);
		if (hasRowNames)
			for (int i = 0, j = hasColNames ? 1 : 0; i < nrow; i++, j++)
				rowNames[i] = tbl[j][0];
		if (isAllNumeric)
		{
			for (int i = 0, ii = hasColNames ? 1 : 0; i < nrow; i++, ii++)
				for (int j = 0, jj = hasRowNames ? 1 : 0; j < ncol; j++, jj++)
					data[i][j] = "".equals(tbl[ii][jj]) ? kUndefinedValue : Double.parseDouble(tbl[ii][jj]);
		} else
		{
			Set<String>[] colcats = new Set[colNames.length];
			colcatsArray = new String[colNames.length][];
			for (int j = 0, jj = hasRowNames ? 1 : 0; j < ncol; j++, jj++)
			{
				try {
					for (int i = 0, ii = hasColNames ? 1 : 0; i < nrow; i++, ii++)
						data[i][j] = "".equals(tbl[ii][jj]) ? kUndefinedValue : Double.parseDouble(tbl[ii][jj]);
				} catch (NumberFormatException e) {
					colcats[j] = new TreeSet<String>();
				}
			}
			for (int j = 0, jj = hasRowNames ? 1 : 0; j < ncol; j++, jj++)
			{
				if (colcats[j] == null)
					continue;
				for (int i = 0, ii = hasColNames ? 1 : 0; i < nrow; i++, ii++)
					colcats[j].add(tbl[ii][jj]);
				colcatsArray[j] = colcats[j].toArray(new String[colcats[j].size()]);
				List<String> cats = new ArrayList<String>();
				cats.addAll(colcats[j]);

				for (int i = 0, ii = hasColNames ? 1 : 0; i < nrow; i++, ii++)
					data[i][j] = "".equals(tbl[ii][jj]) ? kUndefinedValue : cats.indexOf(tbl[ii][jj]);
			}
		}

		QTableData tblData = new QTableData(data, rowNames, colNames);
		if (colcatsArray != null)
			tblData.setColumnCategories(colcatsArray);
		return tblData;
	}

	public static final String[][] readCSV(String filename) throws IOException
	{	return readDelimitedFile(filename, cCommaDelimiter, "#"); }

	public static final String[][] readCSV(File file) throws IOException
	{	return readDelimitedFile(file, cCommaDelimiter, "#"); }

	/**
	 * Read a comma-separated file. Same as <tt>readDelimitedFile(r, ",")</tt>.
	 * @param reader
	 * @return
	 * @throws IOException
	 */
	public static final String[][] readCSV(Reader reader) throws IOException
	{	return readDelimitedFile(reader, cCommaDelimiter, "#"); }

	public static final String[][] readTabFile(String filename) throws IOException
	{	return readDelimitedFile(filename, cTabDelimiter, "#"); }

	public static final String[][] readTabFile(File file) throws IOException
	{	return readDelimitedFile(file, cTabDelimiter, "#"); }

	/**
	 * Read a tab-separated file. Same as <tt>readDelimitedFile(r, "\t")</tt>.
	 * @param reader
	 * @return
	 * @throws IOException
	 */
	public static final String[][] readTabFile(Reader reader) throws IOException
	{	return readDelimitedFile(reader, cTabDelimiter, "#"); }

	/**
	 * Parse a delimited string and return tokens. Handles quotes.<br>
	 * Original author: Ian Darwin<br>
	 * Taken from: http://www.java2s.com/Code/Java/Development-Class/SimpledemoofCSVparserclass.htm<br>
	 * Optimized and fixed by: Roby Joehanes
	 * @param line
	 * @param sdelim
	 * @return
	 */
	public static final String[] tokenize(String line, String sdelim)
	{
		StringBuffer buf = new StringBuffer();
		List<String> tokenList = new ArrayList<String>();
		char delim = sdelim.charAt(0);
		int len = line.length();
		if (len > 0)
		{
			int i = 0;
			do {
				buf.setLength(0);
				if (i < len && line.charAt(i) == '"')
				{
					// skip quote
					int j;
					for (j = i + 1; j<len; j++)
					{
						if (line.charAt(j) == '"' && j+1 < len)
						{
							if (line.charAt(j+1) == '"')
							{
								j++; // skip escape char
							} else if (line.charAt(j+1) == delim) { //next delimeter
								i = j + 1; // skip end quotes
								break;
							}
						} else if (line.charAt(j) == '"' && j+1 == len) { // end quotes at end of line
							i = j;
							break; //done
						}
						buf.append(line.charAt(j));  // regular character.
					}
				}
				else
				{
					int j = line.indexOf(delim, i);
					if (j == -1)
					{
						buf.append(line.substring(i));
						i = len;
					} else {
						buf.append(line.substring(i, j));
						i = j;
					}
				}
				tokenList.add(buf.toString().intern());
			} while (i++ < len);
		}
		int numTokens = tokenList.size();
		return numTokens == 0 ? null : tokenList.toArray(new String[numTokens]);
	}

	public static final int countNumLines(String filename) throws IOException {
	    InputStream is = new BufferedInputStream(new FileInputStream(filename));
	    try {
	        byte[] c = new byte[1024];
	        int count = 0;
	        int readChars = 0;
	        boolean empty = true;
	        while ((readChars = is.read(c)) != -1) {
	            empty = false;
	            for (int i = 0; i < readChars; ++i) {
	                if (c[i] == '\n') {
	                    ++count;
	                }
	            }
	        }
	        return (count == 0 && !empty) ? 1 : count;
	    } finally {
	        is.close();
	    }
	}

	public static final Map<String, Map<String, String>> readIniFile(String filename, boolean caseSensitiveKey) throws IOException {
		return readIniFile(new BufferedReader(new FileReader(filename)), caseSensitiveKey);
	}

	public static final Map<String, Map<String, String>> readIniFile(File f, boolean caseSensitiveKey) throws IOException {
		return readIniFile(new BufferedReader(new FileReader(f)), caseSensitiveKey);
	}

	public static final Map<String, Map<String, String>> readIniFile(Reader r, boolean caseSensitiveKey) throws IOException {
		String[] lines = QFileUtils.readTextFileAsArray(r, new String[] { ";", "#", "//"});
		if (lines == null || lines.length == 0)
			return null;
		Map<String, Map<String, String>> table = new HashMap<String, Map<String,String>>();
		String curSectionName = "";
		for (String curLine: lines) {
			if (curLine.startsWith("[") && curLine.endsWith("]")) {
				curSectionName = curLine.substring(1, curLine.length() - 1);
				continue;
			}
			int eqPos = curLine.indexOf('=');
			if (eqPos < 0)
				throw new RuntimeException("There is no equal sign in the ini file!");
			String
				key = curLine.substring(0, eqPos).trim(),
				value = curLine.substring(eqPos+1).trim();
			if (key.length() == 0)
				throw new RuntimeException("The key is empty!");
			if (value.length() == 0)
				throw new RuntimeException("The value is empty!");
			Map<String, String> curSection = table.get(curSectionName);
			if (curSection == null) {
				curSection = new HashMap<String, String>();
				table.put(curSectionName, curSection);
			}
			curSection.put(caseSensitiveKey ? key : key.toLowerCase(Locale.ENGLISH), value);
		}
		return table;
	}

	public static void main(String[] args) {
		try {
			long time1 = System.currentTimeMillis();
			QTableData tbl = QTableData.load(args[0], cCommaDelimiter, "#", true, true, false);
			long time2 = System.currentTimeMillis();
			System.out.println(tbl.getNumberOfRows() + " x " + tbl.getNumberOfColumns());
			System.out.println("Time = " + (time2 - time1));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
