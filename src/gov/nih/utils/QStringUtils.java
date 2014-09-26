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

import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.lang.reflect.Array;
import java.text.Collator;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.StringTokenizer;
import java.util.regex.Pattern;

import static gov.nih.utils.QDataUtils.kUndefinedValue;

/**
 * A utility class to process strings.
 * 
 * @author Roby Joehanes and others
 */
public class QStringUtils
{
	public static final String
		sWordDelimiter = ", ", //$NON-NLS-1$
		sLn = System.getProperty("line.separator"), //$NON-NLS-1$
		sEq = " = ", //$NON-NLS-1$
		sEmpty = "", //$NON-NLS-1$
		sCm = sWordDelimiter,
		sDoubleQuote = "\"", //$NON-NLS-1$
		sSingleQuote = "\'", //$NON-NLS-1$
		sComma = ",", //$NON-NLS-1$
		sTab = "\t", //$NON-NLS-1$
		sSp = " ", //$NON-NLS-1$
		sPoundSign = "#", //$NON-NLS-1$
		sDoubleSlash = "//", //$NON-NLS-1$
		sPoundOnlyComment[] = new String[] { sPoundSign },
		sDoubleSlashOnlyComment[] = new String[] { sDoubleSlash },
		sPoundOrDoubleSlashComment[] = new String[] { sPoundSign, sDoubleSlash };

	public static final char[]
		cCommaDelimiter = new char[] {','},
		cTabDelimiter = new char[] {'\t'},
		cSpaceOrTabDelimiter = new char[] {' ', '\t'},
		cCommonDelimiter = new char[] { ' ', '\t', ','};

	/**
	 *	Formats a float to a nice string for printing
	 */
	public static final String floatToString(float floatValue, int numDecimalPlaces)
	{
		NumberFormat numFormat = NumberFormat.getInstance();
		numFormat.setMaximumFractionDigits(numDecimalPlaces);
		return numFormat.format(floatValue);	// (double) removed CN 2/29/04
	}

	/**
	 *	Formats a float to a nice string for printing
	 */
	public static final String doubleToString(double value, int numDecimalPlaces)
	{
		NumberFormat numFormat = NumberFormat.getInstance();
		numFormat.setMaximumFractionDigits(numDecimalPlaces);
		return numFormat.format(value);	// (double) removed CN 2/29/04
	}

	// For debugging purposes only -- RJ 2004/11/30
	public static final String toString(float[][][] f)
	{
		StringBuilder buffer = new StringBuilder();
		int dim1 = f.length;
		for (int i = 0; i < dim1; i++)
		{
			float[][] f0 = f[i];
			int dim2 = f0.length;
			for (int j = 0; j < dim2; j++)
			{
				float[] f1 = f0[i];
				int dim3 = f1.length;
				buffer.append(i+","+j+": "); //$NON-NLS-1$ //$NON-NLS-2$
				for (int k = 0; k < dim3; k++)
					buffer.append(f[i][j][k]+","); //$NON-NLS-1$
				buffer.append(sLn);
			}
		}
		return buffer.toString();
	}

	/**
	* Checks if the input string is a number.
	* Taken out of QSwingFactory.QNumDocument
	* Simplified by RJ 2006/01/04
	*/
	public static final boolean isNumericalString(String inputStr)
	{
		try {
			// If Double can parse it, then it's a number
			Double.parseDouble(inputStr);
			return true;
		} catch (Exception e) {
			return false;
		}
	}

	/**
	 * Unescape the string s. The string s is expected to have Java's escape string
	 * Adapted from Apache's org.apache.commons.lang.StringEscapeUtils
	 * RJ 01/02/2006
	 * 
	 * @param s
	 * @return
	 */
	public static final String unescape(String str)
	{
		if (str == null)
			return null;

		int strLength = str.length();
		StringBuffer
			unicode = new StringBuffer(4),
			output = new StringBuffer(strLength);
		boolean
			hadSlash = false,
			inUnicode = false;
		for (int i = 0; i < strLength; i++)
		{
			char ch = str.charAt(i);

			// If we previously have a /u instead, then we need to get the next four
			// digits to translate to unicode string.
			if (inUnicode)
			{
				// if in unicode, then we're reading unicode values in somehow
				unicode.append(ch);
				if (unicode.length() == 4)
				{
					// unicode now contains the four hex digits
					// which represents our unicode character
					try
					{
						int value = Integer.parseInt(unicode.toString(), 16);
						output.append((char) value);
					}
					catch (NumberFormatException nfe)
					{
						// If the unicode is invalid, then append the code as it is
						output.append(unicode.toString());
	                }
					finally
					{
						unicode.setLength(0);
						inUnicode = hadSlash = false;
					}
	                continue;
	            }
			}

			if (hadSlash)
			{
				// handle an escaped value
				hadSlash = false;
				switch (ch)
				{
					case '\\':
						output.append('\\');
						break;
					case '\'':
						output.append('\'');
						break;
					case '\"':
						output.append('"');
						break;
					case 'r':
						output.append('\r');
						break;
					case 'f':
						output.append('\f');
						break;
					case 't':
						output.append('\t');
						break;
					case 'n':
						output.append('\n');
						break;
					case 'b':
						output.append('\b');
						break;
					case 'u':
						inUnicode = true;
						break;
					default :
						if (!inUnicode)
							output.append(ch);
						break;
				}
				continue;
			}
			else if (ch == '\\')
			{
				hadSlash = true;
				continue;
			}
			if (!inUnicode)
				output.append(ch);
		}

		if (hadSlash)
		{
			// then we're in the weird case of a \ at the end of the
			// string, let's output it anyway.
			output.append('\\');
		}
		return output.toString();
	}

	/**
	 * Escape the string s. Special characters such as \n \r in str will
	 * be replaced with appropriate escape code
	 * Adapted from Apache's org.apache.commons.lang.StringEscapeUtils
	 * RJ 01/02/2006
	 * 
	 * @param str
	 * @return
	 */
	public static final String escape(String str)
	{
		if (str == null)
			return null;

		int strLength = str.length();
		StringBuffer output = new StringBuffer(strLength);
		for (int i = 0; i < strLength; i++)
		{
			char ch = str.charAt(i);

			// handle unicode
			if (ch > 0xfff)
			{
				output.append("\\u" + hex(ch)); //$NON-NLS-1$
			}
			else if (ch > 0xff)
			{
				output.append("\\u0" + hex(ch)); //$NON-NLS-1$
			}
			else if (ch > 0x7f)
			{
				output.append("\\u00" + hex(ch)); //$NON-NLS-1$
			}
			else if (ch < 32)
			{
				switch (ch)
				{
					case '\b':
						output.append("\\b"); //$NON-NLS-1$
						break;
					case '\n':
						output.append("\\n"); //$NON-NLS-1$
						break;
					case '\t':
						output.append("\\t"); //$NON-NLS-1$
						break;
					case '\f':
						output.append("\\f"); //$NON-NLS-1$
						break;
					case '\r':
						output.append("\\r"); //$NON-NLS-1$
						break;
                    default :
                    	if (ch > 0xf)
                    		output.append("\\u00" + hex(ch)); //$NON-NLS-1$
                    	else
                    		output.append("\\u000" + hex(ch)); //$NON-NLS-1$
                    	break;
				}
            }
			else
			{
				switch (ch)
				{
					case '"':
						output.append("\\\'"); //$NON-NLS-1$
						break;
					case '\\':
						output.append("\\\\"); //$NON-NLS-1$
						break;
					default :
						output.append(ch);
						break;
				}
			}
		}
		return output.toString();
    }

	/**
	 * <p>Returns an upper case hexadecimal <code>String</code> for the given
	 * character.</p>
	 * 
	 * @param ch The character to convert.
	 * @return An upper case hexadecimal <code>String</code>
	 */
	private static final String hex(char ch)
	{
		return Integer.toHexString(ch).toUpperCase(Locale.ENGLISH);
	}

	/**
	 * Split the string s on whitespaces and commas. For example:
	 * The string "a, b  ,c ,   d e" (without the quotes)
	 * shall return String[] { "a", "b", "c", "d", "e" }
	 * @param s
	 * @return array of strings
	 */
	public static final String[] toStringArray(String s)
	{	return toStringArray(s, sWordDelimiter); }

	public static final String[] toStringArray(String s, String delimiter)
	{
		StringTokenizer tokenizer = new StringTokenizer(s, delimiter);
		List<String> tokens = new ArrayList<String>();
		while (tokenizer.hasMoreTokens())
			tokens.add(tokenizer.nextToken().trim());
		return tokens.toArray(new String[tokens.size()]);
	}

	/**
	 * Convert an array of int into an array of string
	 * @param values
	 * @return
	 */
	public static final String[] toStringArray(int[] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		String[] array = new String[numValues];
		for (int i = 0; i < numValues; i++)
			array[i] = String.valueOf(values[i]);
		return array;
	}

	/**
	 * Convert an array of double into an array of string
	 * @param values
	 * @return
	 */
	public static final String[] toStringArray(double[] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		String[] array = new String[numValues];
		for (int i = 0; i < numValues; i++)
			array[i] = String.valueOf(values[i]);
		return array;
	}

	/**
	 * Convert an array of int into an array of string
	 * @param values
	 * @return
	 */
	public static final String[] toStringArray(float[] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		String[] array = new String[numValues];
		for (int i = 0; i < numValues; i++)
			array[i] = String.valueOf(values[i]);
		return array;
	}

	/**
	 * Convert an array of double into an array of string
	 * @param values
	 * @return
	 */
	public static final String[][] toStringArray(double[][] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		String[][] array = new String[numValues][];
		for (int i = 0; i < numValues; i++)
		{
			double[] valueRow = values[i];
			int numRows = valueRow.length;
			String[] stringRow = array[i] = new String[numRows];
			for (int j = 0; j < numRows; j++)
				stringRow[j] = String.valueOf(valueRow[j]);
		}
		return array;
	}

	/**
	 * Convert an array of float into an array of string
	 * @param values
	 * @return
	 */
	public static final String[][] toStringArray(float[][] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		String[][] array = new String[numValues][];
		for (int i = 0; i < numValues; i++)
		{
			float[] valueRow = values[i];
			int numRows = valueRow.length;
			String[] stringRow = array[i] = new String[numRows];
			for (int j = 0; j < numRows; j++)
				stringRow[j] = String.valueOf(valueRow[j]);
		}
		return array;
	}

	/**
	 * Convert an array of int into an array of string
	 * @param values
	 * @return
	 */
	public static final String[][] toStringArray(int[][] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		String[][] array = new String[numValues][];
		for (int i = 0; i < numValues; i++)
		{
			int[] valueRow = values[i];
			int numRows = valueRow.length;
			String[] stringRow = array[i] = new String[numRows];
			for (int j = 0; j < numRows; j++)
				stringRow[j] = String.valueOf(valueRow[j]);
		}
		return array;
	}

	/**
	 * Return the string as int arrays. The value is assumed to
	 * be delimited with spaces and/or commas
	 * @param s
	 * @return
	 */
	public static final int[] toIntArray(String s)
	{	return toIntArray(s, sWordDelimiter); }

	/**
	 * Return the string as int arrays. The value is assumed to
	 * be delimited with spaces and/or commas
	 * @param s
	 * @return
	 */
	public static final int[] toIntArray(String s, String delimiter)
	{	return toIntArray(toStringArray(s, delimiter), kUndefinedValue); }

	public static final int[] toIntArray(String[] values)
	{	return toIntArray(values, kUndefinedValue); }

	public static final int[] toIntArray(String[] values, int defaultValue)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		int[] result = new int[numValues];
		for (int i = 0; i < numValues; i++)
			try {
				result[i] = Integer.parseInt(values[i]);
			} catch (Exception e)
			{	result[i] = defaultValue; }
		return result;
	}

	/**
	 * Return the string as long arrays. The value is assumed to
	 * be delimited with spaces and/or commas
	 * @param s
	 * @return
	 */
	public static final long[] toLongArray(String s)
	{
		String[] values = toStringArray(s);
		if (values == null)
			return null;
		int numValues = values.length;
		long[] result = new long[numValues];
		for (int i = 0; i < numValues; i++)
			result[i] = Long.parseLong(values[i]);
		return result;
	}

	/**
	 * Return the string as float arrays. The value is assumed to
	 * be delimited with spaces and/or commas
	 * @param s
	 * @return
	 */
	public static final float[] toFloatArray(String s)
	{
		String[] values = toStringArray(s);
		if (values == null)
			return null;
		int numValues = values.length;
		float[] result = new float[numValues];
		for (int i = 0; i < numValues; i++)
			result[i] = Float.parseFloat(values[i]);
		return result;
	}

	/**
	 * Return the string as double arrays. The value is assumed to
	 * be delimited with spaces and/or commas
	 * @param s
	 * @return
	 */
	public static final double[] toDoubleArray(String s)
	{	return toDoubleArray(toStringArray(s), kUndefinedValue); }

	public static final double[] toDoubleArray(String[] values)
	{	return toDoubleArray(values, kUndefinedValue); }

	public static final double[] toDoubleArray(String[] values, double defaultValue)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		double[] result = new double[numValues];
		for (int i = 0; i < numValues; i++)
			try {
				result[i] = Double.parseDouble(values[i]);
			} catch (Exception e)
			{	result[i] = defaultValue; }
		return result;
	}

	/**
	 * Return the string as boolean arrays. The value is assumed to
	 * be delimited with spaces and/or commas
	 * @param s
	 * @return
	 */
	public static final boolean[] toBooleanArray(String s)
	{
		String[] values = toStringArray(s);
		if (values == null)
			return null;
		int numValues = values.length;
		boolean[] result = new boolean[numValues];
		for (int i = 0; i < numValues; i++)
			result[i] = Boolean.parseBoolean(values[i]);
		return result;
	}

	public static final <T> String toString(T[] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert an array of string values into a single,
	 * delimited by delimiter
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final <T> String toString(T[] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		if (numValues > 0)
		{
			buffer.append(values[0].toString());
			for (int i = 1; i < numValues; i++)
				buffer.append(delimiter + values[i].toString());
		}
		return buffer.toString();
	}

	public static final String toString(int[] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert an array of int values into a single string,
	 * delimited by delimiter
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(int[] values, String delimiter)
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

	public static final String toString(byte[] values)
	{	return toString(values, sWordDelimiter); }

	public static final String toString(byte[] values, String delimiter)
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

	public static final String toIndexString(int[] values)
	{	return toIndexString(values, sWordDelimiter); }

	/**
	 * Convert an array of int values into a single string,
	 * delimited by delimiter, showing indices
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toIndexString(int[] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int
			sum = 0,
			numValues = values.length;
		for (int i = 0; i < numValues; i++)
		{
			int val = values[i];
			sum += val;
			buffer.append(i + sEq + val + delimiter);
		}
		buffer.append(" Total = " + sum); //$NON-NLS-1$
		return buffer.toString();
	}

	public static final String toString(int[][] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert a double dimension array of int values into string.
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(int[][] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		for (int i = 0; i < numValues; i++)
			buffer.append(toString(values[i], delimiter) + sLn);
		return buffer.toString();
	}

	public static final String toString(long[] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert an array of int values into a single,
	 * delimited by delimiter
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(long[] values, String delimiter)
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

	public static final String toString(long[][] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert a double dimension array of long values into string.
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(long[][] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		for (int i = 0; i < numValues; i++)
			buffer.append(toString(values[i], delimiter) + sLn);
		return buffer.toString();
	}

	public static final String toString(short[] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert an array of int values into a single,
	 * delimited by delimiter
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(short[] values, String delimiter)
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

	public static final String toString(float[] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert an array of float values into a single,
	 * comma delimited string
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
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

	public static final String toString(float[][] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert a double dimension array of float values into string.
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(float[][] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		for (int i = 0; i < numValues; i++)
			buffer.append(toString(values[i], delimiter) + sLn);
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

	public static final String toString(double[][] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert a double dimension array of double values into string.
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toString(double[][] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		for (int i = 0; i < numValues; i++)
			buffer.append(toString(values[i], delimiter) + sLn);
		return buffer.toString();
	}

	public static final String toStringAsMatrix(double[][] values)
	{	return toStringAsMatrix(values, sWordDelimiter); }

	/**
	 * Convert a double dimension array of double values into string as matrix.
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final String toStringAsMatrix(double[][] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		String close = "]" + sLn; //$NON-NLS-1$
		for (int i = 0; i < numValues; i++)
			buffer.append("[" + toString(values[i], delimiter) + close); //$NON-NLS-1$
		return buffer.toString();
	}

	public static final <T> String toString(T[][] values)
	{	return toString(values, sWordDelimiter); }

	/**
	 * Convert a double dimension array of Object values into string.
	 * 
	 * @param values
	 * @param delimiter
	 * @return
	 */
	public static final <T> String toString(T[][] values, String delimiter)
	{
		if (values == null)
			return null;
		StringBuilder buffer = new StringBuilder();
		int numValues = values.length;
		for (int i = 0; i < numValues; i++)
			buffer.append(toString(values[i], delimiter) + sLn);
		return buffer.toString();
	}

	public static final String toDelimitedString(double[][] data, String[] rowNames, String[] colNames, String delimiter)
	{
		StringBuilder buf = new StringBuilder();
		int
			numCols = data[0].length,
			numRows = data.length;

		if (colNames != null)
		{
			if (colNames.length != numCols)
				throw new RuntimeException();
			if (rowNames != null)
				buf.append(delimiter);
			buf.append(colNames[0]);
			for (int colNo = 1; colNo < numCols; colNo++)
			{
				buf.append(delimiter);
				buf.append(colNames[colNo]);
			}
			buf.append(sLn);
		}

		if (rowNames != null && rowNames.length != numRows)
			throw new RuntimeException();

		if (rowNames != null)
		{
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				buf.append(rowNames[rowNo]);
				double[] curResult = data[rowNo];
				for (int colNo = 0; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		} else {
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				double[] curResult = data[rowNo];
				buf.append(curResult[0]);
				for (int colNo = 1; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		}
		return buf.toString();
	}

	public static final String toDelimitedString(double[][] data, int[] rowNames, String[] colNames, String delimiter)
	{
		StringBuilder buf = new StringBuilder();
		int
			numCols = data[0].length,
			numRows = data.length;

		if (colNames != null)
		{
			if (colNames.length != numCols)
				throw new RuntimeException();
			if (rowNames != null)
				buf.append(delimiter);
			buf.append(colNames[0]);
			for (int colNo = 1; colNo < numCols; colNo++)
			{
				buf.append(delimiter);
				buf.append(colNames[colNo]);
			}
			buf.append(sLn);
		}

		if (rowNames != null && rowNames.length != numRows)
			throw new RuntimeException();

		if (rowNames != null)
		{
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				buf.append(rowNames[rowNo]);
				double[] curResult = data[rowNo];
				for (int colNo = 0; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		} else {
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				double[] curResult = data[rowNo];
				buf.append(curResult[0]);
				for (int colNo = 1; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		}
		return buf.toString();
	}

	public static final String toDelimitedString(int[][] data, String[] rowNames, String[] colNames, String delimiter)
	{
		StringBuilder buf = new StringBuilder();
		int
			numCols = data[0].length,
			numRows = data.length;

		if (colNames != null)
		{
			if (colNames.length != numCols)
				throw new RuntimeException();
			if (rowNames != null)
				buf.append(delimiter);
			buf.append(colNames[0]);
			for (int colNo = 1; colNo < numCols; colNo++)
			{
				buf.append(delimiter);
				buf.append(colNames[colNo]);
			}
			buf.append(sLn);
		}

		if (rowNames != null && rowNames.length != numRows)
			throw new RuntimeException();

		if (rowNames != null)
		{
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				buf.append(rowNames[rowNo]);
				int[] curResult = data[rowNo];
				for (int colNo = 0; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		} else {
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				int[] curResult = data[rowNo];
				buf.append(curResult[0]);
				for (int colNo = 1; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		}
		return buf.toString();
	}

	public static final String toDelimitedString(int[][] data, int[] rowNames, String[] colNames, String delimiter)
	{
		StringBuilder buf = new StringBuilder();
		int
			numCols = data[0].length,
			numRows = data.length;

		if (colNames != null)
		{
			if (colNames.length != numCols)
				throw new RuntimeException();
			if (rowNames != null)
				buf.append(delimiter);
			buf.append(colNames[0]);
			for (int colNo = 1; colNo < numCols; colNo++)
			{
				buf.append(delimiter);
				buf.append(colNames[colNo]);
			}
			buf.append(sLn);
		}

		if (rowNames != null && rowNames.length != numRows)
			throw new RuntimeException();

		if (rowNames != null)
		{
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				buf.append(rowNames[rowNo]);
				int[] curResult = data[rowNo];
				for (int colNo = 0; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		} else {
			for (int rowNo = 0; rowNo < numRows; rowNo++)
			{
				int[] curResult = data[rowNo];
				buf.append(curResult[0]);
				for (int colNo = 1; colNo < numCols; colNo++)
				{
					buf.append(delimiter);
					buf.append(curResult[colNo]);
				}
				buf.append(sLn);
			}
		}
		return buf.toString();
	}

	/**
	 * Search a set of strings in a collection that matched the specified
	 * regular expression
	 * @param collection
	 * @param regex
	 * @return The set of matched values
	 */
	public static final Collection<String> searchCollection(Collection<String> collection, String regex)
	{
		Pattern pattern = Pattern.compile(regex);
		List<String> values = new ArrayList<String>();
		for (String value : collection)
		{
			if (pattern.matcher(value).matches())
				values.add(value);
		}
		return values;
	}

	/**
	 * Divide f by 10 and return the result in string. Needed for pretty print locus positions.
	 * @param f
	 * @return
	 */
	public static final String divideByTen(int f)
	{
		String str = String.valueOf(f);
		int strlen = str.length();
		char lastChar = str.charAt(strlen - 1);
		str = str.substring(0, strlen-1) + "." + lastChar; //$NON-NLS-1$
		if (strlen < 2)
			str = '0' + str;
		return str;
	}

	/**
	 * Save text to system clipboard.
	 * Taken from: http://www.javapractices.com/Topic82.cjp
	 * @param text
	 */
	public static final void saveTextToClipboard(String text)
	{
	    StringSelection stringSelection = new StringSelection(text);
	    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
	    clipboard.setContents(stringSelection, null);
	}

	/**
	 * Get text from system clipboard.
	 */
	public static final String fetchTextFromClipboard()
	{
		Transferable t = Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
		DataFlavor textFlavor = DataFlavor.getTextPlainUnicodeFlavor();

		try
		{
			if (t != null && t.isDataFlavorSupported(textFlavor))
				return (String) t.getTransferData(textFlavor);
		} catch (Exception e) {}
		return null;
	}

	/**
	 * Stringize generic types
	 * @param <T>
	 * @param data
	 * @return
	 */
	public static final <T> String toString(T data)
	{
		Class<?> dataClass = data.getClass();
		if (!dataClass.isArray())
			return data.toString();
		int length = Array.getLength(data);
		StringBuilder buffer = new StringBuilder();
		buffer.append("{" + sLn); //$NON-NLS-1$
		for (int i = 0; i < length; i++)
		{
			Object element = Array.get(data, i);
			buffer.append(toString(element) + sLn); //$NON-NLS-1$
		}
		buffer.append("}" + sLn); //$NON-NLS-1$
		return buffer.toString();
	}

	/**
	 * Replicates the functionality of strdup from C
	 * (i.e. repeat string s, t times) 
	 * @param s
	 * @param t
	 * @return the repeated string
	 */
	public static final String strdup(String s, int t)
	{
		StringBuilder buf = new StringBuilder(s);
		for (int i = 0; i < t; i++)
			buf.append(s);
		return buf.toString();
	}

	public static final String strdup(char c, int t)
	{	return strdup(String.valueOf(c), t); }

	public static final int[] parseVersionString(String version)
	{
		String[] ver = version.split("-")[0].split("[\\._]"); //$NON-NLS-1$ //$NON-NLS-2$
		int
			numTokens = ver.length,
			verNo[] = new int[numTokens];
		for (int i = 0; i < numTokens; i++)
			verNo[i] = Integer.parseInt(ver[i]);
		return verNo;
	}

	/**
	 * Unquote an array of strings. Only works with matched quotes (e.g., "word" or 'word', but not "word')
	 * @param strs
	 * @param trim, if true, remove any leading or trailing spaces (after quote removal).
	 * @return
	 */
	public static final String[] unquoteStrings(String[] strs, boolean trim)
	{
		int n = strs.length;
		String[] result = new String[n];
		for (int i = 0; i < n; i++)
			result[i] = unquoteString(strs[i], trim);
		return result;
	}

	public static final String unquoteString(String str, boolean trim)
	{
		if (str != null)
		{
			str = str.trim();
			if (str.length() > 1) {
				if (str.startsWith(sDoubleQuote) && str.endsWith(sDoubleQuote))
					str = str.substring(1, str.length() - 1);
				else if (str.startsWith(sSingleQuote) && str.endsWith(sSingleQuote))
					str = str.substring(1, str.length() - 1);
			}
			if (trim)
				str = str.trim();
		}
		return str;
	}

	// This is from Stephen Friedrich's blog:
	// http://weblogs.java.net/blog/skelvin/archive/2006/01/natural_string.html
	// RJ's modification: Make the method's modifiers final.
	// STEPHEN FRIEDRICH'S CODE BEGINS
	/*
	 * Copyright (c) 2006, Stephen Kelvin Friedrich,  All rights reserved.
	 *
	 * This a BSD license. If you use or enhance the code, I'd be pleased if you sent a mail to s.friedrich@eekboom.com
	 *
	 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
	 * following conditions are met:
	 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
	 *       following disclaimer.
	 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
	         following disclaimer in the documentation and/or other materials provided with the distribution.
	 *     * Neither the name of the "Stephen Kelvin Friedrich" nor the names of its contributors may be used to endorse
	 *       or promote products derived from this software without specific prior written permission.
	 *
	 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
	 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
	 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	 */

	/**
	 * <p>A string comparator that does case sensitive comparisons and handles embedded numbers correctly.</p>
	 * <p><b>Do not use</b> if your app might ever run on any locale that uses more than 7-bit ascii characters.</p>
	 */
	private static final Comparator<String> NATURAL_COMPARATOR_ASCII = new Comparator<String>() {
		public int compare(String o1, String o2) {
			return compareNaturalAscii(o1, o2);
		}
	};

	/**
	 * <p>A string comparator that does case insensitive comparisons and handles embedded numbers correctly.</p>
	 * <p><b>Do not use</b> if your app might ever run on any locale that uses more than 7-bit ascii characters.</p>
	 */
	private static final Comparator<String> IGNORE_CASE_NATURAL_COMPARATOR_ASCII = new Comparator<String>() {
		public int compare(String o1, String o2) {
			return compareNaturalIgnoreCaseAscii(o1, o2);
		}
	};

	/**
	 * Returns a comparator that compares contained numbers based on their numeric values and compares other parts
	 * using the current locale's order rules.
	 * <p>For example in German locale this will be a comparator that handles umlauts correctly and ignores
	 * upper/lower case differences.</p>
	 *
	 * @return <p>A string comparator that uses the current locale's order rules and handles embedded numbers
	 *         correctly.</p>
	 * @see #getNaturalComparator(java.text.Collator)
	 */
	public static final Comparator<String> getNaturalComparator() {
		Collator collator = Collator.getInstance();
		return getNaturalComparator(collator);
	}

	/**
	 * Returns a comparator that compares contained numbers based on their numeric values and compares other parts
	 * using the given collator.
	 *
	 * @param collator used for locale specific comparison of text (non-number) subwords - must not be null
	 * @return <p>A string comparator that uses the given Collator to compare subwords and handles embedded numbers
	 *         correctly.</p>
	 * @see #getNaturalComparator()
	 */
	public static final Comparator<String> getNaturalComparator(final Collator collator) {
		if(collator == null) {
			// it's important to explicitly handle this here - else the bug will manifest anytime later in possibly
			// unrelated code that tries to use the comparator
			throw new NullPointerException("collator must not be null");
		}
		return new Comparator<String>() {
			public int compare(String o1, String o2) {
				return compareNatural(collator, o1, o2);
			}
		};
	}

	/**
	 * Returns a comparator that compares contained numbers based on their numeric values and compares other parts
	 * based on each character's Unicode value.
	 *
	 * @return <p>a string comparator that does case sensitive comparisons on pure ascii strings and handles embedded
	 *         numbers correctly.</p>
	 *         <b>Do not use</b> if your app might ever run on any locale that uses more than 7-bit ascii characters.
	 * @see #getNaturalComparator()
	 * @see #getNaturalComparator(java.text.Collator)
	 */
	public static final Comparator<String> getNaturalComparatorAscii() {
		return NATURAL_COMPARATOR_ASCII;
	}

	/**
	 * Returns a comparator that compares contained numbers based on their numeric values and compares other parts
	 * based on each character's Unicode value while ignore upper/lower case differences.
	 * <b>Do not use</b> if your app might ever run on any locale that uses more than 7-bit ascii characters.
	 *
	 * @return <p>a string comparator that does case insensitive comparisons on pure ascii strings and handles embedded
	 *         numbers correctly.</p>
	 * @see #getNaturalComparator()
	 * @see #getNaturalComparator(java.text.Collator)
	 */
	public static final Comparator<String> getNaturalComparatorIgnoreCaseAscii() {
		return IGNORE_CASE_NATURAL_COMPARATOR_ASCII;
	}

	/**
	 * <p>Compares two strings using the current locale's rules and comparing contained numbers based on their numeric
	 * values.</p>
	 * <p>This is probably the best default comparison to use.</p>
	 * <p>If you know that the texts to be compared are in a certain language that differs from the default locale's
	 * langage, then get a collator for the desired locale ({@link java.text.Collator#getInstance(java.util.Locale)})
	 * and pass it to {@link #compareNatural(java.text.Collator, String, String)}</p>
	 *
	 * @param s first string
	 * @param t second string
	 * @return zero iff <code>s</code> and <code>t</code> are equal,
	 *         a value less than zero iff <code>s</code> lexicographically precedes <code>t</code>
	 *         and a value larger than zero iff <code>s</code> lexicographically follows <code>t</code>
	 */
	public static final int compareNatural(String s, String t) {
		return compareNatural(s, t, false, Collator.getInstance());
	}

	/**
	 * <p>Compares two strings using the given collator and comparing contained numbers based on their numeric
	 * values.</p>
	 *
	 * @param s first string
	 * @param t second string
	 * @return zero iff <code>s</code> and <code>t</code> are equal,
	 *         a value less than zero iff <code>s</code> lexicographically precedes <code>t</code>
	 *         and a value larger than zero iff <code>s</code> lexicographically follows <code>t</code>
	 */
	public static final int compareNatural(Collator collator, String s, String t) {
		return compareNatural(s, t, true, collator);
	}

	/**
	 * <p>Compares two strings using each character's Unicode value for non-digit characters and the numeric values off
	 * any contained numbers.</p>
	 * <p>(This will probably make sense only for strings containing 7-bit ascii characters only.)</p>
	 *
	 * @return zero iff <code>s</code> and <code>t</code> are equal,
	 *         a value less than zero iff <code>s</code> lexicographically precedes <code>t</code>
	 *         and a value larger than zero iff <code>s</code> lexicographically follows <code>t</code>
	 */
	public static final int compareNaturalAscii(String s, String t) {
		return compareNatural(s, t, true, null);
	}

	/**
	 * <p>Compares two strings using each character's Unicode value - ignoring upper/lower case - for non-digit
	 * characters and the numeric values of any contained numbers.</p>
	 * <p>(This will probably make sense only for strings containing 7-bit ascii characters only.)</p>
	 *
	 * @return zero iff <code>s</code> and <code>t</code> are equal,
	 *         a value less than zero iff <code>s</code> lexicographically precedes <code>t</code>
	 *         and a value larger than zero iff <code>s</code> lexicographically follows <code>t</code>
	 */
	public static final int compareNaturalIgnoreCaseAscii(String s, String t) {
		return compareNatural(s, t, false, null);
	}

	/**
	 * @param s             first string
	 * @param t             second string
	 * @param caseSensitive treat characters differing in case only as equal - will be ignored if a collator is given
	 * @param collator      used to compare subwords that aren't numbers - if null, characters will be compared
	 *                      individually based on their Unicode value
	 * @return zero iff <code>s</code> and <code>t</code> are equal,
	 *         a value less than zero iff <code>s</code> lexicographically precedes <code>t</code>
	 *         and a value larger than zero iff <code>s</code> lexicographically follows <code>t</code>
	 */
	private static final int compareNatural(String s, String t, boolean caseSensitive, Collator collator) {
		int sIndex = 0;
		int tIndex = 0;

		int sLength = s.length();
		int tLength = t.length();

		while(true) {
			// both character indices are after a subword (or at zero)

			// Check if one string is at end
			if(sIndex == sLength && tIndex == tLength) {
				return 0;
			}
			if(sIndex == sLength) {
				return -1;
			}
			if(tIndex == tLength) {
				return 1;
			}

			// Compare sub word
			char sChar = s.charAt(sIndex);
			char tChar = t.charAt(tIndex);

			boolean sCharIsDigit = Character.isDigit(sChar);
			boolean tCharIsDigit = Character.isDigit(tChar);

			if(sCharIsDigit && tCharIsDigit) {
				// Compare numbers

				// skip leading 0s
				int sLeadingZeroCount = 0;
				while(sChar == '0') {
					++sLeadingZeroCount;
					++sIndex;
					if(sIndex == sLength) {
						break;
					}
					sChar = s.charAt(sIndex);
				}
				int tLeadingZeroCount = 0;
				while(tChar == '0') {
					++tLeadingZeroCount;
					++tIndex;
					if(tIndex == tLength) {
						break;
					}
					tChar = t.charAt(tIndex);
				}
				boolean sAllZero = sIndex == sLength || !Character.isDigit(sChar);
				boolean tAllZero = tIndex == tLength || !Character.isDigit(tChar);
				if(sAllZero && tAllZero) {
					continue;
				}
				if(sAllZero && !tAllZero) {
					return -1;
				}
				if(tAllZero) {
					return 1;
				}

				int diff = 0;
				do {
					if(diff == 0) {
						diff = sChar - tChar;
					}
					++sIndex;
					++tIndex;
					if(sIndex == sLength && tIndex == tLength) {
						return diff != 0 ? diff : sLeadingZeroCount - tLeadingZeroCount;
					}
					if(sIndex == sLength) {
						if(diff == 0) {
							return -1;
						}
						return Character.isDigit(t.charAt(tIndex)) ? -1 : diff;
					}
					if(tIndex == tLength) {
						if(diff == 0) {
							return 1;
						}
						return Character.isDigit(s.charAt(sIndex)) ? 1 : diff;
					}
					sChar = s.charAt(sIndex);
					tChar = t.charAt(tIndex);
					sCharIsDigit = Character.isDigit(sChar);
					tCharIsDigit = Character.isDigit(tChar);
					if(!sCharIsDigit && !tCharIsDigit) {
						// both number sub words have the same length
						if(diff != 0) {
							return diff;
						}
						break;
					}
					if(!sCharIsDigit) {
						return -1;
					}
					if(!tCharIsDigit) {
						return 1;
					}
				} while(true);
			}
			else {
				// Compare words
				if(collator != null) {
					// To use the collator the whole subwords have to be compared - character-by-character comparision
					// is not possible. So find the two subwords first
					int aw = sIndex;
					int bw = tIndex;
					do {
						++sIndex;
					} while(sIndex < sLength && !Character.isDigit(s.charAt(sIndex)));
					do {
						++tIndex;
					} while(tIndex < tLength && !Character.isDigit(t.charAt(tIndex)));

					String as = s.substring(aw, sIndex);
					String bs = t.substring(bw, tIndex);
					int subwordResult = collator.compare(as, bs);
					if(subwordResult != 0) {
						return subwordResult;
					}
				}
				else {
					// No collator specified. All characters should be ascii only. Compare character-by-character.
					do {
						if(sChar != tChar) {
							if(caseSensitive) {
								return sChar - tChar;
							}
							sChar = Character.toUpperCase(sChar);
							tChar = Character.toUpperCase(tChar);
							if(sChar != tChar) {
								sChar = Character.toLowerCase(sChar);
								tChar = Character.toLowerCase(tChar);
								if(sChar != tChar) {
									return sChar - tChar;
								}
							}
						}
						++sIndex;
						++tIndex;
						if(sIndex == sLength && tIndex == tLength) {
							return 0;
						}
						if(sIndex == sLength) {
							return -1;
						}
						if(tIndex == tLength) {
							return 1;
						}
						sChar = s.charAt(sIndex);
						tChar = t.charAt(tIndex);
						sCharIsDigit = Character.isDigit(sChar);
						tCharIsDigit = Character.isDigit(tChar);
					} while(!sCharIsDigit && !tCharIsDigit);
				}
			}
		}
	}
	// STEPHEN FRIEDRICH'S CODE ENDS

}
