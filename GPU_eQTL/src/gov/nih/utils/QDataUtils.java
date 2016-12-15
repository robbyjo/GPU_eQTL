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

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

/**
 * Utilities for manipulating data especially in array
 * @author Roby Joehanes
 *
 */
public class QDataUtils
{
	public static final int kUndefinedValue = Integer.MIN_VALUE;
	private static final int kInsertionTreshold = 4;

	/**
	 * Utility to shorten an array to newLength
	 * @param source
	 * @param newLength
	 * @return
	 */
	public static final int[] shortenArray(int[] source, int newLength)
	{
		int[] newArray = new int[newLength];
		System.arraycopy(source, 0, newArray, 0, newLength);
		return newArray;
	}

	public static final float[] shortenArray(float[] source, int newLength)
	{
		float[] newArray = new float[newLength];
		System.arraycopy(source, 0, newArray, 0, newLength);
		return newArray;
	}

	public static final double[] shortenArray(double[] source, int newLength)
	{
		double[] newArray = new double[newLength];
		System.arraycopy(source, 0, newArray, 0, newLength);
		return newArray;
	}

	public static final <T> T[] shortenArray(T[] source, int newLength)
	{
		T[] newArray = (T[]) Array.newInstance(source.getClass().getComponentType(), newLength);
		source.clone();
		System.arraycopy(source, 0, newArray, 0, newLength);
		return newArray;
	}

	/**
	 * Concatenate two arrays
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double[] concatArray(double[] a, double[] b)
	{
		int
			lenA = a.length,
			lenB = b.length;
		double[] result = new double[lenA + lenB];
		System.arraycopy(a, 0, result, 0, lenA);
		System.arraycopy(b, 0, result, lenA, lenB);
		return result;
	}

	/**
	 * Concatenate two arrays
	 * @param a
	 * @param b
	 * @return
	 */
	public static final float[] concatArray(float[] a, float[] b)
	{
		int
			lenA = a.length,
			lenB = b.length;
		float[] result = new float[lenA + lenB];
		System.arraycopy(a, 0, result, 0, lenA);
		System.arraycopy(b, 0, result, lenA, lenB);
		return result;
	}

	/**
	 * Flatten 2D array into a unidimensional array.
	 * Array is assumed to be of the same lengths.
	 * @param arr
	 * @return
	 */
	public static final double[] flattenArray(double[][] arr)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		double[] flatArr = new double[rows * cols];
		for (int i = 0; i < rows; i++, idx += cols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final float[] flattenArray(float[][] arr)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		float[] flatArr = new float[rows * cols];
		for (int i = 0; i < rows; i++, idx += cols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final int[] flattenArray(int[][] arr)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		int[] flatArr = new int[rows * cols];
		for (int i = 0; i < rows; i++, idx += cols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final long[] flattenArray(long[][] arr)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		long[] flatArr = new long[rows * cols];
		for (int i = 0; i < rows; i++, idx += cols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final double[] flattenArray(double[][] arr, int assumedNRows, int assumedNCols)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		assert(assumedNRows >= rows && assumedNCols >= cols);
		double[] flatArr = new double[assumedNRows * assumedNCols];
		for (int i = 0; i < rows; i++, idx += assumedNCols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final float[] flattenArray(float[][] arr, int assumedNRows, int assumedNCols)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		assert(assumedNRows >= rows && assumedNCols >= cols);
		float[] flatArr = new float[assumedNRows * assumedNCols];
		for (int i = 0; i < rows; i++, idx += assumedNCols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final int[] flattenArray(int[][] arr, int assumedNRows, int assumedNCols)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		assert(assumedNRows >= rows && assumedNCols >= cols);
		int[] flatArr = new int[assumedNRows * assumedNCols];
		for (int i = 0; i < rows; i++, idx += assumedNCols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final long[] flattenArray(long[][] arr, int assumedNRows, int assumedNCols)
	{
		int
			rows = arr.length,
			cols = arr[0].length,
			idx = 0;
		assert(assumedNRows >= rows && assumedNCols >= cols);
		long[] flatArr = new long[assumedNRows * assumedNCols];
		for (int i = 0; i < rows; i++, idx += assumedNCols)
			System.arraycopy(arr[i], 0, flatArr, idx, cols);
		return flatArr;
	}

	public static final double[][][][] copyArray(double[][][][] array)
	{
		int n1 = array.length;
		double[][][][] result = new double[n1][][][];
		for (int idx1 = 0; idx1 < n1; idx1++)
		{
			double[][][] arr1 = array[idx1];
			int n2 = arr1.length;
			double[][][] result1 = result[idx1] = new double[n2][][];
			for (int idx2 = 0; idx2 < n2; idx2++)
			{
				double[][] arr2 = arr1[idx2];
				int n3 = arr2.length;
				double[][] result2 = result1[idx2] = new double[n3][];
				for (int idx3 = 0; idx3 < n3; idx3++)
				{
					double[] arr3 = arr2[idx3];
					int n4 = arr3.length;
					double[] result3 = result2[idx3] = new double[n4];
					System.arraycopy(arr3, 0, result3, 0, n4);
				}
			}
		}
		return result;
	}

	/**
	 * Make a clone of an array and then reverse the contents.
	 * That is, the last element becomes the first element and
	 * the first element becomes the last element
	 * @param array
	 * @return
	 */
	public static final int[] reverseArray(int[] array)
	{
		int
			n = array.length,
			nMin1 = n - 1,
			result[] = new int[n];
		for (int i = 0; i < n; i++)
			result[i] = array[nMin1 - i];
		return result;
	}

	public static final double[] reverseArray(double[] array)
	{
		int
			n = array.length,
			nMin1 = n - 1;
		double result[] = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = array[nMin1 - i];
		return result;
	}

	public static final <T> T[] reverseArray(T[] array)
	{
		int
			n = array.length,
			nMin1 = n - 1;
		T result[] = (T[]) Array.newInstance(array.getClass().getComponentType(), n);
		for (int i = 0; i < n; i++)
			result[i] = array[nMin1 - i];
		return result;
	}

	public static final void reverseArrayInPlace(double[] array)
	{
		int
			i = 0,
			j = array.length - 1;
		double temp;
		while (i < j)
		{
			temp = array[i]; array[i] = array[j]; array[j] = temp;
			i++; j--;
		}
	}

	public static final void reverseArrayInPlace(int[] array)
	{
		int
			i = 0,
			j = array.length - 1;
		int temp;
		while (i < j)
		{
			temp = array[i]; array[i] = array[j]; array[j] = temp;
			i++; j--;
		}
	}

	public static final <T> void reverseArrayInPlace(T[] array)
	{
		int
			i = 0,
			j = array.length - 1;
		T temp;
		while (i < j)
		{
			temp = array[i]; array[i] = array[j]; array[j] = temp;
			i++; j--;
		}
	}

	/**
	 * Create interval array
	 * @param length
	 * @param interval
	 * @return
	 */
	public static final double[] createIntervalArray(int length, int interval)
	{
		int numIntervals = length / interval + 1;
		double[] xData = new double[numIntervals];
		for (int i = 0; i < numIntervals; i++)
			xData[i] = i * interval;
		return xData;
	}

	// --- Merging QEmbedding.QDataUtils Begin ---
	// That class is ambiguous, so I'll merge it here.
	// I threw out reading / writing to file functions -- RJ 2005/06/04
	public static final int[][] makeSubArray( int fromArr[][], int fromRow, int toRow, int fromCol, int toCol)
	{	  
	  int
	  	numRows = toRow - fromRow + 1,
	  	numCols = toCol - fromCol + 1;
	  
		int subArr[][] = new int[numRows][numCols];
		for (int i = fromRow; i <= toRow; i++)
			System.arraycopy(fromArr[i], fromCol, subArr[i - fromRow], 0, numCols);	// replaced 0, subArr[i - 1], CN 3/26/03
	  // System.arraycopy(array, 1, neighborsArray, 0, rowLength);  
		return subArr;
	}	
	
	public static final float[] intToFloatArray(int[] intValues)
	{
		int numValues = intValues.length;
		float convertArray[] = new float[numValues];
		for (int i = 0; i < intValues.length; i++)
			convertArray[i] = intValues[i];
		return convertArray;
	}

	public static final double[] intToDoubleArray(int[] intValues)
	{
		int numValues = intValues.length;
		double convertArray[] = new double[numValues];
		for (int i = 0; i < intValues.length; i++)
			convertArray[i] = intValues[i];
		return convertArray;
	}

	public static final int[] toIntArray(Collection<Integer> list)
	{
		int[] array = new int[list.size()];
		int idx = 0;
		for (int i: list)
			array[idx++] = i;
		return array;
	}

	public static final int[] toIntArray(double[] arr)
	{
		int[] array = new int[arr.length];
		int idx = 0;
		for (double i: arr)
			array[idx++] = (int) i;
		return array;
	}

	public static final long[] toLongArray(Collection<Long> list)
	{
		long[] array = new long[list.size()];
		int idx = 0;
		for (long i: list)
			array[idx++] = i;
		return array;
	}

	public static final float[] toFloatArray(Collection<Float> list)
	{
		float[] array = new float[list.size()];
		int idx = 0;
		for (float i: list)
			array[idx++] = i;
		return array;
	}

	public static final double[] toArray(Collection<Double> list)
	{
		double[] array = new double[list.size()];
		int idx = 0;
		for (double i: list)
			array[idx++] = i;
		return array;
	}

	public static final double[][] toArray1D(Collection<double[]> list)
	{
		double[][] array = new double[list.size()][];
		int idx = 0;
		for (double[] i: list)
			array[idx++] = i;
		return array;
	}

	public static final double[][][] toArray2D(Collection<double[][]> list)
	{
		double[][][] array = new double[list.size()][][];
		int idx = 0;
		for (double[][] i: list)
			array[idx++] = i;
		return array;
	}

	public static final double[][] toArrayTranspose(Collection<double[]> dataList)
	{
		int
			numData = dataList.size(),
			numColumns = dataList.iterator().next().length,
			i = 0;
		double[][] data = new double[numColumns][numData];
		for (double[] curData: dataList)
		{
			for (int j = 0; j < numColumns; j++)
				data[j][i] = curData[j];
			i++;
		}
		return data;
	}

	public static final Double[][] toObjectArray(double[][] values)
	{
		if (values == null)
			return null;
		int numValues = values.length;
		Double[][] array = new Double[numValues][];
		for (int i = 0; i < numValues; i++)
		{
			double[] valueRow = values[i];
			int numRows = valueRow.length;
			Double[] stringRow = array[i] = new Double[numRows];
			for (int j = 0; j < numRows; j++)
				stringRow[j] = new Double(valueRow[j]);
		}
		return array;
	}

	public static final double[][] toStaggeredArray(double[] arr)
	{
		int
			numData = arr.length,
			n = (int) Math.floor((1 + Math.sqrt(1 + 8 * numData)) / 2);
		double[][] result = new double[n][];
		for (int i = 0; i < n; i++) {
			int idx = i * (i - 1) / 2;
			double[] curResult = result[i] = new double[i];
			for (int j = 0; j < i; j++)
				curResult[j] = arr[idx + j];
		}
		return result;
	}

	public static final float[][] toStaggeredArray(float[] arr)
	{
		int
			numData = arr.length,
			n = (int) Math.floor((1 + Math.sqrt(1 + 8 * numData)) / 2);
		float[][] result = new float[n][];
		for (int i = 0; i < n; i++) {
			int idx = i * (i - 1) / 2;
			float[] curResult = result[i] = new float[i];
			for (int j = 0; j < i; j++)
				curResult[j] = arr[idx + j];
		}
		return result;
	}

	public static final int[][] toStaggeredArray(int[] arr)
	{
		int
			numData = arr.length,
			n = (int) Math.floor((1 + Math.sqrt(1 + 8 * numData)) / 2);
		int[][] result = new int[n][];
		for (int i = 0; i < n; i++) {
			int idx = i * (i - 1) / 2;
			int[] curResult = result[i] = new int[i];
			for (int j = 0; j < i; j++)
				curResult[j] = arr[idx + j];
		}
		return result;
	}

	public static final float[][] to2DArray(float[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		float result[][] = new float[nrow][ncol];
		for (int i = 0; i < nElems; i += ncol)
			System.arraycopy(array, i, result[i / ncol], 0, ncol);
		return result;
	}

	public static final double[][] to2DArray(double[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		double result[][] = new double[nrow][ncol];
		for (int i = 0; i < nElems; i += ncol)
			System.arraycopy(array, i, result[i / ncol], 0, ncol);
		return result;
	}

	public static final int[][] to2DArray(int[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		int result[][] = new int[nrow][ncol];
		for (int i = 0; i < nElems; i += ncol)
			System.arraycopy(array, i, result[i / ncol], 0, ncol);
		return result;
	}

	public static final long[][] to2DArray(long[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		long result[][] = new long[nrow][ncol];
		for (int i = 0; i < nElems; i += ncol)
			System.arraycopy(array, i, result[i / ncol], 0, ncol);
		return result;
	}

	public static final float[][] to2DArrayTranspose(float[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		float result[][] = new float[nrow][ncol];
		for (int i = 0; i < nElems; i++) {
			int x = i % ncol, y = i / ncol;
			result[x][y] = array[i];
		}
		return result;
	}

	public static final double[][] to2DArrayTranspose(double[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		double result[][] = new double[nrow][ncol];
		for (int i = 0; i < nElems; i++) {
			int x = i % ncol, y = i / ncol;
			result[x][y] = array[i];
		}
		return result;
	}

	public static final int[][] to2DArrayTranspose(int[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		int result[][] = new int[nrow][ncol];
		for (int i = 0; i < nElems; i++) {
			int x = i % ncol, y = i / ncol;
			result[x][y] = array[i];
		}
		return result;
	}

	public static final long[][] to2DArrayTranspose(long[] array, int nrow, int ncol)
	{
		int nElems = array.length;
		long result[][] = new long[nrow][ncol];
		for (int i = 0; i < nElems; i++) {
			int x = i % ncol, y = i / ncol;
			result[x][y] = array[i];
		}
		return result;
	}

	public static final float[][] to2DArray(float[] array, int nrow, int ncol, int nrowTotal, int ncolTotal)
	{
		assert(nrow <= nrowTotal && ncol <= ncolTotal);
		int nElems = ncolTotal * nrow;
		float result[][] = new float[nrow][ncol];
		for (int i = 0; i < nElems; i += ncolTotal)
			System.arraycopy(array, i, result[i / ncolTotal], 0, ncol);
		return result;
	}

	public static final double[][] to2DArray(double[] array, int nrow, int ncol, int nrowTotal, int ncolTotal)
	{
		assert(nrow <= nrowTotal && ncol <= ncolTotal);
		int nElems = ncolTotal * nrow;
		double result[][] = new double[nrow][ncol];
		for (int i = 0; i < nElems; i += ncolTotal)
			System.arraycopy(array, i, result[i / ncolTotal], 0, ncol);
		return result;
	}

	public static final int[][] to2DArray(int[] array, int nrow, int ncol, int nrowTotal, int ncolTotal)
	{
		assert(nrow <= nrowTotal && ncol <= ncolTotal);
		int nElems = ncolTotal * nrow;
		int result[][] = new int[nrow][ncol];
		for (int i = 0; i < nElems; i += ncolTotal)
			System.arraycopy(array, i, result[i / ncolTotal], 0, ncol);
		return result;
	}

	public static final long[][] to2DArray(long[] array, int nrow, int ncol, int nrowTotal, int ncolTotal)
	{
		assert(nrow <= nrowTotal && ncol <= ncolTotal);
		int nElems = ncolTotal * nrow;
		long result[][] = new long[nrow][ncol];
		for (int i = 0; i < nElems; i += ncolTotal)
			System.arraycopy(array, i, result[i / ncolTotal], 0, ncol);
		return result;
	}

	/**
	 * Mimicking R's match function. This function is to return an array with the same length
	 * as that of array x, each element of which is the corresponding index in array y. That is,
	 * for each element in array x, find whether it is also found in array y and if so, return
	 * the corresponding index in array y. If not, return -1.
	 * <P>Example: x = {2, 9, 3, 7, 8, 10, 6, 4, 5, 1, 20, 21, 22}<br>
	 * y = {6, 2, 1, 5, 7, 4, 9, 8, 3, 10}<br>
	 * Output: {1, 6, 8, 4, 7, 9, 0, 5, 3, 2, -1, -1, -1}
	 * @param x
	 * @param y
	 * @return
	 */
	public static final int[] match(int[] x, int[] y) {
		int nx = x.length, ny = y.length, idx[] = new int[ny], matchIdx[] = new int[nx];
		int[] y_copy = new int[ny];
		for (int i = 0; i < ny; i++)
			idx[i] = i;
		System.arraycopy(y, 0, y_copy, 0, ny);
		QDataUtils.sort2DByRow(new int[][] {y_copy, idx}, 0);
		for (int i = 0; i < nx; i++) {
			int idx_y = Arrays.binarySearch(y_copy, x[i]);
			matchIdx[i] = idx_y < 0 ? -1 : idx[idx_y];
		}
		return matchIdx;
	}

	/**
	 * See match(int[], int[])
	 * @param x
	 * @param y
	 * @return
	 */
	public static final int[] match(String[] x, String[] y) {
		int nx = x.length, ny = y.length, idx[] = new int[ny], matchIdx[] = new int[nx];
		String[] y_copy = new String[ny];
		for (int i = 0; i < ny; i++)
			idx[i] = i;
		System.arraycopy(y, 0, y_copy, 0, ny);
		QDataUtils.sort(y_copy, idx);
		for (int i = 0; i < nx; i++) {
			int idx_y = Arrays.binarySearch(y_copy, x[i]);
			matchIdx[i] = idx_y < 0 ? -1 : idx[idx_y];
		}
		return matchIdx;
	}

	/**
	 * An implementation of Radix sort. A modified source from http://en.literateprograms.org/Radix_sort_(Java)?oldid=12461.
	 * <pre>
	 *  // These numbers are for integer arrays, in ms. Apple Intel core i7 2.8GHz, 16GB RAM.
	 * 	// N = 10^3: Radix = 1, Quick = 1
	 * 	// N = 10^4: Radix = 6, Quick = 14
	 * 	// N = 10^5: Radix = 29, Quick = 43
	 * 	// N = 10^6: Radix = 60, Quick = 116
	 * 	// N = 10^7: Radix = 441, Quick = 1203
	 * 	// N = 10^8: Radix = 4090, Quick = 13379, My 2D sort = 14550
	 * 	// N = 5*10^8: Radix = 23416, Quick = 72447, My 2D sort = 78962
	 * </pre>
	 * @param a
	 */
	public static final void radixSort(int[] a)
	{
		int
			n = a.length,
			// Roby's modification: Optimum number of bins
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2))))),
			b[] = new int[n],
			b_orig[] = b;
		for (int mask = ~(-1 << numBins), rshift = 0; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(a[p] & mask) >>> rshift]; // RJ's fix to handle negatives
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (a[p] & mask) >>> rshift; // RJ's fix to handle negatives
				--cntarray[key];
				b[cntarray[key]] = a[p];
			}
			int[] temp = b; b = a; a = temp;
		}
		if (a == b_orig)
			System.arraycopy(a, 0, b, 0, n);
		// Roby's modification to handle negatives -- begin
		int numNegs = 0;
		// Negatives are always placed at the end of the array in the correct order.
		// So, scan to find the point
		for (int i = n - 1; a[i] < 0; i--, numNegs++) ;
		if (numNegs > 0)
		{
			System.arraycopy(a, n - numNegs, b, 0, numNegs);
			System.arraycopy(a, 0, a, numNegs, n - numNegs);
			System.arraycopy(b, 0, a, 0, numNegs);
		}
		// Roby's modification to handle negatives -- end
	}

	public static final void radixSort(long[] a)
	{
		int
			n = a.length,
			// Roby's modification: Optimum number of bins
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2)))));
		long
			b[] = new long[n],
			b_orig[] = b;
		for (long mask = ~(-1L << numBins), rshift = 0; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(int) ((a[p] & mask) >>> rshift)]; // RJ's fix to handle negatives
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (int) ((a[p] & mask) >>> rshift); // RJ's fix to handle negatives
				--cntarray[key];
				b[cntarray[key]] = a[p];
			}
			long[] temp = b; b = a; a = temp;
		}
		if (a == b_orig)
			System.arraycopy(a, 0, b, 0, n);
		// Roby's modification to handle negatives -- begin
		int numNegs = 0;
		// Negatives are always placed at the end of the array in the correct order.
		// So, scan to find the point
		for (int i = n - 1; a[i] < 0; i--, numNegs++) ;
		if (numNegs > 0)
		{
			System.arraycopy(a, n - numNegs, b, 0, numNegs);
			System.arraycopy(a, 0, a, numNegs, n - numNegs);
			System.arraycopy(b, 0, a, 0, numNegs);
		}
		// Roby's modification to handle negatives -- end
	}

	public static final void radixSort(float[] data)
	{
		int
			n = data.length,
			a[] = new int[n],
			b[] = new int[n],
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2)))));
		for (int i = 0; i < n; i++)
			a[i] = Float.floatToIntBits(data[i]);
		for (int mask = ~(-1 << numBins), rshift = 0; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(a[p] & mask) >>> rshift];
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (a[p] & mask) >>> rshift;
				--cntarray[key];
				b[cntarray[key]] = a[p];
			}
			int[] temp = b; b = a; a = temp;
		}

		int numNegs = 0;
		// Negatives are always placed at the end of the array in the reverse order.
		// So, scan to find the point
		if (a[n - 1] < 0)
			for (int i = n - 1; a[i] < 0; i--)
				data[numNegs++] = Float.intBitsToFloat(a[i]);
		if (numNegs > 0)
			for (int i = numNegs; i < n; i++)
				data[i] = Float.intBitsToFloat(a[i - numNegs]);
		else
			for (int i = 0; i < n; i++)
				data[i] = Float.intBitsToFloat(a[i]);
	}

	public static final void radixSort(double[] data)
	{
		// N=5*10^8: Radix = 9919, Quick = 16971 (ms).
		int
			n = data.length,
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2)))));
		long[]
			a = new long[n],
			b = new long[n];
		for (int i = 0; i < n; i++)
			a[i] = Double.doubleToLongBits(data[i]);
		for (long mask = ~(-1L << numBins), rshift = 0L; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(int) ((a[p] & mask) >>> rshift)];
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (int) ((a[p] & mask) >>> rshift);
				--cntarray[key];
				b[cntarray[key]] = a[p];
			}
			long[] temp = b; b = a; a = temp;
		}

		int numNegs = 0;
		// Negatives are always placed at the end of the array in the reverse order.
		// So, scan to find the point
		if (a[n - 1] < 0)
			for (int i = n - 1; a[i] < 0; i--)
				data[numNegs++] = Double.longBitsToDouble(a[i]);
		if (numNegs > 0)
			for (int i = numNegs; i < n; i++)
				data[i] = Double.longBitsToDouble(a[i - numNegs]);
		else
			for (int i = 0; i < n; i++)
				data[i] = Double.longBitsToDouble(a[i]);
	}

	public static final void radixSort(double[] data, int[] idx)
	{
		// N=5*10^8: Radix = 9919, Quick = 16971 (ms).
		int
			n = data.length,
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2)))));
		long[]
			a = new long[n],
			b = new long[n];
		int[]
			idx_a = new int[n],
			idx_b = new int[n];
		for (int i = 0; i < n; i++)
		{
			a[i] = Double.doubleToLongBits(data[i]);
			idx_a[i] = idx[i];
		}
		for (long mask = ~(-1L << numBins), rshift = 0L; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(int) ((a[p] & mask) >>> rshift)];
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (int) ((a[p] & mask) >>> rshift);
				--cntarray[key];
				b[cntarray[key]] = a[p];
				idx_b[cntarray[key]] = idx_a[p];
			}
			long[] temp = b; b = a; a = temp;
			int[] temp_idx = idx_b; idx_b = idx_a; idx_a = temp_idx;
		}

		int numNegs = 0;
		// Negatives are always placed at the end of the array in the reverse order.
		// So, scan to find the point
		if (a[n - 1] < 0)
			for (int i = n - 1; a[i] < 0; i--) {
				data[numNegs] = Double.longBitsToDouble(a[i]);
				idx[numNegs++] = idx_a[i];
			}
		if (numNegs > 0)
			for (int i = numNegs; i < n; i++) {
				data[i] = Double.longBitsToDouble(a[i - numNegs]);
				idx[i] = idx_a[i - numNegs];
			}
		else
			for (int i = 0; i < n; i++) {
				data[i] = Double.longBitsToDouble(a[i]);
				idx[i] = idx_a[i];
			}
	}

	/**
	 * Sort two dimensional array of data
	 * @param data the 2D data arranged row x column
	 * @param columnNo the column number the data is sorted upon
	 */
	public static final void sort2D(double[][] data, int columnNo)
	{
		int maxData = data.length - 1;
		quickSort(data, columnNo, 0, maxData);
		insertionSort(data, columnNo, 0, maxData);
	}

	/**
	 * Sort two dimensional array of data
	 * @param data the 2D data arranged row x column
	 * @param columnNo the column number the data is sorted upon
	 */
	public static final void sort2D(int[][] data, int columnNo)
	{
		int maxData = data.length - 1;
		quickSort(data, columnNo, 0, maxData);
		insertionSort(data, columnNo, 0, maxData);
	}

	public static final void sort2DByRow(double[][] data, int rowNo)
	{
		int maxData = data[0].length - 1;
		quickSortByRow(data, rowNo, 0, maxData);
		insertionSortByRow(data, rowNo, 0, maxData);
	}

	public static final void sort2DByRow(int[][] data, int rowNo)
	{
		int maxData = data[0].length - 1;
		quickSortByRow(data, rowNo, 0, maxData);
		insertionSortByRow(data, rowNo, 0, maxData);
	}

	/**
	 * Sort by the value array. All arrays must have the same length.
	 * @param xData
	 * @param yData
	 * @param value
	 */
	public static final void sortByValue(int[] xData, int[] yData, double[] value)
	{
		int length = value.length;
		assert (xData.length == length && yData.length == length);
		quickSort(xData, yData, value, 0, length - 1);
		insertionSort(xData, yData, value, 0, length - 1);
	}

	public static final void sort(String[] data, int[] idx)
	{
		int maxData = data.length - 1;
		quickSort(data, idx, 0, maxData);
		insertionSort(data, idx, 0, maxData);
	}

	public static final <T> void sort(int[] data, T[] idx)
	{
		int maxData = data.length - 1;
		quickSort(data, idx, 0, maxData);
		insertionSort(data, idx, 0, maxData);
	}

	private static final void swap(double[][] data, int i, int j)
	{
		double[] rowData = data[i];
		data[i] = data[j];
		data[j] = rowData;
	}

	private static final void swap(int[][] data, int i, int j)
	{
		int[] rowData = data[i];
		data[i] = data[j];
		data[j] = rowData;
	}

	private static final void swap(String[] data, int[] idx, int i, int j)
	{
		String s = data[i]; data[i] = data[j]; data[j] = s;
		int p = idx[i]; idx[i] = idx[j]; idx[j] = p;
	}

	private static final <T> void swap(int[] data, T[] idx, int i, int j)
	{
		int s = data[i]; data[i] = data[j]; data[j] = s;
		T p = idx[i]; idx[i] = idx[j]; idx[j] = p;
	}

	private static final void swapRow(double[][] data, int i, int j)
	{
		int numRows = data.length;
		for (int rowNo = 0; rowNo < numRows; rowNo++)
		{
			double
				curRow[] = data[rowNo],
				temp = curRow[i];
			curRow[i] = curRow[j];
			curRow[j] = temp;
		}
	}

	private static final void swapRow(int[][] data, int i, int j)
	{
		int numRows = data.length;
		for (int rowNo = 0; rowNo < numRows; rowNo++)
		{
			int
				curRow[] = data[rowNo],
				temp = curRow[i];
			curRow[i] = curRow[j];
			curRow[j] = temp;
		}
	}

	private static final void swap(int[] xData, int[] yData, double[] value, int i, int j)
	{
		int temp = xData[i];
		xData[i] = xData[j];
		xData[j] = temp;
		temp = yData[i];
		yData[i] = yData[j];
		yData[j] = temp;
		double temp2 = value[i];
		value[i] = value[j];
		value[j] = temp2;
	}

	/**
	 * Three median quick sort
	 *
	 * @param data       the data
	 * @param rowNo   the column index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void quickSortByRow(double[][] data, int rowNo, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (data[rowNo][lo] > data[rowNo][midpoint])
			swapRow(data, lo, midpoint);
		if (data[rowNo][lo] > data[rowNo][hi])
			swapRow(data, lo, hi);
		if (data[rowNo][midpoint] > data[rowNo][hi])
			swapRow(data, midpoint, hi);

		right = hi - 1;
		swapRow(data, midpoint, right);
		left = lo;
		double
			dataRow[] = data[rowNo],
			pivotElement = dataRow[right];
		for(;;)
		{
			while(dataRow[++left] < pivotElement){}
			while(dataRow[--right] > pivotElement){}
			if (right < left)
				break;
			swapRow(data, left, right);
		}
		swapRow(data, left, hi-1);
		quickSortByRow(data, rowNo, lo, right);
		quickSortByRow(data, rowNo, left+1, hi);
	}

	/**
	 * Three median quick sort
	 *
	 * @param data       the data
	 * @param rowNo   the column index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void quickSortByRow(int[][] data, int rowNo, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (data[rowNo][lo] > data[rowNo][midpoint])
			swapRow(data, lo, midpoint);
		if (data[rowNo][lo] > data[rowNo][hi])
			swapRow(data, lo, hi);
		if (data[rowNo][midpoint] > data[rowNo][hi])
			swapRow(data, midpoint, hi);

		right = hi - 1;
		swapRow(data, midpoint, right);
		left = lo;
		int
			dataRow[] = data[rowNo],
			pivotElement = dataRow[right];
		for(;;)
		{
			while(dataRow[++left] < pivotElement){}
			while(dataRow[--right] > pivotElement){}
			if (right < left)
				break;
			swapRow(data, left, right);
		}
		swapRow(data, left, hi-1);
		quickSortByRow(data, rowNo, lo, right);
		quickSortByRow(data, rowNo, left+1, hi);
	}

	/**
	 * Three median quick sort
	 *
	 * @param data       the data
	 * @param columnNo   the column index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void quickSort(double[][] data, int columnNo, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (data[lo][columnNo] > data[midpoint][columnNo])
			swap(data, lo, midpoint);
		if (data[lo][columnNo] > data[hi][columnNo])
			swap(data, lo, hi);
		if (data[midpoint][columnNo] > data[hi][columnNo])
			swap(data, midpoint, hi);

		right = hi - 1;
		swap(data, midpoint, right);
		left = lo;
		double pivotElement = data[right][columnNo];
		for(;;)
		{
			while(data[++left][columnNo] < pivotElement){}
			while(data[--right][columnNo] > pivotElement){}
			if (right < left)
				break;
			swap(data, left, right);
		}
		swap(data, left, hi-1);
		quickSort(data, columnNo, lo, right);
		quickSort(data, columnNo, left+1, hi);
	}

	/**
	 * Three median quick sort
	 *
	 * @param data       the data
	 * @param columnNo   the column index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void quickSort(int[][] data, int columnNo, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (data[lo][columnNo] > data[midpoint][columnNo])
			swap(data, lo, midpoint);
		if (data[lo][columnNo] > data[hi][columnNo])
			swap(data, lo, hi);
		if (data[midpoint][columnNo] > data[hi][columnNo])
			swap(data, midpoint, hi);

		right = hi - 1;
		swap(data, midpoint, right);
		left = lo;
		int pivotElement = data[right][columnNo];
		for(;;)
		{
			while(data[++left][columnNo] < pivotElement){}
			while(data[--right][columnNo] > pivotElement){}
			if (right < left)
				break;
			swap(data, left, right);
		}
		swap(data, left, hi-1);
	    quickSort(data, columnNo, lo, right);
	    quickSort(data, columnNo, left+1, hi);
	}

	private static final void quickSort(String[] data, int[] idx, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (data[lo].compareTo(data[midpoint]) > 0)
			swap(data, idx, lo, midpoint);
		if (data[lo].compareTo(data[hi]) > 0)
			swap(data, idx, lo, hi);
		if (data[midpoint].compareTo(data[hi]) > 0)
			swap(data, idx, midpoint, hi);

		right = hi - 1;
		swap(data, idx, midpoint, right);
		left = lo;
		String pivotElement = data[right];
		for(;;)
		{
			while(data[++left].compareTo(pivotElement) < 0){}
			while(data[--right].compareTo(pivotElement) > 0){}
			if (right < left)
				break;
			swap(data, idx, left, right);
		}
		swap(data, idx, left, hi-1);
	    quickSort(data, idx, lo, right);
	    quickSort(data, idx, left+1, hi);
	}

	private static final <T> void quickSort(int[] data, T[] idx, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (data[lo] > data[midpoint])
			swap(data, idx, lo, midpoint);
		if (data[lo] > data[hi])
			swap(data, idx, lo, hi);
		if (data[midpoint] > data[hi])
			swap(data, idx, midpoint, hi);

		right = hi - 1;
		swap(data, idx, midpoint, right);
		left = lo;
		int pivotElement = data[right];
		for(;;)
		{
			while(data[++left] < pivotElement){}
			while(data[--right] > pivotElement){}
			if (right < left)
				break;
			swap(data, idx, left, right);
		}
		swap(data, idx, left, hi-1);
	    quickSort(data, idx, lo, right);
	    quickSort(data, idx, left+1, hi);
	}

	private static final void quickSort(int[] xData, int[] yData, double[] value, int lo, int hi)
	{
		int left, right, midpoint;

		if ((hi - lo) <= kInsertionTreshold)
			return;

		midpoint = (hi + lo)/2;
		if (value[lo] > value[midpoint])
			swap(xData, yData, value, lo, midpoint);
		if (value[lo] > value[hi])
			swap(xData, yData, value, lo, hi);
		if (value[midpoint] > value[hi])
			swap(xData, yData, value, midpoint, hi);

		right = hi - 1;
		swap(xData, yData, value, midpoint, right);
		left = lo;
		double pivotElement = value[right];
		for(;;)
		{
			while(value[++left] < pivotElement){}
			while(value[--right] > pivotElement){}
			if (right < left)
				break;
			swap(xData, yData, value, left, right);
		}
		swap(xData, yData, value, left, hi-1);
	    quickSort(xData, yData, value, lo, right);
	    quickSort(xData, yData, value, left+1, hi);
	}

	/**
	 * Insertion Sort
	 * 
	 * @param data       the data
	 * @param columnNo   the column index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void insertionSort(double[][] data, int columnNo, int lo, int hi)
	{
		double[] row_i, row_j;
		double element_i, element_j;
		int j;
		for (int i = 1; i <= hi; i++)
		{
			row_i = data[i];
			element_i = row_i[columnNo];
			for (j = i - 1; j >= lo; j--)
			{
				row_j = data[j];
				element_j = row_j[columnNo];
				if (element_j <= element_i)
					break;
				data[j+1] = row_j; // a[j+1] = a[j]
			}
			data[j+1] = row_i;
		}
	}

	/**
	 * Insertion Sort
	 * 
	 * @param data       the data
	 * @param columnNo   the column index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void insertionSort(int[][] data, int columnNo, int lo, int hi)
	{
		int[] row_i, row_j;
		int element_i, element_j;
		int j;
		for (int i = 1; i <= hi; i++)
		{
			row_i = data[i];
			element_i = row_i[columnNo];
			for (j = i - 1; j >= lo; j--)
			{
				row_j = data[j];
				element_j = row_j[columnNo];
				if (element_j <= element_i)
					break;
				data[j+1] = row_j; // a[j+1] = a[j]
			}
			data[j+1] = row_i;
		}
	}

	private static final void insertionSort(String[] data, int[] idx, int lo, int hi)
	{
		int idx_i, idx_j;
		String element_i, element_j;
		int j;
		for (int i = 1; i <= hi; i++)
		{
			element_i = data[i]; idx_i = idx[i];
			for (j = i - 1; j >= lo; j--)
			{
				element_j = data[j]; idx_j = idx[j];
				if (element_j.compareTo(element_i) <= 0)
					break;
				int jp1 = j+1;
				data[jp1] = element_j; // a[j+1] = a[j]
				idx[jp1] = idx_j;
			}
			j++;
			data[j] = element_i;
			idx[j] = idx_i;
		}
	}

	private static final <T> void insertionSort(int[] data, T[] idx, int lo, int hi)
	{
		T idx_i, idx_j;
		int element_i, element_j;
		int j;
		for (int i = 1; i <= hi; i++)
		{
			element_i = data[i]; idx_i = idx[i];
			for (j = i - 1; j >= lo; j--)
			{
				element_j = data[j]; idx_j = idx[j];
				if (element_j <= element_i)
					break;
				int jp1 = j+1;
				data[jp1] = element_j; // a[j+1] = a[j]
				idx[jp1] = idx_j;
			}
			j++;
			data[j] = element_i;
			idx[j] = idx_i;
		}
	}

	private static final void insertionSort(int[] xData, int[] yData, double[] value, int lo, int hi)
	{
		double value_i, value_j;
		int x_i, y_i, x_j, y_j;
		int j;
		for (int i = 1; i <= hi; i++)
		{
			value_i = value[i];
			x_i = xData[i];
			y_i = yData[i];
			for (j = i - 1; j >= lo; j--)
			{
				value_j = value[j];
				x_j = xData[j];
				y_j = yData[j];
				if (value_j <= value_i)
					break;
				int jp1 = j+1;
				value[jp1] = value_j; // a[j+1] = a[j]
				xData[jp1] = x_j;
				yData[jp1] = y_j;
			}
			j++;
			value[j] = value_i;
			xData[j] = x_i;
			yData[j] = y_i;
		}
	}

	/**
	 * Insertion Sort
	 * 
	 * @param data       the data
	 * @param rowNo   the row index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void insertionSortByRow(double[][] data, int rowNo, int lo, int hi)
	{
		double element_i, element_j;
		int j, jp1, numRows = data.length;
		double[] col_i = new double[numRows];
		for (int i = 1; i <= hi; i++)
		{
			for (int r = 0; r < numRows; r++)
				col_i[r] = data[r][i];
			element_i = data[rowNo][i];
			for (j = i - 1; j >= lo; j--)
			{
				element_j = data[rowNo][j];
				if (element_j <= element_i)
					break;
				jp1 = j+1;
				for (int r = 0; r < numRows; r++)
					data[r][jp1] = data[r][j]; // a[j+1] = a[j]
			}
			jp1 = j+1;
			for (int r = 0; r < numRows; r++)
				data[r][jp1] = col_i[r];
		}
	}

	/**
	 * Insertion Sort
	 * 
	 * @param data       the data
	 * @param rowNo   the row index of which the table is compared against
	 * @param lo         left boundary of array partition
	 * @param hi         right boundary of array partition
	 */
	private static final void insertionSortByRow(int[][] data, int rowNo, int lo, int hi)
	{
		int element_i, element_j;
		int j, jp1, numRows = data.length;
		int[] col_i = new int[numRows];
		for (int i = 1; i <= hi; i++)
		{
			for (int r = 0; r < numRows; r++)
				col_i[r] = data[r][i];
			element_i = data[rowNo][i];
			for (j = i - 1; j >= lo; j--)
			{
				element_j = data[rowNo][j];
				if (element_j <= element_i)
					break;
				jp1 = j+1;
				for (int r = 0; r < numRows; r++)
					data[r][jp1] = data[r][j]; // a[j+1] = a[j]
			}
			jp1 = j+1;
			for (int r = 0; r < numRows; r++)
				data[r][jp1] = col_i[r];
		}
	}

	/**
	 * Find peaks in one dimensional data that exceed threshold. However, between
	 * peaks must be separated by at least a specified distance
	 * @return The indices that denote the peak
	 */
	public static final int[] findPeaks(double[] xData, double[] yData,
		double threshold, double minDistance)
	{
		int
			dataLength = yData.length,
			dataLengthMin1 = dataLength - 1;
		assert dataLength == xData.length;
		// This is quite a bad work around that we cannot use inner class :(
		// What we want is a sortable set (that's what TreeSet does), but
		// TreeSet wants the object to be Comparable or we have to put our own
		// comparator. So, here we go. -- RJ 2005/06/30
		Set<double[]> peaks = new TreeSet<double[]>(new Comparator<double[]>()
			{
				// The sort is descending, so, we put t2 - t1 instead of t1 - t2
				public int compare(double[] t1, double[] t2)
				{	return (int) (t2[1] - t1[1]); }
				@Override
				public boolean equals(Object o)
				{	return this == o; }
			});
	
		// Check boundary cases
		// Exactly one data
		double val = yData[0];
		if (dataLength == 1)
		{
			// If it's above threshold, just return the index
			// Otherwise, say it's not found
			return (val >= threshold) ? new int[] { 0 } : null;
		}
	
		// Check the far left value
		if (val >= threshold && val > yData[1])
			peaks.add(new double[] {xData[0], val, 0});
		// Check the far right value
		val = yData[dataLengthMin1];
		if (val >= threshold && val > yData[dataLengthMin1-1])
			peaks.add(new double[] {xData[dataLengthMin1], val, dataLengthMin1});
		// Iterate linearly to check whether one point actually is
		// bigger than the threshold and its left and right neighbors
		// If it is so, consider it as a peak, this is the first pass
		for (int i = 1; i < dataLengthMin1; i++)
		{
			val = yData[i];
			if (val >= threshold && val >= yData[i - 1] && val >= yData[i + 1])
				peaks.add(new double[] {xData[i], val, i});
		}
	
		// If nothing's found inside the threshold, say null
		if (peaks.size() == 0)
			return null;
	
		// N^2 clustering routine work is here
		Set<Integer> truePeaks = new TreeSet<Integer>();
		for(double[] peak: peaks)
		{
			boolean isFeasiblePeak = true;
			double distance = peak[0];
			for (int truePeak: truePeaks)
			{
				// Is there any peaks within +/- minDistance already chosen?
				// If so, consider it as a duplicate
				if (Math.abs(xData[truePeak] - distance) < minDistance)
				{
					isFeasiblePeak = false;
					break;
				}
			}
			if (isFeasiblePeak)
				truePeaks.add((int) peak[2]);
		}
		return toIntArray(truePeaks);
	}

	/**
	 * <P>Find duplicate indices, assuming arr is sorted.
	 * It's similar to R's <tt>match</tt> function, but it
	 * further assumes that <tt>arr</tt> is sorted.<br>
	 * Example:
	 * <br>Input: {1, 2, 2, 2, 3, 3, 4, 5, 6, 6}
	 * <br>Output: {0, 1, 1, 1, 2, 2, 3, 4, 5, 5}
	 * 
	 * <P>So, the number of unique elements can be obtained
	 * by (result[arr.length-1] + 1).
	 * @param arr
	 * @return
	 */
	public static final int[] findDuplicateIndices(double[] arr)
	{
		int
			n = arr.length,
			idx = 0,
			result[] = new int[n];
		double lastX = arr[0];
		for (int i = 1; i < n; i++)
		{
			double curX = arr[i];
			if (curX != lastX)
			{
				idx++;
				lastX = curX;
			}
			result[i] = idx;
		}
		return result;
	}

	public static final float[] convertArray(double[] arr)
	{
		int length = arr.length;
		float[] result = new float[length];
		for (int i = 0; i < length; i++)
			result[i] = (float) arr[i];
		return result;
	}

	public static final double[] convertArray(float[] arr)
	{
		int length = arr.length;
		double[] result = new double[length];
		for (int i = 0; i < length; i++)
			result[i] = arr[i];
		return result;
	}

	/**
	 * Create a double dimensional array clone safely.
	 * @return
	 */
	public static final double[][] safeClone(double[][] data)
	{
		if (data == null)
			return null;
		int dim = data.length;
		double[][] newData = new double[dim][];
		for (int i = 0; i < dim; i++)
		{
			if (data[i] == null)
				continue;
			int d = data[i].length;
			newData[i] = new double[d];
			System.arraycopy(data[i], 0, newData[i], 0, d);
		}
		return newData;
	}

	/**
	 * Create a triple dimensional array clone safely.
	 * @return
	 */
	public static final double[][][] safeClone(double[][][] data)
	{
		if (data == null)
			return null;
		int dim = data.length;
		double[][][] newData = new double[dim][][];
		for (int i = 0; i < dim; i++)
		{
			if (data[i] == null)
				continue;
			int d = data[i].length;
			newData[i] = new double[d][];
			double[][]
				di = data[i],
				ndi = newData[i];
			for (int j = 0; j < d; j++)
			{
				if (di[j] == null)
					continue;
				int dd = di[j].length;
				ndi[j] = new double[dd];
 				System.arraycopy(di[j], 0, ndi[j], 0, dd);
			}
		}
		return newData;
	}

	/**
	 * Cloning using arraycopy. Up to 5x faster
	 * @param data
	 * @return
	 */
	public static final double[] safeClone(double[] data)
	{
		if (data == null)
			return null;
		int n = data.length;
		double[] d = new double[n];
		System.arraycopy(data, 0, d, 0, n);
		return d;
	}

	public static final boolean[] safeClone(boolean[] data)
	{
		if (data == null)
			return null;
		int n = data.length;
		boolean[] d = new boolean[n];
		System.arraycopy(data, 0, d, 0, n);
		return d;
	}

	/**
	 * Cloning using arraycopy. Up to 5x faster
	 * @param data
	 * @return
	 */
	public static final int[] safeClone(int[] data)
	{
		if (data == null)
			return null;
		int
			n = data.length,
			d[] = new int[n];
		System.arraycopy(data, 0, d, 0, n);
		return d;
	}

	public static final int[][] safeClone(int[][] data)
	{
		if (data == null)
			return null;
		int dim = data.length;
		int[][] newData = new int[dim][];
		for (int i = 0; i < dim; i++)
		{
			if (data[i] == null)
				continue;
			int d = data[i].length;
			newData[i] = new int[d];
			System.arraycopy(data[i], 0, newData[i], 0, d);
		}
		return newData;
	}

	/**
	 * Add a value to a table with key type K and List of type V
	 * @param <K, V> (key, value) pair
	 * @param table
	 * @param key
	 * @param value
	 */
	public static final <K, V> void addToTableList(Map<K, List<V>> table, K key, V value)
	{
		List<V> resultList = table.get(key);
		if (resultList == null)
		{
			resultList = new ArrayList<V>();
			table.put(key, resultList);
		}
		resultList.add(value);
	}

	/**
	 * Add a value to a table with key type K and Set of type V
	 * @param <K, V> (key, value) pair
	 * @param table
	 * @param key
	 * @param value
	 */
	public static final <K, V> void addToTableSet(Map<K, Set<V>> table, K key, V value)
	{
		Set<V> resultList = table.get(key);
		if (resultList == null)
		{
			resultList = new TreeSet<V>();
			table.put(key, resultList);
		}
		resultList.add(value);
	}

	/**
	 * Uses .equals to check for equality
	 * @param <T>
	 * @param s
	 * @param arrays
	 * @return
	 */
	public static final <T> int findIndex(T s, T[] arrays)
	{
		int n = arrays.length;
		for (int i = 0; i < n; i++)
			if (arrays[i].equals(s))
				return i;
		return -1;
	}

	/**
	 * Uses == to check for equality
	 * @param <T>
	 * @param s
	 * @param arrays
	 * @return
	 */
	public static final <T> int findIndexSimple(T s, T[] arrays)
	{
		int n = arrays.length;
		for (int i = 0; i < n; i++)
			if (arrays[i] == s)
				return i;
		return -1;
	}

	public static final int findIndex(int s, int[] arrays)
	{
		int n = arrays.length;
		for (int i = 0; i < n; i++)
			if (arrays[i] == s)
				return i;
		return -1;
	}

	public static final int findIndex(double s, double[] arrays)
	{
		int n = arrays.length;
		for (int i = 0; i < n; i++)
			if (arrays[i] == s)
				return i;
		return -1;
	}

	public static final <T> void appendToList(T[] array, List<T> list)
	{
		if (array != null)
			for (T t : array)
				list.add(t);
	}

	public static final <T> Set<T> convertToSet(T[] array)
	{
		Set<T> set = new HashSet<T>();
		if (array != null)
			for (T data: array)
				set.add(data);
		return set;
	}

	public static final Set<Integer> convertToSet(int[] array)
	{
		Set<Integer> set = new HashSet<Integer>();
		if (array != null)
			for (int data: array)
				set.add(new Integer(data));
		return set;
	}

	public static final Set<Long> convertToSet(long[] array)
	{
		Set<Long> set = new HashSet<Long>();
		if (array != null)
			for (long data: array)
				set.add(new Long(data));
		return set;
	}

	public static final Set<Double> convertToSet(double[] array)
	{
		Set<Double> set = new HashSet<Double>();
		if (array != null)
			for (double data: array)
				set.add(new Double(data));
		return set;
	}

	public static final Set<Float> convertToSet(float[] array)
	{
		Set<Float> set = new HashSet<Float>();
		if (array != null)
			for (float data: array)
				set.add(new Float(data));
		return set;
	}

	public static final <T> List<T> convertToList(T[] array)
	{
		List<T> list = new ArrayList<T>();
		if (array != null)
			for (T data: array)
				list.add(data);
		return list;
	}

	public static final List<Integer> convertToList(int[] array)
	{
		List<Integer> list = new ArrayList<Integer>();
		if (array != null)
			for (int data: array)
				list.add(data);
		return list;
	}

	public static final List<Double> convertToList(double[] array)
	{
		List<Double> list = new ArrayList<Double>();
		if (array != null)
			for (double data: array)
				list.add(data);
		return list;
	}

	public static final List<Float> convertToList(float[] array)
	{
		List<Float> list = new ArrayList<Float>();
		if (array != null)
			for (float data: array)
				list.add(data);
		return list;
	}

	public static final List<Long> convertToList(long[] array)
	{
		List<Long> list = new ArrayList<Long>();
		if (array != null)
			for (long data: array)
				list.add(data);
		return list;
	}

	public static final double[][] filterValuesByColumns(double[][] values, BitSet missingBits)
	{
		int
			numRows = values.length,
			numCols = values[0].length;

		if (missingBits == null)
		{
			double[][] result = new double[numRows][numCols];
			for (int i = 0; i < numRows; i++)
				System.arraycopy(values[i], 0, result[i], 0, numCols);
			return result;
		}

		int numMissing = missingBits.cardinality();
		double[][] result = new double[numRows][numCols - numMissing];

		for (int i = 0; i < numRows; i++)
		{
			int idx = 0;
			for (int j = 0; j < numCols; j++)
			{
				if (missingBits.get(j))
					continue;
				result[i][idx] = values[i][j];
				idx++;
			}
		}
		return result;
	}

	/**
	 * Filter a matrix by rows given a missing bits
	 * @param values
	 * @param missingBits
	 * @return
	 */
	public static final double[][] filterValuesByRows(double[][] values, BitSet missingBits)
	{
		int
			numRows = values.length,
			numCols = values[0].length;

		if (missingBits == null)
		{
			double[][] result = new double[numRows][numCols];
			for (int i = 0; i < numRows; i++)
				System.arraycopy(values[i], 0, result[i], 0, numCols);
			return result;
		}

		int numMissing = missingBits.cardinality();
		double[][] result = new double[numRows - numMissing][numCols];

		for (int i = 0, idx = 0; i < numRows; i++)
		{
			if (missingBits.get(i))
				continue;
			System.arraycopy(values[i], 0, result[idx], 0, numCols);
			idx++;
		}
		return result;
	}

	public static final double[] detectAndFilterMissingValues(double[] values)
	{	return detectAndFilterMissingValues(values, null); }

	/**
	 * Detect and filter out missing values.
	 * 
	 * @param values unfiltered values
	 * @param bitset An empty bitset, which upon the return of this method,
	 * indicate with which individuals is missing. If bitset is null, then
	 * the information of missing individuals are not returned.
	 * 
	 * @return Filtered trait values.
	 */
	public static final double[] detectAndFilterMissingValues(double[] values, BitSet bitset)
	{
		int length = values.length;
		if (bitset == null)
			bitset = new BitSet(length);
		//else
		//	bitset.clear();
		// If bitset is not null, don't clear it. Just in case that the caller wants
		// to "OR" it with the previous values.

		for (int i = 0; i < length; i++)
			if (values[i] == kUndefinedValue)
				bitset.set(i);

		double[] trimmedValues = new double[length - bitset.cardinality()];
		int index = 0;

		for (int i = 0; i < length; i++)
		{
			if (values[i] == kUndefinedValue)
				continue;
			trimmedValues[index] = values[i];
			index++;
		}
		return trimmedValues;
	}

	public static final double[][] detectAndFilterMissingValues(double[][] values)
	{	return detectAndFilterMissingValues(values, null, false); }

	/**
	 * Detect and filter out missing multivariate values. It assumes that the array of values has
	 * common length (numTraits x numIndividuals). The missing data is OR-ed.
	 * 
	 * @param values
	 * @param bitset
	 */
	public static final double[][] detectAndFilterMissingValues(double[][] values, BitSet bitset)
	{	return detectAndFilterMissingValues(values, bitset, false); }

	public static final double[][] detectAndFilterMissingValues(double[][] values, BitSet bitset, boolean transposeResult)
	{
		int
			numTraits = values.length,
			numIndividuals = values[0].length;
		if (bitset == null)
			bitset = new BitSet(numIndividuals);
		//else
		//	bitset.clear();
		// If bitset is not null, don't clear it. Just in case that the caller wants
		// to "OR" it with the previous values.

		for (int i = 0; i < numTraits; i++)
			for (int j = 0; j < numIndividuals; j++)
				if (values[i][j] == kUndefinedValue)
					bitset.set(j);

		double[][] trimmedTraitValues =
			transposeResult ?
				new double[numIndividuals - bitset.cardinality()][numTraits] :
				new double[numTraits][numIndividuals - bitset.cardinality()];
		int index = 0;

		if (transposeResult)
		{
			for (int j = 0; j < numIndividuals; j++)
				if (!bitset.get(j))
				{
					for (int i = 0; i < numTraits; i++)
						trimmedTraitValues[index][i] = values[i][j];
					index++;
				}
		}
		else
		{
			for (int j = 0; j < numIndividuals; j++)
				if (!bitset.get(j))
				{
					for (int i = 0; i < numTraits; i++)
						trimmedTraitValues[i][index] = values[i][j];
					index++;
				}
		}
		return trimmedTraitValues;
	}

	public static final double[][] detectAndFilterMissingValuesTranspose(double[][] values, BitSet bitset)
	{
		int
			numTraits = values[0].length,
			numIndividuals = values.length;
		if (bitset == null)
			bitset = new BitSet(numIndividuals);
		//else
		//	bitset.clear();
		// If bitset is not null, don't clear it. Just in case that the caller wants
		// to "OR" it with the previous values.

		for (int i = 0; i < numTraits; i++)
			for (int j = 0; j < numIndividuals; j++)
				if (values[j][i] == kUndefinedValue)
					bitset.set(j);

		double[][] trimmedTraitValues = new double[numIndividuals - bitset.cardinality()][numTraits];
		for (int j = 0, index = 0; j < numIndividuals; j++)
			if (!bitset.get(j))
			{
				double[] v = values[j];
				for (int i = 0; i < numTraits; i++)
					trimmedTraitValues[index][i] = v[i];
				index++;
			}
		return trimmedTraitValues;
	}

	/**
	 * Translate IDs into values given a translation table
	 * @param ids array of IDs
	 * @param table Translation table
	 * @param keyColumn The key column of the translation table against which the IDs are matched
	 * @param valueColumn The value column of the translation table
	 * @return Some array elements may be null upon return due to mismatch
	 */
	public static final String[] translate(String[] ids, String[][] table, int keyColumn, int valueColumn)
	{
		int
			numAnnot = table.length,
			numIDs = ids.length;
		Map<String, String> map = new HashMap<String, String>(numAnnot);
		String[] xlat = new String[numIDs];
		for (int i = 0; i < numAnnot; i++)
			map.put(table[i][keyColumn], table[i][valueColumn]);
		for (int i = 0; i < numIDs; i++)
			xlat[i] = map.get(ids[i]);
		return xlat;
	}

	public static final void main(String[] args)
	{
		/*
		// Test match
		int[] x1 = new int[] {2, 9, 3, 7, 8, 10, 6, 4, 5, 1, 20, 21, 22},
			x2 = new int[] {6, 2, 1, 5, 7, 4, 9, 8, 3, 10};
		System.out.println(QStringUtils.toString(match(x1, x2)));
		String[] example1 = new String[] { "Egyptian", "Indian", "American", "Chinese", "Filipino", "Indonesian",
			"Korean", "Japanese", "English", "German", "French", "Italian", "Irish", "Vietnamese", "Australian",
			"Mongolian", "Russian", "Austrian", "Brazilian"};
		String[] example2 = new String[] {
			"Korean", "Filipino", "Irish", "Indian", "Egyptian", "Australian", "Italian", "Japanese", "Brazilian", "Chinese",
			"English", "Indonesian", "German", "Austrian", "Vietnamese", "American", "Russian", "French", "Mongolian", "Dumb"
		};
		System.out.println(QStringUtils.toString(match(example1, example2)));
		//*/
	}
}
