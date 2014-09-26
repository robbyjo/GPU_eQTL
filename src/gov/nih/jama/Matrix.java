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
package gov.nih.jama;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StreamTokenizer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Locale;

import gov.nih.utils.matrix.QMatrixUtils;

import static java.lang.Math.hypot;

/**
 * jama = Java Matrix class.
 * 
 * <P>The Java Matrix Class provides the fundamental operations of numerical
 * linear algebra.  Various constructors create Matrices from two dimensional
 * arrays of double precision floating point numbers.  Various "gets" and
 * "sets" provide access to submatrices and matrix elements.  Several methods 
 * implement basic matrix arithmetic, including matrix addition and
 * multiplication, matrix norms, and element-by-element array operations.
 * Methods for reading and printing matrices are also included.  All the
 * operations in this version of the Matrix Class involve real matrices.
 * Complex matrices may be handled in a future version.
 * 
 * <P>Five fundamental matrix decompositions, which consist of pairs or triples
 * of matrices, permutation vectors, and the like, produce results in five
 * decomposition classes.  These decompositions are accessed by the Matrix
 * class to compute solutions of simultaneous linear equations, determinants,
 * inverses and other matrix functions.  The five decompositions are:
 * 
 * <P><UL>
 * <LI>Cholesky Decomposition of symmetric, positive definite matrices.
 * <LI>LU Decomposition of rectangular matrices.
 * <LI>QR Decomposition of rectangular matrices.
 * <LI>Singular Value Decomposition of rectangular matrices.
 * <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
 * </UL>
 * 
 * <DL>
 * <DT><B>Example of use:</B></DT>
 * <P>
 * <DD>Solve a linear system matrix x = b and compute the residual norm, ||b - matrix x||.
 * <P><PRE>
      double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
      Matrix mMatrix = new Matrix(vals);
      Matrix b = Matrix.random(3,1);
      Matrix x = mMatrix.solve(b);
      Matrix r = mMatrix.times(x).minus(b);
      double rnorm = r.normInf();
   </PRE></DD>
 * </DL>
 * @author The MathWorks, Inc. and the National Institute of Standards and Technology.
 * 
 * @author Roby Joehanes
 * 
 * <P>Roby Joehanes has modified this file for speed and consistency
 */
public class Matrix implements Cloneable {

	public static final String
		sSubmatrixIndices = "Submatrix indices", // $NON-NLS-1$
		sArrayMultipleOfM = "Array length must be a multiple of m.", // $NON-NLS-1$
		sRowIsTooShort = "Row %d is too short", // $NON-NLS-1$
		sRowIsTooLong = "Row %d is too long", // $NON-NLS-1$
		sMatrixRowMustAgree = "Matrix row dimensions must agree.", // $NON-NLS-1$
		sMatrixNotSymPosDef = "Matrix is not symmetric positive definite.", // $NON-NLS-1$
		sMatrixNotSquare = "Matrix must be square."; // $NON-NLS-1$

	/* ------------------------
	 * Class variables
	 * ------------------------ */

	/**
	 * Array for internal storage of elements.
	 * @serial internal array storage.
	 */
	private double[][] mMatrix;

	/**
	 * Row and column dimensions.
	 * @serial row dimension.
	 * @serial column dimension.
	 */
	private int mRows, mColumns;

	/* ------------------------
	 * Constructors
	 * ------------------------ */
	/**
	 * Construct an m-by-n matrix of zeros. 
	 * @param m    Number of rows.
	 * @param n    Number of colums.
	 */
	public Matrix (int m, int n) {
		mRows = m;
		mColumns = n;
		mMatrix = new double[m][n];
	}

	/**
	 * Construct an mRow-by-mColumn constant matrix.
	 * @param m    Number of rows.
	 * @param n    Number of colums.
	 * @param s    Fill the matrix with this scalar value.
	 */
	public Matrix (int m, int n, double s) {
		mRows = m;
		mColumns = n;
		mMatrix = new double[m][n];
		for (int i = 0; i < m; i++)
			Arrays.fill(mMatrix[i], s);
	}

	/**
	 * Construct a matrix from a 2-D array. (RJ's patch: Removed dimension checking)
	 * @param A    Two-dimensional array of doubles.
	 * @exception  IllegalArgumentException All rows must have the same length
	 * @see        #constructWithCopy
	 */
	public Matrix (double[][] A) {
		mRows = A.length;
		mColumns = A[0].length;
		// RJ-patch: Eliminate checking to speed up (2004/11/30)
		//for (int i = 0; i < mRow; i++) {
		//   if (mMatrix[i].length != mColumn) {
		//      throw new IllegalArgumentException("All rows must have the same length.");
		//   }
		//}
		mMatrix = A;
	}

	/**
	 * Construct a matrix quickly without checking arguments.
	 * @param A    Two-dimensional array of doubles.
	 * @param m    Number of rows.
	 * @param n    Number of colums.
	 */
	public Matrix (double[][] A, int m, int n) {
		mMatrix = A;
		mRows = m;
		mColumns = n;
	}

	/**
	 * Construct a matrix from a one-dimensional packed array
	 * @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
	 * @param m    Number of rows.
	 * @exception  IllegalArgumentException Array length must be a multiple of m.
	 */
	public Matrix (double vals[], int m) {
		mRows = m;
		mColumns = (m != 0 ? vals.length/m : 0);
		if (m*mColumns != vals.length) {
			throw new IllegalArgumentException(sArrayMultipleOfM);
		}
		mMatrix = new double[m][mColumns];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < mColumns; j++) {
				mMatrix[i][j] = vals[i+j*m];
			}
		}
	}

	/* ------------------------
	 * Public Methods
	 * ------------------------ */

	/**
	 * Construct a matrix from a copy of a 2-D array.
	 * @param A    Two-dimensional array of doubles.
	 * @exception  IllegalArgumentException All rows must have the same length
	 */
	public static final Matrix constructWithCopy(double[][] A) {
		int m = A.length;
		int n = A[0].length;
		Matrix X = new Matrix(m,n);
		double[][] C = X.mMatrix;
		for (int i = 0; i < m; i++)
			System.arraycopy(A[i], 0, C[i], 0, n);

		return X;
	}

	/**
	 * Make a deep copy of a matrix
	 */
	public Matrix copy () {
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.mMatrix;
		for (int i = 0; i < mRows; i++)
			System.arraycopy(mMatrix[i], 0, C[i], 0, mColumns);
		return X;
	}

	/**
	 * Clone the Matrix object.
	 */
	@Override
	public Object clone ()
	{	return copy();}

	/**
	 * Access the internal two-dimensional array.
	 * @return     Pointer to the two-dimensional array of matrix elements.
	 */
	public double[][] getArray ()
	{	return mMatrix; }

	/**
	 * Copy the internal two-dimensional array.
	 * @return     Two-dimensional array copy of matrix elements.
	 */
	public double[][] getArrayCopy () {
		double[][] C = new double[mRows][mColumns];
		for (int i = 0; i < mRows; i++)
			System.arraycopy(mMatrix[i], 0, C[i], 0, mColumns);
		return C;
	}

	/**
	 * Make a one-dimensional column packed copy of the internal array.
	 * @return     Matrix elements packed in a one-dimensional array by columns.
	 */
	public double[] getColumnPackedCopy () {
		double[] vals = new double[mRows*mColumns];
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				vals[i+j*mRows] = mMatrix[i][j];
		return vals;
	}

	/**
	 * Make a one-dimensional row packed copy of the internal array.
	 * @return     Matrix elements packed in a one-dimensional array by rows.
	 */
	public double[] getRowPackedCopy () {
		double[] vals = new double[mRows*mColumns];
		for (int i = 0; i < mRows; i++)
			System.arraycopy(mMatrix[i], 0, vals, i*mColumns, mColumns);
		return vals;
	}

	/**
	 * Get row dimension.
	 * @return     mRow, the number of rows.
	 */
	public int getRowDimension ()
	{	return mRows; }

	/**
	 * Get column dimension.
	 * @return     mColumn, the number of columns.
	 */
	public int getColumnDimension ()
	{	return mColumns; }

	/**
	 * Get a single element.
	 * @param i    Row index.
	 * @param j    Column index.
	 * @return     Matrix(i,j)
	 * @exception  ArrayIndexOutOfBoundsException
	 */
	public double get (int i, int j)
	{	return mMatrix[i][j]; }

	/**
	 * Get a submatrix.
	 * @param i0   Initial row index
	 * @param i1   Final row index
	 * @param j0   Initial column index
	 * @param j1   Final column index
	 * @return     Matrix(i0:i1,j0:j1)
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public Matrix getMatrix (int i0, int i1, int j0, int j1) {
		int colLen = j1-j0+1;
		Matrix X = new Matrix(i1-i0+1, colLen);
		double[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++)
				System.arraycopy(mMatrix[i], j0, B[i-i0], 0, colLen); // RJ's speedup
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
		return X;
	}

	/**
	 * Get a submatrix.
	 * @param r    Array of row indices.
	 * @param c    Array of column indices.
	 * @return     Matrix(r(:),c(:))
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public Matrix getMatrix (int[] r, int[] c) {
		int
			rl = r.length,
			cl = c.length;  // RJ's speedup
		Matrix X = new Matrix(rl,cl);
		double[][] B = X.getArray();
		try {
			for (int i = 0; i < rl; i++) {
				double[] Bi = B[i]; // RJ's speedup
				for (int j = 0; j < cl; j++)
					Bi[j] = mMatrix[r[i]][c[j]];
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
		return X;
	}

	/**
	 * Get a submatrix.
	 * @param i0   Initial row index
	 * @param i1   Final row index
	 * @param c    Array of column indices.
	 * @return     Matrix(i0:i1,c(:))
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public Matrix getMatrix (int i0, int i1, int[] c) {
		int cl = c.length;
		Matrix X = new Matrix(i1-i0+1,cl);
		double[][] B = X.getArray();
		try {
			for (int i = i0; i <= i1; i++) {
				double[]
					Mi = mMatrix[i], // RJ's speedup
					Bi = B[i - i0];
				for (int j = 0; j < cl; j++) {
					Bi[j] = Mi[c[j]];
				}
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
		return X;
	}

	/**
	 * Get a submatrix.
	 * @param r    Array of row indices.
	 * @param i0   Initial column index
	 * @param i1   Final column index
	 * @return     Matrix(r(:),j0:j1)
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public Matrix getMatrix (int[] r, int j0, int j1) {
		int colLen = j1 - j0 + 1; // RJ's speedup
		Matrix X = new Matrix(r.length, colLen);
		double[][] B = X.getArray();
		try {
			for (int i = 0; i < r.length; i++)
				System.arraycopy(mMatrix[r[i]], j0, B[i], 0, colLen); // RJ's speedup
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
		return X;
	}

	/**
	 * Set a single element.
	 * @param i    Row index.
	 * @param j    Column index.
	 * @param s    Matrix(i,j).
	 * @exception  ArrayIndexOutOfBoundsException
	 */
	public void set (int i, int j, double s)
	{	mMatrix[i][j] = s; }

	/**
	 * Set a submatrix.
	 * @param i0   Initial row index
	 * @param i1   Final row index
	 * @param j0   Initial column index
	 * @param j1   Final column index
	 * @param X    Matrix(i0:i1,j0:j1)
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public void setMatrix (int i0, int i1, int j0, int j1, Matrix X) {
		try {
			int colLen = j1 - j0 + 1;
			double[][] x = X.getArray();
			for (int i = i0; i <= i1; i++)
				System.arraycopy(x[i-i0], 0, mMatrix[i], j0, colLen); // RJ's speed up
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
	}

	/**
	 * Set a submatrix.
	 * @param r    Array of row indices.
	 * @param c    Array of column indices.
	 * @param X    Matrix(r(:),c(:))
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public void setMatrix (int[] r, int[] c, Matrix X) {
		try {
			double[][] x = X.getArray();
			for (int i = 0; i < r.length; i++) {
				double[]
					xi = x[i], // RJ's speed up
					mi = mMatrix[r[i]];
				for (int j = 0; j < c.length; j++)
					mi[c[j]] = xi[j];
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
	}

	/**
	 * Set a submatrix.
	 * @param r    Array of row indices.
	 * @param j0   Initial column index
	 * @param j1   Final column index
	 * @param X    Matrix(r(:),j0:j1)
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public void setMatrix (int[] r, int j0, int j1, Matrix X) {
		try {
			int colLen = j1 - j0 + 1; // RJ's speed up
			double[][] x = X.getArray();
			for (int i = 0; i < r.length; i++)
				System.arraycopy(x[i], 0, mMatrix[r[i]], j0, colLen); // RJ's speed up
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
	}

	/**
	 * Set a submatrix.
	 * @param i0   Initial row index
	 * @param i1   Final row index
	 * @param c    Array of column indices.
	 * @param X    Matrix(i0:i1,c(:))
	 * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	 */
	public void setMatrix (int i0, int i1, int[] c, Matrix X) {
		try {
			int cl = c.length;
			double[][] x = X.getArray();
			for (int i = i0; i <= i1; i++)
			{
				double[]
					mi = mMatrix[i], // RJ's speed up
					xi = x[i-i0];
				for (int j = 0; j < cl; j++)
					mi[c[j]] = xi[j];
			}
		} catch(ArrayIndexOutOfBoundsException e) {
			throw new ArrayIndexOutOfBoundsException(sSubmatrixIndices);
		}
	}

	/**
	 * Matrix transpose.
	 * @return    Matrix'
	 */
	public Matrix transpose () {
		Matrix X = new Matrix(mColumns,mRows);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[j][i] = mMatrix[i][j];
		return X;
	}

	/**
	 * One norm
	 * @return    maximum column sum.
	 */
	public double norm1 () {
		double f = 0;
		for (int j = 0; j < mColumns; j++) {
			double s = 0;
			for (int i = 0; i < mRows; i++)
				s += Math.abs(mMatrix[i][j]);
			f = Math.max(f,s);
		}
		return f;
	}

	/**
	 * Two norm
	 * @return    maximum singular value.
	 */
	public double norm2 ()
	{	return (new SingularValueDecomposition(this, false, false).norm2()); }

	/**
	 * Infinity norm
	 * @return    maximum row sum.
	 */
	public double normInf () {
		double f = 0;
		for (int i = 0; i < mRows; i++) {
			double s = 0;
			double[] m_i = mMatrix[i];
			for (int j = 0; j < mColumns; j++)
				s += Math.abs(m_i[j]);
			f = Math.max(f,s);
		}
		return f;
	}

	/**
	 * Frobenius norm
	 * @return    sqrt of sum of squares of all elements.
	 */
	public double normF () {
		double f = 0;
		for (int i = 0; i < mRows; i++) {
			double[] m_i = mMatrix[i];
			for (int j = 0; j < mColumns; j++)
				f = hypot(f,m_i[j]);
		}
		return f;
	}

	/**
	 * Unary minus
	 * @return    -Matrix
	 */
	public Matrix uminus () {
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = -mMatrix[i][j];
		return X;
	}

	/**
	 * C = Matrix + B
	 * @param B    another matrix
	 * @return     Matrix + B
	 */
	public Matrix plus (Matrix B) {
		//checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = mMatrix[i][j] + B.mMatrix[i][j];
		return X;
	}

	/**
	 * C = Matrix + B'
	 * @param B    another matrix
	 * @return     Matrix + B'
	 */
	public Matrix plusTranspose (Matrix B) {
		//checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = mMatrix[i][j] + B.mMatrix[j][i];
		return X;
	}

	/**
	 * Matrix = Matrix + B
	 * @param B    another matrix
	 * @return     Matrix + B
	 */
	public Matrix plusEquals (Matrix B) {
		//checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		double[][] b = B.mMatrix;
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = mMatrix[i][j] + b[i][j];
		return this;
	}

	/**
	 * Matrix = Matrix + B'
	 * @param B    another matrix
	 * @return     Matrix + B'
	 */
	public Matrix plusTransposeEquals (Matrix B) {
		//checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		double[][] b = B.mMatrix;
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = mMatrix[i][j] + b[j][i];
		return this;
	}

	/**
	 * C = Matrix - B
	 * @param B    another matrix
	 * @return     Matrix - B
	 */
	public Matrix minus (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][]
			b = B.mMatrix,
			C = X.mMatrix;
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = mMatrix[i][j] - b[i][j];
		return X;
	}

	/**
	 * C = Matrix - B'
	 * @param B    another matrix
	 * @return     Matrix - B'
	 */
	public Matrix minusTranspose (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][]
			b = B.mMatrix,
			C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = mMatrix[i][j] - b[j][i];
		return X;
	}

	/** Matrix = Matrix - B
   @param B    another matrix
   @return     Matrix - B
	 */
	public Matrix minusEquals (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		double[][] b = B.mMatrix;
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = mMatrix[i][j] - b[i][j];
		return this;
	}

	/**
	 * Matrix = Matrix - B'
	 * @param B    another matrix
	 * @return     Matrix - B'
	 */
	public Matrix minusTransposeEquals (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = mMatrix[i][j] - B.mMatrix[j][i];
		return this;
	}

	/**
	 * Element-by-element multiplication, C = Matrix.*B
	 * @param B    another matrix
	 * @return     Matrix.*B
	 */
	public Matrix arrayTimes (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = mMatrix[i][j] * B.mMatrix[i][j];
		return X;
	}

	/**
	 * Element-by-element multiplication in place, Matrix = Matrix.*B
	 * @param B    another matrix
	 * @return     Matrix.*B
	 */
	public Matrix arrayTimesEquals (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = mMatrix[i][j] * B.mMatrix[i][j];
		return this;
	}

	/**
	 * Element-by-element right division, C = Matrix./B
	 * @param B    another matrix
	 * @return     Matrix./B
	 */
	public Matrix arrayRightDivide (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = mMatrix[i][j] / B.mMatrix[i][j];
		return X;
	}

	/**
	 * Element-by-element right division in place, Matrix = Matrix./B
	 * @param B    another matrix
	 * @return     Matrix./B
	 */
	public Matrix arrayRightDivideEquals (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = mMatrix[i][j] / B.mMatrix[i][j];
		return this;
	}

	/**
	 * Element-by-element left division, C = Matrix.\B
	 * @param B    another matrix
	 * @return     Matrix.\B
	 */
	public Matrix arrayLeftDivide (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = B.mMatrix[i][j] / mMatrix[i][j];
		return X;
	}

	/**
	 * Element-by-element left division in place, Matrix = Matrix.\B
	 * @param B    another matrix
	 * @return     Matrix.\B
	 */
	public Matrix arrayLeftDivideEquals (Matrix B) {
		// checkMatrixDimensions(B); // RJ-patch: Eliminate checking -- 2004/11/30
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = B.mMatrix[i][j] / mMatrix[i][j];
		return this;
	}

	/**
	 * Multiply a matrix by a scalar, C = s*Matrix
	 * @param s    scalar
	 * @return     s*Matrix
	 */
	public Matrix times (double s) {
		Matrix X = new Matrix(mRows,mColumns);
		double[][] C = X.getArray();
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				C[i][j] = s*mMatrix[i][j];
		return X;
	}

	/**
	 * Multiply a matrix by a scalar in place, Matrix = s*Matrix
	 * @param s    scalar
	 * @return     replace Matrix by s*Matrix
	 */
	public Matrix timesEquals (double s) {
		for (int i = 0; i < mRows; i++)
			for (int j = 0; j < mColumns; j++)
				mMatrix[i][j] = s*mMatrix[i][j];
		return this;
	}

	/**
	 * Linear algebraic matrix multiplication, Matrix * B
	 * @param B    another matrix
	 * @return     Matrix product, Matrix * B
	 * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */
	public Matrix times (Matrix B) {
		//if (B.mRow != mColumn) { // RJ-patch: Eliminate checking -- 2004/11/30
		//   throw new IllegalArgumentException("Matrix inner dimensions must agree.");
		//}
		Matrix X = new Matrix(mRows,B.mColumns);
		double[][] C = X.getArray();
		double[] Bcolj = new double[mColumns];
		for (int j = 0; j < B.mColumns; j++) {
			for (int k = 0; k < mColumns; k++) {
				Bcolj[k] = B.mMatrix[k][j];
			}
			for (int i = 0; i < mRows; i++) {
				double[] Arowi = mMatrix[i];
				double s = 0;
				for (int k = 0; k < mColumns; k++) {
					s += Arowi[k]*Bcolj[k];
				}
				C[i][j] = s;
			}
		}
		return X;
	}

	/**
	 * Linear algebraic matrix multiplication, Matrix * B, in parallel fashion.
	 * @author Roby Joehanes
	 * @param B    another matrix
	 * @return     Matrix product, Matrix * B
	 * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */
	public Matrix parallelTimes (Matrix B) {
		return new Matrix(QMatrixUtils.calculateXYParallel(mMatrix, B.getArray()));
	}

	/**
	 * Linear algebraic matrix multiplication, Matrix * B'
	 * @author Roby Joehanes
	 * @param B    another matrix
	 * @return     Matrix product, Matrix * B'
	 * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */
	public Matrix timesTranspose (Matrix B) {
		return new Matrix(QMatrixUtils.calculateXYt(mMatrix, B.getArray()));
	}

	/**
	 * Linear algebraic matrix multiplication, Matrix' * B
	 * @author Roby Joehanes
	 * @param B    another matrix
	 * @return     Matrix product, Matrix' * B
	 * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */
	public Matrix transposeTimes (Matrix B) {
		return new Matrix(QMatrixUtils.calculateXtY(mMatrix, B.getArray()));
	}

	/**
	 * Linear algebraic matrix multiplication, Matrix * B', in parallel fashion
	 * @author Roby Joehanes
	 * @param B    another matrix
	 * @return     Matrix product, Matrix * B'
	 * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */
	public Matrix parallelTimesTranspose(Matrix B) {
		return new Matrix(QMatrixUtils.calculateXYtParallel(mMatrix, B.getArray()));
	}

	/**
	 * Linear algebraic matrix multiplication, Matrix' * B, in parallel fashion
	 * @author Roby Joehanes
	 * @param B    another matrix
	 * @return     Matrix product, Matrix' * B
	 * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	 */
	public Matrix parallelTransposeTimes (Matrix B) {
		return new Matrix(QMatrixUtils.calculateXtYParallel(mMatrix, B.getArray()));
	}

	/**
	 * LU Decomposition
	 * @return     LUDecomposition
	 * @see LUDecomposition
	 */
	public LUDecomposition lu ()
	{	return new LUDecomposition(this); }

	/**
	 * QR Decomposition
	 * @return     QRDecomposition
	 * @see QRDecomposition
	 */
	public QRDecomposition qr ()
	{	return new QRDecomposition(this); }

	/**
	 * Cholesky Decomposition
	 * @return     CholeskyDecomposition
	 * @see CholeskyDecomposition
	 */
	public CholeskyDecomposition chol ()
	{	return new CholeskyDecomposition(this); }

	/**
	 * Singular Value Decomposition
	 * @return     SingularValueDecomposition
	 * @see SingularValueDecomposition
	 */
	public SingularValueDecomposition svd ()
	{	return new SingularValueDecomposition(this); }

	/**
	 * Eigenvalue Decomposition
	 * @return     EigenvalueDecomposition
	 * @see EigenvalueDecomposition
	 */
	public EigenvalueDecomposition eig ()
	{	return new EigenvalueDecomposition(this); }

	/**
	 * Solve Matrix*X = B
	 * @param B    right hand side
	 * @return     solution if mMatrix is square, least squares solution otherwise
	 */
	public Matrix solve (Matrix B) {
		return mRows == mColumns ? new LUDecomposition(this).solve(B) :
			new QRDecomposition(this).solve(B);
	}

	/**
	 * Solve X*mMatrix = B, which is also mMatrix'*X' = B'
	 * @param B    right hand side
	 * @return     solution if mMatrix is square, least squares solution otherwise.
	 */
	public Matrix solveTranspose (Matrix B)
	{	return transpose().solve(B.transpose()); }

	/**
	 * Matrix inverse or pseudoinverse
	 * @return     inverse(mMatrix) if mMatrix is square, pseudoinverse otherwise.
	 */
	public Matrix inverse () {
		return mRows == mColumns ? new LUDecomposition(this).inverse() :
			new QRDecomposition(this).inverse();
	}

	/**
	 * Matrix determinant
	 * @return     determinant
	 */
	public double det ()
	{	return new LUDecomposition(this).det(); }

	/**
	 * Matrix rank
	 * @return     effective numerical rank, obtained from SVD.
	 */
	public int rank ()
	{	return new SingularValueDecomposition(this, false, false).rank(); }

	/**
	 * Matrix condition (2 norm)
	 * @return     ratio of largest to smallest singular value.
	 */
	public double cond ()
	{	return new SingularValueDecomposition(this, false, false).cond(); }

	/**
	 * Matrix trace.
	 * @return     sum of the diagonal elements.
	 */
	public double trace()
	{
		double t = 0;
		int len = Math.min(mRows,mColumns);
		for (int i = 0; i < len; i++)
			t += mMatrix[i][i];
		return t;
	}

	/**
	 * Generate identity matrix
	 * @param m    Number of rows.
	 * @param n    Number of columns.
	 * @return     An mRow-by-mColumn matrix with ones on the diagonal and zeros elsewhere.
	 */
	public static final Matrix identity (int m, int n)
	{
		Matrix A = new Matrix(m,n);
		double[][] X = A.getArray();
		for (int i = 0; i < m; i++)
			X[i][i] = 1.0;
		return A;
	}


	/**
	 * Print the matrix to stdout.   Line the elements up in columns
	 * with a Fortran-like 'Fw.d' style format.
	 * @param w    Column width.
	 * @param d    Number of digits after the decimal.
	 */
	public void print (int w, int d)
	{	print(new PrintWriter(System.out,true),w,d); }

	// Random method is struck out by RJ

	/**
	 * Print the matrix to the output stream.   Line the elements up in
	 * columns with a Fortran-like 'Fw.d' style format.
	 * @param output Output stream.
	 * @param w      Column width.
	 * @param d      Number of digits after the decimal.
	 */
	public void print (PrintWriter output, int w, int d) {
		DecimalFormat format = new DecimalFormat();
		format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
		format.setMinimumIntegerDigits(1);
		format.setMaximumFractionDigits(d);
		format.setMinimumFractionDigits(d);
		format.setGroupingUsed(false);
		print(output,format,w+2);
	}

	/**
	 * Print the matrix to stdout.  Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters.
	 * Note that is the matrix is to be read back in, you probably will want
	 * to use a NumberFormat that is set to US Locale.
	 * @param format mMatrix  Formatting object for individual elements.
	 * @param width     Field width for each column.
	 * @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */
	public void print (NumberFormat format, int width)
	{	print(new PrintWriter(System.out,true),format,width); }

	// DecimalFormat is a little disappointing coming from Fortran or C's printf.
	// Since it doesn't pad on the left, the elements will come out different
	// widths.  Consequently, we'll pass the desired column width in as an
	// argument and do the extra padding ourselves.

	/**
	 * Print the matrix to the output stream.  Line the elements up in columns.
	 * Use the format object, and right justify within columns of width
	 * characters.
	 * Note that is the matrix is to be read back in, you probably will want
	 * to use a NumberFormat that is set to US Locale.
	 * @param output the output stream.
	 * @param format mMatrix formatting object to format the matrix elements 
	 * @param width  Column width.
	 * @see java.text.DecimalFormat#setDecimalFormatSymbols
	 */
	public void print (PrintWriter output, NumberFormat format, int width) {
		output.println();  // start on new line.
		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < mColumns; j++) {
				String s = format.format(mMatrix[i][j]); // format the number
				int padding = Math.max(1,width-s.length()); // At _least_ 1 space
				for (int k = 0; k < padding; k++)
					output.print(' ');
				output.print(s);
			}
			output.println();
		}
		output.println();   // end with blank line.
	}

	// Added some toString method -- RJ 2005/06/28
	@Override
	public String toString()
	{
		java.io.StringWriter wr = new java.io.StringWriter();
		print(new java.io.PrintWriter(wr), 6, 4);
		return wr.toString();
	}

	/**
	 * Read a matrix from a stream.  The format is the same the print method,
	 * so printed matrices can be read back in (provided they were printed using
	 * US Locale).  Elements are separated by whitespace, all the elements for
	 * each row appear on a single line, the last row is followed by a blank line.
	 * @param input the input stream.
	 */
	public static final Matrix read (BufferedReader input) throws java.io.IOException {
		StreamTokenizer tokenizer= new StreamTokenizer(input);

		// Although StreamTokenizer will parse numbers, it doesn't recognize
		// scientific notation (E or D); however, Double.valueOf does.
		// The strategy here is to disable StreamTokenizer's number parsing.
		// We'll only get whitespace delimited words, EOL's and EOF's.
		// These words should all be numbers, for Double.valueOf to parse.

		tokenizer.resetSyntax();
		tokenizer.wordChars(0,255);
		tokenizer.whitespaceChars(0, ' ');
		tokenizer.eolIsSignificant(true);
		java.util.Vector<Double> v = new java.util.Vector<Double>();

		// Ignore initial empty lines
		while (tokenizer.nextToken() == StreamTokenizer.TT_EOL){}
		if (tokenizer.ttype == StreamTokenizer.TT_EOF)
			throw new java.io.IOException("Unexpected EOF on matrix read."); //$NON-NLS-1$
		do {
			v.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
		} while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);

		int n = v.size();  // Now we've got the number of columns!
		double row[] = new double[n];
		for (int j = 0; j < n; j++)  // extract the elements of the 1st row.
			row[j] = v.elementAt(j).doubleValue();
		v.removeAllElements();
		java.util.Vector<double[]> vv = new java.util.Vector<double[]>();
		vv.addElement(row);  // Start storing rows instead of columns.
		while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
			// While non-empty lines
			vv.addElement(row = new double[n]);
			int j = 0;
			do {
				if (j >= n) throw new IOException(String.format(sRowIsTooLong, vv.size()));
				row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
			} while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
			if (j < n) throw new IOException(String.format(sRowIsTooShort, vv.size()));
		}
		int m = vv.size();  // Now we've got the number of rows.
		double[][] A = new double[m][];
		vv.copyInto(A);  // copy the rows out of the vector
		return new Matrix(A);
	}

	/**
	 * Equals method. RJ patch 2004/12/01
	 */
	@Override
	public boolean equals(Object o)
	{
		if (!(o instanceof Matrix))
			return false;
		Matrix m = (Matrix) o;
		if (m.mRows != mRows || m.mColumns != mColumns)
			return false;
		double[][] mm = m.mMatrix;
		for (int i = 0; i < mRows; i++) {
			for (int j = 0; j < mColumns; j++) {
				if (mm[i][j] != mMatrix[i][j])
					return false;
			}
		}
		return true;
	}
}
