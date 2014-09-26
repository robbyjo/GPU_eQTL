/*
 * Roby Joehanes
 * 
 * Copyright 2007 Roby Joehanes.
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

import static gov.nih.jama.Matrix.sMatrixNotSquare;
import static gov.nih.jama.Matrix.sMatrixRowMustAgree;

/**
 * LU Decomposition.
 * 
 * <P>For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
 * unit lower triangular matrix L, an n-by-n upper triangular matrix U, and
 * a permutation vector piv of length m so that A(piv,:) = L*U.
 * If m < n, then L is m-by-m and U is m-by-n.
 * 
 * <P>The LU decompostion with pivoting always exists, even if the matrix is
 * singular, so the constructor will never fail.  The primary use of the
 * LU decomposition is in the solution of square systems of simultaneous
 * linear equations.  This will fail if isNonsingular() returns false.
 */

public class LUDecomposition
{
	/* ------------------------
	 * Class variables
	 * ------------------------ */

	/**
	 * Array for internal storage of decomposition.
	 * @serial internal array storage.
	 */
	private double[][] LU;

	/**
	 * Row and column dimensions, and pivot sign.
	 * @serial column dimension.
	 * @serial row dimension.
	 * @serial pivot sign.
	 */
	private int m, n, pivsign; 

	/**
	 * Internal storage of pivot vector.
	 * @serial pivot vector.
	 */
	private int[] piv;

	private boolean isNonSingular;

	/* ------------------------
	 * Constructor
	 * ------------------------ */

	/**
	 * LU Decomposition
	 * @param  A   Rectangular matrix
	 * @return     Structure to access L, U and piv.
	 */
	public LUDecomposition (Matrix A)
	{
		// Use a "left-looking", dot-product, Crout/Doolittle algorithm.
		LU = A.getArrayCopy();
		m = A.getRowDimension();
		n = A.getColumnDimension();
		isNonSingular = true;
		piv = new int[m];
		for (int i = 0; i < m; i++)
			piv[i] = i;

		pivsign = 1;
		double[] LUrowi;
		double[] LUcolj = new double[m];

		// Outer loop.
		for (int j = 0; j < n; j++)
		{
			// Make a copy of the j-th column to localize references.
			for (int i = 0; i < m; i++)
				LUcolj[i] = LU[i][j];

			// Apply previous transformations.
			for (int i = 0; i < m; i++) {
				LUrowi = LU[i];

				// Most of the time is spent in the following dot product.
				int kmax = Math.min(i,j);
				double s = 0.0;
				for (int k = 0; k < kmax; k++)
					s += LUrowi[k]*LUcolj[k];
				LUrowi[j] = LUcolj[i] -= s;
			}

			// Find pivot and exchange if necessary.
			int p = j;
			for (int i = j+1; i < m; i++)
				if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p]))
					p = i;

			if (p != j) {
				for (int k = 0; k < n; k++) {
					double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
				}
				int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
				pivsign = -pivsign;
			}

			// Compute multipliers.
			if (j < m)
			{
				double nrm = LU[j][j];
				if (nrm != 0.0)
					for (int i = j+1; i < m; i++)
						LU[i][j] /= nrm;
				else
					isNonSingular = false;
			}
		}
	}

	/* ------------------------
	 * Public Methods
	 * ------------------------ */

	/**
	 * Is the matrix nonsingular?
	 * @return     true if U, and hence A, is nonsingular.
	 */
	public boolean isNonsingular()
	{	return isNonSingular; }

	/**
	 * Return lower triangular factor
	 * @return     L
	 */
	public Matrix getL () {
		Matrix X = new Matrix(m,n);
		double[][] L = X.getArray();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (i > j)
					L[i][j] = LU[i][j];
				else if (i == j)
					L[i][j] = 1.0;
				else
					L[i][j] = 0.0;
			}
		}
		return X;
	}

	/**
	 * Return upper triangular factor
	 * @return     U
	 */
	public Matrix getU () {
		Matrix X = new Matrix(n,n);
		double[][] U = X.getArray();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i <= j) {
					U[i][j] = LU[i][j];
				} else {
					U[i][j] = 0.0;
				}
			}
		}
		return X;
	}

	/**
	 * Return pivot permutation vector
	 * @return     piv
	 */
	public int[] getPivot () {
		int[] p = new int[m];
		System.arraycopy(piv, 0, p, 0, m); // RJ's modification
		return p;
	}

	/**
	 * Return pivot permutation vector as a one-dimensional double array
	 * @return     (double) piv
	 */
	public double[] getDoublePivot () {
		double[] vals = new double[m];
		for (int i = 0; i < m; i++)
			vals[i] = piv[i];
		return vals;
	}

	/** Determinant
	 * @return     det(A)
	 * @exception  IllegalArgumentException  Matrix must be square
	 */
	public double det () {
		if (m != n)
			throw new IllegalArgumentException(sMatrixNotSquare);

		double d = pivsign;
		for (int j = 0; j < n; j++)
			d *= LU[j][j];
		return d;
	}

	/**
	 * Solve A*X = B
	 * @param  B   A Matrix with as many rows as A and any number of columns.
	 * @return     X so that L*U*X = B(piv,:)
	 * @exception  IllegalArgumentException Matrix row dimensions must agree.
	 * @exception  RuntimeException  Matrix is singular.
	 */
	public Matrix solve (Matrix B) {
		if (B.getRowDimension() != m)
			throw new IllegalArgumentException(sMatrixRowMustAgree);
		//if (!this.isNonsingular())
		//	throw new RuntimeException("Matrix is singular.");

		// Copy right hand side with pivoting
		int nx = B.getColumnDimension();
		double[][]
			b = B.getArray(),
			X = new double[m][nx];
		for (int i = 0; i < m; i++)
			System.arraycopy(b[piv[i]], 0, X[i], 0, nx);

		// Solve L*Y = B(piv,:)
		for (int k = 0; k < n; k++)
			for (int i = k+1; i < n; i++)
				for (int j = 0; j < nx; j++)
					X[i][j] -= X[k][j]*LU[i][k];

		// Solve U*X = Y;
		for (int k = n-1; k >= 0; k--) {
			for (int j = 0; j < nx; j++)
				X[k][j] /= LU[k][k];

			for (int i = 0; i < k; i++)
				for (int j = 0; j < nx; j++)
					X[i][j] -= X[k][j]*LU[i][k];
		}
		return new Matrix(X);
	}

	/**
	 * Invert a matrix using LU decomposition. Added by Roby Joehanes
	 * @return
	 */
	public Matrix inverse() {
		//if (!this.isNonsingular())
		//	throw new RuntimeException("Matrix is singular.");

		double[][] X = new double[m][m];
		for (int i = 0; i < m; i++)
			X[piv[i]][i] = 1.0;

		// Solve L*Y = B(piv,:)
		for (int k = 0; k < n; k++) {
			double[] Xk = X[k];
			for (int i = k+1; i < n; i++) {
				double[] Xi = X[i];
				double luik = LU[i][k];
				for (int j = 0; j < m; j++)
					Xi[j] -= Xk[j]*luik;
			}
		}

		// Solve U*X = Y;
		for (int k = n-1; k >= 0; k--) {
			double lukk = LU[k][k];
			double[] Xk = X[k];
			for (int j = 0; j < m; j++)
				Xk[j] /= lukk;

			for (int i = 0; i < k; i++) {
				double luik = LU[i][k];
				double[] Xi = X[i];
				for (int j = 0; j < m; j++)
					Xi[j] -= Xk[j]*luik;
			}
		}
		return new Matrix(X);
	}
}
