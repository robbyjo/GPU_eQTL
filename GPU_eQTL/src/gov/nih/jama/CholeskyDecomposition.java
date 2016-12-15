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

import static gov.nih.jama.Matrix.sMatrixNotSymPosDef;
import static gov.nih.jama.Matrix.sMatrixRowMustAgree;


/**
 * Cholesky Decomposition.
 * <P>For a symmetric, positive definite matrix A, the Cholesky decomposition
 * is an lower triangular matrix L so that A = L*L'.
 * 
 * <P>If the matrix is not symmetric or positive definite, the constructor
 * returns a partial decomposition and sets an internal flag that may
 * be queried by the isSPD() method.
 */
public class CholeskyDecomposition {

	/* ------------------------
	 * Class variables
	 * ------------------------ */

	/**
	 * Array for internal storage of decomposition.
	 * @serial internal array storage.
	 */
	private double[][] L;

	/**
	 * Row and column dimension (square matrix).
	 * @serial matrix dimension.
	 */
	private int n;

	/**
	 * Symmetric and positive definite flag.
	 * @serial is symmetric and positive definite flag.
	 */
	private boolean isspd;

	/* ------------------------
	 * Constructor
	 * ------------------------ */

	/**
	 * Cholesky algorithm for symmetric and positive definite matrix.
	 * @param  A   Square, symmetric matrix.
	 * @return     Structure to access L and isspd flag.
	 */
	public CholeskyDecomposition (Matrix Arg) {
		// Initialize.
		double[][] A = Arg.getArray();
		n = A.length;
		L = new double[n][n];
		isspd = (A[0].length == n);
		// Main loop.
		for (int j = 0; j < n; j++) {
			double[] Lrowj = L[j];
			double d = 0.0;
			for (int k = 0; k < j; k++) {
				double[] Lrowk = L[k];
				double s = 0.0;
				for (int i = 0; i < k; i++)
					s += Lrowk[i]*Lrowj[i];
				Lrowj[k] = s = (A[j][k] - s)/L[k][k];
				d = d + s*s;
				isspd = isspd & (A[k][j] == A[j][k]); 
			}
			d = A[j][j] - d;
			isspd = isspd & (d > 0.0);
			L[j][j] = Math.sqrt(Math.max(d,0.0));
			for (int k = j+1; k < n; k++)
				L[j][k] = 0.0;
		}
	}

	/* ------------------------
	 * Public Methods
	 * ------------------------ */

	/**
	 * Is the matrix symmetric and positive definite?
	 * @return     true if A is symmetric and positive definite.
	 */
	public boolean isSPD ()
	{	return isspd; }

	/**
	 * Return triangular factor.
	 * @return     L
	 */
	public Matrix getL ()
	{	return new Matrix(L); }

	/**
	 * Solve A*X = B
	 * @param  B   A Matrix with as many rows as A and any number of columns.
	 * @return     X so that L*L'*X = B
	 * @exception  IllegalArgumentException  Matrix row dimensions must agree.
	 * @exception  RuntimeException  Matrix is not symmetric positive definite.
	 */
	public Matrix solve (Matrix B) {
		if (B.getRowDimension() != n)
			throw new IllegalArgumentException(sMatrixRowMustAgree);
		if (!isspd)
			throw new RuntimeException(sMatrixNotSymPosDef);

		// Copy right hand side.
		double[][] X = B.getArrayCopy();
		int nx = X[0].length;

		// Solve L*Y = B;
		for (int k = 0; k < n; k++) {
			for (int i = k+1; i < n; i++) {
				for (int j = 0; j < nx; j++) {
					X[i][j] -= X[k][j]*L[i][k];
				}
			}
			for (int j = 0; j < nx; j++) {
				X[k][j] /= L[k][k];
			}
		}

		// Solve L'*X = Y;
		for (int k = n-1; k >= 0; k--) {
			for (int j = 0; j < nx; j++) {
				X[k][j] /= L[k][k];
			}
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < nx; j++) {
					X[i][j] -= X[k][j]*L[k][i];
				}
			}
		}
		return new Matrix(X,n,nx);
	}

	/**
	 * Invert a matrix using Cholesky decomposition. Added by Roby Joehanes
	 * @return
	 */
	public Matrix inverse() {
		if (!isspd)
			throw new RuntimeException(sMatrixNotSymPosDef);

		// Copy right hand side.
		double[][] X = new double[n][n];
		for (int i = 0; i < n; i++)
			X[i][i] = 1;

		// Solve L*Y = B;
		for (int k = 0; k < n; k++) {
			for (int i = k+1; i < n; i++)
				for (int j = 0; j < n; j++)
					X[i][j] -= X[k][j]*L[i][k];
			for (int j = 0; j < n; j++)
				X[k][j] /= L[k][k];
		}

		// Solve L'*X = Y;
		for (int k = n-1; k >= 0; k--) {
			for (int j = 0; j < n; j++)
				X[k][j] /= L[k][k];
			for (int i = 0; i < k; i++)
				for (int j = 0; j < n; j++)
					X[i][j] -= X[k][j]*L[k][i];
		}
		return new Matrix(X);
	}
}
