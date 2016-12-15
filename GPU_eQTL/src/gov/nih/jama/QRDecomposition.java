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

import java.util.Set;
import java.util.TreeSet;

import static gov.nih.jama.Matrix.sMatrixRowMustAgree;
import static gov.nih.utils.matrix.QMatrixUtils.createSubMatrixByColumns;
import static gov.nih.utils.matrix.QMatrixUtils.createSubVector;

/**
 * QR Decomposition.
 * 
 * <P>For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
 * orthogonal matrix Q and an n-by-n upper triangular matrix R so that
 * A = Q*R.
 * 
 * <P>The QR decompostion always exists, even if the matrix does not have
 * full rank, so the constructor will never fail.  The primary use of the
 * QR decomposition is in the least squares solution of nonsquare systems
 * of simultaneous linear equations.  <s>This will fail if isFullRank()
 * returns false.</s>
 * 
 * <P>Roby Joehanes' modification:
 * <ul>
 * <li>Optimized for speed</li>
 * <li>Added partial Householder transformation</li>
 * <li>Added variable deletion</li>
 * <li>Added pivoting support</li>
 * <li>Handle rank deficient matrices gracefully, instead of failing</li>
 * <li>Bugfixes</li>
 * </ul>
 */
public class QRDecomposition implements java.io.Serializable
{
	private static final double kDefaultQRTolerance = 1e-6;
	private static final int kInvalidColumnIndex = -1;

	/* ------------------------
	 * Class variables
	 * ------------------------ */

	/**
	 * Array for internal storage of decomposition.
	 * @serial internal array storage.
	 */
	private double[][] mQR;

	/**
	 * Row and column dimensions.
	 * @serial column dimension.
	 * @serial row dimension.
	 */
	private int mNumRows, mNumColumns;

	/**
	 * Array for internal storage of diagonal of R.
	 * @serial diagonal of R.
	 */
	private double[] mRDiag;

	private int mPivot[], mRank;

	/* ------------------------
	 * Constructor
	 * ------------------------ */

	private QRDecomposition() {}

	/**
	 * QR Decomposition, computed by Householder reflections, with default tolerance (1e-6).
	 * @param A    Rectangular matrix
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QRDecomposition (Matrix A) {
		this(A.getArray(), kDefaultQRTolerance, false, false);
	}

	/**
	 * QR Decomposition, computed by Householder reflections.
	 * @param A    Rectangular matrix
	 * @param tolerance Small integer below which a number is considered zero
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QRDecomposition (Matrix A, double tolerance) {
		this(A.getArray(), tolerance, false, false);
	}

	public QRDecomposition (Matrix A, boolean transpose) {
		this(A.getArray(), kDefaultQRTolerance, transpose, false);
	}

	/**
	 * QR Decomposition, computed by Householder reflections, with default tolerance (1e-6).
	 * @param A    Rectangular array
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QRDecomposition (double[][] A) {
		this(A, kDefaultQRTolerance, false, false);
	}

	public QRDecomposition (double[][] A, boolean transpose) {
		this(A, kDefaultQRTolerance, transpose, false);
	}

	/**
	 * QR Decomposition, computed by Householder reflections.
	 * @param A    Rectangular array
	 * @param tolerance Small integer below which a number is considered zero
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QRDecomposition (double[][] A, double tolerance) {
		this(A, tolerance, false, false);
	}

	/**
	 * QR Decomposition, computed by Householder reflections.
	 * @param A    Rectangular array
	 * @param tolerance Small integer below which a number is considered zero
	 * @param transpose If true, A is in transpose position
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QRDecomposition (double[][] A, double tolerance, boolean transpose) {
		this(A, tolerance, transpose, false);
	}

	public QRDecomposition (double[][] A, double tolerance, boolean transpose, boolean addMean) {
		// Initialize
		int k = addMean ? 1 : 0;
		if (!transpose) {
			mNumRows = A.length;
			int col = A[0].length;
			mNumColumns = col + k;
			mQR = new double[mNumRows][mNumColumns];
			for (int i = 0; i < mNumRows; i++) {
				mQR[i][0] = 1;
				System.arraycopy(A[i], 0, mQR[i], k, col);
			}
		} else {
			mNumRows = A[0].length;
			int col = A.length;
			mNumColumns = col + k;
			mQR = new double[mNumRows][mNumColumns];
			for (int i = 0; i < mNumRows; i++) {
				double[] curQR = mQR[i];
				curQR[0] = 1;
				for (int j = 0; j < col; j++)
					curQR[j + k] = A[j][i];
			}
		}
		mRDiag = new double[mNumColumns];
		mPivot = new int[mNumColumns];
		for (int i = 0; i < mNumColumns; i++)
			mPivot[i] = i;
		mRank = doHouseholder(mQR, mRDiag, mPivot, 0, mNumColumns, tolerance);
	}

	/**
	 * Main loop of the Householder algorithm
	 * @param qr QR matrix
	 * @param rDiag R diagonal array
	 * @param pivot pivot indices
	 * @param startCol Starting column
	 * @param oldRank previous rank
	 * @param tol Tolerance
	 * @return new rank
	 */
	private static final int doHouseholder(double[][] qr, double[] rDiag, int[] pivot, int startCol, int oldRank, double tol)
	{
		int
			numRows = qr.length,
			numCols = qr[0].length,
			rank = oldRank;
		// double curExtremeNorm = 0;
		// Main loop.
		for (int colNo = startCol; colNo < rank; colNo++) {
			// Compute 2-norm of k-th column without under/overflow.
			double nrm = 0;
			// The old norm code did not handle ill-conditioned data.
			// This one is much slower but safer
			for (int rowNo = colNo; rowNo < numRows; rowNo++)
				nrm = Math.hypot(nrm,qr[rowNo][colNo]);
			double absNrm = Math.abs(nrm);
			// if (absNrm > curExtremeNorm)
			//	curExtremeNorm = absNrm;

			// if (curExtremeNorm > 0 && absNrm / curExtremeNorm > tol) { // This is WRONG!
			if (absNrm > tol) {
				// Form k-th Householder vector.
				if (qr[colNo][colNo] < 0)
					nrm = -nrm;
				for (int rowNo = colNo; rowNo < numRows; rowNo++)
					qr[rowNo][colNo] /= nrm;
				qr[colNo][colNo] += 1.0;

				// Apply transformation to remaining columns.
				for (int col = colNo+1; col < rank; col++) {
					double s = 0.0; 
					for (int rowNo = colNo; rowNo < numRows; rowNo++)
						s += qr[rowNo][colNo]*qr[rowNo][col];
					s = -s/qr[colNo][colNo];
					for (int rowNo = colNo; rowNo < numRows; rowNo++)
						qr[rowNo][col] += s*qr[rowNo][colNo];
				}
				rDiag[colNo] = -nrm;
			} else {
				rank--;
				if (colNo < numCols - 1) {
					int
						cp1 = colNo + 1,
						colMin1 = numCols - 1;
					for (int rowNo = 0; rowNo < numRows; rowNo++) {
						double
							curQR[] = qr[rowNo],
							temp = curQR[colNo];
						for (int col = cp1; col < numCols; col++)
							curQR[col - 1] = curQR[col];
						curQR[colMin1] = temp;
					}
					int temp = pivot[colNo];
					for (int col = cp1; col < numCols; col++) {
						pivot[col - 1] = pivot[col];
						rDiag[col - 1] = rDiag[col];
					}
					pivot[colMin1] = temp;
					rDiag[colMin1] = -nrm;
					colNo--;
				}
			}
		}
		return rank;
	}

	protected int findQRColumnIndex(int colNo)
	{
		int n = mPivot.length;
		if (colNo >= n)
			return kInvalidColumnIndex;
		if (colNo == mPivot[colNo]) // Easy, non-singular case
			return colNo;
		for (int i = 0; i < n; i++)
			if (mPivot[i] == colNo)
				return i;
		return kInvalidColumnIndex;
	}

	/**
	 * This is buggy. You can only delete the last column without any impact
	 * @param idx
	 * @return
	 */
	private QRDecomposition deleteColumn(int idx)
	{
		int colNo = findQRColumnIndex(idx);
		if (colNo == kInvalidColumnIndex)
			return null;
		int numNewColumns = mNumColumns - 1;
		double[][] newQR = new double[mNumRows][numNewColumns];
		double[] newRDiag = new double[numNewColumns];
		int[] newPivot = new int[numNewColumns];
		if (colNo == 0)
		{
			System.arraycopy(mRDiag, 1, newRDiag, 0, numNewColumns);
			for (int i = 0; i < mNumRows; i++)
				System.arraycopy(mQR[i], 1, newQR[i], 0, numNewColumns);
			for (int i = 1; i < mNumColumns; i++)
				newPivot[i - 1] = mPivot[i] - 1;
		}
		else if (colNo == numNewColumns)
		{
			System.arraycopy(mRDiag, 0, newRDiag, 0, numNewColumns);
			for (int i = 0; i < mNumRows; i++)
				System.arraycopy(mQR[i], 0, newQR[i], 0, numNewColumns);
			System.arraycopy(mPivot, 0, newPivot, 0, numNewColumns);
		}
		else
		{
			int numRemainingCol = numNewColumns - colNo;
			System.arraycopy(mRDiag, 0, newRDiag, 0, colNo);
			System.arraycopy(mRDiag, colNo + 1, newRDiag, colNo, numRemainingCol);
			for (int i = 0; i < mNumRows; i++)
			{
				System.arraycopy(mQR[i], 0, newQR[i], 0, colNo);
				System.arraycopy(mQR[i], colNo + 1, newQR[i], colNo, numRemainingCol);
			}
			System.arraycopy(mPivot, 0, newPivot, 0, colNo);
			for (int i = colNo + 1; i < mNumColumns; i++)
				newPivot[i - 1] = mPivot[i] - 1;
		}
		if (colNo < mRank)
			doHouseholder(newQR, newRDiag, newPivot, colNo, mRank - 1, kDefaultQRTolerance);

		QRDecomposition qr = new QRDecomposition();
		qr.mQR = newQR;
		qr.mRDiag = newRDiag;
		qr.mPivot = newPivot;
		qr.mNumColumns = numNewColumns;
		qr.mNumRows = mNumRows;
		qr.mRank = colNo < mRank ? mRank - 1: mRank;
		return qr;
	}

	QRDecomposition deleteColumn(int[] idx)
	{
		if (idx.length == 1)
			return deleteColumn(idx[0]);
		Set<Integer> deletedIdxSet = new TreeSet<Integer>();
		for (int i: idx)
		{
			int colNo = findQRColumnIndex(i);
			if (colNo != kInvalidColumnIndex)
				deletedIdxSet.add(i);
		}
		if (deletedIdxSet.size() == 0) // Nothing to delete
			return null;
		int
			numColsToDelete = idx.length,
			numNewColumns = mNumColumns - numColsToDelete;
		int[] idxToKeep = new int[numNewColumns];
		for (int i = 0, j = 0; i < mNumColumns; i++)
		{
			if (!deletedIdxSet.contains(i))
				idxToKeep[j++] = i;
		}
		double[][] newQR = createSubMatrixByColumns(idxToKeep, mQR);
		double[] newRDiag = createSubVector(idxToKeep, mRDiag);
		int
			newPivot[] = new int[numNewColumns],
			dec = 0;
		for (int i = 0, j = 0; i < mNumColumns; i++)
		{
			if (!deletedIdxSet.contains(i))
				newPivot[j++] = mPivot[i] - dec;
			else
				dec++;
		}
		for (int i: deletedIdxSet)
			if (i >= mRank)
				dec--;
		if (dec > 0)
			doHouseholder(newQR, newRDiag, newPivot, deletedIdxSet.iterator().next(), mRank - dec, kDefaultQRTolerance);

		QRDecomposition qr = new QRDecomposition();
		qr.mQR = newQR;
		qr.mRDiag = newRDiag;
		qr.mPivot = newPivot;
		qr.mNumColumns = numNewColumns;
		qr.mNumRows = mNumRows;
		qr.mRank = mRank - dec;
		return qr;
	}

	public QRDecomposition addMoreColumns(double[][] A)
	{	return addMoreColumns(A, kDefaultQRTolerance);}

	/**
	 * Partial Householder transformation.
	 * Added by Roby Joehanes
	 * @param A numRows x numAddlCols
	 * @return
	 */
	public QRDecomposition addMoreColumns(double[][] A, double tolerance)
	{
		if (A.length != mNumRows)
			throw new RuntimeException(); // Has to be of the same row
		int
			numAddlCol = A[0].length,
			numTotalColumns = mNumColumns + numAddlCol,
			numSingularVars = mNumColumns - mRank,
			singularIdx = mRank + numAddlCol,
			newPivot[] = new int[numTotalColumns];
		double
			newRDiag[] = new double[numTotalColumns],
			newQR[][] = new double[mNumRows][numTotalColumns];
		if (mRank == mNumColumns)
		{
			for (int i = 0; i < mNumRows; i++)
			{
				System.arraycopy(mQR[i], 0, newQR[i], 0, mNumColumns);
				System.arraycopy(A[i], 0, newQR[i], mNumColumns, numAddlCol);
			}
			System.arraycopy(mRDiag, 0, newRDiag, 0, mNumColumns);
			System.arraycopy(mPivot, 0, newPivot, 0, mNumColumns);
		} else
		{
			int nextCol = mRank + numAddlCol;
			for (int i = 0; i < mNumRows; i++)
			{
				System.arraycopy(mQR[i], 0, newQR[i], 0, mRank);
				System.arraycopy(A[i], 0, newQR[i], mRank, numAddlCol);
				System.arraycopy(mQR[i], mRank, newQR[i], nextCol, numSingularVars);
			}
			System.arraycopy(mRDiag, 0, newRDiag, 0, mRank);
			System.arraycopy(mRDiag, mRank, newRDiag, singularIdx, numSingularVars);
			System.arraycopy(mPivot, 0, newPivot, 0, mRank);
			System.arraycopy(mPivot, mRank, newPivot, singularIdx, numSingularVars);
			for (int i = mNumColumns; i < numTotalColumns; i++)
				newPivot[mRank + i - mNumColumns] = i;
		}

		// Apply transformation to additional columns.
		for (int colNo = 0; colNo < mRank; colNo++)
		{
			if (Math.abs(newRDiag[colNo]) > tolerance)
			{
				for (int col = mRank; col < singularIdx; col++)
				{
					double s = 0.0; 
					for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
						s += newQR[rowNo][colNo]*newQR[rowNo][col];
					s = -s/newQR[colNo][colNo];
					for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
						newQR[rowNo][col] += s*newQR[rowNo][colNo];
				}
			}
		}

		int newRank = doHouseholder(newQR, newRDiag, newPivot, mRank, singularIdx, tolerance);

		QRDecomposition qr = new QRDecomposition();
		qr.mQR = newQR;
		qr.mRDiag = newRDiag;
		qr.mNumRows = mNumRows;
		qr.mNumColumns = numTotalColumns;
		qr.mPivot = newPivot;
		qr.mRank = newRank;
		return qr;
	}

	/**
	 * Same as <tt>addMoreColumns</tt>, but A is transposed
	 * @param A numAddlCols x numRows
	 * @return
	 */
	public QRDecomposition addMoreColumnsTransposed(double[][] A)
	{	return addMoreColumnsTransposed(A, kDefaultQRTolerance); }

	/**
	 * Same as <tt>addMoreColumns</tt>, but A is transposed
	 * @param A numAddlCols x numRows
	 * @param tolerance
	 * @return
	 */
	public QRDecomposition addMoreColumnsTransposed(double[][] A, double tolerance)
	{
		if (A[0].length != mNumRows)
			throw new RuntimeException(); // Has to be of the same row
		int
			numAddlCol = A.length,
			numTotalColumns = mNumColumns + numAddlCol,
			numSingularVars = mNumColumns - mRank,
			singularIdx = mRank + numAddlCol,
			newPivot[] = new int[numTotalColumns];
		double
			newRDiag[] = new double[numTotalColumns],
			newQR[][] = new double[mNumRows][numTotalColumns];
		if (mRank == mNumColumns)
		{
			for (int i = 0; i < mNumRows; i++)
			{
				double[] curNewQR = newQR[i];
				System.arraycopy(mQR[i], 0, curNewQR, 0, mNumColumns);
				for (int j = 0; j < numAddlCol; j++)
					curNewQR[j + mNumColumns] = A[j][i];
			}
			System.arraycopy(mRDiag, 0, newRDiag, 0, mNumColumns);
			System.arraycopy(mPivot, 0, newPivot, 0, mNumColumns);
			for (int i = mNumColumns; i < numTotalColumns; i++)
				newPivot[i] = i;
		} else
		{
			int nextCol = mRank + numAddlCol;
			for (int i = 0; i < mNumRows; i++)
			{
				double[] curNewQR = newQR[i];
				System.arraycopy(mQR[i], 0, curNewQR, 0, mRank);
				for (int j = 0; j < numAddlCol; j++)
					curNewQR[j + mRank] = A[j][i];
				System.arraycopy(mQR[i], mRank, newQR[i], nextCol, numSingularVars);
			}
			System.arraycopy(mRDiag, 0, newRDiag, 0, mRank);
			System.arraycopy(mRDiag, mRank, newRDiag, singularIdx, numSingularVars);
			System.arraycopy(mPivot, 0, newPivot, 0, mRank);
			System.arraycopy(mPivot, mRank, newPivot, singularIdx, numSingularVars);
			int deltaRank = mRank - mNumColumns;
			for (int i = mNumColumns; i < numTotalColumns; i++)
				newPivot[deltaRank + i] = i;
		}

		// Apply transformation to additional columns.
		for (int colNo = 0; colNo < mRank; colNo++)
		{
			if (Math.abs(newRDiag[colNo]) > tolerance)
			{
				for (int col = mRank; col < singularIdx; col++)
				{
					double s = 0.0; 
					for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
						s += newQR[rowNo][colNo]*newQR[rowNo][col];
					s = -s/newQR[colNo][colNo];
					for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
						newQR[rowNo][col] += s*newQR[rowNo][colNo];
				}
			}
		}

		int newRank = doHouseholder(newQR, newRDiag, newPivot, mRank, singularIdx, tolerance);

		QRDecomposition qr = new QRDecomposition();
		qr.mQR = newQR;
		qr.mRDiag = newRDiag;
		qr.mNumRows = mNumRows;
		qr.mNumColumns = numTotalColumns;
		qr.mPivot = newPivot;
		qr.mRank = newRank;
		return qr;
	}

	/* ------------------------
	 * Public Methods
	 * ------------------------ */

	/**
	 * Is the matrix full rank?
	 * @return     true if R, and hence A, has full rank.
	 */
	public boolean isFullRank ()
	{	return mRank == mNumColumns; }

	public int getNumColumns()
	{	return mNumColumns; }

	public int getNumRows()
	{	return mNumRows; }

	/**
	 * Convenience method to get the diagonal elements of the R matrix
	 * @return
	 */
	public double[] getRDiag()
	{	return mRDiag; }

	/**
	 * Convenience method to get the compact QR matrix
	 * @return
	 */
	public double[][] getQRArray()
	{	return mQR; }

	public int[] getPivot()
	{	return mPivot; }

	public int getRank()
	{	return mRank; }

	/**
	 * Return the Householder vectors
	 * @return     Lower trapezoidal matrix whose columns define the reflections
	 */
	public Matrix getH () {
		Matrix X = new Matrix(mNumRows,mNumColumns);
		double[][] H = X.getArray();
		for (int rowNo = 0; rowNo < mNumRows; rowNo++)
			System.arraycopy(mQR[rowNo], rowNo, H[rowNo], rowNo, mNumColumns - rowNo);
		return X;
	}

	/**
	 * Return the upper triangular factor
	 * @return     R
	 */
	public Matrix getR () {
		return new Matrix(getRArray());
	}

	public double[][] getRArray () {
		double[][] R = new double[mRank][mRank];
		for (int columnNo = 0; columnNo < mRank; columnNo++) {
			R[columnNo][columnNo] = mRDiag[columnNo];
			int idx = columnNo + 1;
			System.arraycopy(mQR[columnNo], idx, R[columnNo], idx, mRank - idx);
		}
		return R;
	}

	/**
	 * Generate and return the (economy-sized) orthogonal factor
	 * @return     Q
	 */
	public Matrix getQ () {
		return new Matrix(getQArray());
	}

	public double[][] getQArray () {
		double[][] Q = new double[mNumRows][mNumColumns];
		for (int colNo = mNumColumns-1; colNo >= 0; colNo--)
		{
			Q[colNo][colNo] = 1.0;
			for (int col = colNo; col < mNumColumns; col++)
			{
				if (mQR[colNo][colNo] != 0)
				{
					double s = 0.0;
					for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
						s += mQR[rowNo][colNo]*Q[rowNo][col];
					s = -s/mQR[colNo][colNo];
					for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
						Q[rowNo][col] += s*mQR[rowNo][colNo];
				}
			}
		}
		return Q;
	}

	/**
	 * Least squares solution of A*X = B
	 * @param B    A Matrix with as many rows as A and any number of columns.
	 * @return     Beta that minimizes the two norm of Q*R*X-B.
	 * @exception  IllegalArgumentException  Matrix row dimensions must agree.
	 */
	public Matrix solve (Matrix B)
	{
		double[][] barr = B.getArray();
		int
			nrow = barr.length,
			nx = barr[0].length;
		double[][]
			X = new double[nx][nrow],
			beta = new double[mNumColumns][nx];
		return new Matrix(solve(barr, X, beta));
	}

	/**
	 * Least squares solution of A*X = B
	 * @param B    An array with as many rows as A and any number of columns.
	 * @return     Beta that minimizes the two norm of Q*R*X-B.
	 * @exception  IllegalArgumentException  Matrix row dimensions must agree.
	 */
	public double[][] solve (double[][] B)
	{
		int
			nrow = B.length,
			nx = B[0].length;
		double[][]
			X = new double[nx][nrow],
			beta = new double[mNumColumns][nx];
		return solve(B, X, beta);
	}

	/**
	 * Least squares solution of A*X = B
	 * @param B    An array with as many rows as A and any number of columns.
	 * @param X    A temporary array with the same dimension as B transpose
	 * @param beta An array to hold the result. Must have the same number of rows as <tt>getNumColumns()</tt>
	 *             and as many columns as B's.
	 * @return     Beta that minimizes the two norm of Q*R*X-B.
	 * @exception  IllegalArgumentException  Matrix row dimensions must agree.
	 */
	public double[][] solve (double[][] B, double[][] X, double[][] beta)
	{
		// Copy right hand side
		int
			nrow = B.length,
			nx = B[0].length;
		if (nrow != mNumRows)
			throw new IllegalArgumentException(sMatrixRowMustAgree);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < nrow; j++)
				X[i][j] = B[j][i];

		// Compute Y = transpose(Q)*B
		for (int j = 0; j < nx; j++)
		{
			double[] Xj = X[j];
			for (int colNo = 0; colNo < mRank; colNo++)
			{
				double s = 0.0; 
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					s += mQR[rowNo][colNo]*Xj[rowNo];
				s = -s/mQR[colNo][colNo];
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					Xj[rowNo] += s*mQR[rowNo][colNo];
			}
			for (int colNo = 0; colNo < mRank; colNo++)
				beta[colNo][j] = Xj[colNo];
		}
		// Solve R*X = Y;
		for (int k = mRank-1; k >= 0; k--)
		{
			double rd = mRDiag[k];
			for (int j = 0; j < nx; j++)
				beta[k][j] /= rd;
			for (int i = 0; i < k; i++)
			{
				double[]
					Xi = beta[i],
					QRi = mQR[i];
				for (int j = 0; j < nx; j++)
					Xi[j] -= beta[k][j]*QRi[k];
			}
		}

		for (int i = mRank - 1; i >= 0; i--)
			for (int j = 0; j < nx; j++)
				beta[mPivot[i]][j] = beta[i][j];

		for (int i = mRank; i < mNumColumns; i++)
			for (int j = 0; j < nx; j++)
				beta[mPivot[i]][j] = 0;
		return beta;
	}

	/**
	 * Same as <tt>solve</tt>, but matrix B is in transpose position
	 * @param B numDim x numRows
	 * @return
	 */
	public Matrix solveTranspose (Matrix B)
	{
		double[][] barr = B.getArray();
		int
			nrow = barr[0].length,
			nx = barr.length;
		double[][]
			X = new double[nx][nrow],
			beta = new double[nx][mNumColumns];
		return new Matrix(solveTranspose(barr, X, beta));
	}

	/**
	 * Same as <tt>solve</tt>, but array B is in transpose position
	 * @param B numDim x numRows
	 * @return
	 */
	public double[][] solveTranspose (double[][] B)
	{
		int
			nrow = B[0].length,
			nx = B.length;
		double[][]
			X = new double[nx][nrow],
			beta = new double[nx][mNumColumns];
		return solveTranspose(B, X, beta);
	}

	/**
	 * Same as <tt>solve</tt>, but array B is in transpose position
	 * @param B numDim x numRows
	 * @param X    A temporary array with the same dimension as B
	 * @param beta An array to hold the result. Must have the same number of columns as <tt>getNumColumns()</tt>
	 *             and as many rows as B's.
	 * @return
	 */
	public double[][] solveTranspose (double[][] B, double[][] X, double[][] beta)
	{
		// Copy right hand side
		int
			nrow = B[0].length,
			nx = B.length;
		if (nrow != mNumRows)
			throw new IllegalArgumentException(sMatrixRowMustAgree);
		for (int i = 0; i < nx; i++)
			System.arraycopy(B[i], 0, X[i], 0, nrow);

		// Compute Y = transpose(Q)*B
		for (int j = 0; j < nx; j++)
		{
			double[] Xj = X[j];
			for (int colNo = 0; colNo < mRank; colNo++)
			{
				double s = 0.0; 
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					s += mQR[rowNo][colNo]*Xj[rowNo];
				s = -s/mQR[colNo][colNo];
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					Xj[rowNo] += s*mQR[rowNo][colNo];
			}
			System.arraycopy(Xj, 0, beta[j], 0, mRank);
		}
		// Solve R*X = Y;
		for (int k = mRank-1; k >= 0; k--)
		{
			double rd = mRDiag[k];
			for (int j = 0; j < nx; j++)
				beta[j][k] /= rd;
			for (int i = 0; i < k; i++)
			{
				double[] QRi = mQR[i];
				for (int j = 0; j < nx; j++)
					beta[j][i] -= beta[j][k]*QRi[k];
			}
		}
		for (int j = 0; j < nx; j++)
		{
			double[] betaj = beta[j];
			for (int i = mRank - 1; i >= 0; i--)
				betaj[mPivot[i]] = betaj[i];
			for (int i = mRank; i < mNumColumns; i++)
				betaj[mPivot[i]] = 0;
		}
		return beta;
	}

	public double[] calculateQtY(double[] Y)
	{
		return calculateQtYTranspose(new double[][] {Y}, new double[1][Y.length])[0];
	}

	public double[][] calculateQtY(double[][] Y)
	{
		return calculateQtY(Y, new double[Y.length][Y[0].length]);
	}

	/**
	 * Calculate Q'Y
	 * @param Y numRows x numDim
	 * @param X Temporary array with the same dimension as Y
	 * @return
	 */
	public double[][] calculateQtY(double[][] Y, double[][] X)
	{
		int
			nrow = Y.length,
			nx = Y[0].length;
		if (nrow != mNumRows)
			throw new IllegalArgumentException(sMatrixRowMustAgree);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < nrow; j++)
				X[i][j] = Y[j][i];
	
		// Compute Y = transpose(Q)*B
		for (int j = 0; j < nx; j++)
		{
			double[] Xj = X[j];
			for (int colNo = 0; colNo < mRank; colNo++)
			{
				double s = 0.0; 
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					s += mQR[rowNo][colNo]*Xj[rowNo];
				s = -s/mQR[colNo][colNo];
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					Xj[rowNo] += s*mQR[rowNo][colNo];
			}
		}
		return X;
	}

	public double[][] calculateQtYTranspose(double[][] Y)
	{
		return calculateQtYTranspose(Y, new double[Y.length][Y[0].length]);
	}

	/**
	 * Calculate Q'Y'
	 * @param Y numDim x numRows
	 * @param X Temporary array with the same dimension as Y
	 * @return
	 */
	public double[][] calculateQtYTranspose(double[][] Y, double[][] X)
	{
		int
			nrow = Y[0].length,
			nx = Y.length;
		if (nrow != mNumRows)
			throw new IllegalArgumentException(sMatrixRowMustAgree);
		for (int i = 0; i < nx; i++)
			System.arraycopy(Y[i], 0, X[i], 0, nrow);
	
		// Compute Y = transpose(Q)*B
		for (int j = 0; j < nx; j++)
		{
			double[] Xj = X[j];
			for (int colNo = 0; colNo < mRank; colNo++)
			{
				double s = 0.0; 
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					s += mQR[rowNo][colNo]*Xj[rowNo];
				s = -s/mQR[colNo][colNo];
				for (int rowNo = colNo; rowNo < mNumRows; rowNo++)
					Xj[rowNo] += s*mQR[rowNo][colNo];
			}
		}
		return X;
	}

	/**
	 * Invert a matrix using QR decomposition. Added by Roby Joehanes
	 * @return
	 */
	public Matrix inverse()
	{
		// Copy right hand side
		int nx = mNumRows;
		double[][] X = new double[mNumRows][mNumRows];
		for (int i = 0; i < mNumRows; i++)
			X[i][i] = 1.0;

		// Compute Y = transpose(Q)*B
		for (int k = 0; k < mNumColumns; k++)
		{
			for (int j = 0; j < nx; j++)
			{
				double s = 0.0; 
				for (int i = k; i < mNumRows; i++)
					s += mQR[i][k]*X[i][j];
				s = -s/mQR[k][k];
				for (int i = k; i < mNumRows; i++)
					X[i][j] += s*mQR[i][k];
			}
		}
		// Solve R*X = Y;
		for (int k = mNumColumns-1; k >= 0; k--)
		{
			double
				rd = mRDiag[k],
				Xk[] = X[k];
			for (int j = 0; j < nx; j++)
				Xk[j] /= rd;
			for (int i = 0; i < k; i++)
			{
				double[]
					Xi = X[i],
					QRi = mQR[i];
				for (int j = 0; j < nx; j++)
					Xi[j] -= Xk[j]*QRi[k];
			}
		}
		return new Matrix(X);
	}
}
