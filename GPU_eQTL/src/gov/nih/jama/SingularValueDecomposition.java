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

import gov.nih.utils.matrix.QMatrixUtils;

import static java.lang.Math.hypot;

/**
 * Singular Value Decomposition.
 * <P>For an m-by-n matrix A with m >= n, the singular value decomposition is
 * an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
 * an n-by-n orthogonal matrix V so that A = U*S*V'.
 * 
 * <P>The singular values, sigma[k] = S[k][k], are ordered so that
 * sigma[0] >= sigma[1] >= ... >= sigma[n-1].
 * 
 * <P>The singular value decomposition always exists, so the constructor will
 * never fail.  The matrix condition number and the effective numerical
 * rank can be computed from this decomposition.
 * 
 * <P>Roby's note: This is Jama's SingularValueDecomposition class that I modified
 * to fix the bug that occurs when n > m. I also optimized it for speed.
 */
public class SingularValueDecomposition {
	static final double eps = Math.pow(2.0,-52.0);
	/* ------------------------
	 * Class variables
	 * ------------------------ */

	/** Arrays for internal storage of U and V.
	 @serial internal storage of U.
	 @serial internal storage of V.
	*/
	private double[][] U, V;

	/** Array for internal storage of singular values.
	 @serial internal storage of singular values.
	*/
	private double[] s;

	/** Row and column dimensions.
	 @serial row dimension.
	 @serial column dimension.
	 */
	private int m, n;

	/* ------------------------
	 * Constructor
	 * ------------------------ */

	public SingularValueDecomposition (double[][] arg)
	{	this (arg, true, true); }

	public SingularValueDecomposition (Matrix arg)
	{	this (arg.getArray(), true, true); }

	public SingularValueDecomposition (Matrix arg, boolean wantu, boolean wantv)
	{	this (arg.getArray(), wantu, wantv); }

	/**
	 * <P>Construct the singular value decomposition. Optimized by Roby Joehanes.
	 * 
	 * <P>Rationale:<br>
	 * SVD decomposes n x p matrix X into USV', where U is an n x m matrix,
	 * V' is an m x p matrix, S is an m x m diagonal matrix, and m = min(n, p),
	 * while SVD of X' equals VSU'.
	 * 
	 * <P>Notice that both U and V are orthonormal matrices, i.e. U'U = I and V'V = I.
	 * 
	 * <P>So, the SVD of X'X = VSU'USV' = V(S^2)V' and SVD of XX' = USV'VSU' = U(S^2)U'.
	 * If n << p, it is advantageous to do SVD on XX' matrix because the SVD is
	 * performed in an n x n matrix. Conversely, when n >> p, we do SVD on X'X. Thus,
	 * saving the computational effort.
	 * 
	 * <P>In X'X case, we need to obtain U = XVS^(-1), and in XX' case, we obtain
	 * V = X'US^(-1). Since S is a diagonal matrix, inverting S takes only O(m) time.
	 * Thus, the decision whether we should obtain the SVD on the quadratic form hinges
	 * upon how expensive to compute U or V (i.e. the missing component when performing
	 * SVD on the quadratic case).
	 * 
	 * <P>Funnily, it's numerically more stable than straight SVD when the ratio between
	 * the highest and the lowest eigenvalues is high. 
	 * 
	 * @author Roby Joehanes
	 * @param A    Rectangular matrix
	 * @return     Structure to access U, S and V.
	 */
	public SingularValueDecomposition (double[][] arg, boolean wantu, boolean wantv)
	{
		// Detect matrix structure and optimize SVD
		int
			row = arg.length,
			col = arg[0].length;
		if (row >= 2 * col) // col is small
		{
			double[][] quad = QMatrixUtils.calculateXtYParallel(arg, arg);
			boolean ww = wantu | wantv;
			computeSVD(quad, ww, ww);
			double tol = Math.max(m,n)*s[0]*eps;
			for (int i = 0; i < col; i++)
				s[i] = s[i] < tol ? 0 : Math.sqrt(s[i]);
			if (wantu)
			{
				U = QMatrixUtils.calculateXYParallel(arg, V);
				for (int i = 0; i < row; i++)
				{
					double[] Ui = U[i];
					for (int j = 0; j < col; j++)
					{
						double sj = s[j];
						Ui[j] = sj == 0 ? 0 : Ui[j] / sj;
					}
				}
				if (!wantv)
					V = null;
			}
			else
				U = null;
		}
		else if (col >= 2 * row) // row is small
		{
			double[][] quad = QMatrixUtils.calculateXYtParallel(arg, arg);
			boolean ww = wantu | wantv;
			computeSVD(quad, ww, ww);
			double tol = Math.max(m,n)*s[0]*eps;
			for (int i = 0; i < row; i++)
				s[i] = s[i] < tol ? 0 : Math.sqrt(s[i]);
			if (wantv)
			{
				V = QMatrixUtils.calculateXtYParallel(arg, U);
				for (int i = 0; i < col; i++)
				{
					double[] Vi = V[i];
					for (int j = 0; j < row; j++)
					{
						double sj = s[j];
						Vi[j] = sj == 0 ? 0 : Vi[j] / sj;
					}
				}
				if (!wantu)
					U = null;
			}
			else
				V = null;
		}
		else
			computeSVD(arg, wantu, wantv);
		
	}

	private void computeSVD(double[][] arg, boolean wantu, boolean wantv)
	{
		// Derived from LINPACK code.
		// Initialize.
		m = arg.length;
		n = arg[0].length;
		int nu = Math.min(m,n);
		s = new double [nu]; // RJ's bugfix
		U = new double [m][nu];
		V = new double [n][nu]; // RJ's bugfix
		double
			e[] = new double [n],
			work[] = new double[m],
			A[][] = new double[m][n];
		for (int i = 0; i < m; i++)
			System.arraycopy(arg[i], 0, A[i], 0, n);

		// Reduce A to bidiagonal form, storing the diagonal elements
		// in s and the super-diagonal elements in e.
		int
			nct = Math.min(m-1,n),
			nrt = Math.max(0,Math.min(n-2,m)),
			mcrtnrt = Math.max(nct,nrt);
		for (int k = 0; k < mcrtnrt; k++)
		{
			if (k < nct)
			{
				// Compute the transformation for the k-th column and
				// place the k-th diagonal in s[k].
				// Compute 2-norm of k-th column without under/overflow.
				double sk = 0;  // RJ's optimization
				for (int i = k; i < m; i++)
					sk = hypot(sk,A[i][k]);

				if (sk != 0.0)
				{
					if (A[k][k] < 0.0)
						sk = -sk;
					for (int i = k; i < m; i++)
						A[i][k] /= sk;
					A[k][k] += 1.0;
				}
				s[k] = -sk;
			}

			for (int j = k+1; j < n; j++)
			{
				if (k < nct && s[k] != 0.0) // RJ's bugfix
				{
					// Apply the transformation.
					double t = 0;
					for (int i = k; i < m; i++)
						t += A[i][k] * A[i][j];
					t = -t / A[k][k];
					for (int i = k; i < m; i++)
						A[i][j] += t * A[i][k];
				}

				// Place the k-th row of A into e for the
				// subsequent calculation of the row transformation.
				e[j] = A[k][j];
			}

			if (wantu && k < nct) // RJ's bugfix
			{
				// Place the transformation in U for subsequent back multiplication.
				for (int i = k; i < m; i++)
					U[i][k] = A[i][k];
			}

			if (k < nrt)
			{
				// Compute the k-th row transformation and place the
				// k-th super-diagonal in e[k].
				// Compute 2-norm without under/overflow.
				double ek = 0; // RJ's optimization
				for (int i = k+1; i < n; i++)
					ek = hypot(ek,e[i]);

				if (ek != 0.0)
				{
					if (e[k+1] < 0.0)
						ek = -ek;

					for (int i = k+1; i < n; i++)
						e[i] /= ek;

					e[k+1] += 1.0;
				}
				e[k] = -ek;

				if (k+1 < m && ek != 0.0) // RJ's bugfix
				{
					// Apply the transformation.
					for (int i = k+1; i < m; i++)
						work[i] = 0.0f;

					for (int j = k+1; j < n; j++)
						for (int i = k+1; i < m; i++)
							work[i] += e[j]*A[i][j];

					for (int j = k+1; j < n; j++)
					{
						double t = -e[j]/e[k+1];
						for (int i = k+1; i < m; i++)
							A[i][j] += t*work[i];
					}
				}

				if (wantv)
				{
					// Place the transformation in V for subsequent
					// back multiplication.
					for (int i = k+1; i < n; i++)
						V[i][k] = e[i];
				}
			}
		}

		// Set up the final bidiagonal matrix or order p.

		int p = nu; // RJ's bugfix
		if (nct < n)
			s[nct] = A[nct][nct];

		if (m < p)
			s[p-1] = 0.0f;

		if (nrt+1 < p)
			e[nrt] = A[nrt][p-1];

		e[p-1] = 0.0f;
		A = null; // Dispose unused buffer -- RJ

		// If required, generate U.
		if (wantu)
		{
			for (int j = nct; j < nu; j++)
			{
				for (int i = 0; i < m; i++)
					U[i][j] = 0.0f;
				U[j][j] = 1.0f;
			}

			for (int k = nct-1; k >= 0; k--)
			{
				if (s[k] != 0.0)
				{
					for (int j = k+1; j < nu; j++)
					{
						double t = 0;
						for (int i = k; i < m; i++)
							t += U[i][k]*U[i][j];
						t = -t/U[k][k];
						for (int i = k; i < m; i++)
							U[i][j] += t*U[i][k];
					}

					for (int i = k; i < m; i++ )
						U[i][k] = -U[i][k];

					U[k][k] = 1.0f + U[k][k];
					for (int i = 0; i < k-1; i++)
						U[i][k] = 0.0f;
				} else {
					for (int i = 0; i < m; i++)
						U[i][k] = 0.0f;
					U[k][k] = 1.0f;
				}
			}
		}

		// If required, generate V.
		if (wantv)
		{
			for (int k = nu-1; k >= 0; k--) // RJ's bugfix
			{
				if (k < nrt && e[k] != 0.0) // RJ's bugfix
				{
					for (int j = k+1; j < nu; j++)
					{
						double t = 0;
						for (int i = k+1; i < n; i++)
							t += V[i][k]*V[i][j];

						t = -t/V[k+1][k];
						for (int i = k+1; i < n; i++)
							V[i][j] += t*V[i][k];
					}
				}
				for (int i = 0; i < n; i++)
					V[i][k] = 0.0f;

				V[k][k] = 1.0f;
			}
		}

		// Main iteration loop for the singular values.

		int
			pp = p-1;
		while (p > 0)
		{
			int k,kase;

			// Here is where a test for too many iterations would go.

			// This section of the program inspects for
			// negligible elements in the s and e arrays.  On
			// completion the variables kase and k are set as follows.

			// kase = 1     if s(p) and e[k-1] are negligible and k<p
			// kase = 2     if s(k) is negligible and k<p
			// kase = 3     if e[k-1] is negligible, k<p, and
			//              s(k), ..., s(p) are not negligible (qr step).
			// kase = 4     if e(p-1) is negligible (convergence).

			for (k = p-2; k >= -1; k--)
			{
				if (k == -1)
					break;

				if (Math.abs(e[k]) <= eps*(Math.abs(s[k]) + Math.abs(s[k+1])))
				{
					e[k] = 0.0f;
					break;
				}
			}

			if (k == p-2)
				kase = 4;
			else
			{
				int ks;
				for (ks = p-1; ks >= k; ks--)
				{
					if (ks == k)
						break;

					double t = (ks != p ? Math.abs(e[ks]) : 0.) + 
								(ks != k+1 ? Math.abs(e[ks-1]) : 0.);
					if (Math.abs(s[ks]) <= eps*t)
					{
						s[ks] = 0.0f;
						break;
					}
				}
				if (ks == k)
					kase = 3;
				else if (ks == p-1)
					kase = 1;
				else {
					kase = 2;
					k = ks;
				}
			}
			k++;

			// Perform the task indicated by kase.

			switch (kase)
			{
				case 1: { // Deflate negligible s(p).
					double f = e[p-2];
					e[p-2] = 0.0f;
					for (int j = p-2; j >= k; j--)
					{
						double
							t = hypot(s[j],f),
							cs = s[j]/t,
							sn = f/t;
						s[j] = t;
						if (j != k)
						{
							f = -sn*e[j-1];
							e[j-1] = cs*e[j-1];
						}
						if (wantv)
						{
							for (int i = 0; i < n; i++)
							{
								t = cs*V[i][j] + sn*V[i][p-1];
								V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
								V[i][j] = t;
							}
						}
					}
				}
				break;

				case 2: { // Split at negligible s(k).
					double f = e[k-1];
					e[k-1] = 0.0f;
					for (int j = k; j < p; j++)
					{
						double
							t = hypot(s[j],f),
							cs = s[j]/t,
							sn = f/t;
						s[j] = t;
						f = -sn*e[j];
						e[j] = cs*e[j];
						if (wantu)
						{
							for (int i = 0; i < m; i++)
							{
								t = cs*U[i][j] + sn*U[i][k-1];
								U[i][k-1] = (-sn*U[i][j] + cs*U[i][k-1]);
								U[i][j] = t;
							}
						}
					}
				}
				break;

				case 3: { // Perform one qr step.
					// Calculate the shift.
					double
						scale = Math.max(Math.max(Math.max(Math.max(
							Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])), 
							Math.abs(s[k])),Math.abs(e[k])),
						sp = s[p-1]/scale,
						spm1 = s[p-2]/scale,
						epm1 = e[p-2]/scale,
						sk = s[k]/scale,
						ek = e[k]/scale,
						b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0f,
						c = (sp*epm1)*(sp*epm1),
						shift = 0.0f;

					if (b != 0.0 || c != 0.0)
					{
						shift = Math.sqrt(b*b + c);
						if (b < 0.0)
							shift = -shift;
						shift = c/(b + shift);
					}
					double
						f = (sk + sp)*(sk - sp) + shift,
						g = sk*ek;
   
					// Chase zeroes.
					for (int j = k; j < p-1; j++)
					{
						double
							t = hypot(f,g),
							cs = f/t,
							sn = g/t;
						if (j != k)
							e[j-1] = t;

						f = cs*s[j] + sn*e[j];
						e[j] = cs*e[j] - sn*s[j];
						g = sn*s[j+1];
						s[j+1] = cs*s[j+1];
						if (wantv)
						{
							for (int i = 0; i < n; i++)
							{
								t = cs*V[i][j] + sn*V[i][j+1];
								V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
								V[i][j] = t;
							}
						}
						t = hypot(f,g);
						cs = f/t;
						sn = g/t;
						s[j] = t;
						f = cs*e[j] + sn*s[j+1];
						s[j+1] = -sn*e[j] + cs*s[j+1];
						g = sn*e[j+1];
						e[j+1] = cs*e[j+1];
						if (wantu && j < m-1) // RJ's bugfix
						{
							for (int i = 0; i < m; i++)
							{
								t = cs*U[i][j] + sn*U[i][j+1];
								U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
								U[i][j] = t;
							}
						}
					}
					e[p-2] = f;
				}
				break;

				case 4: { // Convergence.
					// Make the singular values positive.
					if (s[k] <= 0.0)
					{
						s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
						if (wantv) {
							for (int i = 0; i <= pp; i++)
								V[i][k] = -V[i][k];
						}
					}

					// Order the singular values.
					while (k < pp)
					{
						if (s[k] >= s[k+1])
							break;

						double t = s[k];
						s[k] = s[k+1];
						s[k+1] = t;
						if (wantv && k < n-1) // RJ's bugfix
						{
							for (int i = 0; i < n; i++)
							{
								t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
							}
						}
						if (wantu && k < m-1) // RJ's bugfix
						{
							for (int i = 0; i < m; i++)
							{
								t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
							}
						}
						k++;
					}
					p--;
				}
				break;
			}
		}
	}

	/**
	 * Return the left singular vectors
	 * @return     U
	 */
	public Matrix getU ()
	{	return new Matrix(U); }

	public double[][] getUArray()
	{	return U; }

	/**
	 * Return the right singular vectors
	 * @return     V
	 */
	public Matrix getV ()
	{	return new Matrix(V); }

	public double[][] getVArray()
	{	return V; }

	/**
	 * Return the one-dimensional array of singular values
	 * @return     diagonal of S.
	 */
	public double[] getSingularValues ()
	{	return s; }

	public Matrix getS ()
	{
		int nu = m < n ? m : n;
		double[][] S = new double[nu][nu];
		for (int i = 0; i < nu; i++)
			S[i][i] = s[i];
		return new Matrix(S);
	}

	/** Two norm
	@return     max(S)
	 */
	public double norm2 ()
	{	return s[0]; }
	
	/**
	 * Two norm condition number
	 * @return     max(S)/min(S)
	 */
	public double cond ()
	{	return s[0]/s[Math.min(m,n)-1]; }

	/**
	 * Effective numerical matrix rank
	 * @return     Number of nonnegligible singular values.
	 */
	public int rank ()
	{
		double tol = Math.max(m,n)*s[0]*eps;
		// Exploit the fact that the eigenvalues are sorted in descending order -- RJ
		int r = s.length;
		for (int i = 0; i < r; i++)
			if (s[i] <= tol)
				return i;
		return 0;
	}

	public void purge()
	{
		U = V = null; s = null;
	}

	public static void main(String[] args)
	{
		// When row = 4, col = 8, standard computeSVD breaks down (i.e. V'V != I,
		// especially at the last row/column). But in many other occassions, it's okay.
		int
			numRow = 4,
			numCol = 8;
		double[][] arr = new double[numRow][numCol];
		for (int i = 0; i < numRow; i++)
			for (int j = 0; j < numCol; j++)
				arr[i][j] = 1.0 / (i + j + 1);
		SingularValueDecomposition svd = new SingularValueDecomposition(arr);
		System.out.println(svd.getU());
		System.out.println(svd.getV());
		System.out.println(svd.getU().transpose().times(svd.getU())); // Should be == I
		System.out.println(svd.getV().transpose().times(svd.getV())); // Should be == I
	}
}
