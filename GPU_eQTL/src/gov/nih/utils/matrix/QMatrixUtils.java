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
/**
 * 
 */
package gov.nih.utils.matrix;

import gov.nih.jama.Matrix;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import jdistlib.rng.RandomEngine;
import gov.nih.parallel.QSynchronizedCounter;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

/**
 * Utilities to perform matrix computation. Refactored from QMathUtils.
 * 
 * @author Roby Joehanes
 *
 */
public class QMatrixUtils
{
	public static final double[][] transpose(double[][] mtx)
	{
		int
			ncol = mtx[0].length,
			nrow = mtx.length;
		double[][] mat = new double[ncol][nrow];
		for (int i = 0; i < nrow; i++)
			for (int j = 0; j < ncol; j++)
				mat[j][i] = mtx[i][j];
		return(mat);
	}

	/**
	 * Calculate X - Y'
	 * @param X (m x n)
	 * @param Y (n x m)
	 * @return (m x n)
	 */
	public static final double[][] calculateXMinusYt(double[][] X, double[][] Y)
	{
		int
			m = X.length,
			n = X[0].length;
		double[][] result = new double[m][n];
	
		for (int i = 0; i < m; i++)
		{
			double[]
				Xrow = X[i],
				resultRow = result[i];
			for (int j = 0; j < n; j++)
				resultRow[j] = Xrow[j] - Y[j][i];
		}
		return result;
	}

	public static final double[][] calculateXMinusY(double[][] X, double[][] Y)
	{
		int
			m = X.length,
			n = X[0].length;
		double[][] result = new double[m][n];
	
		for (int i = 0; i < m; i++)
		{
			double[]
				Xrow = X[i],
				Yrow = Y[i],
				resultRow = result[i];
			for (int j = 0; j < n; j++)
				resultRow[j] = Xrow[j] - Yrow[j];
		}
		return result;
	}

	/**
	 * Calculate XY. Eliminate the need to instantiate a matrix.
	 * @param X (n x m)
	 * @param Y (m x p)
	 * @return
	 */
	public static final double[][] calculateXY(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Y[k][j];
				result[i][j] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate XY'. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXYt(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y[0].length,
			cols = Y.length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = Y[j];
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Yrow[k];
				result[i][j] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate X'Y. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXtY(double[][] X, double[][] Y)
	{
		int
			rows = X[0].length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += X[k][i] * Y[k][j];
				result[i][j] = sum;
			}
		return result;
	}

	/**
	 * Calculate I - XY. Eliminate the need to instantiate a matrix.
	 * @param X (n x m)
	 * @param Y (m x p)
	 * @return
	 */
	public static final double[][] calculateIMinusXY(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Y[k][j];
				result[i][j] = i == j ? 1 - sum : -sum;
			}
		}
		return result;
	}

	/**
	 * Calculate I - XY'. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateIMinusXYt(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y[0].length,
			cols = Y.length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = Y[j];
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Yrow[k];
				result[i][j] = i == j ? 1 - sum : -sum;
			}
		}
		return result;
	}

	/**
	 * Calculate I - X'Y. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateIMinusXtY(double[][] X, double[][] Y)
	{
		int
			rows = X[0].length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += X[k][i] * Y[k][j];
				result[i][j] = i == j ? 1 - sum : -sum;
			}
		return result;
	}

	/**
	 * Calculate c(I - XY). Eliminate the need to instantiate a matrix.
	 * @param c
	 * @param X (n x m)
	 * @param Y (m x p)
	 * @return
	 */
	public static final double[][] calculatecIMinusXY(double c, double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Y[k][j];
				result[i][j] = c * (i == j ? 1 - sum : -sum);
			}
		}
		return result;
	}

	/**
	 * Calculate I - XY'. Eliminate the need to transpose.
	 * @param c
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculatecIMinusXYt(double c, double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			m = Y[0].length,
			cols = Y.length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
		{
			double[] Xrow = X[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = Y[j];
				for (int k = 0; k < m; k++)
					sum += Xrow[k] * Yrow[k];
				result[i][j] = c * (i == j ? 1 - sum : -sum);
			}
		}
		return result;
	}

	/**
	 * Calculate I - X'Y. Eliminate the need to transpose.
	 * @param c
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculatecIMinusXtY(double c, double[][] X, double[][] Y)
	{
		int
			rows = X[0].length,
			m = Y.length,
			cols = Y[0].length;
		assert (m == X.length);
		double[][] result = new double[rows][cols];
	
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < m; k++)
					sum += X[k][i] * Y[k][j];
				result[i][j] = c * (i == j ? 1 - sum : -sum);
			}
		return result;
	}

	/**
	 * Calculate X * Y in parallel fashion. Only useful for large matrices (row * col > 90000).
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXYParallel(double[][] X, double[][] Y)
	{
		return QMatrixUtils.parallelMatrixMultiplication(X, Y, null, 1, X.length, Y[0].length, EMultiplicationMode.XY);
	}

	/**
	 * Calculate X * Y' in parallel fashion. Only useful for large matrices (row * col > 90000).
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXYtParallel(double[][] X, double[][] Y)
	{
		return QMatrixUtils.parallelMatrixMultiplication(X, Y, null, 1, X.length, Y.length, EMultiplicationMode.XYt);
	}

	/**
	 * Calculate X' * Y in parallel fashion. Only useful for large matrices (row * col > 90000).
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXtYParallel(double[][] X, double[][] Y)
	{
		return QMatrixUtils.parallelMatrixMultiplication(X, Y, null, 1, X[0].length, Y[0].length, EMultiplicationMode.XtY);
	}

	/**
	 * <P>Return Kronecker product of vectors v and w.
	 * 
	 * <P>Let V = [v1, v2, ..., vm] and W = [w1, w2, ..., wn].
	 * Kronecker product of V and W is denoted by V \otimes W
	 * (the symbol \otimes is an encircled X) and defined as:
	 * 
	 * <p>V \otimes W = [v1*w1, v1*w2, ..., v1*wn, v2*w1, ..., v2*wn, ..., vm*w1, ..., vm*wn]
	 * 
	 * <P>Thus, the resulting vector is of dimension m * n.
	 * Note: the Kronecker product is NOT symmetric!
	 * 
	 * <P>This routine is useful to construct orthogonal contrasts for epistatic effects.
	 * Taking the Kronecker product is in itself an ordinary mathematical routine.
	 * 
	 * @param v
	 * @param w
	 * @return Kronecker product of v and w.
	 */
	public static final double[] kroneckerize(double[] v, double[] w)
	{
		int
			m = v.length,
			n = w.length,
			idx = 0;
		double[] result = new double[m*n];
	
		for (int i = 0; i < m; i++)
		{
			double value = v[i];
			for (int j = 0; j < n; j++, idx++)
				result[idx] = value * w[j];
		}
		return result;
	}

	/**
	 * Invert an upper triangular matrix. What we do here is a series of back substitutions.
	 * We take advantage of the special structure. Instead of the usual O(n^3),
	 * we have O(n^3/4).
	 * 
	 * @param R matrix to invert
	 * @param The inverse
	 */
	public static final double[][] invertUpperTriangularMatrix(double[][] R)
	{
		int rLength = R.length;
		double[][] Rinv = new double[rLength][rLength];
	
		for (int k = rLength - 1; k >= 0; k--)
		{
			Rinv[k][k] = 1;
			for (int i = k; i >= 0; i--)
			{
				for (int j = rLength - 1; j > i; j--)
					Rinv[i][k] -= R[i][j] * Rinv[j][k];
				Rinv[i][k] /= R[i][i];
			}
		}
		return Rinv;
	}

	public static final double[][] invertLowerTriangularMatrix(double[][] R)
	{
		int rLength = R.length;
		double[][] Rinv = new double[rLength][rLength];
	
		for (int k = 0; k < rLength; k++)
		{
			Rinv[k][k] = 1;
			for (int i = k; i < rLength; i++)
			{
				for (int j = 0; j < i; j++)
					Rinv[i][k] -= R[i][j] * Rinv[j][k];
				Rinv[i][k] /= R[i][i];
			}
		}
		return Rinv;
	}

	/**
	 * Parallel matrix multiplication kernel.
	 * 
	 * @see EMultiplicationMode
	 * @see QMultiplicationJob
	 * 
	 * @param X
	 * @param Y
	 * @param Z
	 * @param cons
	 * @param rows
	 * @param cols
	 * @param mode
	 * @return
	 */
	public static final double[][] parallelMatrixMultiplication(double[][] X, double[][] Y, double[][] Z, double cons,
		int rows, int cols, EMultiplicationMode mode)
	{
		int
			numElements = rows * cols,
			numThreads = QSystemUtils.kNumCPUCores > numElements ? numElements : QSystemUtils.kNumCPUCores;
		ExecutorService threadPool = Executors.newFixedThreadPool(numThreads + 1);
		QSynchronizedCounter counter = new QSynchronizedCounter(0, rows);
		double[][] result = new double[rows][cols];
		assert (X[0].length == Y.length);

		if (cons == 1)
			for (int i = 0; i < numThreads; i++)
				threadPool.execute(new QMultiplicationJob(X, Y, Z, result, counter, mode));
		else
			for (int i = 0; i < numThreads; i++)
				threadPool.execute(new QMultiplicationJob(X, Y, Z, result, cons, counter, mode));
	
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException exx) {}
		return result;
	}

	public static final double[][] parallelMatrixMultiplication(double[][] X, double[][] Y, double[][] Z, double cons,
		double[][] result, EMultiplicationMode mode)
	{
		int
			rows = result.length,
			cols = result[0].length,
			numElements = rows * cols,
			numThreads = QSystemUtils.kNumCPUCores > numElements ? numElements : QSystemUtils.kNumCPUCores;
		ExecutorService threadPool = Executors.newFixedThreadPool(numThreads + 1);
		QSynchronizedCounter counter = new QSynchronizedCounter(0, rows);
		assert (X[0].length == Y.length);

		if (cons == 1)
			for (int i = 0; i < numThreads; i++)
				threadPool.execute(new QMultiplicationJob(X, Y, Z, result, counter, mode));
		else
			for (int i = 0; i < numThreads; i++)
				threadPool.execute(new QMultiplicationJob(X, Y, Z, result, cons, counter, mode));

		threadPool.shutdown();
		try {
			threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException exx) {}
		return result;
	}

	/**
	 * Compute X * beta
	 * @param X an n x p matrix
	 * @param beta a vector of length p
	 * @return
	 */
	public static final double[] calculateXBeta(double[][] X, double[] beta)
	{
		int
			n = X.length,
			p = beta.length;
		assert(p == X[0].length);
		double[] yhat = new double[n];

		if (p == 0)
			return yhat;
		for (int i = 0; i < n; i++)
		{
			double
				Xi[] = X[i],
				sum = Xi[0] * beta[0];
			for (int j = 1; j < p; j++)
				sum += Xi[j] * beta[j];
			yhat[i] = sum;
		}
		return yhat;
	}

	/**
	 * Calculate XVX'.
	 * 
	 * <P>Ordinary matrix multiplication will need O(np(n + p)) time.
	 * This routine exploits symmetry and needs only half as much.
	 * If V is not symmetric, it WON'T WORK!
	 * 
	 * @param X an n x p matrix
	 * @param V a p x p SYMMETRIC matrix
	 * @return
	 */
	public static final double[][] calculateQuadraticForm(double[][] X, double[][] V)
	{
		int
			n = X.length,
			p = X[0].length;
		assert (p == V.length && p == V[0].length);
		double
			sum,
			Vk[],
			XVi[] = new double[p],
			result[][] = new double[n][n];
	
		for (int i = 0; i < n; i++)
		{
			double[] Xi = X[i];
			for (int k = 0; k < p; k++)
			{
				sum = 0;
				Vk = V[k];
				for (int l = 0; l < p; l++)
					sum += Xi[l] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = X[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += XVi[k] * Xj[k];
				result[i][j] = result[j][i] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate X'VX.
	 * 
	 * <P>Ordinary matrix multiplication will need O(np(n + p)) time.
	 * This routine exploits symmetry and needs only half as much.
	 * If V is not symmetric, it WON'T WORK!
	 * 
	 * @param X an p x n matrix
	 * @param V a p x p SYMMETRIC matrix.
	 * @return
	 */
	public static final double[][] calculateQuadraticFormTranspose(double[][] X, double[][] V)
	{
		int
			n = X[0].length,
			p = X.length;
		assert (p == V.length && p == V[0].length);
		double
			sum,
			Vk[],
			XVi[] = new double[p],
			result[][] = new double[n][n];
	
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < p; k++)
			{
				sum = 0;
				Vk = V[k];
				for (int l = 0; l < p; l++)
					sum += X[l][i] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += XVi[k] * X[k][j];
				result[i][j] = result[j][i] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate XVX'.
	 * 
	 * <P>Ordinary matrix multiplication will need O(np(n + p)) time.
	 * This routine exploits symmetry and needs only half as much.
	 * V is a diagonal matrix
	 * 
	 * @param X an n x p matrix
	 * @param V an array of p elements that represents the diagonal elements of the matrix
	 * @return
	 */
	public static final double[][] calculateQuadraticFormDiagV(double[][] X, double[] V)
	{
		int
			n = X.length,
			p = X[0].length;
		assert (p == V.length && p == V.length);
		double
			sum,
			XVi[] = new double[p],
			result[][] = new double[n][n];
	
		for (int i = 0; i < n; i++)
		{
			double[] Xi = X[i];
			for (int k = 0; k < p; k++)
				XVi[k] = Xi[k] * V[k];
	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = X[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += XVi[k] * Xj[k];
				result[i][j] = result[j][i] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate X'VX.
	 * 
	 * <P>Ordinary matrix multiplication will need O(np(n + p)) time.
	 * This routine exploits symmetry and needs only half as much.
	 * If V is a diagonal matrix
	 * 
	 * @param X an p x n matrix
	 * @param V an array of p elements that represents the diagonal elements of the matrix
	 * @return
	 */
	public static final double[][] calculateQuadraticFormTransposeDiagV(double[][] X, double[] V)
	{
		int
			n = X[0].length,
			p = X.length;
		assert (p == V.length && p == V.length);
		double
			sum,
			XVi[] = new double[p],
			result[][] = new double[n][n];
	
		for (int i = 0; i < n; i++)
		{
			for (int k = 0; k < p; k++)
				XVi[k] = X[k][i] * V[k];
	
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += XVi[k] * X[k][j];
				result[i][j] = result[j][i] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate XX'.
	 * 
	 * <P>Ordinary matrix multiplication will need O(np(n + p)) time.
	 * This routine exploits symmetry and needs only half as much.
	 * 
	 * @param X an n x p matrix
	 * @return n x n matrix
	 */
	public static final double[][] calculateQuadraticForm(double[][] X)
	{
		int
			n = X.length,
			p = X[0].length;
		double
			sum,
			result[][] = new double[n][n];
	
		for (int i = 0; i < n; i++)
		{
			double[] Xi = X[i];	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = X[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += Xi[k] * Xj[k];
				result[i][j] = result[j][i] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate X'X.
	 * 
	 * <P>Ordinary matrix multiplication will need O(np(n + p)) time.
	 * This routine exploits symmetry and needs only half as much.
	 * 
	 * @param X an n x p matrix
	 * @return p x p SYMMETRIC matrix.
	 */
	public static final double[][] calculateQuadraticFormTranspose(double[][] X)
	{
		int
			n = X[0].length,
			p = X.length;
		double
			sum,
			result[][] = new double[n][n];
	
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += X[k][i] * X[k][j];
				result[i][j] = result[j][i] = sum;
			}
		}
		return result;
	}

	/**
	 * Calculate Y-bar, in case there's no cofactors selected.
	 * So, Y-bar is used in place of mu.
	 * Actually, it's sufficient to return only a single-dimension
	 * array of numTraits. But for the sake of conformity with
	 * computeMu (and ease the computation afterwards), I duplicated
	 * the ybar values across the individuals.
	 * 
	 * @param traitData numTraits x numNonMissingIndividuals
	 * @return an array of numTraits x numIndividuals
	 */
	public static final double[][] calculateYBar(double[][] traitData)
	{
		int
			numTraits = traitData.length,
			numIndividuals = traitData[0].length;
		double result[][] = new double[numTraits][numIndividuals];
	
		for (int traitNo = 0; traitNo < numTraits; traitNo++)
		{
			double
				y[] = traitData[traitNo],
				sum = 0;
			for (int indNo = 0; indNo < numIndividuals; indNo++)
				sum += y[indNo];
			Arrays.fill(result[traitNo], sum / numIndividuals);
		}
		return result;
	}

	/**
	 * Calculate Y-Bar single trait version.
	 * @param traitData
	 * @return
	 */
	public static final double[] calculateYBar(double[] traitData)
	{
		int	numIndividuals = traitData.length;
		double
			sum = 0,
			result[] = new double[numIndividuals];
	
		for (int indNo = 0; indNo < numIndividuals; indNo++)
			sum += traitData[indNo];
		Arrays.fill(result, sum / numIndividuals);
		return result;
	}

	/**
	 * <P>Calculate trace(V<sup>-1</sup>ZZ').
	 * 
	 * <P>Let Z = [z<sub>ij</sub>], where i = 1,...,n, and j = 1,...,r
	 * 
	 * <P>Note that [ZZ']<sub>ij</sub> = \sum_{k=1}^{r} z<sub>ik</sub> z<sub>jk</sub>
	 * 
	 * <P>Let V<sup>-1</sup> = [v<sub>ij</sub>], where i, j = 1,..., n
	 * 
	 * <P>Note that [V<sup>-1</sup>ZZ']<sub>ij</sub> = \sum_{l=1}^{n} v<sub>il</sub> \sum_{k=1}^{r} z<sub>lk</sub> z<sub>jk</sub>
	 * 
	 * <P>Thus, trace(V<sup>-1</sup>ZZ') = \sum_{i=1}^{n} [V<sup>-1</sup>ZZ']<sub>ii</sub> =
	 * \sum_{i=1}^{n} \sum_{l=1}^{n} v<sub>il</sub> \sum_{k=1}^{r} z<sub>lk</sub> z<sub>ik</sub>
	 * 
	 * @param VInv is V<sup>-1</sup> matrix of n x n dimension
	 * @param Z is the matrix of random covariates of n x r dimension
	 * @return
	 */
	public static final double calculateTraceVInvZZTranspose(double[][] VInv, double[][] Z)
	{
		int
			n = VInv.length,
			r = Z[0].length;
		double sum = 0;
		assert (n == VInv[0].length && n == Z.length);
		for (int i = 0; i < n; i++)
		{
			double[]
				Zi = Z[i],
				VInvi = VInv[i];
			for (int l = 0; l < n; l++)
			{
				double[] Zl = Z[l];
				int sub = 0;
				for (int k = 0; k < r; k++)
					sub += Zi[k] * Zl[k];
				sum += VInvi[l] * sub;
			}
		}
		return sum;
	}

	/**
	 * Return the sum of all elements in the matrix.
	 * Required to compute weight in V matrix for QxE experiment
	 * 
	 * @param vector
	 * @return
	 */
	public static final double sum(double[] vector)
	{
		double sum = 0;
		int dim = vector.length;
		for (int i = 0; i < dim; i++)
			sum += vector[i];
		return sum;
	}

	/**
	 * Return the sum of the matrix, column wise.
	 * 
	 * @param matrix
	 * @return
	 */
	public static final double[] sumColumns(double[][] matrix)
	{
		int
			numRows = matrix.length,
			numCols = matrix[0].length;
		double[] result = new double[numCols];
		for (int i = 0; i < numCols; i++)
		{
			double sum = 0;
			for (int j = 0; j < numRows; j++)
				sum += matrix[j][i];
			result[i] = sum;
		}
		return result;
	}

	/**
	 * Return the mean of each column of matrix
	 * @param matrix
	 * @return
	 */
	public static final double[] meanColumns(double[][] matrix)
	{
		int
			numRows = matrix.length,
			numCols = matrix[0].length;
		double[] result = new double[numCols];
		for (int i = 0; i < numCols; i++)
		{
			double sum = 0;
			for (int j = 0; j < numRows; j++)
				sum += matrix[j][i];
			result[i] = sum / numRows;
		}
		return result;
	}

	/**
	 * Create identity matrix of n x n
	 * @param n
	 * @return
	 */
	public static final double[][] createIdentityMatrix(int n)
	{
		double[][] m = new double[n][n];
		for (int i = 0; i < n; i++)
			m[i][i] = 1.0;
		return m;
	}

	/**
	 * Joins matrices by rows
	 * WARNING: This routine DOES NOT check the dimensions!
	 * @param A
	 * @param B
	 * @return
	 */
	public static final double[][] joinMatricesByRows(double[][] A, double[][]... B)
	{
		if (B == null)
			return A;
		int numBMatrices = B.length;
		if (A == null)
		{
			if (numBMatrices == 1)
				return B[0];
			double[][][] newB = new double[numBMatrices-1][][];
			System.arraycopy(B, 1, newB, 0, numBMatrices-1);
			return joinMatricesByRows(B[0], newB);
		}
	
		int
			lenA = A.length,
			lenB = 0;
		for (int i = 0; i < numBMatrices; i++)
			if (B[i] != null)
				lenB += B[i].length; 
		double[][] result = new double[lenA + lenB][];
	
		for (int i = 0; i < lenA; i++)
			result[i] = A[i];
		for (int i = 0, idx = lenA; i < numBMatrices; i++)
		{
			double[][] Bi = B[i];
			if (Bi == null)
				continue;
			int lenBi = Bi.length;
			for (int j = 0; j < lenBi; j++, idx++)
				result[idx] = Bi[j];
		}
	
		return result;
	}

	/**
	 * Join matrices by columns.
	 * WARNING: This routine DOES NOT check the dimensions!
	 * @param A
	 * @param B
	 * @return
	 */
	public static final double[][] joinMatricesByColumns(double[][] A, double[][]... B)
	{
		if (B == null)
			return A;
		int numBMatrices = B.length;
		if (A == null)
		{
			if (numBMatrices == 1)
				return B[0];
			double[][][] newB = new double[numBMatrices-1][][];
			System.arraycopy(B, 1, newB, 0, numBMatrices-1);
			return joinMatricesByColumns(B[0], newB);
		}
	
		int
			n = A.length,
			lenA = A[0].length,
			lenB = 0;
		for (int i = 0; i < numBMatrices; i++)
			if (B[i] != null)
				lenB += B[i][0].length; 
		double[][] result = new double[n][lenA + lenB];
	
		for (int i = 0; i < n; i++)
			System.arraycopy(A[i], 0, result[i], 0, lenA);
		for (int i = 0, idx = lenA; i < numBMatrices; i++)
		{
			double[][] Bi = B[i];
			if (Bi == null)
				continue;
			int lenBi = Bi[0].length;
			if (lenBi > 0)
				for (int j = 0; j < n; j++)
					System.arraycopy(Bi[j], 0, result[j], idx, lenBi);
			idx += lenBi;
		}
	
		return result;
	}

	/**
	 * Accumulate values in the array
	 * @param values
	 * @param reverse
	 * @return
	 */
	public static final double[][] accumulate(double[][] values, boolean reverse)
	{
		int
			dim = values.length,
			n = values[0].length;
		double[][] cumulative = new double[dim][n];
		if (reverse)
		{
			System.arraycopy(values[dim - 1], 0, cumulative[dim - 1], 0, n);
			for (int i = dim - 2; i >= 0; i--)
			{
				double[]
					values_i = values[i],
					cum_prev = cumulative[i + 1],
					cum_i = cumulative[i];
				for (int j = 0; j < n; j++)
					cum_i[j] = values_i[j] + cum_prev[j];
			}
		}
		else
		{
			System.arraycopy(values[0], 0, cumulative[0], 0, n);
			for (int i = 1; i < dim; i++)
			{
				double[]
					values_i = values[i],
					cum_prev = cumulative[i - 1],
					cum_i = cumulative[i];
				for (int j = 0; j < n; j++)
					cum_i[j] = values_i[j] + cum_prev[j];
			}
		}
		return cumulative;
	}

	/**
	 * Take differences of the values in the array
	 * @param values
	 * @param reverse
	 * @return
	 */
	public static final double[][] difference(double[][] values, boolean reverse)
	{
		int
			dim = values.length,
			n = values[0].length;
		double[][] diff = new double[dim][n];
		if (reverse)
		{
			System.arraycopy(values[dim - 1], 0, diff[dim - 1], 0, n);
			for (int i = dim - 2; i >= 0; i--)
			{
				double[]
					diff_i = diff[i],
					res_i = values[i],
					res_next = values[i + 1];
				for (int j = 0; j < n; j++)
					diff_i[j] = res_i[j] - res_next[j];
			}
		}
		else
		{
			System.arraycopy(values[0], 0, diff[0], 0, n);
			for (int i = 1; i < dim; i++)
			{
				double[]
					diff_i = diff[i],
					res_i = values[i],
					res_next = values[i - 1];
				for (int j = 0; j < n; j++)
					diff_i[j] = res_i[j] - res_next[j];
			}
		}
		return diff;
	}

	public static final double[][] product(double[][] values, boolean reverse)
	{
		int
			dim = values.length,
			n = values[0].length;
		double[][] prods = new double[dim][n];
		if (reverse)
		{
			System.arraycopy(values[dim - 1], 0, prods[dim - 1], 0, n);
			for (int i = dim - 2; i >= 0; i--)
			{
				double[]
					values_i = values[i],
					cum_prev = prods[i + 1],
					cum_i = prods[i];
				for (int j = 0; j < n; j++)
					cum_i[j] = values_i[j] * cum_prev[j];
			}
		}
		else
		{
			System.arraycopy(values[0], 0, prods[0], 0, n);
			for (int i = 1; i < dim; i++)
			{
				double[]
					values_i = values[i],
					cum_prev = prods[i - 1],
					cum_i = prods[i];
				for (int j = 0; j < n; j++)
					cum_i[j] = values_i[j] * cum_prev[j];
			}
		}
		return prods;
	}

	/**
	 * Given a 1D array of double, creates a column vector of class Matrix
	 */
	public static final Matrix array2ColumnVector(double[] data)
	{
		int
			len = data.length,
			rows = len,
			cols = 1,
			offs = 0;
		
		double[][] d2 = new double[rows][cols];
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				d2[i][j] = data[offs++];
		return new Matrix(d2);
	}

	/**
	 * <P>Return pairwise row-based Hadamard product. Let D_i denote a
	 * vector from row i of matrix D of (m x n) dimension.
	 * Hadamard product is defined as: D_i # D_j = { d_ik * d_jk }
	 * with k = 1, ..., n
	 * 
	 * <P>This routine returns the pairwise Hadamard products arranged as
	 * follows: { D_1 # D_1, D_1 # D_2, ..., D_1 # D_m, D_2 # D_2, D_2 # D_3, ... ,
	 * D_2 # D_m, D_3 # D_3, ... , D_m # D_m }
	 * 
	 * <P>So the resulting matrix will have m(m + 1)/2 rows.
	 * 
	 * @param matrix
	 * @return
	 */
	public static final double[][] hadamardize(double[][] matrix)
	{
		int
			numRows = matrix.length,
			numCols = matrix[0].length,
			rowIndex = 0;
		double[][] result = new double[numRows*(numRows + 1)/2][numCols];
	
		for (int i = 0; i < numRows; i++)
		{
			double[] d_i = matrix[i];
			for (int j = i; j < numRows; j++, rowIndex++)
			{
				double[]
					d_j = matrix[j],
					vector = result[rowIndex];
				for (int k = 0; k < numCols; k++)
					vector[k] = d_i[k] * d_j[k];
			}
		}
		return result;
	}

	/**
	 * <P>Calculate the array index of D_i # D_j from the Hadamard matrix
	 * produced by hadamardize routine.
	 * 
	 * <P>Since the pairwise Hadamard products arranged like this:
	 * { D_1 # D_1, D_1 # D_2, ..., D_1 # D_m, D_2 # D_2, D_2 # D_3, ... ,
	 * D_2 # D_m, D_3 # D_3, ... , D_m # D_m }
	 * 
	 * <P>Thus D_1 # D_1 will be stored at index 0, D_2 # D_2 will be
	 * stored at index m, D_3 # D_3 will be stored at index m + m - 1,
	 * and so forth. Thus, D_i # D_i will be stored at index:<br>
	 * sum_{k = m - i + 2}^{m} k = sum_{k = 1}^{m} k - sum_{k = 1}^{m - i + 1} = m(m+1)/2 - (m-i+1)(m-i+2)/2 = (m(m+1) - (m-i+1)(m-i+2))/2
	 * 
	 * <P>It follows that D_i # D_j will be stored at index:<br>
	 * (m(m + 1) - (m - i + 1)( m - i + 2))/2 + (j - i) + 1
	 * 
	 * <P>After some algebra, this expression simplifies to: (m - i/2)(i - 1) + j
	 * 
	 * <P>But since we have a zero-based array index, and i and j are expressed
	 * in zero-based index as well, we substitute i for i + 1 and j for j + 1 and then
	 * subtract the entire expression from 1. [DO YOU MEAN REDUCE THE ENTIRE EXPRESSION
	 * BY 1? CN 10.17.07] The formula then becomes:<br>
	 * (m - (i + 1)/2)i + j = (2m - i - 1)* i/2 + j
	 * 
	 * @param numRows
	 * @param i
	 * @param j
	 * @return
	 */
	public static final int calcHadamardIndex(int numRows, int i, int j)
	{
		assert (i < numRows && j < numRows);
		if (j < i) // Swap it
		{
			int temp = i; i = j; j = temp;
		}
		return j + (2 * numRows - i - 1) * i / 2;
	}

	/**
	 * <P>Cholesky decomposition decomposes matrix A into LL'. Where L is a lower
	 * triangular matrix. Cholesky decomposition requires A to be symmetric
	 * (and hence square matrix) and positive definite.
	 * 
	 * <P>This decomposition is very useful for linear models (Y = Xb), especially when
	 * we want to solve multiple Y values at once for a given X, which is great
	 * for multiple-trait analysis.
	 * 
	 * <P>In single column of Y, the solution of the linear model, beta (or just b), is defined as:<br>
	 * beta = (X' X)^-1 X'Y
	 * 
	 * <P>However, with multiple columns of Y values, beta becomes:<br>
	 * beta = (X'V^-1 X)^-1 X'V^-1 Y
	 * 
	 * <P>with V equal to the variance-covariance matrix of the Y columns.
	 * 
	 * <P>We certainly can solve beta simply using the formula above, which is
	 * commonly referred as the normal equation. However, the result of normal
	 * equation is numerically unstable (and hence inaccurate) due to rounding
	 * errors.
	 * 
	 * <P>To overcome the numerical stability problem of the normal equation, we
	 * use what is commonly use QR decomposition. See <tt>doGramSchmidt</tt> comment
	 * for details.
	 * 
	 * <P>Having beta solution defined as above in multiple columns of Y values,
	 * we cannot use QR decomposition directly on X.
	 * 
	 * <P>It turns out that the solution is pretty simple. If the variance covariance
	 * matrix V is decomposed into Q and Q' (i.e. V = QQ'), the usual (Y = Xb)
	 * linear model can be rewritten as: Q^-1 Y = Q^-1 X b
	 * 
	 * <P>Thus, the solution is:<br>
	 * Q^-1 Y = Q^-1 X b<br>
	 * X' Q^-1' Q^-1 Y = X' Q^-1' Q^-1 X b (multiply each side by X' Q^-1')<br>
	 * X' V^-1 Y = X' V^-1 X b (because V = QQ', V^-1 = Q^-1'Q^-1)<br>
	 * b = (X'V^-1 X)^-1 X'V^-1 Y
	 * 
	 * <P>Thus, it's evident that we must do the QR decomposition on (Q^-1 X) instead
	 * of just X. In addition, we must do the computation on (Q^-1 Y) instead of just Y.
	 * 
	 * <P>Cholesky decomposition is a method to calculate the Q matrix needed here.
	 * Since V is definitely symmetric and positive definite, Cholesky decomposition
	 * should work perfectly for this purpose.
	 * 
	 * <P>The result of Cholesky decomposition of a matrix X is a lower triangular matrix
	 * L, whose entry l_ij is defined as follows:
	 * 
	 * <ul>
	 * <li>l_ii = sqrt[ x_ii - sum_{k=1}^{i-1} (l_ik)^2 ]</li>
	 * <li>l_ij = (1/l_jj) * [ x_ij - sum_{k=1}^{j-1} (l_ik * l_jk) ], for i > j</li>
	 * </ul>
	 * 
	 * <P>For details, see: http://en.wikipedia.org/wiki/Cholesky_decomposition
	 * 
	 * <P>We do Cholesky-Banachiewicz routine here, which starts from the upper left and
	 * move from row to row, top to bottom.
	 * 
	 * @param X
	 * @return The L matrix.
	 */
	public static final double[][] computeCholeskyDecomposition(double[][] X)
	{
		int dim = X.length;
		assert (dim == X[0].length);
		double[][] L = new double[dim][dim];
		double sum;
	
		for (int i = 0; i < dim; i++)
		{
			double[]
				Li = L[i],
				Xi = X[i];
			for (int j = 0; j < i; j++)
			{
				double[] Lj = L[j];
				sum = 0;
				for (int k = 0; k < j; k++)
					sum += Li[k] * Lj[k];
				Li[j] = (Xi[j] - sum) / Lj[j];
			}
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += Li[k] * Li[k];
			Li[i] = sqrt(Xi[i] - sum);
		}
		return L;
	}

	/**
	 * Vector V times matrix M
	 * @param V (n elements)
	 * @param M (n x m)
	 * @return
	 */
	public static final double[] vectorTimes(double[] V, double[][] M)
	{
		int
			n = V.length,
			m = M[0].length;
		assert (n == M.length);
		double[] result = new double[m];
		for (int i = 0; i < m; i++)
		{
			double sum = 0;
			for (int j = 0; j < n; j++)
				sum += M[j][i] * V[j];
			result[i] = sum;
		}
		return result;
	}

	/**
	 * Compute X - Y
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[] vectorMinus(double[] X, double[] Y)
	{
		int n = X.length;
		assert (n == Y.length);
		double[] result = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = X[i] - Y[i];
		return result;
	}

	/**
	 * X + Y
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[] vectorPlus(double[] X, double[] Y)
	{
		int n = X.length;
		assert (n == Y.length);
		double[] result = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = X[i] + Y[i];
		return result;
	}

	/**
	 * X + Y
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final int[] vectorPlus(int[] X, int[] Y)
	{
		int n = X.length;
		assert (n == Y.length);
		int[] result = new int[n];
		for (int i = 0; i < n; i++)
			result[i] = X[i] + Y[i];
		return result;
	}

	/**
	 * Compute XX'
	 * @param X
	 * @return
	 */
	public static final double[][] vectorSquare(double[] X)
	{
		int n = X.length;
		double[][] result = new double[n][n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				result[i][j] = X[i] * X[j];
		return result;
	}

	/**
	 * Permute vector vec
	 * @param vec
	 * @param idx
	 * @return
	 */
	public static final double[] vectorPermute(double[] vec, int[] idx)
	{
		int n = idx.length;
		double[] newX = new double[n];
		for (int i = 0; i < n; i++)
			newX[i] = vec[idx[i]];
		return newX;
	}

	/**
	 * Like vectorPermute, but vec is clobbered
	 * @param vec
	 * @param idx
	 */
	public static final void vectorPermuteInPlace(double[] vec, int[] idx)
	{
		int n = idx.length;
		for (int i = 0; i < n; i++)
		{
			int j = idx[i];
			if (j == i)
				continue;
			if (j < 0)
			{
				idx[i] = -j;
				continue;
			}
			double t = vec[i];
			int k = i;
			do {
				vec[k] = vec[j];
				k = j;
				j = idx[k];
				idx[k] = -j;
			} while (j > i);
			vec[k] = t;
		}
	}

	/**
	 * Invert permutation vector (I7PNVR)
	 * @param vec
	 * @return
	 */
	public static final double[] invertVectorPermute(int[] vec)
	{
		int n = vec.length;
		double[] newX = new double[n];
		for (int i = 0; i < n; i++)
			newX[vec[i]] = i;
		return newX;
	}

	public static final void shiftVector(int[] vec, int k)
	{
		int
			n = vec.length,
			temp, nmk, t;
		if (abs(k) >= n)
			return;
		if (k < 0)
		{
			k = -k;
			temp = vec[n-1];
			nmk = n - k;
			for (int i = 0; i < nmk; i++)
			{
				t = n - i;
				vec[t + 1] = vec[t];
			}
			vec[k] = temp;
		}
		else
		{
			nmk = n - 1;
			temp = vec[k];
			for (int i = k; i < nmk; i++)
				vec[i] = vec[i+1];
			vec[nmk] = temp;
		}
	}
	/**
	 * Component-wise multiplication of vectors X and Y.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[] vectorTimes(double[] X, double[] Y)
	{
		int n = X.length;
		double[] result = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = X[i] * Y[i];
		return result;
	}

	/**
	 * Component-wise division of vectors X and Y, i.e. Z[i] = X[i] / Y[i];
	 * 
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[] vectorDivide(double[] X, double[] Y)
	{
		int n = X.length;
		double[] result = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = X[i] / Y[i];
		return result;
	}

	/**
	 * BLAS routine aX + Y
	 * @param a
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[] vectorDaXpY(double a, double[] X, double[] Y)
	{
		int n = X.length;
		double[] result = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = a * X[i] + Y[i];
		return result;
	}

	/**
	 * Calculate the dot product of vectors X and Y.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double calculateDotProduct(double[] X, double[] Y)
	{
		int n = X.length;
		double result = 0;
		for (int i = 0; i < n; i++)
			result += X[i] * Y[i];
		return result;
	}

	/**
	 * Convert a single-column matrix into an array
	 * @param matrix
	 * @return
	 */
	public static final double[] extractToVector(Matrix matrix)
	{
		double[][] array = matrix.getArray();
		int n = array.length;
		double[] result = new double[n];
	
		for (int i = 0; i < n; i++)
			result[i] = array[i][0];
		return result;
	}

	/**
	 * Calculate norm (i.e. X'X)
	 * @param X
	 * @return
	 */
	public static final double calculateNorm2(double[] X)
	{
		int n = X.length;
		double
			nrm = 0,
			scale = 0,
			ssq = 1;
		for (int i = 0; i < n; i++) {
			double absxi = Math.abs(X[i]);
			if (scale < absxi) {
				nrm = scale / absxi;
				ssq = 1 + ssq * nrm * nrm;
				scale = absxi;
			} else {
				nrm = absxi / scale;
				ssq += nrm * nrm;
			}
		}
		return scale * Math.sqrt(ssq);
	}

	/**
	 * Extract diagonal elements out of square array V
	 * @param V
	 * @return
	 */
	public static final double[] extractDiagonalElements(double[][] V)
	{
		int n = V.length;
		assert(n == V[0].length);
		double[] result = new double[n];
		for (int i = 0; i < n; i++)
			result[i] = V[i][i];
		return result;
	}

	/**
	 * Convert an array of values into a square matrix with the values as its diagonal elements
	 * @param V
	 * @return
	 */
	public static final double[][] diagonalize(double[] V)
	{
		int n = V.length;
		double[][] result = new double[n][n];
		for (int i = 0; i < n; i++)
			result[i][i] = V[i];
		return result;
	}

	/**
	 * Emulating R's <tt>seq</tt>.
	 * @param from
	 * @param to
	 * @param inc
	 * @return
	 */
	public static final double[] seq(double from, double to, double inc)
	{
		if ((from < to && inc < 0) || (from > to && inc > 0) || (from == to))
			return new double[] { from };
		int numElements = 1 + (int) Math.floor(Math.abs(to - from) / inc + 1e-7);
		double[] result = new double[numElements];
		if (from > to && inc > 0)
			inc = -inc;
		for (int i = 0; i < numElements; i++)
			result[i] = from + inc * i;
		return result;
	}

	/**
	 * Create a sub matrix out of X
	 * @param r row indices of X to include
	 * @param c column indices of X to include
	 * @param X
	 * @return
	 */
	public static final double[][] createSubMatrix(int[] r, int[] c, double[][] X)
	{
		int
			nrows = r.length,
			ncols = c.length;
		double[][] m = new double[nrows][ncols];
		for (int i = 0; i < nrows; i++) {
			double[]
				xi = X[r[i]], // RJ's speed up
				mi = m[i];
			for (int j = 0; j < ncols; j++)
				mi[j] = xi[c[j]];
		}
		return m;
	}

	/**
	 * Create sub matrix of X by column selection
	 * @param c column indices of X to include
	 * @param X
	 * @return
	 */
	public static final double[][] createSubMatrixByColumns(int[] c, double[][] X)
	{
		int
			nrows = X.length,
			ncols = c.length;
		double[][] m = new double[nrows][ncols];
		for (int i = 0; i < nrows; i++) {
			double[]
				xi = X[i], // RJ's speed up
				mi = m[i];
			for (int j = 0; j < ncols; j++)
				mi[j] = xi[c[j]];
		}
		return m;
	}

	/**
	 * Create sub matrix of X by row selection
	 * @param r row indices of X to include
	 * @param X
	 * @return
	 */
	public static final double[][] createSubMatrixByRows(int[] r, double[][] X)
	{
		int
			nrows = r.length,
			ncols = X[0].length;
		double[][] m = new double[nrows][ncols];
		for (int i = 0; i < nrows; i++)
			System.arraycopy(X[r[i]], 0, m[i], 0, ncols);
		return m;
	}

	/**
	 * Create a sub vector out of v
	 * @param c indices to include
	 * @param v
	 * @return
	 */
	public static final double[] createSubVector(int[] c, double[] v)
	{
		int ncols = c.length;
		double[] mi = new double[ncols];
		for (int j = 0; j < ncols; j++)
			mi[j] = v[c[j]];
		return mi;
	}

	/**
	 * Delete matrix columns
	 * @param X
	 * @param deletedIdx
	 * @return
	 */
	public static final double[][] deleteMatrixColumn(double[][] X, int[] deletedIdx)
	{
		int[] idx = calcIdxToKeep(deletedIdx, X[0].length);
		if (idx == null)
			throw new RuntimeException();
		return createSubMatrixByColumns(idx, X);
	}

	/**
	 * Delete single column from X
	 * @param X
	 * @param colNo must be a valid index. Else, there will be an exception
	 * @return
	 */
	public static final double[][] deleteMatrixColumn(double[][] X, int colNo)
	{
		int
			nrows = X.length,
			ncols = X[0].length,
			newncols = ncols - 1;
		if (colNo < 0 || colNo >= ncols)
			throw new RuntimeException(); // I told you
		double[][] result = new double[nrows][newncols];
		if (colNo == 0)
			for (int i = 0; i < nrows; i++)
				System.arraycopy(X[i], 1, result[i], 0, newncols);
		else if (colNo == newncols)
			for (int i = 0; i < nrows; i++)
				System.arraycopy(X[i], 0, result[i], 0, newncols);
		else
		{
			int numRemainingCol = newncols - colNo;
			for (int i = 0; i < nrows; i++)
			{
				System.arraycopy(X[i], 0, result[i], 0, colNo);
				System.arraycopy(X[i], colNo + 1, result[i], colNo, numRemainingCol);
			}
		}
		return result;
	}

	/**
	 * Utility function to determine which column indices will be kept
	 * @param deletedIdx indices of deleted columns
	 * @param numColumns the number of columns
	 * @return null if there is nothing to delete
	 */
	public static final int[] calcIdxToKeep(int[] deletedIdx, int numColumns)
	{
		if (deletedIdx == null || deletedIdx.length == 0)
			return null;
		Set<Integer> deletedIdxSet = new TreeSet<Integer>();
		for (int i: deletedIdx)
		{
			if (i >= 0 && i < numColumns)
				deletedIdxSet.add(i);
		}
		if (deletedIdxSet.size() == 0) // Nothing to delete
			return null;
		int
			numColsToDelete = deletedIdx.length,
			numNewColumns = numColumns - numColsToDelete;
		int[] idxToKeep = new int[numNewColumns];
		for (int i = 0, j = 0; i < numColumns; i++)
		{
			if (!deletedIdxSet.contains(i))
				idxToKeep[j++] = i;
		}
		return idxToKeep;
	}

	/*
	 * A method needed for evaluating polynomials. CN 1.27.06
	 */
	public static final void fillOutPowersVector(double d, double[] powers)
	{
		double term;
		powers[0] = term = 1.0;
		for (int i = 1; i < powers.length; i++)
			powers[i] = term *= d;
	}

	/**
	 * Calculate X - XYY'. Eliminate the need to transpose.
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXMinusXYYt(double[][] X, double[][] Y)
	{
		int
			rows = X.length,
			n = X[0].length,
			p = Y[0].length;
		assert (n == Y.length);
		double[][] result = new double[rows][n];
		double[] temp = new double[p];
		double sum;
	
		for (int i = 0; i < rows; i++)
		{
			double[] xi = X[i];
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += xi[k] * Y[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = Y[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				result[i][j] = xi[j] - sum;
			}
		}
		return result;
	}

	public static final double[] calculateXMinusXYYt(double[] X, double[][] Y)
	{
		int
			n = Y.length,
			rows = X.length / n,
			p = Y[0].length;
		double[] result = new double[rows*n];
		double[] temp = new double[p];
		double sum;
	
		for (int i = 0; i < rows; i++)
		{
			int i_n = i * n;
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += X[i_n + k] * Y[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = Y[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				result[i_n + j] = X[i_n + j] - sum;
			}
		}
		return result;
	}

	public static final void calculateXMinusXYYtInPlace(double[] X, double[][] Y)
	{
		int
			n = Y.length,
			rows = X.length / n,
			p = Y[0].length;
		double[] temp = new double[p];
		double sum;

		for (int i = 0; i < rows; i++)
		{
			int i_n = i * n;
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += X[i_n + k] * Y[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = Y[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				X[i_n + j] = X[i_n + j] - sum;
			}
		}
	}

	/**
	 * Calculate X - YY'X. Eliminate the need to transpose.
	 * FIXME: SLOW!
	 * @param X
	 * @param Y
	 * @return
	 */
	public static final double[][] calculateXMinusYYtX(double[][] X, double[][] Y)
	{
		int
			n = X.length,
			m = X[0].length,
			p = Y[0].length;
		assert (n == Y.length);
		double[][] result = new double[n][m];
		double[] temp = new double[p];
		double sum;

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += X[k][i] * Y[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = Y[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				result[j][i] = X[j][i] - sum;
			}
		}
		return result;
	}

	public static final double[] calculateXMinusYYtX(double[] X, double[][] Y)
	{
		int
			n = Y.length,
			p = Y[0].length,
			m = X.length / n;
		assert (m*n == X.length);
		double[] result = new double[n*m];
		double[] temp = new double[p];
		double sum;
	
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += X[k*m+i] * Y[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = Y[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				result[j*m + i] = X[j*m+i] - sum;
			}
		}
		return result;
	}

	public static final void calculateXMinusYYtXInPlace(double[] X, double[][] Y)
	{
		int
			n = Y.length,
			p = Y[0].length,
			m = X.length / n;
		assert (m*n == X.length);
		double[] temp = new double[p];
		double sum;
	
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += X[k*m+i] * Y[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = Y[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				X[j*m + i] = X[j*m+i] - sum;
			}
		}
	}

	public static final double[][] standardizeByRowKeepSD(double[][] X, double[] sd) {
		int
			nrow = X.length,
			ncol = X[0].length;
		assert(nrow == sd.length);
		double[][] result = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++) {
			double sum = 0, sumsq = 0;
			double[] xi = X[i], resulti = result[i];
			for (int j = 0; j < ncol; j++) {
				double v = xi[j];
				sum += v; sumsq += v*v;
			}
			sumsq = sqrt((sumsq - sum * sum / ncol) / (ncol-1));
			sum /= ncol;
			sd[i] = sumsq;
			for (int j = 0; j < ncol; j++)
				resulti[j] = (xi[j] - sum) / sumsq;
		}
		return result;
	}

	public static final void standardizeByRowKeepSDInPlace(double[][] X, double[] sd) {
		int
			nrow = X.length,
			ncol = X[0].length;
		assert(nrow == sd.length);
		for (int i = 0; i < nrow; i++) {
			double sum = 0, sumsq = 0;
			double[] xi = X[i];
			for (int j = 0; j < ncol; j++) {
				double v = xi[j];
				sum += v; sumsq += v*v;
			}
			sumsq = sqrt((sumsq - sum * sum / ncol) / (ncol-1));
			sum /= ncol;
			sd[i] = sumsq;
			for (int j = 0; j < ncol; j++)
				xi[j] = (xi[j] - sum) / sumsq;
		}
	}

	public static final double[][] standardizeByColKeepSD(double[][] X, double[] sd) {
		int
			nrow = X.length,
			ncol = X[0].length;
		assert(ncol == sd.length);
		double[][] result = new double[nrow][ncol];
		for (int i = 0; i < ncol; i++) {
			double sum = 0, sumsq = 0;
			for (int j = 0; j < nrow; j++) {
				double v = X[j][i];
				sum += v; sumsq += v*v;
			}
			sumsq = sqrt((sumsq - sum * sum / nrow) / (nrow-1));
			sum /= nrow;
			sd[i] = sumsq;
			for (int j = 0; j < nrow; j++)
				result[j][i] = (X[j][i] - sum) / sumsq;
		}
		return result;
	}

	public static final void standardizeByColKeepSDInPlace(double[][] X, double[] sd) {
		int
			nrow = X.length,
			ncol = X[0].length;
		assert(ncol == sd.length);
		for (int i = 0; i < ncol; i++) {
			double sum = 0, sumsq = 0;
			for (int j = 0; j < nrow; j++) {
				double v = X[j][i];
				sum += v; sumsq += v*v;
			}
			sumsq = sqrt((sumsq - sum * sum / nrow) / (nrow-1));
			sum /= nrow;
			sd[i] = sumsq;
			for (int j = 0; j < nrow; j++)
				X[j][i] = (X[j][i] - sum) / sumsq;
		}
	}

	public static final double[][] createRandomDoubleMatrix(int nrow, int ncol, RandomEngine random) {
		double[][] mat = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
			for (int j = 0; j < ncol; j++)
				mat[i][j] = random.nextDouble();
		return mat;
	}

	public static final double[][] createRandomBinaryMatrix(int nrow, int ncol, RandomEngine random) {
		double[][] mat = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
			for (int j = 0; j < ncol; j++)
				mat[i][j] = random.nextDouble() >= 0.5 ? 1 : 0;
		return mat;
	}

	public static final double calculateSSE(float[][] mat1, float[][] mat2) {
		double sum = 0;
		int
		nrow = mat1.length,
		ncol = mat1[0].length;
		for (int i = 0; i < nrow; i++)
			for (int j = 0; j < ncol; j++) {
				double v = (mat1[i][j] - mat2[i][j]);
				sum += v * v;
			}
		return sum;
	}

	public static final double calculateSSE(double[][] mat1, double[][] mat2) {
		double sum = 0;
		int
			nrow = mat1.length,
			ncol = mat1[0].length;
		for (int i = 0; i < nrow; i++)
			for (int j = 0; j < ncol; j++) {
				double v = (mat1[i][j] - mat2[i][j]);
				sum += v * v;
			}
		return sum;
	}

	static void testHadamardize()
	{
		double[][]
			matrix = new double[][] { {1, 0, -1}, {-0.5, 0.5, -0.5} },
			result = hadamardize(matrix);
		System.out.println(QStringUtils.toStringAsMatrix(result));
	}

	static void testCholesky()
	{
		double[][]
			X = new double[][] {{25, -5, 10}, {-5, 17, 10}, {10, 10, 62}},
			L = computeCholeskyDecomposition(X);
		System.out.println(QStringUtils.toStringAsMatrix(L));
		System.out.println(QStringUtils.toStringAsMatrix(calculateXtY(L, L)));
		System.out.println(QStringUtils.toStringAsMatrix(calculateXYt(L, L)));
		System.out.println(QStringUtils.toStringAsMatrix(calculateXY(L, L)));
	}

	static void testQuadraticForm()
	{
		// Oxygen uptake example from Dr. Paul Nelson's STATS 705 handout 12.1
		double
			X[][] = new double[][]
			{
				{8.4, 132.0, 29.1, 14.4},
				{8.7, 135.5, 29.7, 14.5},
				{8.9, 127.7, 28.4, 14.0},
				{9.9, 131.1, 28.8, 14.2},
				{9.0, 130.0, 25.9, 13.6},
				{7.7, 127.6, 27.6, 13.9},
				{7.3, 129.9, 29.0, 14.0},
				{9.9, 138.1, 33.6, 14.6},
				{9.3, 126.6, 27.7, 13.9},
				{8.1, 131.8, 30.8, 14.5}
			},
			V[][];
		Matrix
			matX = new Matrix(X),
			matXt = matX.transpose(),
			XX = matXt.times(matX);

		V = XX.inverse().getArray();
		double[] vv = extractDiagonalElements(V);
		System.out.println(QStringUtils.toStringAsMatrix(V));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(calculateQuadraticForm(X, V)));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(calculateQuadraticFormTranspose(matXt.getArray(), V)));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(calculateQuadraticFormDiagV(X, vv)));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(calculateQuadraticForm(X, diagonalize(vv))));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(XX.getArray()));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(calculateQuadraticFormTranspose(X)));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(calculateQuadraticForm(X)));
		System.out.println();
		System.out.println(QStringUtils.toStringAsMatrix(matX.times(matXt).getArray()));
		System.out.println();
	}

	public static void main(String[] args)
	{
		/*
		int
			numRow = 1500,
			numCol = 1500;
		double[][] arr = new double[numRow][numCol];
		for (int i = 0; i < numRow; i++)
			for (int j = 0; j < numCol; j++)
				arr[i][j] = 1.0 / (i + j + 1);
		long time1, time2;
		System.out.println("Benchmark 2 started"); //$NON-NLS-1$
		time1 = System.currentTimeMillis();
		QMatrixUtils.calculateXYtParallel(arr, arr);
		time2 = System.currentTimeMillis();
		System.out.println("Time = " + (time2 - time1)); //$NON-NLS-1$

		System.out.println("Running garbage collector..."); //$NON-NLS-1$
		QSystemUtils.runGCAggressively();
		System.out.println("Benchmark 1 started"); //$NON-NLS-1$
		time1 = System.currentTimeMillis();
		QMatrixUtils.calculateXYt(arr, arr);
		time2 = System.currentTimeMillis();
		System.out.println("Time = " + (time2 - time1)); //$NON-NLS-1$
		//*/

		//System.out.println(qutils.QStringUtils.toString(seq(0, 1.05, 0.1)));
		//testQuadraticForm();
		System.exit(0);
	}
}
