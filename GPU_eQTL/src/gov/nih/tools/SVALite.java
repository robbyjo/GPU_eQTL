/*
 * Roby Joehanes
 * 
 * Copyright 2013 Roby Joehanes
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
package gov.nih.tools;

import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import com.csvreader.CsvReader;

import jdistlib.F;
import jdistlib.Normal;
import jdistlib.math.density.Density;
import jdistlib.math.spline.SmoothSpline;
import jdistlib.math.spline.SmoothSplineResult;
import jdistlib.rng.MersenneTwister;
import jdistlib.util.Utilities;
import static java.lang.Math.*;

/**
 * @author Roby Joehanes
 */
public class SVALite {
	private static final int kNumCPUCores = Runtime.getRuntime().availableProcessors();
	private static final double[][] parallelResidualize(double[][] X, double[][] Y) {
		int
			rows = X.length,
			cols = X[0].length,
			numThreads = kNumCPUCores > rows ? rows : kNumCPUCores;
		ExecutorService threadPool = Executors.newFixedThreadPool(numThreads + 1);
		SynchronizedCounter counter = new SynchronizedCounter(0, rows);
		double[][] result = new double[rows][cols];
		assert (X[0].length == Y.length);
		for (int i = 0; i < numThreads; i++)
			threadPool.execute(new ResidualizeJob(X, Y, result, counter));
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException exx) {}
		return result;
	}

	private static final double[][] parallelCrossProduct(double[][] X) {
		int
			cols = X[0].length,
			numThreads = kNumCPUCores > cols ? cols : kNumCPUCores;
		ExecutorService threadPool = Executors.newFixedThreadPool(numThreads + 1);
		SynchronizedCounter counter = new SynchronizedCounter(0, cols);
		double[][] result = new double[cols][cols];
		for (int i = 0; i < numThreads; i++)
			threadPool.execute(new CrossProductJob(X, result, counter));
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException exx) {}
		return result;
	}
	private static final double[] calcRSS(double[][] X, double[][] Y) {
		double[][] mat = parallelResidualize(X, new QR(Y).getQ());
		int n = mat.length, m = mat[0].length;
		double[] rss = new double[n];
		for (int i = 0; i < n; i++) {
			double sum = 0, v[] = mat[i];
			for (int j = 0; j < m; j++)
				sum += v[j] * v[j];
			rss[i] = sum;
		}
		mat = null;
		return rss;
	}

	private static double[] edge_lfdr(double[] p) {
		return edge_lfdr(p, 1.5, 0.8, 1e-8);
	}

	private static double[] edge_lfdr(double[] p, double adj, double lambda, double eps) {
		int n = p.length;
		double pi0 = 0, eps0 = 1-eps, x[] = new double[n];
		for (int i = 0; i < n; i++) {
			if (p[i] >= lambda) pi0++;
			if (p[i] < eps) p[i] = eps; else if (p[i] > eps0) p[i] = eps0;
			x[i] = Normal.quantile(p[i], 0, 1, true, false);
		}
		pi0 = min(1, (pi0/n)/(1-lambda));

		Density myd = Density.density(x, adj);
		SmoothSplineResult mys = SmoothSpline.fit(myd.x, myd.y);
		double[] lfdr = new double[n];
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			double val = pi0 * Normal.density(x[i], 0, 1, false) / SmoothSpline.predict(mys, x[i], 0);
			lfdr[i] = val > 1 ? 1 : val;
			order[i] = i;
		}
		double[] lfdr_copy = new double[n];
		System.arraycopy(lfdr, 0, lfdr_copy, 0, n);
		Utilities.sort(lfdr_copy, order);
		for (int i = 1; i < n; i++)
			lfdr_copy[i] = lfdr_copy[i] < lfdr[i - 1] ? lfdr[i - 1] : lfdr[i];
		for (int i = 0; i < n; i++)
			lfdr[order[i]] = lfdr_copy[i];
		return lfdr;
	}

	private static void shuffle(double[][] src, double[][] dest) {
		/*
		int
			nrow = src.length,
			numThreads = kNumCPUCores > nrow ? nrow : kNumCPUCores;
		ExecutorService threadPool = Executors.newFixedThreadPool(numThreads + 1);
		SynchronizedCounter counter = new SynchronizedCounter(0, nrow);
		for (int i = 0; i < numThreads; i++)
			threadPool.execute(new ShuffleJob(src, dest, counter));
		threadPool.shutdown();
		try {
			threadPool.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException exx) {}
		/*/
		MersenneTwister rng = new MersenneTwister();
		int nrow = src.length, ncol = src[0].length, k = ncol / 2;
		for (int i = 0; i < nrow; i++) {
			System.arraycopy(src[i], 0, dest[i], 0, ncol);
			double[] v = dest[i];
			for (int j = 0; j < k; j++) {
				int a = rng.nextInt(ncol), b = rng.nextInt(ncol);
				double temp = v[a]; v[a] = v[b]; v[b] = temp;
			}
		}
		//*/
	}

	private static double[] calcProportion(double[] d, int ndf) {
		double sum = 0;
		int n = d.length;
		for (int i = 0; i < n; i++) sum += d[i];
		double[] prop = new double[ndf];
		for (int i = 0; i < ndf; i++) prop[i] = d[i] / sum;
		return prop;
	}

	private static double[][] makeColumnVector(int n, double val) {
		double[][] mat = new double[n][1];
		for (int i = 0; i < n; i++) mat[i][0] = val;
		return mat;
	}

	static double[][] cbind(double[][] orig, double[][] uu) {
		int n = orig.length, p1 = orig[0].length, p2 = uu[0].length;
		double[][] result = new double[n][p1+p2];
		for (int i = 0; i < n; i++) {
			System.arraycopy(orig[i], 0, result[i], 0, p1);
			System.arraycopy(uu[i], 0, result[i], p1, p2);
		}
		return result;
	}

	private static void copyAt(double[][] orig, double[][] uu, int pos) {
		int n = orig.length, p = uu[0].length;
		for (int i = 0; i < n; i++) {
			System.arraycopy(uu[i], 0, orig[i], pos, p);
		}
	}

	/**
	 * Bootstrap version of num.sv
	 * @param data
	 * @param altModel
	 * @param svsig
	 * @param rng
	 * @return
	 */
	public static int calcNumSV(double[][] data, double[][] altModel, double svsig, int B) {
		if (svsig <= 0 || svsig > 1) throw new RuntimeException();
		QR qr = new QR(altModel);
		int
			nrow = data.length,
			ncol = data[0].length,
			ndf = ncol - qr.getRank();
		double[][]
			qq = qr.getQ(),
			residual = parallelResidualize(data, qq),
			res0 = new double[nrow][ncol];
		SVD svd = new SVD(parallelCrossProduct(residual), false, false);
		double[]
			dstat = calcProportion(svd.getSingularValues(), ndf),
			dstat0 = new double[ndf];

		for (int i = 0; i < B; i++) {
			System.out.print("Iteration " + i + ": Shuffle");
			shuffle(residual, res0);
			System.out.print(", residualize");
			res0 = parallelResidualize(res0, qq);
			System.out.print(", cross product");
			double[][] temp = parallelCrossProduct(res0);
			System.out.print(", SVD");
			svd = new SVD(temp, false, false);
			System.out.print(", counting");
			double[] newd = calcProportion(svd.getSingularValues(), ndf);
			for (int j = 0; j < ndf; j++)
				if (newd[j] >= dstat[j]) dstat0[j]++;
			System.out.println();
		}

		int numSVs = 0;
		for (int j = 0; j < ndf; j++) {
			dstat0[j] /= B;
			if (j > 0) dstat0[j] = max(dstat0[j-1], dstat0[j]);
			if (dstat0[j] <= svsig) numSVs++;
		}
		return numSVs;
	}

	public static double[][] irwSVA(double[][] data, double[][] altModel, double[][] nullModel, int numSVs, int B) {
		if (numSVs <= 0) throw new RuntimeException();
		int
			n = data[0].length,
			m = data.length,
			pa = altModel[0].length,
			p0 = nullModel == null ? 1: nullModel[0].length,
			pa_sv = pa + numSVs,
			p0_sv = p0 + numSVs;
		if (nullModel == null) nullModel = makeColumnVector(n, 1);
		double[][] resid = parallelResidualize(data, new QR(altModel).getQ());
		SVD svd = new SVD(parallelCrossProduct(resid), true, false);
		double[][]
			uu = svd.getV(numSVs),
			mod_b = new double[m][pa_sv],
			mod0_b = new double[m][p0_sv];
		double[]
			pprob_b_comp = new double[m],
			pprob_gam_comp = new double[m],
			rss_null = calcRSS(data, nullModel);
		copyAt(mod_b, altModel, 0);
		copyAt(mod0_b, nullModel, 0);

		for (int i = 0; i < B; i++) {
			copyAt(mod_b, uu, pa);
			copyAt(mod0_b, uu, p0);
			double[] rss1 = calcRSS(data, mod_b), rss0 = calcRSS(data, mod0_b);
			int df1 = pa_sv, df0 = p0_sv, dfr = df1 - df0, dfe = n - df1, dfr0 = df0 - p0, dfe0 = n - df0;
			for (int j = 0; j < m; j++) {
				pprob_b_comp[j] = F.cumulative(((rss0[j] - rss1[j]) / dfr) / (rss1[j] / dfe), dfr, dfe, false, false);
				pprob_gam_comp[j] = F.cumulative(((rss_null[j] - rss0[j]) / dfr0) / (rss0[j] / dfe0), dfr0, dfe0, false, false);
			}
			pprob_b_comp = edge_lfdr(pprob_b_comp);
			pprob_gam_comp = edge_lfdr(pprob_gam_comp);
			for (int j = 0; j < m; j++) {
				double v[] = data[j], w[] = resid[j], wt = pprob_b_comp[j] * (1 - pprob_gam_comp[j]), sum = 0;
				for (int k = 0; k < n; k++) sum += v[k];
				sum /= n;
				for (int k = 0; k < n; k++) w[k] = (v[k] - sum) * wt;
			}
			svd = new SVD(parallelCrossProduct(resid), true, false);
			uu = svd.getV(numSVs);
		}
		return uu;
	}

	static final void usage() {
		System.exit(0);
	}
	static final Map<String, String> parseArgs(String[] args) {
		if (args == null || args.length == 0) usage();
		Map<String, String> argMap = new HashMap<String, String>();
		int n = args.length;
		for (int i = 0; i < n; i++) {
			String tok = args[i];
			if (tok.startsWith("--")) {
				tok = tok.substring(2);
				int eqPos = tok.indexOf('=');
				if (eqPos < 0) {
					argMap.put(tok, "true");
				} else {
					String val = tok.substring(eqPos + 1);
					argMap.put(tok, val);
				}
			} else
				throw new RuntimeException("Unknown option "+tok);
		}
		return argMap;
	}

	public static final void main(String[] args) {
		Map<String, String> argMap = parseArgs(args);

		String
			dataFile = argMap.get("input"),
			altModelFile = argMap.get("alt"),
			nullModelFile = argMap.get("null"),
			outputFile = argMap.get("output"),
			numIterNumSVString = argMap.get("num-iter-numsv"),
			numIterMainString = argMap.get("num-iter-main"),
			svSigString = argMap.get("sv-sig"),
			numSVString = argMap.get("numsv");
		if (dataFile == null || altModelFile == null || outputFile == null) usage();
		if ((numIterNumSVString != null || svSigString != null) && numSVString != null) usage();

		int numIterMain = numIterMainString == null ? 5 : Integer.parseInt(numIterMainString);
		double svsig = svSigString == null ? 0.1 : Double.parseDouble(svSigString);
		int numIterNumSV = numIterNumSVString == null ? (int) ceil(2.0 / svsig) : Integer.parseInt(numIterNumSVString);
		int numSVs = numSVString == null ? -1 : Integer.parseInt(numSVString);
		if (svsig <= 0 || svsig >= 1 || numIterMain <= 0 || numIterNumSV < (int) ceil(2.0 / svsig) || numSVs <= 0) usage();

		Table altModelTable, nullModelTable = null, dataTable;
		try {
			System.out.println("Loading...");
			altModelTable = Table.load(altModelFile);
			if (nullModelFile != null) {
				nullModelTable = Table.load(nullModelFile);
				if (altModelTable.matrix.length != nullModelTable.matrix.length)
					throw new RuntimeException("Number of individuals does not match between alternative and null matrix tables!");
			}
			dataTable = Table.load(dataFile);
			System.out.println("Loading is done!");
			if (altModelTable.matrix.length != dataTable.matrix[0].length)
				throw new RuntimeException("Number of individuals does not match between alternative matrix table and input data table!");
			if (numSVs <= 0)
				numSVs = calcNumSV(dataTable.matrix, altModelTable.matrix, svsig, numIterNumSV);
			System.out.println("Num SVs = " + numSVs);

			double[][] sv = irwSVA(dataTable.matrix, altModelTable.matrix, nullModelTable == null ? null : nullModelTable.matrix, numSVs, numIterMain);
			int nrow = sv.length, ncol = numSVs;
			PrintWriter writer = new PrintWriter(outputFile);
			writer.print("ID");
			for (int i = 0; i < ncol; i++)
				writer.print(",SV" + (i+1));
			String[] rn = altModelTable.rowNames;
			for (int i = 0; i < nrow; i++) {
				if (rn != null) writer.print(rn[i]);
				double[] v = sv[i];
				for (int j = 0; j < ncol; j++)
					writer.print("," + v[j]);
			}
			writer.close(); writer = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class ResidualizeJob implements Runnable
{
	protected double[][] mX, mY, mResult;
	protected SynchronizedCounter mCounter;
	public ResidualizeJob(double[][] x, double[][] y, double[][] result, SynchronizedCounter counter) {
		mX = x; mY = y; mResult = result; mCounter = counter;
	}

	public void run() {
		int
			i = mCounter.next(),
			n = mX[0].length,
			p = mY[0].length;
		double sum;
		double[] temp = new double[p];
		if (i < 0) return;
		do {
			double[] xi = mX[i];
			for (int j = 0; j < p; j++) {
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += xi[k] * mY[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++) {
				double[] yj = mY[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				mResult[i][j] = xi[j] - sum;
			}
			i = mCounter.next();
		} while (i >= 0);
		mX = mY = mResult = null;
		mCounter = null;
	}
}

class CrossProductJob implements Runnable
{
	protected double[][] mX, mResult;
	protected SynchronizedCounter mCounter;
	public CrossProductJob(double[][] x, double[][] result, SynchronizedCounter counter) {
		mX = x; mResult = result; mCounter = counter;
	}

	public void run() {
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0) return;
		do {
			for (int j = 0; j <= i; j++) {
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += mX[k][i] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}
}

class ShuffleJob implements Runnable
{
	protected double[][] src, dest;
	protected SynchronizedCounter mCounter;
	protected MersenneTwister rng = new MersenneTwister();

	public ShuffleJob(double[][] x, double[][] result, SynchronizedCounter counter) {
		src = x; dest = result; mCounter = counter;
	}

	public void run() {
		int
			i = mCounter.next(),
			ncol = src[0].length,
			k = ncol / 2;
		if (i < 0) return;
		do {
			System.arraycopy(src[i], 0, dest[i], 0, ncol);
			double[] v = dest[i];
			for (int j = 0; j < k; j++) {
				int a = rng.nextInt(ncol), b = rng.nextInt(ncol);
				double temp = v[a]; v[a] = v[b]; v[b] = temp;
			}

			i = mCounter.next();
		} while (i >= 0);
	}
}

class SynchronizedCounter {
	protected int mStart, mStop, mCurrent;

	public SynchronizedCounter(int start, int stop) {
		mStart = mCurrent = start; mStop = stop;
	}

	public synchronized int next()
	{	return mCurrent < mStop ? mCurrent++ : -1; }

	public synchronized boolean hasNext()
	{	return mCurrent < mStop; }

	public int getStart()
	{	return mStart; }

	public int getStop()
	{	return mStop; }

	public void reset()
	{	mCurrent = mStart; }

	@Override
	public String toString()
	{	return mCurrent + "/" + mStop; } //$NON-NLS-1$
}

class Table {
	public static final int kUndefinedValue = Integer.MIN_VALUE;
	private static final String sLineFormat = "%d vs %d, at line %d"; //$NON-NLS-1$

	public String[] colNames, rowNames;
	public double[][] matrix;
	public Table() {}

	public Table(double[][] mat, String[] rn, String[] cn) {
		matrix = mat; rowNames = rn; colNames = cn;
	}

	/**
	 * Read a comma delimited table
	 * @param fn File name
	 * @return the table
	 * @throws IOException
	 */
	public static final Table load(String fn) throws IOException {
		return load(new FileReader(fn), new char[] { ','}, "#", true, true);
	}

	public static final Table load(String fn, char[] delim, String comments, boolean hasRowNames, boolean hasColNames) throws IOException {
		return load(new FileReader(fn), delim, comments, hasRowNames, hasColNames);
	}

	/**
	 * Read a delimited file (tab or comma or whatever)
	 * @param r Reader
	 * @param delim The delimiter
	 * @param comments Comment tokens
	 * @param hasRowNames set to true if the first column is the row names
	 * @param hasColNames set to true if the first row is the column names
	 * @return the table
	 * @throws IOException
	 */
	public static final Table load(Reader r, char[] delim, String comments, boolean hasRowNames, boolean hasColNames) throws IOException {
		CsvReader reader = null;
		List<double[]> list = new ArrayList<double[]>();
		List<String> rowNames = new ArrayList<String>();
		int lineNo = 0, numTokens = 0;
		String[] tokens = null, colNames = null;
		try {
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
				if (hasRowNames) {
					curData = new double[numTokens - 1];
					tok = tokens[0].trim();
					rowNames.add(tok);
					for (int i = 1; i < curNumTokens; i++) {
						tok = tokens[i].trim();
						curData[i - 1] = "".equals(tok) ? kUndefinedValue : Double.parseDouble(tok);
					}
				} else {
					curData = new double[numTokens];
					for (int i = 0; i < curNumTokens; i++) {
						tok = tokens[i].trim();
						curData[i] = "".equals(tok) ? kUndefinedValue : Double.parseDouble(tok);
					}
				}
				list.add(curData);
			}
		} finally {
			if (reader != null) reader.close();
		}
		int size = list.size();
		if (size == 0) return null;
		return new Table(list.toArray(new double[size][]),
			hasRowNames ? rowNames.toArray(new String[rowNames.size()]) : null,
			hasColNames ? colNames : null);
	}
}

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
class QR {
	private static final double kDefaultQRTolerance = 1e-6;
	private static final int kInvalidColumnIndex = -1;
	private static final String sMatrixRowMustAgree = "Matrix row dimensions must agree."; // $NON-NLS-1$

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

	private QR() {}

	/**
	 * QR Decomposition, computed by Householder reflections, with default tolerance (1e-6).
	 * @param A    Rectangular array
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QR (double[][] A) {
		this(A, kDefaultQRTolerance, false, false);
	}

	public QR (double[][] A, boolean transpose) {
		this(A, kDefaultQRTolerance, transpose, false);
	}

	/**
	 * QR Decomposition, computed by Householder reflections.
	 * @param A    Rectangular array
	 * @param tolerance Small integer below which a number is considered zero
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QR (double[][] A, double tolerance) {
		this(A, tolerance, false, false);
	}

	/**
	 * QR Decomposition, computed by Householder reflections.
	 * @param A    Rectangular array
	 * @param tolerance Small integer below which a number is considered zero
	 * @param transpose If true, A is in transpose position
	 * @return     Structure to access R and the Householder vectors and compute Q.
	 */
	public QR (double[][] A, double tolerance, boolean transpose) {
		this(A, tolerance, transpose, false);
	}

	public QR (double[][] A, double tolerance, boolean transpose, boolean addMean) {
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
				nrm = hypot(nrm,qr[rowNo][colNo]);
			double absNrm = abs(nrm);
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
	 * Partial Householder transformation.
	 * Added by Roby Joehanes
	 * @param A numRows x numAddlCols
	 * @return
	 */
	public QR addMoreColumns(double[][] A, double tolerance)
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
			if (abs(newRDiag[colNo]) > tolerance)
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

		QR qr = new QR();
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
	public QR addMoreColumnsTransposed(double[][] A)
	{	return addMoreColumnsTransposed(A, kDefaultQRTolerance); }

	/**
	 * Same as <tt>addMoreColumns</tt>, but A is transposed
	 * @param A numAddlCols x numRows
	 * @param tolerance
	 * @return
	 */
	public QR addMoreColumnsTransposed(double[][] A, double tolerance)
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
			if (abs(newRDiag[colNo]) > tolerance)
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

		QR qr = new QR();
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
	public double[][] getH() {
		double[][] H = new double[mNumRows][mNumColumns];
		for (int rowNo = 0; rowNo < mNumRows; rowNo++)
			System.arraycopy(mQR[rowNo], rowNo, H[rowNo], rowNo, mNumColumns - rowNo);
		return H;
	}

	public double[][] getR() {
		double[][] R = new double[mRank][mRank];
		for (int columnNo = 0; columnNo < mRank; columnNo++) {
			R[columnNo][columnNo] = mRDiag[columnNo];
			int idx = columnNo + 1;
			System.arraycopy(mQR[columnNo], idx, R[columnNo], idx, mRank - idx);
		}
		return R;
	}

	public double[][] getQ() {
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
	public double[][] inverse()
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
		return X;
	}
}

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
class SVD {
	static final double eps = pow(2.0,-52.0);
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

	public SVD (double[][] arg)
	{	this (arg, true, true); }

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
	public SVD (double[][] arg, boolean wantu, boolean wantv) {
		computeSVD(arg, wantu, wantv);
	}

	private void computeSVD(double[][] arg, boolean wantu, boolean wantv)
	{
		// Derived from LINPACK code.
		// Initialize.
		m = arg.length;
		n = arg[0].length;
		int nu = min(m,n);
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
			nct = min(m-1,n),
			nrt = max(0,min(n-2,m)),
			mcrtnrt = max(nct,nrt);
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

				if (abs(e[k]) <= eps*(abs(s[k]) + abs(s[k+1])))
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

					double t = (ks != p ? abs(e[ks]) : 0.) + 
								(ks != k+1 ? abs(e[ks-1]) : 0.);
					if (abs(s[ks]) <= eps*t)
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
						scale = max(max(max(max(
							abs(s[p-1]),abs(s[p-2])),abs(e[p-2])), 
							abs(s[k])),abs(e[k])),
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
						shift = sqrt(b*b + c);
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
	public double[][] getU()
	{	return U; }

	/**
	 * Return the right singular vectors
	 * @return     V
	 */
	public double[][] getV()
	{	return V; }

	private static double[][] copyNColumns(double[][] mat, int ndf) {
		int n = mat.length;
		if (ndf > mat[0].length) throw new RuntimeException();
		double[][] newMat = new double[n][ndf];
		for (int i = 0; i < n; i++)
			System.arraycopy(mat[i], 0, newMat, 0, ndf);
		return newMat;
	}

	/**
	 * Return the left singular vectors
	 * @return     U
	 */
	public double[][] getU(int ndf)
	{	return copyNColumns(U, ndf); }

	/**
	 * Return the right singular vectors
	 * @return     V
	 */
	public double[][] getV(int ndf)
	{	return copyNColumns(V, ndf); }

	/**
	 * Return the one-dimensional array of singular values
	 * @return     diagonal of S.
	 */
	public double[] getSingularValues ()
	{	return s; }

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
	{	return s[0]/s[min(m,n)-1]; }

	/**
	 * Effective numerical matrix rank
	 * @return     Number of nonnegligible singular values.
	 */
	public int rank () {
		double tol = max(m,n)*s[0]*eps;
		// Exploit the fact that the eigenvalues are sorted in descending order -- RJ
		int r = s.length;
		for (int i = 0; i < r; i++)
			if (s[i] <= tol)
				return i;
		return 0;
	}

	public void purge()
	{	U = V = null; s = null; }
}
