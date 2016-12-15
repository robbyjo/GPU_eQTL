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

import static java.lang.Math.sqrt;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QStatsUtils {
	public static final double calcStdDev(double[] data) {
		int n = data.length;
		double sum = 0, sumSq = 0;

		for (int i = 0; i < n; i++) {
			double val = data[i];
			sum += val;
			sumSq += val * val;
		}
		return sqrt((sumSq - (sum / n) * sum) / (n - 1));
	}

	public static final double calcStdDevAndStandardize(double[] data) {
		int n = data.length;
		double sum = 0, sumSq = 0;

		for (int i = 0; i < n; i++) {
			double val = data[i];
			sum += val;
			sumSq += val * val;
		}
		sumSq = sqrt((sumSq - (sum / n) * sum) / (n - 1));
		sum = sum / n;

		for (int i = 0; i < n; i++)
			data[i] = (data[i] - sum) / sumSq;
		return sumSq;
	}

	public static final double[] calcStdDevAndStandardizeByColumn(double[] data, int nrow, int ncol) {
		int n = data.length;
		assert (n == nrow * ncol);
		double[] sum = new double[ncol], sumSq = new double[ncol];

		for (int i = 0; i < n; i++) {
			double val = data[i];
			int j = i % nrow;
			sum[j] = sum[j] + val;
			sumSq[j] = sumSq[j] + val * val;
		}
		for (int i = 0; i < ncol; i++) {
			sumSq[i] = sqrt((sumSq[i] - (sum[i] / nrow) * sum[i]) / (nrow - 1));
			sum[i] = sum[i] / nrow;
		}

		for (int i = 0; i < n; i++) {
			int j = i % nrow;
			data[i] = (data[i] - sum[j]) / sumSq[j];
		}
		return sumSq;
	}
}
