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
/**
 * 
 */
package gov.nih.utils.matrix;

import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.IJobOwner;
import gov.nih.parallel.QSynchronizedCounter;

/**
 * <P>Parallel matrix multiplication thread, depending on the mode.
 * See EMultiplicationMode for details.
 * 
 * <P>The result matrix MUST be initialized and ALL matrices MUST be in
 * conformable dimensions! NO dimension checks will be performed!
 * 
 * <P>Note: For XYXt and XtYX modes, Y MUST BE SYMMETRIC!<BR>
 * For XDiagYXt, XtDiagYX, cXDiagYXt, cXtDiagYX, and similar variants,
 * matrix Y MUST HAVE ONLY ONE ROW and Y[0] MUST contain all diagonal elements.
 * This encoding is to conserve space.
 * 
 * <P>Note that this multiplication thread is optimized for each respective operation mode.
 * 
 * @see EMultiplicationMode
 * 
 * @author Roby Joehanes
 */
public class QMultiplicationJob implements IGenericParallelJob, Runnable
{
	protected double[][]
		mX,
		mY,
		mZ,
		mResult;
	protected double mC = 1;
	protected QSynchronizedCounter mCounter;
	protected EMultiplicationMode mMode = EMultiplicationMode.XY;

	/**
	 * 
	 * @param X
	 * @param Y
	 * @param result
	 * @param counter
	 * @param lock
	 */
	public QMultiplicationJob(double[][] X, double[][] Y, double[][] Z, double[][] result,
		QSynchronizedCounter counter)
	{
		mX = X; mY = Y; mZ = Z; mResult = result; mCounter = counter;
	}

	/**
	 * @param X
	 * @param Y
	 * @param result
	 * @param counter
	 * @param lock
	 * @param mode
	 */
	public QMultiplicationJob(double[][] X, double[][] Y, double[][] Z, double[][] result,
		QSynchronizedCounter counter, EMultiplicationMode mode)
	{
		this(X, Y, Z, result, counter);
		mMode = mode;
	}

	public QMultiplicationJob(double[][] X, double[][] Y, double[][] Z, double[][] result, double c,
		QSynchronizedCounter counter, EMultiplicationMode mode)
	{
		this(X, Y, Z, result, counter, mode);
		mC = c;
	}

	/* (non-Javadoc)
	 * @see qparallel.IGenericParallelJob#cancel()
	 */
	public void cancel()
	{}

	public IJobOwner getOwner() // Not important
	{	return null; }

	/* (non-Javadoc)
	 * @see qparallel.IGenericParallelJob#setOwner(qparallel.IJobOwner)
	 */
	public void setOwner(IJobOwner owner) // Not important
	{}

	public void run()
	{	execute(); }

	/* (non-Javadoc)
	 * @see qparallel.IGenericParallelJob#execute()
	 */
	public void execute()
	{
		switch (mMode)
		{
			case XY:
				calcXY();
				break;
			case XtY:
				calcXtY();
				break;
			case XYt:
				calcXYt();
				break;
			case cXY:
				calccXY();
				break;
			case cXtY:
				calccXtY();
				break;
			case cXYt:
				calccXYt();
				break;
			case XXt:
				calcXXt();
				break;
			case XtX:
				calcXtX();
				break;
			case cXXt:
				calccXXt();
				break;
			case cXtX:
				calccXtX();
				break;
			case XYXt:
				calcXYXt();
				break;
			case XtYX:
				calcXtYX();
				break;
			case cXYXt:
				calccXYXt();
				break;
			case cXtYX:
				calccXtYX();
				break;
			case XDiagYXt:
				calcXDiagYXt();
				break;
			case XtDiagYX:
				calcXtDiagYX();
				break;
			case cXDiagYXt:
				calccXDiagYXt();
				break;
			case cXtDiagYX:
				calccXtDiagYX();
				break;
			case XYPlusZ:
				calcXYPlusZ();
				break;
			case XtYPlusZ:
				calcXtYPlusZ();
				break;
			case XYtPlusZ:
				calcXYtPlusZ();
				break;
			case XYMinusZ:
				calcXYMinusZ();
				break;
			case XtYMinusZ:
				calcXtYMinusZ();
				break;
			case XYtMinusZ:
				calcXYtMinusZ();
				break;
			case ZMinusXY:
				calcZMinusXY();
				break;
			case ZMinusXtY:
				calcZMinusXtY();
				break;
			case ZMinusXYt:
				calcZMinusXYt();
				break;
			case cXYPlusZ:
				calccXYPlusZ();
				break;
			case cXtYPlusZ:
				calccXtYPlusZ();
				break;
			case cXYtPlusZ:
				calccXYtPlusZ();
				break;
			case cXYMinusZ:
				calccXYMinusZ();
				break;
			case cXtYMinusZ:
				calccXtYMinusZ();
				break;
			case cXYtMinusZ:
				calccXYtMinusZ();
				break;
			case ZMinuscXY:
				calcZMinuscXY();
				break;
			case ZMinuscXtY:
				calcZMinuscXtY();
				break;
			case ZMinuscXYt:
				calcZMinuscXYt();
				break;
			case IMinusXY:
				calcIMinusXY();
				break;
			case IMinusXtY:
				calcIMinusXtY();
				break;
			case IMinusXYt:
				calcIMinusXYt();
				break;
			case IMinuscXY:
				calcIMinuscXY();
				break;
			case IMinuscXtY:
				calcIMinuscXtY();
				break;
			case IMinuscXYt:
				calcIMinuscXYt();
				break;
			case IMinusXYXt:
				calcIMinusXYXt();
				break;
			case IMinusXtYX:
				calcIMinusXtYX();
				break;
			case IMinuscXYXt:
				calcIMinuscXYXt();
				break;
			case IMinuscXtYX:
				calcIMinuscXtYX();
				break;
			case IMinusXXt:
				calcIMinusXXt();
				break;
			case IMinusXtX:
				calcIMinusXtX();
				break;
			case IMinuscXXt:
				calcIMinuscXXt();
				break;
			case IMinuscXtX:
				calcIMinuscXtX();
				break;
			case IMinusXDiagYXt:
				calcIMinusXDiagYXt();
				break;
			case IMinusXtDiagYX:
				calcIMinusXtDiagYX();
				break;
			case IMinuscXDiagYXt:
				calcIMinuscXDiagYXt();
				break;
			case IMinuscXtDiagYX:
				calcIMinuscXtDiagYX();
				break;
			case XMinusXYYt:
				calcXMinusXYYt();
				break;
		}
		mX = mY = mResult = null;
		mCounter = null;
	}

	private void calcXY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXtY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXYt()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXtY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXYt()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXYPlusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = sum + mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXtYPlusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = sum + mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXYtPlusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = sum + mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXYMinusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = sum - mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXtYMinusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = sum - mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXYtMinusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = sum - mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcZMinusXY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = mZ[i][j] - sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcZMinusXtY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = mZ[i][j] - sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcZMinusXYt()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = mZ[i][j] - sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXYPlusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = mC * sum + mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXtYPlusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = mC * sum + mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXYtPlusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = mC * sum + mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXYMinusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = mC * sum - mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXtYMinusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = mC * sum - mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXYtMinusZ()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = mC * sum - mZ[i][j];
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcZMinuscXY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = mZ[i][j] - mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcZMinuscXtY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = mZ[i][j] - mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcZMinuscXYt()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = mZ[i][j] - mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = i == j ? 1 - sum : -sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXtY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = i == j ? 1 - sum : -sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXYt()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = i == j ? 1 - sum : -sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * mY[k][j];
				mResult[i][j] = i == j ? 1 - mC * sum : -mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXtY()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY.length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < cols; j++)
			{
				double sum = 0;
				for (int k = 0; k < numElements; k++)
					sum += mX[k][i] * mY[k][j];
				mResult[i][j] = i == j ? 1 - mC * sum : -mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXYt()
	{
		int
			i = mCounter.next(),
			cols = mResult[0].length,
			numElements = mY[0].length;
		if (i < 0)
			return;
		do
		{
			double Xrow[] = mX[i];
			for (int j = 0; j < cols; j++)
			{
				double
					sum = 0,
					Yrow[] = mY[j];
				for (int k = 0; k < numElements; k++)
					sum += Xrow[k] * Yrow[k];
				mResult[i][j] = i == j ? 1 - mC * sum : -mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += Xi[l] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXtYX()
	{
		int
			i = mCounter.next(),
			cols = mX.length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += mX[l][i] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += Xi[l] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXtYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += mX[l][i] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXDiagYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
				XVi[k] = Xi[k] * Vk[k];
	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXtDiagYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
				XVi[k] = mX[k][i] * Vk[k];
	
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXDiagYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
				XVi[k] = Xi[k] * Vk[k];
	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXtDiagYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
				XVi[k] = mX[k][i] * Vk[k];
	
			for (int j = 0; j <= i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			double[] Xi = mX[i];	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = mX[j];
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += Xi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXtX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j <= i; j++)
			{
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += mX[k][i] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			double[] Xi = mX[i];	
			for (int j = 0; j <= i; j++)
			{
				double[] Xj = mX[j];
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += Xi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calccXtX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j <= i; j++)
			{
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += mX[k][i] * mX[k][j];
				mResult[i][j] = mResult[j][i] = mC * sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += Xi[l] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j < i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * Xi[k];
			mResult[i][i] = 1 - sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXtYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += mX[l][i] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j < i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * mX[k][i];
			mResult[i][i] = 1 - sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += Xi[l] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j < i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * Xi[k];
			mResult[i][i] = 1 - sum * mC;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXtYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
			{
				sum = 0;
				Vk = mY[k];
				for (int l = 0; l < cols; l++)
					sum += mX[l][i] * Vk[l];
				XVi[k] = sum;
			}
	
			for (int j = 0; j < i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * mX[k][i];
			mResult[i][i] = 1 - sum * mC;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXDiagYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
				XVi[k] = Xi[k] * Vk[k];
	
			for (int j = 0; j < i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * Xi[k];
			mResult[i][i] = 1 - sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXtDiagYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
				XVi[k] = mX[k][i] * Vk[k];
	
			for (int j = 0; j < i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * mX[k][i];
			mResult[i][i] = 1 - sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXDiagYXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			double[] Xi = mX[i];
			for (int k = 0; k < cols; k++)
				XVi[k] = Xi[k] * Vk[k];
	
			for (int j = 0; j < i; j++)
			{
				double[] Xj = mX[j];
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * Xi[k];
			mResult[i][i] = 1 - sum * mC;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXtDiagYX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		double
			sum,
			Vk[] = mY[0],
			XVi[] = new double[cols];
		do
		{
			for (int k = 0; k < cols; k++)
				XVi[k] = mX[k][i] * Vk[k];
	
			for (int j = 0; j < i; j++)
			{
				sum = 0;
				for (int k = 0; k < cols; k++)
					sum += XVi[k] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum * mC;
			}
			sum = 0;
			for (int k = 0; k < cols; k++)
				sum += XVi[k] * mX[k][i];
			mResult[i][i] = 1 - sum * mC;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			double[] Xi = mX[i];	
			for (int j = 0; j < i; j++)
			{
				double[] Xj = mX[j];
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += Xi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = sum;
			}
			double sum = 0;
			for (int k = 0; k < cols; k++)
				sum += Xi[k] * Xi[k];
			mResult[i][i] = 1 - sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinusXtX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < i; j++)
			{
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += mX[k][i] * mX[k][j];
				mResult[i][j] = mResult[j][i] = sum;
			}
			double sum = 0;
			for (int k = 0; k < cols; k++)
				sum += mX[k][i] * mX[k][i];
			mResult[i][i] = 1 - sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXXt()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			double[] Xi = mX[i];	
			for (int j = 0; j < i; j++)
			{
				double[] Xj = mX[j];
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += Xi[k] * Xj[k];
				mResult[i][j] = mResult[j][i] = mC * sum;
			}
			double sum = 0;
			for (int k = 0; k < cols; k++)
				sum += Xi[k] * Xi[k];
			mResult[i][i] = 1 - mC * sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcIMinuscXtX()
	{
		int
			i = mCounter.next(),
			cols = mX[0].length;
		if (i < 0)
			return;
		do
		{
			for (int j = 0; j < i; j++)
			{
				double sum = 0;
				for (int k = 0; k < cols; k++)
					sum += mX[k][i] * mX[k][j];
				mResult[i][j] = mResult[j][i] = mC * sum;
			}
			double sum = 0;
			for (int k = 0; k < cols; k++)
				sum += mX[k][i] * mX[k][i];
			mResult[i][i] = 1 - mC * sum;
			i = mCounter.next();
		} while (i >= 0);
	}

	private void calcXMinusXYYt()
	{
		int
			i = mCounter.next(),
			n = mX[0].length,
			p = mY[0].length;
		double sum;
		double[] temp = new double[p];
		if (i < 0)
			return;
		do
		{
			double[] xi = mX[i];
			for (int j = 0; j < p; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
					sum += xi[k] * mY[k][j];
				temp[j] = sum;
			}
			for (int j = 0; j < n; j++)
			{
				double[] yj = mY[j];
				sum = 0;
				for (int k = 0; k < p; k++)
					sum += temp[k] * yj[k];
				mResult[i][j] = xi[j] - sum;
			}
			i = mCounter.next();
		} while (i >= 0);
	}
}
