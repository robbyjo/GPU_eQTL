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

/**
 * <P>Enumeration that describes modes you can do in multiplication thread (QMultiplicationJob).
 * BLAS-esque naming style is followed here. The multiplication modes are often used in
 * linear models.
 * 
 * <P>Modes explanation / requirements:
 * <ul>
 * <li>XY = XY</li>
 * <li>XtY = X'Y</li>
 * <li>XYt = XY'</li>
 * <li>XXt = XX'</li>
 * <li>XtX = X'X</li>
 * <li>cXXt = cXX'</li>
 * <li>cXtX = cX'X</li>
 * <li>XYXt = XYX' (Y must be symmetric)</li>
 * <li>XtYX = X'YX (Y must be symmetric)</li>
 * <li>cXYXt = cXYX' (Y must be symmetric)</li>
 * <li>cXtYX = cX'YX (Y must be symmetric)</li>
 * <li>XDiagYXt = XDiag(Y)X' (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>XtDiagYX = X'Diag(Y)X (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>cXDiagYXt = cXDiag(Y)X' (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>cXtDiagYX = cX'Diag(Y)X (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>XYPlusZ = XY + Z</li>
 * <li>XtYPlusZ = X'Y + Z</li>
 * <li>XYtPlusZ = XY' + Z</li>
 * <li>cXYPlusZ = cXY + Z</li>
 * <li>cXtYPlusZ = cX'Y + Z</li>
 * <li>cXYtPlusZ = cXY' + Z</li>
 * <li>XYMinusZ = XY - Z</li>
 * <li>XtYMinusZ = X'Y - Z</li>
 * <li>XYtMinusZ = XY' - Z</li>
 * <li>cXYMinusZ = cXY - Z</li>
 * <li>cXtYMinusZ = cX'Y - Z</li>
 * <li>cXYtMinusZ = cXY' - Z</li>
 * <li>ZMinusXY = Z - XY</li>
 * <li>ZMinusXtY = Z - X'Y</li>
 * <li>ZMinusXYt = Z - XY'</li>
 * <li>ZMinuscXY = Z - cXY</li>
 * <li>ZMinuscXtY = Z - cX'Y</li>
 * <li>ZMinuscXYt = Z - cXY'</li>
 * <li>IMinusXY = I - XY</li>
 * <li>IMinusXtY = I - X'Y</li>
 * <li>IMinusXYt = I - XY'</li>
 * <li>IMinuscXY = I - cXY</li>
 * <li>IMinuscXtY = I - cX'Y</li>
 * <li>IMinuscXYt = I - cXY'</li>
 * <li>IMinusXYXt = I - XYX' (Y must be symmetric)</li>
 * <li>IMinusXtYX = I - X'YX (Y must be symmetric)</li>
 * <li>IMinuscXYXt = I - cXYX' (Y must be symmetric)</li>
 * <li>IMinuscXtYX = I - cX'YX (Y must be symmetric)</li>
 * <li>IMinusXXt = I - XX'</li>
 * <li>IMinusXtX = I - X'X</li>
 * <li>IMinuscXXt = I - cXX'</li>
 * <li>IMinuscXtX = I - cX'X</li>
 * <li>IMinusXDiagYXt = I - XDiag(Y)X' (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>IMinusXtDiagYX = I - X'Diag(Y)X (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>IMinuscXDiagYXt = I - cXDiag(Y)X' (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>IMinuscXtDiagYX = I - cX'Diag(Y)X (Y must have only one row. Y[0] contains the diagonal elements)</li>
 * <li>XMinusXYYt = X - XYY'</li>
 * <li>XMinusYYtX = X - YY'X</li>
 * </ul>
 * 
 * <P>Note: I = corresponding identity matrix, c = scalar constant, Diag(X) is
 * a diagonal matrix with X's diagonal elements.
 * 
 * @author Roby Joehanes
 *
 */
public enum EMultiplicationMode
{
	XY, XtY, XYt,
	cXY, cXtY, cXYt,
	XXt, XtX, cXXt, cXtX,
	XYXt, XtYX, cXYXt, cXtYX,
	XDiagYXt, XtDiagYX, cXDiagYXt, cXtDiagYX,
	XYPlusZ, XtYPlusZ, XYtPlusZ,
	XYMinusZ, XtYMinusZ, XYtMinusZ,
	ZMinusXY, ZMinusXtY, ZMinusXYt,
	cXYPlusZ, cXtYPlusZ, cXYtPlusZ,
	cXYMinusZ, cXtYMinusZ, cXYtMinusZ,
	ZMinuscXY, ZMinuscXtY, ZMinuscXYt,
	IMinusXY, IMinusXtY, IMinusXYt,
	IMinuscXY, IMinuscXtY, IMinuscXYt,
	IMinusXYXt, IMinusXtYX, IMinuscXYXt, IMinuscXtYX,
	IMinusXXt, IMinusXtX, IMinuscXXt, IMinuscXtX,
	IMinusXDiagYXt, IMinusXtDiagYX, IMinuscXDiagYXt, IMinuscXtDiagYX,
	XMinusXYYt
}