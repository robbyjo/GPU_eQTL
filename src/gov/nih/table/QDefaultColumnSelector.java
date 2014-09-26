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
package gov.nih.table;

import java.util.BitSet;

/**
 * Simple row data selector. All columns are selected.
 * Trapping for a certain column value
 * 
 * @author Roby Joehanes
 *
 */
public class QDefaultColumnSelector implements IDataFilter
{
	protected BitSet mSelector;

	public QDefaultColumnSelector(BitSet selector)
	{
		mSelector = selector;
	}

	/* (non-Javadoc)
	 * @see qplugin.IDataFilter#isColumnSelected(int, java.lang.String, java.lang.String[])
	 */
	public boolean isColumnSelected(int colNo, Object colName, String[] colCategory)
	{	return mSelector.get(colNo); }

	/* (non-Javadoc)
	 * @see qplugin.IDataFilter#isRowSelected(int, double[], java.lang.String[])
	 */
	public boolean isRowSelected(int rowNo, double[] data, Object[] colNames)
	{	return true; }

	public boolean isAllColumnsSelected()
	{	return false; }
}
