/*
 * Roby Joehanes
 * 
 * Copyright 2009 Roby Joehanes
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
/*
 * Created on Aug 17, 2004
 *
 */
package gov.nih.eqtl.datastructure;

import gov.nih.utils.QSystemUtils;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;

import com.csvreader.CsvReader;

import static gov.nih.utils.QFileUtils.sLineFormat;
import static gov.nih.utils.QStringUtils.cCommonDelimiter;

/**
 * The main data structure class for handling genetic stuffs.
 * Contain map and population object.
 * 
 * @author Roby Joehanes
 */
public class QGeneticSNPData
{
	protected List<QSNPData> mSNPList = new ArrayList<QSNPData>();
	protected QIndividual[] mPedigree;
	protected ESNPDataType mSNPDataType = ESNPDataType.UNKNOWN;

	/**
	 * The number of individuals contained in this population.
	 */
	protected int mNumIndividuals = -1;

	public void addSNP(QSNPData snp) {
		if (mNumIndividuals == -1)
			mNumIndividuals = snp.getNumIndividuals();
		else
			if (mNumIndividuals != snp.getNumIndividuals())
				throw new RuntimeException("Number of individuals mismatch between locus at locus " + snp.mID);
		mSNPList.add(snp);
	}

	public List<QSNPData> getSNPs() {
		return mSNPList;
	}

	public void clearSNPs() {
		mSNPList.clear();
	}

	public void setPedigree(QIndividual[] ped) {
		mPedigree = ped;
	}

	public QIndividual[] getPedigree() {
		return mPedigree;
	}

	public int getNumSNPs() {
		return mSNPList.size();
	}

	public int getNumIndividuals() {
		return mNumIndividuals;
	}

	public void setSNPDataType(ESNPDataType t) {
		mSNPDataType = t;
	}

	public ESNPDataType getSNPDataType() {
		return mSNPDataType;
	}

	public static final QGeneticSNPData load(String filename, char[] delimiters, String commentMarkers, boolean hasRowNames, boolean hasColNames) throws IOException
	{	return load(new FileReader(filename), delimiters, commentMarkers, hasRowNames, hasColNames, false); }

	public static final QGeneticSNPData load(String filename, char[] delimiters, String commentMarkers, boolean hasRowNames, boolean hasColNames, boolean isEncoded) throws IOException
	{	return load(new FileReader(filename), delimiters, commentMarkers, hasRowNames, hasColNames, isEncoded); }

	public static final QGeneticSNPData load(Reader r, char[] delim, String comments,
		boolean hasRowNames, boolean hasColNames, boolean isEncoded) throws IOException
	{
		QGeneticSNPData dat = new QGeneticSNPData();
		CsvReader reader = null;
		int lineNo = 0, numTokens = 0;
		String[] tokens = null, colNames = null;
		try
		{
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
				if (numTokens == 0) numTokens = curNumTokens;
				if (numTokens != curNumTokens)
					throw new RuntimeException(String.format(sLineFormat, curNumTokens, numTokens, lineNo));

				String tok;
				QSNPData snp;
				if (isEncoded) {
					int[] curData = null;
					snp = new QSNPDataReal();
					if (hasRowNames)
					{
						curData = new int[numTokens - 1];
						tok = tokens[0].trim();
						//snp.setID(Integer.parseInt(tok.substring(2)));
						snp.setID(tok);
						for (int i = 1; i < curNumTokens; i++)
						{
							tok = tokens[i].trim();
							curData[i - 1] = "".equals(tok) || "-".equals(tok) || ".".equals(tok) ? QSNPData.kMissingCode : Integer.parseInt(tok);
						}
					}
					else
					{
						curData = new int[numTokens];
						for (int i = 0; i < curNumTokens; i++)
						{
							tok = tokens[i].trim();
							curData[i] = "".equals(tok) || "-".equals(tok) || ".".equals(tok) ? QSNPData.kMissingCode : Integer.parseInt(tok);
						}
					}
					((QSNPDataInt) snp).setSNPValues(curData);
				} else {
					double[] curData = null;
					snp = new QSNPDataReal();
					if (hasRowNames)
					{
						curData = new double[numTokens - 1];
						tok = tokens[0].trim();
						//snp.setID(Integer.parseInt(tok.substring(2)));
						snp.setID(tok);
						for (int i = 1; i < curNumTokens; i++)
						{
							tok = tokens[i].trim();
							curData[i - 1] = Double.parseDouble(tok);
						}
					}
					else
					{
						curData = new double[numTokens];
						for (int i = 0; i < curNumTokens; i++)
						{
							tok = tokens[i].trim();
							curData[i] = Double.parseDouble(tok);
						}
					}
					((QSNPDataReal) snp).setSNPValues(curData);
				}
				dat.addSNP(snp);
			}
		}
		finally
		{
			if (reader != null)
				reader.close();
		}
		int size = dat.getSNPs().size();
		if (size == 0)
			return null;
		return dat;
	}

	public static void main(String[] args) {
		long gb1 = QSystemUtils.usedMemoryAfterGC(), gb2, time1, time2;
		try {
			time1 = System.currentTimeMillis();
			QGeneticSNPData data = QGeneticSNPData.load(args[0], cCommonDelimiter, "#", true, true);
			time2 = System.currentTimeMillis();
			assert (data != null);
			gb2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println(data.getNumSNPs() + "x" + data.getNumIndividuals());
			System.out.println("Memory consumed (MB) = " + String.format("%5.2f",(gb2 - gb1) / (1024.0 * 1024)));
			System.out.println("Load time (ms) = " + (time2 - time1));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}

