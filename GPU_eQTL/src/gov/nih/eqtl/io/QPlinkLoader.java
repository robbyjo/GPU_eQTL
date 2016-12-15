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
package gov.nih.eqtl.io;

import gov.nih.eqtl.datastructure.ESNPDataType;
import gov.nih.eqtl.datastructure.QFamily;
import gov.nih.eqtl.datastructure.QGeneticSNPData;
import gov.nih.eqtl.datastructure.QIndividual;
import gov.nih.eqtl.datastructure.QSNPData;
import gov.nih.eqtl.datastructure.QSNPDataInt;
import gov.nih.eqtl.datastructure.QSNPDataReal;
import gov.nih.utils.QFileUtils;
import gov.nih.utils.QSystemUtils;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringReader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.archivers.ArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

import static gov.nih.eqtl.datastructure.ESex.*;
import static gov.nih.utils.QDefaultLexer.*;
import static gov.nih.utils.QStringUtils.*;

/**
 * <P>Plugin for loading / saving Plink files.
 * 
 * 
 * @author Roby Joehanes
 */
public class QPlinkLoader
{
	protected QGeneticSNPData mData;
	protected int mNumIndividuals = 0;
	protected boolean mAllowComments = true;

	// If the stream is from filename, then it should go here
	protected String[] mFilenames = null;

	protected static final String
		sBZ2 = ".bz2", //$NON-NLS-1$
		sGZ = ".gz", //$NON-NLS-1$
		sTPED = ".tped", //$NON-NLS-1$
		sPED = ".ped", //$NON-NLS-1$
		sFAM = ".fam", //$NON-NLS-1$
		sTFAM = ".tfam", //$NON-NLS-1$
		sDelim = " \t\n\r\f,"; //$NON-NLS-1$
	protected static final int
		kNumColumnsInFamFiles = 6,
		kFamIdxFamilyID = 0,
		kFamIdxSubjectID = 1,
		kFamIdxFatherID = 2,
		kFamIdxMotherID = 3,
		kFamIdxSex = 4,
		kFamIdxPhenotype = 5, // We ignore this column
		kNumHeaderColumnsInTpedFiles = 4,
		kTpedIdxChrom = 0,
		kTpedIdxSNPID = 1,
		kTpedIdxGenDistID = 2,
		kTpedIdxLocationID = 3;

	/**
	 * Wrapper for loading routine if the source is just a file.
	 * @param filename The file name to load
	 * @return
	 * @throws Exception
	 */
	public QGeneticSNPData load(String filename) throws Exception
	{
		FileReader reader = null;
		try
		{
			reader = new FileReader(filename);
			QGeneticSNPData data = load(reader);
			return data;
		}
		catch (Exception e)
		{
			throw new IOException(String.format(
					"Cannot read file %s", filename) +
					sLn + e.getMessage());
		}
		finally
		{
			if (reader != null)
				reader.close();
		}
	}

	/**
	 * Perform loading from multiple files
	 * @param filenames arrays of filenames
	 * @return
	 * @throws Exception
	 */
	public QGeneticSNPData load(String... filenames) throws Exception
	{
		List<String> nonNullFilenames = new ArrayList<String>();
		for (String s: filenames)
			if (s != null)
				nonNullFilenames.add(s);
		int numFiles = nonNullFilenames.size();
		mFilenames = nonNullFilenames.toArray(new String[numFiles]);
		Reader[] readers = new Reader[numFiles];
		for (int i = 0; i < numFiles; i++)
			readers[i] = new FileReader(nonNullFilenames.get(i));
		return load(readers);
	}

	public int getNumberFilesNeededToLoad()
	{	return 2; }

	public QGeneticSNPData loadAsOneFile(String... filenames) throws Exception
	{
		// Force the files to load separately.
		return load(filenames);
	}

	public QIndividual[] readFamFile(Reader r) throws IOException
	{
		String txt = QFileUtils.readTextFile(r);
		String[][] tokens = QFileUtils.readCSV(new StringReader(txt));
		tokens = QFileUtils.readDelimitedFile(new StringReader(txt), cCommonDelimiter, null);
		int
			nrows = tokens.length,
			ncols = tokens[0].length,
			startIdx = 0;
		if (ncols != kNumColumnsInFamFiles)
			throw new IOException("Error: Invalid family file");
		try {
			// As per PLINK manual, the phenotype column has to be encoded as -9, 0, and 1.
			// Anything else means that the first row is a header row.
			// See: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
			Integer.parseInt(tokens[0][kFamIdxPhenotype]);
		} catch (Exception e) {
			startIdx++; // We have a header row
		}
		Map<String, QIndividual> individualCache = new HashMap<String, QIndividual>();
		Map<String, QFamily> familyCache = new HashMap<String, QFamily>();
		Set<String> seenSelfID = new HashSet<String>();
		QIndividual[] inds = new QIndividual[nrows - startIdx];
		for (int i = startIdx; i < nrows; i++)
		{
			String[] toks = tokens[i];
			String
				dadID = toks[kFamIdxFatherID].intern(),
				momID = toks[kFamIdxMotherID].intern(),
				famID = toks[kFamIdxFamilyID].intern(),
				selfID = toks[kFamIdxSubjectID].intern();
			QIndividual
				father = individualCache.get(dadID),
				mother = individualCache.get(momID);
			QFamily family = familyCache.get(famID);
			if (seenSelfID.contains(selfID))
				throw new IOException(String.format("Error: Duplicate subject ID %d in family file", selfID));

			if (family == null)
			{
				family = new QFamily(famID);
				familyCache.put(famID, family);
			}
			if (father == null)
			{
				father = new QIndividual(null, dadID, MALE);
				individualCache.put(dadID, father);
			}
			if (mother == null)
			{
				mother = new QIndividual(null, momID, FEMALE);
				individualCache.put(momID, mother);
			}
			int sexID = -1;
			try {
				sexID = Integer.parseInt(toks[kFamIdxSex]);
			} catch (Exception e)
			{}
			QIndividual subject = inds[i - startIdx] = new QIndividual(family, selfID, sexID == 1 ? MALE : sexID == 2 ? FEMALE : UNKNOWN_SEX);
			subject.setParents(new QIndividual[] {father, mother});
			seenSelfID.add(selfID);
		}
		return inds;
	}

	public QGeneticSNPData readTPEDFile(Reader reader) throws IOException
	{
		int
			numAlleles = 0,
			expectedNumTokens = 0;
		ESNPDataType dataType = ESNPDataType.UNKNOWN;
		LineNumberReader fileReader = new LineNumberReader(reader);
		String curLine;
		boolean encounteredHeader = false;
		QGeneticSNPData allData = new QGeneticSNPData();

		//StringBuilder buf = new StringBuilder();
		//buf.append("Chr, SNP ID, Pos, A, C, G, T, AA, AC, AG, AT, CC, CG, CT, GG, GT, TT, NumDiploidAlleles, Major Homozygote, Minor Homozygote, Heterozygote, Major Homozygote Count, Minor Homozygote Count, Heterozygote Count, Num Missing"+qutils.QStringUtils.sLn);
		//long gb1 = QSystemUtils.usedMemoryAfterGC();
		int lineNo = 0;
		do {
			curLine = fileReader.readLine();
			if (curLine == null)
				break;
			lineNo++;
			//*
			curLine = curLine.trim();
			int commentIndex = curLine.indexOf(sComment1);
			int commentIndex2 = curLine.indexOf(sComment2);
			if (commentIndex != -1)
			{
				if (commentIndex2 > -1 && commentIndex2 < commentIndex)
					commentIndex = commentIndex2;
			}
			else
				commentIndex = commentIndex2;
			if (commentIndex > -1) // Chop the comment off (if any)
				curLine = curLine.substring(0, commentIndex).trim();
			if (curLine.length() == 0)
				continue;
			//*/
			//System.out.println(curLine);
			//*
			StringTokenizer tok = new StringTokenizer(curLine, sDelim);
			int numTokens = tok.countTokens();
			if (expectedNumTokens == 0) {
				expectedNumTokens = numTokens;
				numAlleles = numTokens - kNumHeaderColumnsInTpedFiles;
				if ((numAlleles & 1) == 1 && dataType == ESNPDataType.ALPHANUM)
					throw new IOException("Error: Unexpected number of columns!");
			} else if (numTokens != expectedNumTokens)
				throw new IOException("Error: Unexpected number of columns!");
			String
				chrID = tok.nextToken().intern(),
				snpIDStr = tok.nextToken().intern();
			tok.nextToken(); // We discard the distance info in morgan
			String basePosToken = tok.nextToken();
			int basePos = 0;
			try {
				basePos = Integer.parseInt(basePosToken);
				encounteredHeader = true;
			} catch (Exception e) {
				// This is header line
				if (encounteredHeader) {
					throw new IOException(String.format("Error: Unexpected token %s at line %d", basePosToken, lineNo));
				}
				continue;
			}
			String alleles[] = new String[numAlleles];
			for (int i = 0; i < numAlleles; i++)
				alleles[i] = tok.nextToken();
			if (dataType == ESNPDataType.UNKNOWN)
				try {
					Double.parseDouble(alleles[0]);
					dataType = ESNPDataType.NUMERIC;
				} catch (Exception e) {
					dataType = ESNPDataType.ALPHANUM;
				}

			QSNPData snpData = dataType == ESNPDataType.ALPHANUM ? new QSNPDataInt() : new QSNPDataReal();
			snpData.setID(snpIDStr);
			snpData.setSNPValues(alleles);
			snpData.setChromID(Integer.parseInt(chrID));
			snpData.setLocation(basePos);
			/*
			buf.append(chrID + ", " + snpID + ", " + basePos + ", " + qutils.QStringUtils.toString(snpData.getAlleleCounts())
				+ ", " + qutils.QStringUtils.toString(snpData.getDiploidAlleleCounts())
				+ ", " + snpData.getNumDiploidAlleleTypes()
				+ ", " + snpData.getMajorHomozygoteAllele()
				+ ", " + snpData.getMinorHomozygoteAllele()
				+ ", " + snpData.getHeterozygoteAllele()
				+ ", " + snpData.getMajorHomozygoteCount()
				+ ", " + snpData.getMinorHomozygoteCount()
				+ ", " + snpData.getHeterozygoteCount()
				+ ", " + snpData.getNumMissing()
				+ qutils.QStringUtils.sLn);
			//*/
			allData.addSNP(snpData);
		} while (true);
		allData.setSNPDataType(dataType);
		//long gb2 = QSystemUtils.usedMemoryAfterGC();
		//System.out.println("Memory use for all SNPs = " + (gb2 - gb1));
		//QFileUtils.writeText(buf.toString(), "/Users/joehanesr/Desktop/off-snpcount.txt");
		return allData;
	}

	/**
	 * @see qplugin.io.QIOPlugin#load(java.lang.String)
	 */
	public QGeneticSNPData load(Reader... reader) throws IOException
	{
		if (reader.length > getNumberFilesNeededToLoad())
			throw new IOException("Error: Files need to be read in at most two streams!");
		int
			numFilenames = 0,
			famIdx = 0;
		if (mFilenames != null)
		{
			numFilenames = mFilenames.length;
			if (numFilenames != reader.length)
				throw new IOException("Error: Files need to be read in at most two streams!");
			for (int i = 0; i < numFilenames; i++)
			{
				String fn = mFilenames[i].toLowerCase(Locale.ENGLISH);
				if (fn.endsWith(sBZ2))
				{
					reader[i] = new InputStreamReader(new BZip2CompressorInputStream(new BufferedInputStream(new FileInputStream(fn))));
					fn = fn.substring(0, fn.lastIndexOf(sBZ2));
				} else if (fn.endsWith(sGZ))
				{
					reader[i] = new InputStreamReader(new GZIPInputStream(new BufferedInputStream(new FileInputStream(fn))));
					fn = fn.substring(0, fn.lastIndexOf(sGZ));
				}
				if (fn.endsWith(sPED) || fn.endsWith(sTFAM) || fn.endsWith(sFAM))
				{	famIdx = i; }
			}
		}
		if (reader.length == 1)
		{
			return readTPEDFile(reader[0]);
		} else
		{
			QIndividual inds[] = readFamFile(reader[famIdx]);
			QGeneticSNPData snpData = readTPEDFile(reader[famIdx == 0 ? 1 : 0]);
			snpData.setPedigree(inds);
			return snpData;
		}
	}

	/**
	 * Dump the data into string
	 * @param pluginData
	 * @return
	 */
	public final String toString(QGeneticSNPData pluginData)
	{
		StringBuilder buffer = new StringBuilder();
		// TODO
		return buffer.toString();
	}

	/**
	 * @see qplugin.io.QIOPlugin#save(qplugin.IPluginData, java.lang.String)
	 */
	public void save(QGeneticSNPData pluginData, Writer w) throws IOException
	{
		try
		{
			PrintWriter writer = new PrintWriter(w);
			writer.println(toString(null));
			writer.flush();
		}
		catch (Exception e)
		{
			// Should never reach here
			//e.printStackTrace();
			throw new IOException("Error in writing to files!");
		}
	}

	public static final void saveEncodedFile(QGeneticSNPData popn, Writer pw) throws IOException
	{
		pw.write("Chr SNP_ID Morgan BP");
		for (QIndividual ind: popn.getPedigree())
			pw.write(sSp + ind.getID());
		pw.write(sLn);
		for (QSNPData snp : popn.getSNPs())//popn.selectSNPs(0.05, 0.05))
		{
			pw.write(snp.getChromID() + sSp);
			pw.write("ss"+snp.getID());
			pw.write(" 0 ");
			pw.write(snp.getLocation() + sSp);
			for (int i : snp.getSNPCodes())
				pw.write(i < 0 ? "- " : i+sSp);
			pw.write(sLn);
		}
		pw.flush();
		pw.close();
	}

	public QGeneticSNPData loadAsTarBZip2(String filename) throws Exception
	{
		// TODO This routine will blow up when reading huge files.
		TarArchiveInputStream instr = new TarArchiveInputStream(new BZip2CompressorInputStream(new BufferedInputStream(new FileInputStream(filename))));
		ArchiveEntry entry = null;
		List<Reader> readerList = new ArrayList<Reader>();
		List<String> fileNameList = new ArrayList<String>();
		while ((entry = instr.getNextEntry()) != null) {
			long fileSize = entry.getSize();
			String fileName = entry.getName();
			if (entry.isDirectory())
			{
				// A directory entry shouldn't have a very big entry size. Should be less than a few hundred bytes.
				if (fileSize > Integer.MAX_VALUE)
					throw new RuntimeException();
				instr.read(new byte[(int) fileSize]);
				continue;
			}
			fileNameList.add(fileName);
			Reader reader = new StringReader(QFileUtils.readTextFile(instr));
			readerList.add(reader);
		}
		mFilenames = fileNameList.toArray(new String[fileNameList.size()]);
		return load(readerList.toArray(new Reader[readerList.size()]));
	}

	public QGeneticSNPData loadAsTarGZip(String filename) throws Exception
	{
		// TODO This routine will blow up when reading huge files.
		TarArchiveInputStream instr = new TarArchiveInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(filename))));
		ArchiveEntry entry = null;
		List<Reader> readerList = new ArrayList<Reader>();
		List<String> fileNameList = new ArrayList<String>();
		while ((entry = instr.getNextEntry()) != null) {
			long fileSize = entry.getSize();
			String fileName = entry.getName();
			if (entry.isDirectory())
			{
				// A directory entry shouldn't have a very big entry size. Should be less than a few hundred bytes.
				if (fileSize > Integer.MAX_VALUE)
					throw new RuntimeException();
				instr.read(new byte[(int) fileSize]);
				continue;
			}
			fileNameList.add(fileName);
			Reader reader = new StringReader(QFileUtils.readTextFile(instr));
			readerList.add(reader);
		}
		mFilenames = fileNameList.toArray(new String[fileNameList.size()]);
		return load(readerList.toArray(new Reader[readerList.size()]));
	}

	public static void main(String[] args)
	{
		QPlinkLoader plugin = new QPlinkLoader();
		try
		{
			long time1 = System.currentTimeMillis();
			long gb1 = QSystemUtils.usedMemoryAfterGC();
			//System.out.println(new File(".").getAbsolutePath());
			QGeneticSNPData data = args.length == 1 ? plugin.loadAsTarBZip2(args[0]) : plugin.load(args[0], args[1]);
			long time2 = System.currentTimeMillis();
			long gb2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Total load time (in seconds) = " + (time2 - time1) / 1000.0);
			System.out.println("Memory use for loading (in GB) = " + (gb2 - gb1) * 1e-9);
			System.out.println("Total memory use (in GB) = " + gb2 * 1e-9);
			time1 = System.currentTimeMillis();
			java.io.PrintWriter pw = new java.io.PrintWriter(args[2]);
			saveEncodedFile(data, pw);
			time2 = System.currentTimeMillis();
			System.out.println("Total save time (in seconds) = " + (time2 - time1) / 1000.0);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
}
