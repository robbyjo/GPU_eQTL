package gov.nih.exon.reader;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import qtables.QTableDataExtra;
import qutils.QStringUtils;

import static qmath.QMathConstants.kUndefinedValue;
import static qutils.QFileUtils.readDelimitedFileAsTableDataExtra;
import static qutils.QStringUtils.*;

/**
 *
 * @author Roby Joehanes
 *
 */
public class QExonDataReader {
	private static final int
		kNumHeaderColumns = 2,
		kProbeSetIDColumnNo = 0,
		kGeneIDColumnNo = 1;

	/**
	 * Read exon data file.
	 * Most of the code is copied from {@link QFileUtils#readDelimitedFile(java.io.Reader, String)}.
	 * <br>
	 * NOTE: The first column is assumed to be ProbeSetID and the second is GeneID.
	 * The remaining columns are assumed to be of numerical values. Non-numerical values in
	 * these columns are considered missing (kUndefinedValue).
	 * 
	 * @param filename
	 * @param delim
	 * @param hasHeader if true, the first column is treated as header (i.e. ignored).
	 * @return
	 * @throws IOException
	 */
	public static final Map<String, double[][]> readExonDataFile(String filename, String delim, boolean hasHeader) throws IOException {
		LineNumberReader reader = null;
		HashMap<String, double[][]> map = new HashMap<String, double[][]>();
		Set<String> probeSet = new HashSet<String>();
		int lineNo = 0, numTokens = 0;
		try {
			reader = new LineNumberReader(new BufferedReader(new FileReader(filename), 512 * 1024 * 1024));
			String s = null;
			do {
				s = reader.readLine();
				if (s == null)
					break;
				s = s.trim();
				if (s.length() == 0 || s.startsWith(sPoundSign)) // Empty string or starts with a comment
					continue;
				lineNo++;
				if (lineNo == 1 && hasHeader)
					continue;
				// StringTokenizer doesn't work for white spaces. It gobbles
				// multiple delimiters at once.
				String[] tokens = s.split(delim);
				int curNumTokens = tokens.length;
				if (numTokens == 0)
					numTokens = curNumTokens - kNumHeaderColumns;
				else
					if (curNumTokens - kNumHeaderColumns != numTokens)
						throw new RuntimeException(String.format("%d vs %d, at line %d", curNumTokens, numTokens, lineNo));
				for (int i = 0; i < curNumTokens; i++)
				{
					String tok = tokens[i];
					if (tok.startsWith(sDoubleQuote) && tok.endsWith(sDoubleQuote))
						tokens[i] = tok.substring(1, tok.length() - 1);
				}
				String geneName = tokens[kGeneIDColumnNo];
				probeSet.add(tokens[kProbeSetIDColumnNo]);
				double[] nums = new double[numTokens];
				for (int i = kNumHeaderColumns; i < curNumTokens; i++)
					try {
						nums[i - kNumHeaderColumns] = Double.parseDouble(tokens[i].trim());
					} catch (NumberFormatException e) {
						nums[i - kNumHeaderColumns] = kUndefinedValue;
					}
				double[][] list = map.get(geneName);
				if (list == null) {
					map.put(geneName, list = new double[1][]);
					list[0] = nums;
				}
				else {
					int numList = list.length;
					double[][] newList = new double[numList + 1][];
					System.arraycopy(list, 0, newList, 0, numList);
					newList[numList] = nums;
					map.put(geneName, newList);
				}
			} while (true);
		}
		finally {
			if (reader != null)
				reader.close();
		}
		int size = map.size();
		System.out.println(String.format("%d lines, %d probesets, %d genes.", lineNo, probeSet.size(), size));
		probeSet.clear(); probeSet = null;
		return size == 0 ? null : map;
	}

	public static final QExonData readExonDataFileUngrouped(String filename, char[] delim, boolean hasHeader) throws IOException {
		QTableDataExtra table = readDelimitedFileAsTableDataExtra(filename, delim, "#", 2, hasHeader ? 1 : 0, true);
		String[][] rowNames = table.getAllRowNames();
		QExonData data = new QExonData(QStringUtils.toIntArray(rowNames[0]), rowNames[1], table.getData());
		return data;
	}
}
