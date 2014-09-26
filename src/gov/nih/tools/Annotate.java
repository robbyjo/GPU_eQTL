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
/**
 * 
 */
package gov.nih.tools;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import com.csvreader.CsvReader;

//; Comment
//
//[GeneAnnot]
//file = ../HuEx-1_0-st-v2.na33.1.hg19.probeset-core.txt
//file_format = csv
//columns_kept = probeset_id, transcript_cluster_id, seqname, strand, start, stop, GeneSymbol, EntrezGeneID
//columns_rename = ProbesetID, Transcript_Cluster_ID, Transcript_Chr, Transcript_Strand, Transcript_Start_hg19, Transcript_Stop_hg19, Transcript_GeneSymbol, Transcript_EntrezGeneID
//Transcript_Chr = substring(Transcript_Chr, 3)
//ProbesetID = ProbesetID == "4000000" ? "4e+06" : ProbesetID
//
//[SNPAnnot]
//file = ../1000G-annot.txt
//file_format = csv
//columns_kept = marker, rs, chr, pos, Al1, Al2, Freq1, MAF, Rsq
//columns_rename = Marker, Rs_ID, SNP_Chr, SNP_Pos_hg19, SNP_Fx_Allele, SNP_Non_Fx_Allele, SNP_Fx_Allele_Freq, SNP_MAF, SNP_Imputation_RSq
//
//[Main]
//input = eqtl-1000g-exon-peer-adj7.txt
//output = eqtl-1000g-exon-peer-adj7-annot.txt
//file_format = csv
//matching_rule = GeneAnnot.ProbesetID == Main.ProbesetID, SNPAnnot.Marker == Main.Rs_ID
//Is_Cis = SNP_Chr == Transcript_Chr && abs(Transcript_Start_hg19 - SNP_Pos_hg19) <= 1000000 ? "1" : "0"
//Main.ProbesetID = String.valueOf(round(Main.ProbesetID))
//column_order = Marker, SNPAnnot.Rs_ID, SNP_Chr, SNP_Pos_hg19, SNP_Fx_Allele, SNP_Non_Fx_Allele, SNP_Fx_Allele_Freq, SNP_MAF, SNP_Imputation_RSq, Main.ProbesetID, Transcript_Cluster_ID, Transcript_Chr, Transcript_Strand, Transcript_Start_hg19, Transcript_Stop_hg19, Transcript_GeneSymbol, Transcript_EntrezGeneID, Is_Cis, RSq, Fx, T, log10P
//new_header = Marker, Rs_ID, SNP_Chr, SNP_Pos_hg19, SNP_Fx_Allele, SNP_Non_Fx_Allele, SNP_Fx_Allele_Freq, SNP_MAF, SNP_Imputation_RSq, ProbesetID, Transcript_Cluster_ID, Transcript_Chr, Transcript_Strand, Transcript_Start_hg19, Transcript_Stop_hg19, Transcript_GeneSymbol, Transcript_EntrezGeneID, Is_Cis, RSq, Fx, T, log10P

/**
 * A powerful annotator
 * 
 * @author Roby Joehanes
 */
public class Annotate {
	static File curDir;
	static final String
		ln = System.getProperty("line.separator"), //$NON-NLS-1$
		tab = "    ", //$NON-NLS-1$
		sMain = "Main", //$NON-NLS-1$
		sInput = "input", //$NON-NLS-1$
		sOutput = "output", //$NON-NLS-1$
		sFile = "file", //$NON-NLS-1$
		sFileFormat = "file_format", //$NON-NLS-1$
		sKeyColumn = ":key_column", //$NON-NLS-1$
		sColumnsKept = "columns_kept", //$NON-NLS-1$
		sColumnsRename = "columns_rename", //$NON-NLS-1$
		sMatchingRule = "matching_rule", //$NON-NLS-1$
		sColumnOrder = "column_order", //$NON-NLS-1$
		sNewHeader = "new_header", //$NON-NLS-1$
		sValidIDRegex = "^[A-Za-z_]+[0-9A-Za-z_\\.]*$", //$NON-NLS-1$
		sValidIDWithDotRegex = "^[A-Za-z_]+[0-9A-Za-z_]*\\..+$", //$NON-NLS-1$
		sCommaSeparator = "\\s*,\\s*", //$NON-NLS-1$
		sComments[] = { ";", "#", "//"}, //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
		keywords[] = { sFile, sFileFormat, sInput, sOutput, sKeyColumn, sColumnsKept, sColumnsRename,
		sMatchingRule, sColumnOrder, sNewHeader, sMain},
		sTransformMethodName = "transform"; //$NON-NLS-1$
	static Set<String> defaultSectionKeys = new HashSet<String>(),
		defaultMainSectionKeys = new HashSet<String>();
	static int classCounter = 0;
	static {
		defaultSectionKeys.add(sFile);
		defaultSectionKeys.add(sFileFormat);
		defaultSectionKeys.add(sColumnsKept);
		defaultSectionKeys.add(sColumnsRename);
		defaultSectionKeys.add(sKeyColumn);
		defaultMainSectionKeys.add(sInput);
		defaultMainSectionKeys.add(sOutput);
		defaultMainSectionKeys.add(sFileFormat);
		defaultMainSectionKeys.add(sMatchingRule);
		defaultMainSectionKeys.add(sColumnOrder);
		defaultMainSectionKeys.add(sNewHeader);
	}

	private static final void createStandardMethods(StringBuilder buf) {
		buf.append(tab + "public static final <T> Set<T> makeSet(T[] tokens) {" + ln);
		buf.append(tab + tab + "Set<T> set = new HashSet<T>();" + ln);
		buf.append(tab + tab + "for(T t: tokens) set.add(t);" + ln);
		buf.append(tab + tab + "return set;" + ln);
		buf.append(tab + "}" + ln);

		buf.append(tab + "public static final <T> boolean nonEmptyIntersect(Set<T> s1, Set<T> s2) {" + ln);
		buf.append(tab + tab + "Set<T> set = new HashSet<T>(); set.addAll(s1);" + ln);
		buf.append(tab + tab + "set.retainAll(s2);" + ln);
		buf.append(tab + tab + "return set.size() > 0;" + ln);
		buf.append(tab + "}" + ln);
	}

	public static final Annotation loadAnnotationFile(String filename, char delim,
		String keyColumn, String[] keptColumns, String[] renamedColumns,
		Map<String, String> customRules) throws IOException {
		System.out.println("Reading " + filename);
		if (keyColumn == null || keyColumn.length() == 0)
			throw new RuntimeException("Key column cannot be empty");
		if (renamedColumns != null && keptColumns != null && renamedColumns.length != keptColumns.length)
			throw new RuntimeException("Number of columns in the columns_rename must match with the number of columns in columns_kept");
		if (keptColumns != null && !keyColumn.equals(keptColumns[0]) && renamedColumns != null && !keyColumn.equals(renamedColumns[0])) {
			List<String> strList = new ArrayList<String>();
			strList.add(keyColumn);
			for (String col: keptColumns)
				if (!strList.contains(col))
					strList.add(col);
			keptColumns = strList.toArray(new String[strList.size()]);
		}
		CsvReader reader = null;
		int numLines = countNumLines(filename);
		HashMap<String, String[]> table = new HashMap<String, String[]>(numLines);
		int lineNo = 0, numColumnsKept = keptColumns == null ? 0 : keptColumns.length;
		int[] idx = keptColumns == null ? null : new int[numColumnsKept];
		reader = new CsvReader(filename, new char[] { delim });
		reader.setTrimWhitespace(true);
		reader.setUseTextQualifier(true);
		Method method = null;
		while(reader.readRecord()) {
			String[] tokens = reader.getValues();
			lineNo++;
			if (lineNo == 1) {
				Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
				for (int i = 0; i < tokens.length; i++)
					colnamesToIdx.put(tokens[i], i);
				if (numColumnsKept > 0) {
					for (int i = 0; i < numColumnsKept; i++)
						idx[i] = colnamesToIdx.get(keptColumns[i]);
				} else {
					idx = new int[tokens.length];
					keptColumns = new String[tokens.length];
					idx[0] = colnamesToIdx.get(keyColumn);
					keptColumns[0] = tokens[idx[0]];
					for (int i = 0, j = 1; i < tokens.length; i++) {
						int cur_idx = colnamesToIdx.get(tokens[i]);
						if (cur_idx == idx[0]) continue;
						idx[j] = cur_idx;
						keptColumns[j] = tokens[cur_idx];
						j++;
					}
				}
				if (customRules != null && customRules.size() > 0) {
					colnamesToIdx.clear();
					for (int i = 0; i < idx.length; i++) {
						colnamesToIdx.put(keptColumns[i], i);
						colnamesToIdx.put(renamedColumns[i], i);
					}
					StringBuilder buf = new StringBuilder();
					String className = "Acceptor" + (++classCounter);
					buf.append("import java.util.*;" + ln);
					buf.append("public class " + className + " {" + ln);
					createStandardMethods(buf);
					buf.append(tab + "public static final void "+sTransformMethodName+" (String[] tokens) {" + ln);
					for (String key: customRules.keySet()) {
						String value = customRules.get(key);
						TokenManager tokmgr = TokenManager.tokenize(value);
						ASTNode node = ASTBuilder.parseExpr(tokmgr);
						if (node == null)
							throw new RuntimeException("Parse error!");
						ASTBuilder.inferTypes(node);
						List<ASTNode> termNodes = ASTBuilder.getTerminalNodes(node);
						IdentityHashMap<ASTNode, String> symbolTable = new IdentityHashMap<ASTNode, String>();
						for (ASTNode curNode: termNodes)
							if (curNode.getNodeType() == ASTNodeType.IDENTIFIER) {
								String tok = curNode.getToken();
								Integer iidx = colnamesToIdx.get(tok);
								if (iidx == null)
									throw new RuntimeException("Column " + tok + " is not found!");
								symbolTable.put(curNode, iidx.toString());
							}
						Integer assignee = colnamesToIdx.get(key);
						if (assignee == null)
							throw new RuntimeException("Column " + key + " is not found!");
						JavaProgram prog = new JavaProgram();
						CodeGenerator.generateMethod(node, prog, symbolTable);
						buf.append(tab + tab + "tokens[" + assignee + "] = " + prog.method + ";" + ln);
					}
					buf.append(tab + "}" + ln);
					buf.append("}" + ln);
					System.out.println(buf.toString());
					MemoryClassLoader mcl = new MemoryClassLoader(className, buf.toString());
					try {
						method = mcl.loadClass(className).getMethod(sTransformMethodName, String[].class);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
				continue; // Skip first line
			}
			String[] keptTokens = new String[idx.length];
			for (int i = 0; i < idx.length; i++)
				keptTokens[i] = tokens[idx[i]].trim();
			if (method != null) {
				try {
					method.invoke(null, (Object) keptTokens);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			table.put(keptTokens[idx[0]], keptTokens);
		}
		System.out.println(lineNo + " lines were read.");
		reader.close(); reader = null;
		return new Annotation(null, renamedColumns == null ? keptColumns : renamedColumns, table);
	}

	private static final int countNumLines(String filename) throws IOException {
	    InputStream is = new BufferedInputStream(new FileInputStream(filename));
	    try {
	        byte[] c = new byte[1024];
	        int count = 0;
	        int readChars = 0;
	        boolean empty = true;
	        while ((readChars = is.read(c)) != -1) {
	            empty = false;
	            for (int i = 0; i < readChars; ++i) {
	                if (c[i] == '\n') {
	                    ++count;
	                }
	            }
	        }
	        return (count == 0 && !empty) ? 1 : count;
	    } finally {
	        is.close();
	    }
	}

	private static final Map<String, Map<String, String>> readIniFile(String filename, boolean caseSensitiveKey) throws IOException {
		String[] lines = readTextFileAsArray(new BufferedReader(new FileReader(filename)), sComments);
		if (lines == null || lines.length == 0)
			return null;
		Map<String, Map<String, String>> table = new HashMap<String, Map<String,String>>();
		String curSectionName = "";
		for (String curLine: lines) {
			if (curLine.startsWith("[") && curLine.endsWith("]")) {
				curSectionName = curLine.substring(1, curLine.length() - 1);
				continue;
			}
			int eqPos = curLine.indexOf('=');
			if (eqPos < 0)
				throw new RuntimeException("There is no equal sign in the ini file!");
			String
				key = curLine.substring(0, eqPos).trim(),
				value = curLine.substring(eqPos+1).trim();
			if (key.length() == 0)
				throw new RuntimeException("The key is empty!");
			if (value.length() == 0)
				throw new RuntimeException("The value is empty!");
			if (!key.matches(sValidIDRegex))
				throw new RuntimeException("Malformed key: " + key);
			Map<String, String> curSection = table.get(curSectionName);
			if (curSection == null) {
				curSection = new HashMap<String, String>();
				table.put(curSectionName, curSection);
			}
			curSection.put(caseSensitiveKey ? key : key.toLowerCase(Locale.ENGLISH), value);
		}
		return table;
	}

	private static final String[] readTextFileAsArray(Reader r, String[] comments) throws IOException
	{
		LineNumberReader reader = null;
		List<String> list = new ArrayList<String>();
		try
		{
			reader = new LineNumberReader(r);
			String s = null;
			do {
				s = reader.readLine();
				if (s == null)
					break;
				if (comments != null)
				{
					String temp = s.trim();
					boolean isComments = false;
					if (temp.length() > 0)
					{
						for (int comNo = 0; comNo < comments.length; comNo++)
							if (temp.startsWith(comments[comNo]))
							{
								isComments = true;
								break;
							}
					}
					else
						isComments = true;
					if (isComments)
						continue;
				}
				list.add(s);
			} while (true);
		}
		finally
		{
			if (reader != null)
				reader.close();
		}
		return list.toArray(new String[list.size()]);
	}

	/**
	 * Unescape the string str. The string str is expected to have Java's escape string
	 * Adapted from Apache's org.apache.commons.lang.StringEscapeUtils
	 * RJ 01/02/2006
	 * 
	 * @param str
	 * @return
	 */
	private static final String unescape(String str)
	{
		if (str == null)
			return null;

		int strLength = str.length();
		StringBuffer
			unicode = new StringBuffer(4),
			output = new StringBuffer(strLength);
		boolean
			hadSlash = false,
			inUnicode = false;
		for (int i = 0; i < strLength; i++)
		{
			char ch = str.charAt(i);

			// If we previously have a /u instead, then we need to get the next four
			// digits to translate to unicode string.
			if (inUnicode)
			{
				// if in unicode, then we're reading unicode values in somehow
				unicode.append(ch);
				if (unicode.length() == 4)
				{
					// unicode now contains the four hex digits
					// which represents our unicode character
					try
					{
						int value = Integer.parseInt(unicode.toString(), 16);
						output.append((char) value);
					}
					catch (NumberFormatException nfe)
					{
						// If the unicode is invalid, then append the code as it is
						output.append(unicode.toString());
	                }
					finally
					{
						unicode.setLength(0);
						inUnicode = hadSlash = false;
					}
	                continue;
	            }
			}

			if (hadSlash)
			{
				// handle an escaped value
				hadSlash = false;
				switch (ch)
				{
					case '\\':
						output.append('\\');
						break;
					case '\'':
						output.append('\'');
						break;
					case '\"':
						output.append('"');
						break;
					case 'r':
						output.append('\r');
						break;
					case 'f':
						output.append('\f');
						break;
					case 't':
						output.append('\t');
						break;
					case 'n':
						output.append('\n');
						break;
					case 'b':
						output.append('\b');
						break;
					case 'u':
						inUnicode = true;
						break;
					default :
						if (!inUnicode)
							output.append(ch);
						break;
				}
				continue;
			}
			else if (ch == '\\')
			{
				hadSlash = true;
				continue;
			}
			if (!inUnicode)
				output.append(ch);
		}

		if (hadSlash)
		{
			// then we're in the weird case of a \ at the end of the
			// string, let's output it anyway.
			output.append('\\');
		}
		return output.toString();
	}

	private static final String parseSection(String m1) {
		int dotPos = m1.indexOf('.');
		if (dotPos < 0)
			throw new RuntimeException("Each matching rule must specify section ID! " + m1);
		return m1.substring(0, dotPos);
	}

	private static final void processMatchingRules(Map<String, Map<String, String>> ini) {
		Map<String, String> mainSection = ini.get(sMain);
		String matchRuleStr = mainSection.get(sMatchingRule);
		if (matchRuleStr == null)
			throw new RuntimeException("Unspecified matching rule!");
		String[] matchingRules = matchRuleStr.split(sCommaSeparator);
		List<String> nonEmptyRule = new ArrayList<String>();
		Set<String> seenSections = new HashSet<String>();
		for (String rule: matchingRules) {
			rule = rule.trim();
			if (rule.length() > 0) {
				int eqPos = rule.indexOf("==");
				if (eqPos < 0)
					throw new RuntimeException("Equal sign (==) is needed in each matching rule!");
				String
					m1 = rule.substring(0, eqPos).trim(),
					m2 = rule.substring(eqPos + 2).trim();
				if (m1.length() == 0 || m2.length() == 0)
					throw new RuntimeException("Malformed matching rule " + rule);
				if (!m1.matches(sValidIDWithDotRegex) || !m2.matches(sValidIDWithDotRegex))
					throw new RuntimeException("Matching rule does not conform to ID requirements" + rule);
				String
					v1 = m1.substring(m1.indexOf('.') + 1).trim(),
					v2 = m2.substring(m2.indexOf('.') + 1).trim(),
					section1 = parseSection(m1),
					section2 = parseSection(m2);

				if (ini.get(section1) == null)
					throw new RuntimeException("Matching rule refers to an unknown section ID "+ section1 + " !");
				if (!section1.equals(sMain)) {
					if (seenSections.contains(section1))
						throw new RuntimeException("Section " + section1 + " is mentioned multiple times!");
					seenSections.add(section1);
					Map<String, String> annotSection = ini.get(section1);
					annotSection.put(sKeyColumn, v1);
				} else {
					mainSection.put(":" + v1 , section2);
				}

				if (ini.get(section2) == null)
					throw new RuntimeException("Matching rule refers to an unknown section ID "+ section2 + " !");
				if (!section2.equals(sMain)) {
					if (seenSections.contains(section2))
						throw new RuntimeException("Section " + section2 + " is mentioned multiple times!");
					seenSections.add(section2);
					Map<String, String> annotSection = ini.get(section2);
					annotSection.put(sKeyColumn, v2);
				} else {
					mainSection.put(":" + v2 , section1);
				}
				nonEmptyRule.add(rule);
			}
		}
		if (nonEmptyRule.size() == 0)
			throw new RuntimeException("Matching rule is blank!");
		//return nonEmptyRule.toArray(new String[nonEmptyRule.size()]);
	}

	public static final void annotate(String inputFn, String outputFn, char delim, String[] columnOrder,
		String[] newHeader, Map<String, String> matchingMap, Map<String, String> newColumns,
		Map<String, Annotation> annotTables) throws IOException {

		int
			n1 = columnOrder == null ? 0 : columnOrder.length,
			n2 = newHeader == null ? 0 : newHeader.length;
		if (n1 > 0 && n2 > 0 && n1 != n2)
			throw new RuntimeException("Error: The length of column order does not equal to the length of the new header!");

		List<String> defColOrder = new ArrayList<String>(); // Default column order, in fully qualified name
		Map<String, Set<String>> colNameToAnnot = new HashMap<String, Set<String>>();

		int[] orderIdx = null, revOrderIdx = null;
		String[] lookupCols = null, origCols = null;
		int numOrigCols = 0;
		System.out.println("Reading " + inputFn);
		CsvReader reader = null;
		PrintWriter writer = null;
		String[] tokensInDefaultOrder = null, keptTokens = null;
		Method method = null;

		long lineNo = 0;
		try {
			writer = new PrintWriter(outputFn);
			reader = new CsvReader(inputFn, new char[] { delim });
			reader.setTrimWhitespace(true);
			reader.setUseTextQualifier(true);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					// Construct default column ordering by lining up the referenced annotation columns FIRST!
					// Along the way, set up the referenced annotation in lookupCols.
					// If entry != null, means reference the annotation named in the entry
					numOrigCols = tokens.length;
					origCols = tokens;
					lookupCols = new String[numOrigCols];
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					Set<String> seenColumns = new HashSet<String>();
					for (int i = 0; i < numOrigCols; i++) {
						String lookupSection = matchingMap.get(":" + tokens[i]);
						if (lookupSection != null) {
							lookupCols[i] = lookupSection;
							Annotation annot = annotTables.get(lookupSection);
							int numCols = annot.colNames.length;
							for (int j = 0; j < numCols; j++) {
								String col = annot.colNames[j];
								String fullyQualColName = lookupSection + "." + col;
								// There might be a duplicate fully qualified name.
								// For example. The user may have section A.B with column C
								// and section A with column B.C. Hopefully this is rare.
								if (colnamesToIdx.get(fullyQualColName) != null)
									throw new RuntimeException("Duplicate fully qualified name '" + fullyQualColName + "', in section '"+ lookupSection + "'");
								colnamesToIdx.put(fullyQualColName, defColOrder.size());

								// If the column is duplicate, make sure it is not indexable,
								// forcing the user to use the fully-qualified name.
								if (seenColumns.contains(col)) {
									colnamesToIdx.remove(col);
								} else {
									colnamesToIdx.put(col, defColOrder.size());
								}
								seenColumns.add(col);
								defColOrder.add(fullyQualColName);
								Set<String> keyset = colNameToAnnot.get(col);
								if (keyset == null) {
									keyset = new HashSet<String>();
									colNameToAnnot.put(col, keyset);
								}
								keyset.add(lookupSection);
							}
						}
					}

					// THEN, line up the columns of the main table
					for (int i = 0; i < numOrigCols; i++) {
						String col = tokens[i];
						String fullyQualColName = sMain + "." + col;
						if (colnamesToIdx.get(fullyQualColName) != null)
							throw new RuntimeException("Duplicate fully qualified name '" + fullyQualColName + "', in main section.");
						colnamesToIdx.put(fullyQualColName, defColOrder.size());
						if (seenColumns.contains(col)) {
							colnamesToIdx.remove(col);
						} else {
							colnamesToIdx.put(col, defColOrder.size());
						}
						seenColumns.add(col);
	
						defColOrder.add(fullyQualColName);
						Set<String> keyset = colNameToAnnot.get(col);
						if (keyset == null) {
							keyset = new HashSet<String>();
							colNameToAnnot.put(col, keyset);
						}
						keyset.add(sMain);
					}

					// Is there a specified column ordering?
					int
						numDefaultOrder = defColOrder.size(),
						numNewColumns = 0;
					Map<String, Integer> referencedRules = new HashMap<String, Integer>();
					if (columnOrder != null) {
						revOrderIdx = new int[numDefaultOrder];
						Arrays.fill(revOrderIdx, -1);
						// If so, make sure that the specified ordering is sane (i.e., no duplicates and fully referenced)
						int numOrder = columnOrder.length;
						orderIdx = new int[numOrder];
						for (int i = 0; i < numOrder; i++) {
							String col = columnOrder[i];
							Set<String> keyset = colNameToAnnot.get(col);
							if (keyset == null) {
								int idx = defColOrder.indexOf(col); // Assume fully qualified name
								if (idx == -1) {
									// This column may have rules specified in newColumns
									if (newColumns.get(col) == null)
										throw new RuntimeException("Error: Column '"+col+"' is not found!");
									referencedRules.put(col, i);
									numNewColumns++;
								} else {
									revOrderIdx[idx] = i;
									if (newColumns.get(col) != null) {
										referencedRules.put(col, i);
										numNewColumns++;
									}
								}
								orderIdx[i] = idx;
							} else {
								if (keyset.size() > 1)
									throw new RuntimeException("Error: Duplicate column '"+col+"'. It can be found in these annotations: " + keyset);
								col = keyset.iterator().next() + "." + col; // Convert to fully qualified name
								int idx = defColOrder.indexOf(col);
								if (idx == -1)
									throw new RuntimeException("Should not happen!");
								orderIdx[i] = idx;
								revOrderIdx[idx] = i;
								if (newColumns.get(col) != null) {
									referencedRules.put(col, i);
									numNewColumns++;
								}
							}
						}
						if (numNewColumns < newColumns.size()) {
							System.out.println("Some of the custom columns are not referenced in the column ordering. These columns are ignored!");
						}
					} else {
						// No column ordering? Proceed with the default ordering.
						List<String> newColumnsList = new ArrayList<String>(newColumns.keySet());
						numNewColumns = newColumns.size();
						orderIdx = new int[numDefaultOrder + numNewColumns];
						for (int i = 0; i < numDefaultOrder; i++)
							orderIdx[i] = i;
						for (int i = 0; i < numNewColumns; i++) {
							orderIdx[i + numDefaultOrder] = -1;
							referencedRules.put(newColumnsList.get(i), i + numDefaultOrder);
						}
						columnOrder = defColOrder.toArray(new String[defColOrder.size()]);
					}
					if (newHeader != null) {
						if (newHeader.length != columnOrder.length)
							throw new RuntimeException("The count of the supplied new header does not match the number of columns post-ordering!");
						columnOrder = newHeader; // Override the header with the supplied header
					}

					// Processing custom rules
					if (numNewColumns > 0) {
						StringBuilder buf = new StringBuilder();
						String className = "Acceptor" + (++classCounter);
						buf.append("import java.util.*;" + ln);
						buf.append("import static java.lang.Math.*;" + ln);
						buf.append("public class " + className + " {" + ln);
						createStandardMethods(buf);
						buf.append(tab + "public static final void "+sTransformMethodName+" (String[] tokens) {" + ln);
						for (String ruleName: referencedRules.keySet()) {
							String rule = newColumns.get(ruleName);
							TokenManager tokmgr = TokenManager.tokenize(rule);
							ASTNode node = ASTBuilder.parseExpr(tokmgr);
							if (node == null)
								throw new RuntimeException("Parse error!");
							ASTBuilder.inferTypes(node);
							List<ASTNode> termNodes = ASTBuilder.getTerminalNodes(node);
							IdentityHashMap<ASTNode, String> symbolTable = new IdentityHashMap<ASTNode, String>();
							for (ASTNode curNode: termNodes)
								if (curNode.getNodeType() == ASTNodeType.IDENTIFIER) {
									String tok = curNode.getToken();
									Integer iidx = colnamesToIdx.get(tok);
									if (iidx == null)
										throw new RuntimeException("Column " + tok + " is not found!");
									int val = revOrderIdx[iidx.intValue()];
									if (val != -1) {
										symbolTable.put(curNode, String.valueOf(val));
									} else {
										// TODO
										throw new RuntimeException("Column " + tok + " is not used!");
									}
								}
							Integer assignee = referencedRules.get(ruleName);
							if (assignee == null)
								throw new RuntimeException("Column " + ruleName + " is not found!");
							JavaProgram prog = new JavaProgram();
							CodeGenerator.generateMethod(node, prog, symbolTable);
							buf.append(tab + tab + "tokens[" + assignee + "] = " + prog.method + ";" + ln);
						}
						buf.append(tab + "}" + ln);
						buf.append("}" + ln);
						System.out.println(buf.toString());
						MemoryClassLoader mcl = new MemoryClassLoader(className, buf.toString());
						try {
							method = mcl.loadClass(className).getMethod(sTransformMethodName, String[].class);
						} catch (Exception e) {
							e.printStackTrace();
						}
					}

					// Write reordered header
					writer.print(columnOrder[0]);
					for (int i = 1; i < columnOrder.length; i++)
						writer.print("," + columnOrder[i]);
					writer.println();
					writer.flush();

					// Prepare temporary buffer for reordering / processing
					tokensInDefaultOrder = new String[defColOrder.size()];
					keptTokens = new String[orderIdx.length];
					continue; // Skip first line
				}

				if (tokens.length != numOrigCols)
					throw new RuntimeException("Error: Number of columns detected in header = " +
						numOrigCols + ", number of current column at line " + lineNo + " = " + tokens.length);
				// Line up the referenced annotation columns first
				int j = 0;
				for (int i = 0; i < numOrigCols; i++) {
					if (lookupCols[i] != null) {
						Annotation ann = annotTables.get(lookupCols[i]);
						String[] lookupResult = ann.table.get(tokens[i]);
						if (lookupResult != null)
							System.arraycopy(lookupResult, 0, tokensInDefaultOrder, j, lookupResult.length);
						else {
							// Warn the lookup failure.
							System.out.println("Warning! Lookup failure at line " + lineNo +
								": Column '" + origCols[i] + "' = '" + tokens[i] + "' could not be found");
						}
						j += ann.colNames.length;
					}
				}
				// Then, line up the columns of the main table
				System.arraycopy(tokens, 0, tokensInDefaultOrder, j, tokens.length);

				// Reorder according to the prescribed ordering;
				for (int i = 0; i < orderIdx.length; i++)
					if (orderIdx[i] >= 0)
						keptTokens[i] = tokensInDefaultOrder[orderIdx[i]];

				// Execute any outstanding transformation
				if (method != null) {
					try {
						method.invoke(null, (Object) keptTokens);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}

				// Write the result
				writer.print(keptTokens[0]);
				for (int i = 1; i < keptTokens.length; i++)
					writer.print("," + keptTokens[i]);
				writer.println();
				writer.flush();
				if (lineNo % 1000000 == 0) System.out.println(lineNo);
			}
			System.out.println(lineNo + " lines were read/written.");
		} finally {
			if (reader != null) reader.close(); reader = null;
			if (writer != null) writer.close(); writer = null;
		}
	}

	private static final String convertFormat(String format) {
		if (format == null) {
			format = ",";
			System.out.println("File format missing, assumed csv");
		} else {
			format = format.toLowerCase(Locale.ENGLISH).trim();
			if (format.equals("csv")) {
				format = ",";
			} else if (format.equals("tab")) {
				format = "\t";
			} else if (format.equals("space")) {
				format = " ";
			} else if (format.startsWith("\"") && format.endsWith("\"")) {
				format = format.substring(1, format.length() - 1);
			}
			format = unescape(format);
		}
		return format;
	}

	private static final String getAbsolutePath(String filename) {
		if (!new File(filename).isAbsolute())
			filename = curDir.getAbsolutePath() + File.separator + filename;
		return filename;
	}

	private static final Annotation parseAnnotationSection(String sectionName, Map<String, String> curSection) throws IOException {
		String keyColumn = curSection.get(sKeyColumn);
		if (keyColumn == null) {
			System.out.println("This section is not linked by the matching rule and thereby ignored!");
			return null;
		}
		String
			filename = curSection.get(sFile),
			format = convertFormat(curSection.get(sFileFormat)),
			keptColsStr = curSection.get(sColumnsKept),
			colRenameStr = curSection.get(sColumnsRename);
		if (filename == null)
			throw new RuntimeException("File name missing!");
		System.out.println("File name: "+ filename);
		filename = getAbsolutePath(filename);
		String[] keptCols = null, renamedCols = null;
		if (keptColsStr != null) {
			keptCols = keptColsStr.split(sCommaSeparator);
			if (colRenameStr != null) {
				renamedCols = colRenameStr.split(sCommaSeparator);
				if (keptCols.length != renamedCols.length)
					throw new RuntimeException("Error: The length of columns kept does not equal to the length of renamed columns!");
			}
		}
		System.out.println("Loading...");
		long time1, time2;
		Map<String, String> ruleMap = new HashMap<String, String>(curSection);
		ruleMap.keySet().removeAll(defaultSectionKeys);
		time1 = System.currentTimeMillis();
		Annotation ann = loadAnnotationFile(filename, format.charAt(0), keyColumn, keptCols, renamedCols, ruleMap);
		time2 = System.currentTimeMillis();
		ann.name = sectionName;
		System.out.println("Annotation read " + (time2 - time1) + " ms");
		return ann;
	}

	public static void main(String[] args)
	{
		long time1, time2;
		try {
			curDir = new File(args[0]).getParentFile();
			Map<String, Map<String, String>> ini = readIniFile(args[0], true);
			Map<String, String> mainSection = ini.get(sMain);
			if (mainSection == null)
				throw new RuntimeException("Unspecified main section!");
			processMatchingRules(ini);
			Set<String> keys = new HashSet<String>();
			keys.addAll(ini.keySet());
			keys.remove(sMain);
			Map<String, Annotation> annotTables = new HashMap<String, Annotation>();
			for (String key: keys) {
				System.out.println("\nSection " + key);
				annotTables.put(key, parseAnnotationSection(key, ini.get(key)));
			}

			System.out.println("Processing main section.");
			String
				inputFn = getAbsolutePath(mainSection.get(sInput)),
				outputFn = getAbsolutePath(mainSection.get(sOutput)),
				format = convertFormat(mainSection.get(sFileFormat)),
				colOrderStr = mainSection.get(sColumnOrder),
				newHeaderStr = mainSection.get(sNewHeader);
			if (inputFn == null)
				throw new RuntimeException("Input file name missing!");
			if (outputFn == null)
				throw new RuntimeException("Output file name missing!");
			String[] colOrder = null, newHeader = null;
			if (colOrderStr != null) {
				colOrder = colOrderStr.split(sCommaSeparator);
				if (newHeaderStr != null) {
					newHeader = newHeaderStr.split(sCommaSeparator);
				}
			}
			Map<String, String> ruleMap = new HashMap<String, String>(mainSection);
			Map<String, String> matchingMap = new HashMap<String, String>();
			ruleMap.keySet().removeAll(defaultMainSectionKeys);
			for (String key: ruleMap.keySet()) {
				if (key.startsWith(":"))
					matchingMap.put(key, ruleMap.get(key));
			}
			for (String key: matchingMap.keySet())
				ruleMap.remove(key);

			time1 = System.currentTimeMillis();
			annotate(inputFn, outputFn, format.charAt(0), colOrder, newHeader, matchingMap, ruleMap, annotTables);
			time2 = System.currentTimeMillis();
			System.out.println("Annotation time " + (time2 - time1) + " ms");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class Annotation {
	public String name, colNames[];
	public Map<String, String[]> table;
	private Map<String, Integer> colNameToIdx = null;

	public Annotation(String annName, String[] columnNames, Map<String, String[]> tbl) {
		name = annName;
		colNames = columnNames;
		table = tbl;
	}

	private void buildIndex() {
		colNameToIdx = new HashMap<String, Integer>();
		for (int i = 0; i < colNames.length; i++)
			colNameToIdx.put(colNames[i], i);
	}

	public int getColumnIndex(String str) {
		if (colNameToIdx == null)
			buildIndex();
		Integer idx = colNameToIdx.get(str);
		return idx == null ? -1 : idx.intValue();
	}
}