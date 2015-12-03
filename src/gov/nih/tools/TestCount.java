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

import java.io.FileReader;
import java.io.LineNumberReader;
import java.util.*;

/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class TestCount {
	static HashSet<String>
		markers_cis = new HashSet<String>(20000000),
		markers_trans = new HashSet<String>(20000000),
		probesets_cis = new HashSet<String>(300000),
		probesets_trans = new HashSet<String>(300000);
	static int pairCount_cis = 0, pairCount_trans = 0;

	public static void main(String[] args)
	{
		String pString = args[0];
		int dotPos = pString.indexOf('=');
		boolean pvalMode = true;
		if (dotPos >= 0) {
			String mode = pString.substring(0, dotPos).trim();
			if (mode.equalsIgnoreCase("fdr"))
				pvalMode = false;
			pString = pString.substring(dotPos + 1);
		}
		double pval = Double.parseDouble(pString);
		LineNumberReader reader = null;
		long lineNo = 0;
		String curLine = null;
		int
			pvalCol = -1,
			fdrCol = -1,
			transcriptCol = -1,
			markerCol = -1,
			isCisCol = -1;
		try {
			reader = new LineNumberReader(new FileReader(args[1]));
			do {
				curLine = reader.readLine();
				if (curLine == null)
					break;
				lineNo++;
				String[] tok = curLine.split(",");
				if (lineNo == 1) {
					for (int i = 0; i < tok.length; i++) {
						String t = tok[i];
						if ("log10P".equalsIgnoreCase(t))
							pvalCol = i;
						else if ("log10FDR".equalsIgnoreCase(t))
							fdrCol = i;
						else if ("Marker".equalsIgnoreCase(t))
							markerCol = i;
						else if ("CpG".equalsIgnoreCase(t))
							transcriptCol = i;
						else if ("Is_Cis".equalsIgnoreCase(t))
							isCisCol = i;
					}
					if (!pvalMode)
						pvalCol = fdrCol;
					if (pvalCol == -1 || markerCol == -1 || transcriptCol == -1 || isCisCol == -1) {
						System.err.println("pvalCol = " + pvalCol);
						System.err.println("markerCol = " + markerCol);
						System.err.println("transcriptCol = " + transcriptCol);
						System.err.println("isCisCol = " + isCisCol);
						reader.close(); reader = null;
						throw new RuntimeException();
					}
					continue;
				}
				if (lineNo % 10000000 == 0) System.out.println(lineNo/1000000);
				double curPval = Double.parseDouble(tok[pvalCol]);
				int isCis = Integer.parseInt(tok[isCisCol]);
				if (curPval < pval) {
					if (isCis == 1) {
						markers_cis.add(tok[markerCol]);
						probesets_cis.add(tok[transcriptCol]);
						pairCount_cis++;
					} else {
						markers_trans.add(tok[markerCol]);
						probesets_trans.add(tok[transcriptCol]);
						pairCount_trans++;
					}
				} else break;
			} while (true);
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("There are " + pairCount_cis + " cis pairs and " + pairCount_trans + " trans pairs.");
		System.out.println("There are " + markers_cis.size() + " unique cis SNPs and " + markers_trans.size() + " unique trans SNPs.");
		System.out.println("There are " + probesets_cis.size() + " unique cis probesets and " + probesets_trans.size() + " unique trans probesets.");
		markers_cis.addAll(markers_trans);
		System.out.println("There are " + markers_cis.size() + " unique SNPs total");
		markers_cis.clear(); markers_cis = null;
		markers_trans.clear(); markers_trans = null;
		probesets_cis.addAll(probesets_trans);
		System.out.println("There are " + probesets_cis.size() + " unique probesets total");
	}
}
