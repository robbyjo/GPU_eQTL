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

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.csvreader.CsvReader;

/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class mQTLCis {
	static HashMap<String, String[]>
		probeAnnotTable = new HashMap<String, String[]>(),
		snpAnnotTable = new HashMap<String, String[]>();
	static HashMap<String, IntervalTree<String>>
		probeIDToTree = new HashMap<String, IntervalTree<String>>();
	static final int cisLimit = 1000000;
	static int pairCount = 0;

	public static final void loadProbeAnnot(String filename) {
		CsvReader reader = null;
		HashMap<String, List<IntervalTree.IntervalData<String>>>
			probeIDToNodes = new HashMap<String, List<IntervalTree.IntervalData<String>>>();
		int lineNo = 0;
		try {
			reader = new CsvReader(filename);
			reader.setTrimWhitespace(true);
			reader.setUseTextQualifier(true);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i], i);
					continue; // Skip first line
				}
				String
					chr = tokens[11],
					pos_str = tokens[12];
				if (chr == "") continue;
				//if (chr.startsWith("chr")) chr = chr.substring(3);
				probeAnnotTable.put(tokens[0], new String[] {
					tokens[0],chr,pos_str });
				List<IntervalTree.IntervalData<String>> ivalList = probeIDToNodes.get(chr);
				if (ivalList == null) {
					ivalList = new ArrayList<IntervalTree.IntervalData<String>>();
					probeIDToNodes.put(chr, ivalList);
				}
				try {
				long pos = (long) Double.parseDouble(pos_str);
				IntervalTree.IntervalData<String> node = new IntervalTree.IntervalData<String>(Math.max(0, pos - cisLimit), pos + cisLimit, tokens[0]);
				ivalList.add(node);
				} catch (Exception ee) {
					System.out.println("Error: Line " + lineNo + ": ID = " + tokens[0] + ", Chr: " + chr + ", Pos = " + pos_str);
				}
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
		// Construct tree for each chromosome
		for(String key: probeIDToNodes.keySet())
			probeIDToTree.put(key, new IntervalTree(probeIDToNodes.get(key)));
		probeIDToNodes.clear();
	}

	public static void main(String[] args) {
		CsvReader reader = null;
		PrintWriter writer = null;
		long lineNo = 0, writtenNo = 0;
		long time1, time2;
		try {
			time1 = System.currentTimeMillis();
			loadProbeAnnot(args[0]);
			time2 = System.currentTimeMillis();
			System.out.println("Probe annotation read " + (time2 - time1) + " ms");
			time1 = System.currentTimeMillis();
			reader = new CsvReader(args[1]);
			reader.setTrimWhitespace(true);
			reader.setUseTextQualifier(true);
			writer = new PrintWriter(args[2]);
			writer.println("Marker,Rs_ID,CpG_ID,Delta_Pos");
			writer.flush();
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i], i);
					continue; // Skip first line
				}
				String chr = tokens[1];
				IntervalTree<String> tree = probeIDToTree.get(chr);
				if (tree == null) // Unlikely
					continue;
				IntervalTree.IntervalData<String> data = tree.query((long) Double.parseDouble(tokens[2]));
				if (data == null) // Not found
					continue;
				long pos = (long) Double.parseDouble(tokens[2]);
				for (String id: data.getSet()) {
					String[] toks = probeAnnotTable.get(id);
					long deltapos = (pos - ((long) Double.parseDouble(toks[2])));
					writer.println(tokens[0]+","+tokens[3]+","+toks[0]+","+ deltapos);
					writtenNo++;
				}
				writer.flush();
			}
			System.out.println(lineNo + " lines were read.");
			System.out.println(writtenNo + " lines were written.");
			reader.close(); reader = null;
			time2 = System.currentTimeMillis();
			System.out.println("mQTL cis annotation time " + (time2 - time1) + " ms");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
