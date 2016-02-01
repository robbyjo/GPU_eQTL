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

import gov.nih.tools.IntervalTree;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.csvreader.CsvReader;

// java -cp .:javacsv.jar AnnotateCis HuEx-1_0-st-v2.na33.1.hg19.transcript-core.txt HuEx-1_0-st-v2.na33.1.hg19.probeset-core.txt 1000G-annot-longcis.txt 1000G-annot-longcis-chr.txt
/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class AnnotateCis {
	static final long kLongCis = Long.MAX_VALUE, kShortCis = 1000000;
	static class ParsedAnnotation {
		public Map<String, IntervalTree<String>> probeIDToTree;
		public Map<String, String[]> probeAnnotTable;
		public ParsedAnnotation(Map<String, String[]> annottbl, Map<String, IntervalTree<String>> treeTbl) {
			probeAnnotTable = annottbl;
			probeIDToTree = treeTbl;
		}
	}

	static final Set<String> parse(String token, String sep) {
		Set<String> set = new HashSet<String>();
		for (String curtok : token.split(sep)) {
			set.add(curtok);
		}
		return set;
	}

	static final String setToString(Set<String> set, String sep) {
		StringBuilder buf = new StringBuilder();
		boolean first = true;
		for (String tok : set) {
			if (!first)
				buf.append(sep);
			buf.append(tok);
			first = false;
		}
		return buf.toString();
	}

	public static final ParsedAnnotation loadHuEx(String filename) throws NumberFormatException, IOException {
		Map<String, List<IntervalTree.IntervalData<String>>> probeIDToNodes = new HashMap<String, List<IntervalTree.IntervalData<String>>>();
		Map<String, String[]> probeAnnotTable = new HashMap<String, String[]>();
		Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
		CsvReader reader = null;
		int lineNo = 0, txidx = -1, psidx = -1, strandidx = -1, chridx = -1, startidx = -1, stopidx = -1, symidx = -1, entrezidx = -1;
		boolean encounteredHeader = false;
		reader = new CsvReader(filename);
		reader.setTrimWhitespace(true);
		reader.setUseTextQualifier(true);
		while(reader.readRecord()) {
			String[] tokens = reader.getValues();
			lineNo++;
			if (!encounteredHeader) {
				if (!tokens[0].startsWith("#")) {
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i].toLowerCase(), i);
					txidx = colnamesToIdx.get("transcript_cluster_id");
					psidx = colnamesToIdx.get("probeset_id") == null ? txidx : colnamesToIdx.get("probeset_id");
					chridx = colnamesToIdx.get("seqname");
					strandidx = colnamesToIdx.get("strand");
					startidx = colnamesToIdx.get("start");
					stopidx = colnamesToIdx.get("stop");
					symidx = colnamesToIdx.get("genesymbol");
					entrezidx = colnamesToIdx.get("entrezgeneid");
					encounteredHeader = true;
				}
				continue; // Skip until the right header is met
			}
			String
				txid = tokens[txidx],
				probeid = tokens[psidx],
				chr = tokens[chridx],
				strand = tokens[strandidx],
				genesym = tokens[symidx],
				entrez = tokens[entrezidx],
				start = tokens[startidx],
				stop = tokens[stopidx];
			if (chr.startsWith("chr")) chr = chr.substring(3);
			probeAnnotTable.put(probeid, new String[] {
					probeid, txid, genesym, entrez, strand, chr, start, stop
			});
			List<IntervalTree.IntervalData<String>> ivalList = probeIDToNodes.get(chr);
			if (ivalList == null) {
				ivalList = new ArrayList<IntervalTree.IntervalData<String>>();
				probeIDToNodes.put(chr, ivalList);
			}
			//long startPos = Long.parseLong(start);
			IntervalTree.IntervalData<String> node = new IntervalTree.IntervalData<String>(0, kLongCis, probeid);
			ivalList.add(node);
		}
		System.out.println(lineNo + " lines were read.");
		reader.close(); reader = null;
		Map<String, IntervalTree<String>> treeTbl = new HashMap<String, IntervalTree<String>>();
		for (String chr: probeIDToNodes.keySet()) {
			List<IntervalTree.IntervalData<String>> ivalList = probeIDToNodes.get(chr);
			treeTbl.put(chr, new IntervalTree<String>(ivalList));
		}
		return new ParsedAnnotation(probeAnnotTable, treeTbl);
	}

	/*
	static void example() {
		List<IntervalTree.IntervalData<String>> ll = new ArrayList<IntervalTree.IntervalData<String>>();
		ll.add(new IntervalTree.IntervalData<String>(1, 4, "A"));
		ll.add(new IntervalTree.IntervalData<String>(5, 6, "A"));
		ll.add(new IntervalTree.IntervalData<String>(3, 8, "B"));
		ll.add(new IntervalTree.IntervalData<String>(5, 12, "C"));
		ll.add(new IntervalTree.IntervalData<String>(6, 10, "D"));
		ll.add(new IntervalTree.IntervalData<String>(11, 15, "E"));
		ll.add(new IntervalTree.IntervalData<String>(13, 18, "F"));
		IntervalTree<String> tree = new IntervalTree<String>(ll);
		IntervalTree.IntervalData<String> result = tree.query(6);
		System.out.println(result.getSet());
	}
	//*/

	public static void main(String[] args) {
		ParsedAnnotation parsedAnnotTableGene = null, parsedAnnotTableExon = null;
		CsvReader reader = null;
		PrintWriter writer = null;
		long lineNo = 0;
		long time1, time2;
		try {
			time1 = System.currentTimeMillis();
			parsedAnnotTableGene = loadHuEx(args[0]);
			time2 = System.currentTimeMillis();
			System.out.println(args[0] + " annotation read " + (time2 - time1) + " ms");

			time1 = System.currentTimeMillis();
			parsedAnnotTableExon = loadHuEx(args[1]);
			time2 = System.currentTimeMillis();
			System.out.println(args[1] + " annotation read " + (time2 - time1) + " ms");

			reader = new CsvReader(args[2]);
			reader.setTrimWhitespace(true);
			reader.setUseTextQualifier(true);
			writer = new PrintWriter(args[3]);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i], i);
					writer.println(reader.getRawRecord() + ",LongCisChrGene,LongCisChrExon");
					writer.flush();
					continue; // Skip first line
				}
				String chr = tokens[1];
				long pos = (long) Double.parseDouble(tokens[2]); // Position in SNP annotation
				int countGene = 0, countExon = 0;

				IntervalTree<String> tree = parsedAnnotTableGene.probeIDToTree.get(chr);
				if (tree != null) {
					IntervalTree.IntervalData<String> data = tree.query(pos);
					if (data != null) {
						for (String id: data.getSet()) {
							String[] toks = parsedAnnotTableGene.probeAnnotTable.get(id);
							long deltapos = Math.abs(pos - ((long) Double.parseDouble(toks[6])));
							if (deltapos > kShortCis) countGene++;
						}
					}
				}

				tree = parsedAnnotTableExon.probeIDToTree.get(chr);
				if (tree != null) {
					IntervalTree.IntervalData<String> data = tree.query(pos);
					if (data != null) {
						for (String id: data.getSet()) {
							String[] toks = parsedAnnotTableExon.probeAnnotTable.get(id);
							long deltapos = Math.abs(pos - ((long) Double.parseDouble(toks[6])));
							if (deltapos > kShortCis) countExon++;
						}
					}
				}

				writer.println(reader.getRawRecord()+","+ countGene+","+ countExon);
				writer.flush();
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
			time2 = System.currentTimeMillis();
			System.out.println("Long cis annotation time " + (time2 - time1) + " ms");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
