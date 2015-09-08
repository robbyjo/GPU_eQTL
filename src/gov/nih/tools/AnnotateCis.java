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
import gov.nih.tools.IntervalTree.IntervalData;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.csvreader.CsvReader;

/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class AnnotateCis {
	static final long kLongCis = 5000000, kShortCis = 1000000;
	static class ParsedAnnotation {
		public Map<String, List<IntervalData<String>>> probeIDToNodes;
		public Map<String, String[]> probeAnnotTable;
		public String header;
		public ParsedAnnotation(Map<String, String[]> annottbl, Map<String, List<IntervalData<String>>> nodestbl, String hdr) {
			probeAnnotTable = annottbl;
			probeIDToNodes = nodestbl;
			header = hdr;
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
		Map<String, List<IntervalData<String>>> probeIDToNodes = new HashMap<String, List<IntervalData<String>>>();
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
			List<IntervalData<String>> ivalList = probeIDToNodes.get(chr);
			if (ivalList == null) {
				ivalList = new ArrayList<IntervalData<String>>();
				probeIDToNodes.put(chr, ivalList);
			}
			long startPos = Long.parseLong(start);
			IntervalData<String> node = new IntervalData<String>(startPos > kLongCis ? startPos - kLongCis : 0, startPos + kLongCis, probeid);
			ivalList.add(node);
		}
		System.out.println(lineNo + " lines were read.");
		reader.close(); reader = null;
		return new ParsedAnnotation(probeAnnotTable, probeIDToNodes, "HuEx_ProbesetID,HuEx_Transcript_Cluster_ID,HuEx_GeneSymbol,HuEx_EntrezID,HuEx_Strand,HuEx_Chr,HuEx_Start,HuEx_Stop");
	}

	/*
	static void example() {
		List<IntervalData<String>> ll = new ArrayList<IntervalData<String>>();
		ll.add(new IntervalData<String>(1, 4, "A"));
		ll.add(new IntervalData<String>(5, 6, "A"));
		ll.add(new IntervalData<String>(3, 8, "B"));
		ll.add(new IntervalData<String>(5, 12, "C"));
		ll.add(new IntervalData<String>(6, 10, "D"));
		ll.add(new IntervalData<String>(11, 15, "E"));
		ll.add(new IntervalData<String>(13, 18, "F"));
		IntervalTree<String> tree = new IntervalTree<String>(ll);
		IntervalData<String> result = tree.query(6);
		System.out.println(result.getSet());
	}
	//*/

	public static void main(String[] args) {
		ParsedAnnotation parsedAnnotTable = null;
		Map<String, IntervalTree<String>> probeIDToTree = new HashMap<String, IntervalTree<String>>();
		CsvReader reader = null;
		PrintWriter writer = null;
		long lineNo = 0;
		long time1, time2;
		try {
			time1 = System.currentTimeMillis();
			parsedAnnotTable = loadHuEx(args[0]);
			time2 = System.currentTimeMillis();
			System.out.println(args[0] + " annotation read " + (time2 - time1) + " ms");

			reader = new CsvReader(args[1]);
			reader.setTrimWhitespace(true);
			reader.setUseTextQualifier(true);
			writer = new PrintWriter(args[2]);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i], i);
					writer.println(reader.getRawRecord() + ",LongCis5MGene");
					writer.flush();
					continue; // Skip first line
				}
				String chr = tokens[1];
				IntervalTree<String> tree = probeIDToTree.get(chr);
				if (tree == null) // Unlikely
					continue;
				IntervalTree.IntervalData<String> data = tree.query((long) Double.parseDouble(tokens[2])); // Position in SNP annotation
				if (data == null) // Not found
					continue;
				long pos = (long) Double.parseDouble(tokens[2]); // Position in SNP annotation
				int count = 0;
				for (String id: data.getSet()) {
					String[] toks = parsedAnnotTable.probeAnnotTable.get(id);
					long deltapos = Math.abs(pos - ((long) Double.parseDouble(toks[6])));
					if (deltapos > kShortCis) count++;
				}
				writer.println(reader.getRawRecord()+","+ count);
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
