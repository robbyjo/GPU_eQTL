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
package gov.nih.tools.old;

import gov.nih.tools.IntervalTree.IntervalData;

import java.io.IOException;
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
public class AnnotationCommon {
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

	public static final ParsedAnnotation loadHG133U(String filename) throws NumberFormatException, IOException {
		Map<String, List<IntervalData<String>>> probeIDToNodes = new HashMap<String, List<IntervalData<String>>>();
		Map<String, String[]> probeAnnotTable = new HashMap<String, String[]>();
		CsvReader reader = null;
		boolean encounteredHeader = false;
		int lineNo = 0, pidx = -1, gsymidx = -1, entrezidx = -1, alignidx = -1;
		reader = new CsvReader(filename);
		reader.setTrimWhitespace(true);
		reader.setUseTextQualifier(true);
		while(reader.readRecord()) {
			String[] tokens = reader.getValues();
			lineNo++;
			if (!encounteredHeader) {
				if (!tokens[0].startsWith("#")) {
					Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i].toLowerCase(), i);
					encounteredHeader = true;
					pidx = colnamesToIdx.get("probe set id");
					gsymidx = colnamesToIdx.get("gene symbol");
					entrezidx = colnamesToIdx.get("entrez gene");
					alignidx = colnamesToIdx.get("alignments");
				}
				continue; // Skip until the right header is met
			}
			String
			probeid = tokens[pidx],
			genesym = tokens[gsymidx].replaceAll(" /// ", "|"),
			entrez = tokens[entrezidx].replaceAll(" /// ", "|"),
			align = tokens[alignidx];
			probeAnnotTable.put(probeid, new String[] {
					probeid, genesym, entrez, align
			});

			for (String cur_align: align.split(" /// ")) {
				int pos = cur_align.indexOf("(");
				if (pos == -1)
					continue;
				cur_align = cur_align.substring(0, pos - 1);
				pos = cur_align.indexOf(":");
				String chr = cur_align.substring(0, pos);
				if (chr.startsWith("chr")) chr = chr.substring(3);
				int pos2 = cur_align.indexOf("-");
				String start = cur_align.substring(pos+1, pos2), end = cur_align.substring(pos2+1);

				List<IntervalData<String>> ivalList = probeIDToNodes.get(chr);
				if (ivalList == null) {
					ivalList = new ArrayList<IntervalData<String>>();
					probeIDToNodes.put(chr, ivalList);
				}
				IntervalData<String> node = new IntervalData<String>(Long.parseLong(start), Long.parseLong(end), probeid);
				ivalList.add(node);
			}
		}
		System.out.println(lineNo + " lines were read.");
		reader.close(); reader = null;
		return new ParsedAnnotation(probeAnnotTable, probeIDToNodes, "HG133UP2_ProbesetID,HG133UP2_GeneSymbol,HG133UP2_EntrezID,HG133UP2_Alignments");
	}

	public static final ParsedAnnotation loadIllumina(String filename) throws IOException {
		Map<String, List<IntervalData<String>>> probeIDToNodes = new HashMap<String, List<IntervalData<String>>>();
		Map<String, String[]> probeAnnotTable = new HashMap<String, String[]>();
		Map<String, Integer> colnamesToIdx = new HashMap<String, Integer>();
		CsvReader reader = null;
		int lineNo = 0, pidx = -1, arraddridx = -1, chridx = -1, strandidx = -1, genesymidx = -1, entrezidx = -1, alignidx = -1,
				synidx = -1, obsoleteidx = -1;
		boolean encounteredHeader = false;
		reader = new CsvReader(filename);
		reader.setTrimWhitespace(true);
		reader.setUseTextQualifier(true);
		while(reader.readRecord()) {
			String[] tokens = reader.getValues();
			lineNo++;
			if (!encounteredHeader) {
				if (tokens[0].equalsIgnoreCase("species")) {
					for (int i = 0; i < tokens.length; i++)
						colnamesToIdx.put(tokens[i].toLowerCase(), i);
					pidx = colnamesToIdx.get("probe_id");
					arraddridx = colnamesToIdx.get("array_address_id");
					chridx = colnamesToIdx.get("chromosome");
					strandidx = colnamesToIdx.get("probe_chr_orientation");
					genesymidx = colnamesToIdx.get("symbol");
					entrezidx = colnamesToIdx.get("entrez_gene_id");
					alignidx = colnamesToIdx.get("probe_coordinates");
					synidx = colnamesToIdx.get("synonyms");
					obsoleteidx = colnamesToIdx.get("obsolete_probe_id");
					encounteredHeader = true;
				}
				continue; // Skip until the right header is met
			}
			if (tokens[0].startsWith("["))
				break;
			String
				probeid = tokens[pidx],
				arrayaddr = tokens[arraddridx],
				chr = tokens[chridx],
				strand = tokens[strandidx],
				genesym = tokens[genesymidx],
				entrez = tokens[entrezidx],
				align = tokens[alignidx];
			Set<String> syn = parse(tokens[synidx], "; ");
			syn.addAll(parse(tokens[obsoleteidx], "; "));
			if (!genesym.equals("")) syn.add(genesym);
			genesym = setToString(syn, "|");
			probeAnnotTable.put(probeid, new String[] {
					probeid, arrayaddr, genesym, entrez, strand, chr, align
			});

			for (String cur_align : align.split(":")) {
				int pos = cur_align.indexOf("-");
				if (pos == -1)
					continue;
				String start = cur_align.substring(0, pos), end = cur_align.substring(pos+1);

				List<IntervalData<String>> ivalList = probeIDToNodes.get(chr);
				if (ivalList == null) {
					ivalList = new ArrayList<IntervalData<String>>();
					probeIDToNodes.put(chr, ivalList);
				}
				IntervalData<String> node = new IntervalData<String>(Long.parseLong(start), Long.parseLong(end), probeid);
				ivalList.add(node);
			}
		}
		System.out.println(lineNo + " lines were read.");
		reader.close(); reader = null;
		return new ParsedAnnotation(probeAnnotTable, probeIDToNodes, "HT_ProbeID,HT_ArrayAddress,HT_GeneSymbol,HT_EntrezID,HT_Strand,HT_Chr,HT_Alignment");
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
			IntervalData<String> node = new IntervalData<String>(Long.parseLong(start), Long.parseLong(stop), probeid);
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

	public static final ParsedAnnotation parseAnnotation(String filename, String type) throws NumberFormatException, IOException {
		type = type.toLowerCase();
		if (type.equals("hg133up2"))
			return loadHG133U(filename);
		if (type.equals("ht"))
			return loadIllumina(filename);
		if (type.equals("huex"))
			return loadHuEx(filename);
		throw new RuntimeException("Unknown annotation type " + type);
	}
}
