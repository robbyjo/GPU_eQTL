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
public class CisAnnotation {
	static final int cisLimit = 1000000;
	static String snpAnnotHeader;

	public static final HashMap<String, List<AnnotatedLine>> loadTranscriptAnnot(String filename) {
		HashMap<String, List<AnnotatedLine>> tbl = new HashMap<String, List<AnnotatedLine>>(283805);
		CsvReader reader = null;
		int lineNo = 0;
		int[] idx = new int[7];
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
					if (colnamesToIdx.get("probeset_id") == null)
						colnamesToIdx.put("probeset_id", colnamesToIdx.get("transcript_cluster_id"));
					idx[0] = colnamesToIdx.get("probeset_id");
					idx[1] = colnamesToIdx.get("transcript_cluster_id");
					idx[2] = colnamesToIdx.get("seqname");
					idx[3] = colnamesToIdx.get("strand");
					idx[4] = colnamesToIdx.get("start");
					idx[5] = colnamesToIdx.get("stop");
					idx[6] = colnamesToIdx.get("GeneSymbol");
					continue; // Skip first line
				}
				String chr = tokens[idx[2]];
				if (chr.startsWith("chr")) chr = chr.substring(3);
				AnnotatedLine ann = new AnnotatedLine(null, chr, Long.valueOf(tokens[idx[4]].trim()), Long.valueOf(tokens[idx[5]].trim()));
				List<AnnotatedLine> list = tbl.get(chr);
				if (list == null) {
					list = new ArrayList<AnnotatedLine>();
					tbl.put(chr, list);
				}
				list.add(ann);
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
		//System.out.println("Table:");
		//for(String key: tbl.keySet()) {
		//	System.out.println(key + "="+tbl.get(key).size());
		//}
		return tbl;
	}

	public static final HashMap<String, List<AnnotatedLine>> loadSNPAnnot(String filename) {
		HashMap<String, List<AnnotatedLine>> tbl = new HashMap<String, List<AnnotatedLine>>(39315185);
		CsvReader reader = null;
		int lineNo = 0;
		try {
			reader = new CsvReader(filename);
			reader.setTrimWhitespace(true);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					snpAnnotHeader = reader.getRawRecord();
					continue; // Skip first line
				}
				String chr = tokens[1];
				AnnotatedLine ann = new AnnotatedLine(reader.getRawRecord(), chr, (long) Double.parseDouble(tokens[2]), 0);
				List<AnnotatedLine> list = tbl.get(chr);
				if (list == null) {
					list = new ArrayList<AnnotatedLine>();
					tbl.put(chr, list);
				}
				list.add(ann);
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return tbl;
	}

	public static final void radixSort(AnnotatedLine[] a)
	{
		// The positions are assumed to be always positive!
		int
			n = a.length,
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2)))));
		AnnotatedLine[] b = new AnnotatedLine[n], b_orig = b;
		for (long mask = ~(-1L << numBins), rshift = 0L; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(int) ((a[p].pos_start & mask) >>> rshift)];
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (int) ((a[p].pos_start & mask) >>> rshift);
				--cntarray[key];
				b[cntarray[key]] = a[p];
			}
			AnnotatedLine[] temp = b; b = a; a = temp;
		}

		if (a == b_orig)
			System.arraycopy(a, 0, b, 0, n);
	}

	public static final void computeCis(String outputFn,
			HashMap<String, List<AnnotatedLine>> snpAnnotTable,
			HashMap<String, List<AnnotatedLine>> transcriptAnnotTable,
			HashMap<String, List<AnnotatedLine>> probesetAnnotTable) {
		PrintWriter writer = null;

		try {
			writer = new PrintWriter(outputFn);
			writer.println(snpAnnotHeader+",ProbesetCis,TranscriptCis");
			for (int ch = 1; ch <= 23; ch++) {
				String chr = ch == 23 ? "X" : String.valueOf(ch);
				if (snpAnnotTable.get(chr) == null) continue;
				AnnotatedLine[]
					snpAnnots = snpAnnotTable.get(chr).toArray(new AnnotatedLine[snpAnnotTable.get(chr).size()]),
					geneAnnots = transcriptAnnotTable.get(chr).toArray(new AnnotatedLine[transcriptAnnotTable.get(chr).size()]),
					exonAnnots = probesetAnnotTable.get(chr).toArray(new AnnotatedLine[probesetAnnotTable.get(chr).size()]);
				radixSort(snpAnnots);
				radixSort(geneAnnots);
				radixSort(exonAnnots);
				int
					geneStartIdx = 0,
					exonStartIdx = 0,
					numGenes = geneAnnots.length,
					numExons = exonAnnots.length;
				for (AnnotatedLine snpAnnot: snpAnnots) {
					int geneCisCount = 0, exonCisCount = 0;
					long pos = snpAnnot.pos_start;
					while (geneStartIdx < numGenes && pos - geneAnnots[geneStartIdx].pos_start > cisLimit) {
						geneStartIdx++;
					}
					while((geneStartIdx + geneCisCount) < numGenes && geneAnnots[geneStartIdx + geneCisCount].pos_start - pos <= cisLimit) {
						geneCisCount++;
					}

					while (exonStartIdx < numExons && pos - exonAnnots[exonStartIdx].pos_start > cisLimit) {
						exonStartIdx++;
					}
					while((exonStartIdx + exonCisCount) < numExons && exonAnnots[exonStartIdx + exonCisCount].pos_start - pos <= cisLimit) {
						exonCisCount++;
					}

					//String[] tokens = snpAnnot.line.split(",");
					//if (Integer.parseInt(tokens[tokens.length - 2]) != exonCisCount || Integer.parseInt(tokens[tokens.length - 1]) != geneCisCount)
					//	System.out.println("Here!");
					writer.println(snpAnnot.line+","+exonCisCount+","+geneCisCount);
				}
			}
			writer.close(); writer = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args)
	{
		long time1, time2;
		try {
			time1 = System.currentTimeMillis();
			HashMap<String, List<AnnotatedLine>> transcriptAnnotTable = loadTranscriptAnnot(args[1]);
			time2 = System.currentTimeMillis();
			System.out.println("Gene-level transcript annotation read " + (time2 - time1) + " ms");
			time1 = System.currentTimeMillis();
			HashMap<String, List<AnnotatedLine>> probesetAnnotTable = loadTranscriptAnnot(args[2]);
			time2 = System.currentTimeMillis();
			System.out.println("Exon-level transcript annotation read " + (time2 - time1) + " ms");
			time1 = System.currentTimeMillis();
			HashMap<String, List<AnnotatedLine>> snpAnnotTable = loadSNPAnnot(args[0]);
			time2 = System.currentTimeMillis();
			System.out.println("SNP annotation read " + (time2 - time1) + " ms");
			time1 = System.currentTimeMillis();
			computeCis(args[3], snpAnnotTable, transcriptAnnotTable, probesetAnnotTable);
			time2 = System.currentTimeMillis();
			System.out.println("Cis computation time " + (time2 - time1) + " ms");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class AnnotatedLine {
	String chrom, line;
	long pos_start, pos_end;
	AnnotatedLine() {}
	AnnotatedLine(String l, String c, long s, long e)
	{	chrom = c; line = l; pos_start = s; pos_end = e; }
}
