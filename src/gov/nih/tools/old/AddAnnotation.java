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

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import com.csvreader.CsvReader;

/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class AddAnnotation {
	static HashMap<String, String[]>
		transcriptAnnotTable = new HashMap<String, String[]>(283805),
		snpAnnotTable = new HashMap<String, String[]>(39315185);
	static final int cisLimit = 1000000;
	static int pairCount = 0;

	public static final void loadTranscriptAnnot(String filename) {
		CsvReader reader = null;
		int lineNo = 0;
		int[] idx = new int[9];
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
					idx[7] = colnamesToIdx.get("EntrezGeneID");
					idx[8] = colnamesToIdx.get("AccessionID");
					continue; // Skip first line
				}
				String chr = tokens[idx[2]];
				if (chr.startsWith("chr")) chr = chr.substring(3);
				transcriptAnnotTable.put(tokens[idx[0]], new String[] { tokens[idx[1]], chr, tokens[idx[3]], tokens[idx[4]].trim(), tokens[idx[5]].trim(), tokens[idx[6]],
					tokens[idx[7]], tokens[idx[8]]});
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static final void loadSNPAnnot(String filename) {
		CsvReader reader = null;
		int lineNo = 0;
		try {
			reader = new CsvReader(filename);
			reader.setTrimWhitespace(true);
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) continue; // Skip first line
				String
					marker = tokens[0],
					chr = tokens[1],
					pos = String.valueOf((int) Double.parseDouble(tokens[2])),
					rs = tokens[3],
					//effect = tokens[4],
					//noneffect = tokens[5],
					//ref_eaf = tokens[6],
					//obs_eaf = tokens[7],
					//obs_var = tokens[8],
					//exp_var = tokens[9],
					//ratio = tokens[10],
					Al1 = tokens[11],
					Al2 = tokens[12],
					Freq1 = tokens[13],
					MAF = tokens[14],
					//AvgCall = tokens[15],
					Rsq = tokens[16];
					//Genotyped = tokens[17],
					//LooRsq = tokens[18],
					//EmpR = tokens[19],
					//EmpRsq = tokens[20],
					//Dose1 = tokens[21],
					//Dose2 = tokens[22];

				snpAnnotTable.put(marker, new String[] { rs, chr, pos, Al1, Al2, Freq1, MAF, Rsq });
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static final int countNumLines(String filename) throws IOException {
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

	public static final AnnotatedEQTLResult loadAndSortFile(String filename) {
		int totalNumLines = 0;
		try {
			totalNumLines = countNumLines(filename) - 1;
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		if (totalNumLines <= 0)
			throw new RuntimeException();
		System.out.println("There are "+totalNumLines+" lines.");
		AnnotatedEQTLResult result = new AnnotatedEQTLResult();
		result.markers = new String[totalNumLines];
		result.probeset_ids = new int[totalNumLines];
		result.results = new double[4][totalNumLines];

		CsvReader reader = null;
		int lineNo = 0;
		try {
			reader = new CsvReader(filename);
			reader.setTrimWhitespace(true);
			while(reader.readRecord()) {
				lineNo++;
				if (lineNo == 1) continue; // Skip first line
				String[] tokens = reader.getValues();
				int ln = lineNo-2;
				result.markers[ln] = tokens[0];
				result.probeset_ids[ln] = Integer.parseInt(tokens[1]);
				result.results[0][ln] = Double.parseDouble(tokens[2]);
				result.results[1][ln] = Double.parseDouble(tokens[3]);
				result.results[2][ln] = Double.parseDouble(tokens[4]);
				result.results[3][ln] = Double.parseDouble(tokens[5]);
			}
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
		radixSort(result);
		System.out.println("Sorted.");
		return result;
	}

	public static final void radixSort(AnnotatedEQTLResult result)
	{
		// N=5*10^8: Radix = 9919, Quick = 16971 (ms).
		double[] data = result.results[3]; // Sort by P-value
		int
			n = data.length,
			numBins = (int) Math.max(4,Math.min(10, Math.round(Math.log(n) / (2 * Math.log(2)))));
		long[]
			a = new long[n],
			b = new long[n];
		int[]
			idx_a = new int[n],
			idx_b = new int[n];
		for (int i = 0; i < n; i++)
		{
			a[i] = Double.doubleToLongBits(data[i]);
			idx_a[i] = i;
		}
		for (long mask = ~(-1L << numBins), rshift = 0L; mask != 0; mask <<= numBins, rshift += numBins)
		{
			int[] cntarray = new int[1 << numBins];
			for (int p = 0; p < n; ++p)
				++cntarray[(int) ((a[p] & mask) >>> rshift)];
			for (int i = 1; i < cntarray.length; ++i)
				cntarray[i] += cntarray[i-1];
			for (int p = n-1; p >= 0; --p)
			{
				int key = (int) ((a[p] & mask) >>> rshift);
				--cntarray[key];
				b[cntarray[key]] = a[p];
				idx_b[cntarray[key]] = idx_a[p];
			}
			long[] temp = b; b = a; a = temp;
			int[] temp_idx = idx_b; idx_b = idx_a; idx_a = temp_idx;
		}

		int numNegs = 0;
		// Negatives are always placed at the end of the array in the reverse order.
		// So, scan to find the point
		if (a[n - 1] < 0)
			for (int i = n - 1; i >= 0 && a[i] < 0; i--) {
				data[numNegs] = Double.longBitsToDouble(a[i]);
				idx_b[numNegs++] = idx_a[i];
			}
		if (numNegs > 0)
			for (int i = numNegs; i < n; i++) {
				data[i] = Double.longBitsToDouble(a[i - numNegs]);
				idx_b[i] = idx_a[i - numNegs];
			}
		else
			for (int i = 0; i < n; i++) {
				data[i] = Double.longBitsToDouble(a[i]);
				idx_b[i] = idx_a[i];
			}
		for (int i = 0; i < n; i++) {
			String marker = result.markers[i];
			int psid = result.probeset_ids[i];
			double
				x1 = result.results[0][i],
				x2 = result.results[1][i],
				x3 = result.results[2][i];
			int j = i;
			do {
				int k = idx_b[j];
				idx_b[j] = j;
				if (k == i) break;
				result.markers[j] = result.markers[k];
				result.probeset_ids[j] = result.probeset_ids[k];
				result.results[0][j] = result.results[0][k];
				result.results[1][j] = result.results[1][k];
				result.results[2][j] = result.results[2][k];
				j = k;
			} while (true);
			result.markers[j] = marker;
			result.probeset_ids[j] = psid;
			result.results[0][j] = x1;
			result.results[1][j] = x2;
			result.results[2][j] = x3;
		}
	}

	public static final void loadAndAnnotateFile(String inputFn, String outputFn) {
		CsvReader reader = null;
		long lineNo = 0;
		PrintWriter writer = null;
		String[] tokens = null;

		try {
			writer = new PrintWriter(outputFn);
			reader = new CsvReader(inputFn);
			reader.setTrimWhitespace(true);
			while(reader.readRecord()) {
				tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) {
					writer.print("Marker,Rs_ID,SNP_Chr,SNP_Pos_hg19,SNP_Fx_Allele,SNP_Non_Fx_Allele,SNP_Fx_Allele_Freq,SNP_MAF,SNP_Imputation_RSq,ProbesetID,Transcript_Cluster_ID,Transcript_Chr,Transcript_Strand,Transcript_Start_hg19,Transcript_Stop_hg19,Transcript_GeneSymbol,Transcript_EntrezGeneID,Transcript_AccessionID,Is_Cis");
					for (int i = 2; i < tokens.length; i++)
						writer.print("," + tokens[i]);
					writer.println();
					writer.flush();
					continue; // Skip first line
				}
				String[]
					snpAnnot = snpAnnotTable.get(tokens[0]),
					transcriptAnnot = transcriptAnnotTable.get(tokens[1]);
				writer.print(tokens[0]);
				for (int i = 0; i < snpAnnot.length; i++)
					writer.print("," + snpAnnot[i]);
				writer.print("," + tokens[1]);
				for (int i = 0; i < transcriptAnnot.length; i++)
					writer.print("," + transcriptAnnot[i]);
				// Determine cis-ness
				int isCis = snpAnnot[1].equals(transcriptAnnot[1]) &&
					(Math.abs(Integer.parseInt(snpAnnot[2]) - Integer.parseInt(transcriptAnnot[3])) < cisLimit)
					? 1 : 0;
				writer.print("," + isCis);
				for (int i = 2; i < tokens.length; i++)
					writer.print("," + tokens[i]);
				writer.println();
				writer.flush();
				if (lineNo % 1000000 == 0) System.out.println(lineNo);
			}
			System.out.println(lineNo + " lines were read.");
			writer.close(); writer = null;
			reader.close(); reader = null;
		} catch (Exception e) {
			System.out.println("Line " + lineNo);
			if (tokens != null) {
				for (int i = 0; i < tokens.length; i++) System.out.print(tokens[i] + ",");
				System.out.println();
			}
			e.printStackTrace();
		}
	}

	public static final void loadSortAndAnnotateFile(String inputFn, String outputFn) {
		AnnotatedEQTLResult result = loadAndSortFile(inputFn);
		PrintWriter writer = null;

		try {
			writer = new PrintWriter(outputFn);
			writer.print("Marker,Rs_ID,SNP_Chr,SNP_Pos_hg19,SNP_Fx_Allele,SNP_Non_Fx_Allele,SNP_Fx_Allele_Freq,SNP_MAF,SNP_Imputation_RSq,ProbesetID,Transcript_Cluster_ID,Transcript_Chr,Transcript_Strand,Transcript_Start_hg19,Transcript_Stop_hg19,Transcript_GeneSymbol,Transcript_EntrezGeneID,Transcript_AccessionID,Is_Cis");
			writer.flush();
			for (int i = 0; i < result.markers.length; i++) {
				String
					markerID = result.markers[i],
					transcriptID = String.valueOf(result.probeset_ids[i]),
					snpAnnot[] = snpAnnotTable.get(markerID),
					transcriptAnnot[] = transcriptAnnotTable.get(transcriptID);
				writer.print(result.markers[i]);
				for (int j = 0; j < snpAnnot.length; j++)
					writer.print("," + snpAnnot[j]);
				writer.print("," + transcriptID);
				for (int j = 0; j < transcriptAnnot.length; j++)
					writer.print("," + transcriptAnnot[j]);
				// Determine cis-ness
				int isCis = snpAnnot[1].equals(transcriptAnnot[1]) &&
					(Math.abs(Integer.parseInt(snpAnnot[2]) - Integer.parseInt(transcriptAnnot[3])) < cisLimit)
					? 1 : 0;
				writer.print("," + isCis);
				for (int j = 0; j < 4; j++)
					writer.print("," + result.results[j][i]);
				writer.println();
				writer.flush();
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
			loadTranscriptAnnot(args[1]);
			time2 = System.currentTimeMillis();
			System.out.println("Transcript annotation read " + (time2 - time1) + " ms");
			time1 = System.currentTimeMillis();
			loadSNPAnnot(args[0]);
			time2 = System.currentTimeMillis();
			System.out.println("SNP annotation read " + (time2 - time1) + " ms");
			time1 = System.currentTimeMillis();
			if (args.length < 4) {
				for (int i = 1; i <= 23; i++) {
					String
						fn = "chr" + (i == 23 ? "X" : String.valueOf(i)) + "-exon-result.txt",
						outfn = "chr" + (i == 23 ? "X" : String.valueOf(i)) + "-exon-result-annot.txt";
					loadAndAnnotateFile(fn, outfn);
				}
			} else {
				loadAndAnnotateFile(args[2], args[3]);
			}
			time2 = System.currentTimeMillis();
			System.out.println("Annotation time " + (time2 - time1) + " ms");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

class AnnotatedEQTLResult {
	public String[] markers;
	public int[] probeset_ids;
	public double[][] results;
}
