package gov.nih.tools;

import java.io.FileReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.StringTokenizer;

/**
 * @author Roby Joehanes
 */
public class ProfileLines {
	public static final void main(String[] args) {
		String inputfile = args[0];
		int dotPos = inputfile.lastIndexOf('.');
		String
			prefix = inputfile.substring(0, dotPos),
			suffix = inputfile.substring(dotPos),
			outputfile = prefix + "-profile" + suffix,
			curLine = null;
		System.out.println("Input file: " + inputfile);
		System.out.println("Output file: " + outputfile);
		LineNumberReader reader = null;
		PrintWriter writer = null;
		int lineNo = 0;
		HashMap<String, Integer> genoMap = new HashMap<String, Integer>();
		HashMap<String, String> genoEquivMap = new HashMap<String, String>();
		int[] alleleCounts = new int[5];
		String[] alleleString = new String[] { "0", "A", "C", "G", "T" };

		genoEquivMap.put("AA", "AA");
		genoEquivMap.put("AC", "AC");
		genoEquivMap.put("AG", "AG");
		genoEquivMap.put("AT", "AT");
		genoEquivMap.put("A0", "A0");
		genoEquivMap.put("CA", "AC");
		genoEquivMap.put("CC", "CC");
		genoEquivMap.put("CG", "CG");
		genoEquivMap.put("CT", "CT");
		genoEquivMap.put("C0", "C0");
		genoEquivMap.put("GA", "AG");
		genoEquivMap.put("GC", "CG");
		genoEquivMap.put("GG", "GG");
		genoEquivMap.put("GT", "GT");
		genoEquivMap.put("G0", "G0");
		genoEquivMap.put("TA", "AT");
		genoEquivMap.put("TC", "CT");
		genoEquivMap.put("TG", "GT");
		genoEquivMap.put("TT", "TT");
		genoEquivMap.put("T0", "T0");
		genoEquivMap.put("0A", "A0");
		genoEquivMap.put("0C", "C0");
		genoEquivMap.put("0G", "G0");
		genoEquivMap.put("0T", "T0");
		genoEquivMap.put("00", "00");

		try {
			reader = new LineNumberReader(new FileReader(inputfile));
			writer = new PrintWriter(outputfile);
			writer.print("SNP ID,Major Allele,Minor Allele,Major Allele Count, Minor Allele Count, Missing Allele Count,Major Homozygous Genotype,Minor Homozygous Genotype,Heterozygous Genotype,");
			writer.print("Major Homozygous Genotype Count,Minor Homozygous Genotype Count,Heterozygous Genotype Count,Completely Missing Genotype Count,MAF");
			writer.println();
			do {
				curLine = reader.readLine();
				if (curLine == null)
					break;
				lineNo++;
				if (lineNo == 1) continue; // Skip the first line
				if (lineNo % 10000 == 0) System.out.println(lineNo);
				StringTokenizer tok = new StringTokenizer(curLine, ",");
				tok.nextToken();
				String rsid = tok.nextToken();
				tok.nextToken();
				tok.nextToken();
				genoMap.clear();
				Arrays.fill(alleleCounts, 0);

				do {
					String curToken = genoEquivMap.get(tok.nextToken());
					if (!genoMap.containsKey(curToken)) genoMap.put(curToken, 0);
					genoMap.put(curToken, genoMap.get(curToken) + 1);
					if (!tok.hasMoreTokens()) break;
					for (int i = 0; i < curToken.length(); i++) {
						char ch = curToken.charAt(i);
						switch (ch) {
							case 'A': alleleCounts[1]++; break;
							case 'C': alleleCounts[2]++; break;
							case 'G': alleleCounts[3]++; break;
							case 'T': alleleCounts[4]++; break;
							case '0':
							default:
								alleleCounts[0]++; break;
						}
					}
				} while (true);

				int maj_idx = 0, min_idx = 0, maj_allele_ct = 0, min_allele_ct = Integer.MAX_VALUE;
				for (int i = 1; i < alleleCounts.length; i++)
					if (alleleCounts[i] > 0) {
						if (alleleCounts[i] > maj_allele_ct) { maj_idx = i; maj_allele_ct = alleleCounts[i]; }
						if (alleleCounts[i] <= min_allele_ct) { min_idx = i; min_allele_ct = alleleCounts[i]; }
					}
				if (maj_idx == min_idx) min_idx = maj_idx == 1 ? 2: 1;

				String
					maj_allele = alleleString[maj_idx],
					min_allele = alleleString[min_idx],
					maj_homozygous = maj_allele + maj_allele,
					min_homozygous = min_allele + min_allele,
					heterozygous = genoEquivMap.get(maj_allele + min_allele);
				writer.print(rsid + "," + maj_allele + "," + min_allele + "," + alleleCounts[maj_idx] + "," + alleleCounts[min_idx] + "," + alleleCounts[0] + ",");
				writer.print(maj_homozygous + "," + min_homozygous + "," + heterozygous + ",");
				Integer ct = genoMap.get(maj_homozygous);
				writer.print(ct == null ? "0," : ct + ",");
				ct = genoMap.get(min_homozygous);
				writer.print(ct == null ? "0," : ct + ",");
				ct = genoMap.get(heterozygous);
				writer.print(ct == null ? "0," : ct + ",");
				ct = genoMap.get("00");
				writer.print(ct == null ? "0," : ct + ",");
				/*
				genoMap.remove(maj_homozygous);
				genoMap.remove(min_homozygous);
				genoMap.remove(heterozygous);
				genoMap.remove("00");
				if (genoMap.size() > 0) {
					for (String key: genoMap.keySet()) {
						int count = genoMap.get(key);
						writer.print(key+"="+count+";");
					}
				}
				//*/
				writer.print(alleleCounts[min_idx] *1.0 / (alleleCounts[maj_idx] + alleleCounts[min_idx]));
				writer.println();
			} while (true);
			System.out.println(lineNo + " lines were read.");
			writer.close(); writer = null;
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (reader != null)
				try {
					reader.close();
				} catch (Exception e) {}
			if (writer != null)
				try {
					writer.close();
				} catch (Exception e) {}
		}
	}
}
