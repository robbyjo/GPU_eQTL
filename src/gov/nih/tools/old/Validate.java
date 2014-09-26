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
package gov.nih.tools.old;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import com.csvreader.CsvReader;

/**
 * @author Roby Joehanes
 */
public class Validate {

	static Set<String>
		ludeMarkerPairs = new HashSet<String>(),
		ludeRsIDPairs = new HashSet<String>(),
		mrcaMarkerPairs = new HashSet<String>(),
		mrcaRsIDPairs = new HashSet<String>(),
		mrcaMarkerEntrezPairs = new HashSet<String>(),
		mrcaRsIDEntrezPairs = new HashSet<String>(),
		ncbiEntrezPairs = new HashSet<String>(),
		ncbiPairs = new HashSet<String>();

	static interface TokenProcessor {
		void process(String[] tokens);
	}

	static class NCBIProcessor implements TokenProcessor {
		public void process(String[] tokens) {
			String rsid = tokens[2];
			Set<String> geneSymbol = new HashSet<String>(), entrez = new HashSet<String>();
			for (String geneSym: tokens[30].split("\\|")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				geneSymbol.add(geneSym);
			}
			for (String geneSym: tokens[32].split("\\|")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				geneSymbol.add(geneSym);
			}
			for (String geneSym: tokens[38].split("\\|")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				geneSymbol.add(geneSym);
			}
			for (String geneSym: tokens[29].split("\\|")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				entrez.add(geneSym);
			}
			for (String hugo: geneSymbol)
				ncbiPairs.add(rsid + "|" + hugo);
			for (String hugo: entrez)
				ncbiEntrezPairs.add(rsid + "|" + hugo);
		}
	}

	static class LudeProcessor implements TokenProcessor {
		public void process(String[] tokens) {
			String rsid = tokens[1], marker = tokens[3] + ":" + tokens[4];
			for (String hugo: tokens[13].split(",")) {
				hugo = hugo.trim();
				if (hugo.length() == 0 || hugo.startsWith("-"))
					continue;
				ludeMarkerPairs.add(marker + "|" + hugo);
				ludeRsIDPairs.add(rsid + "|" + hugo);
			}
		}
	}

	static class MRCAProcessor implements TokenProcessor {
		public void process(String[] tokens) {
			String rsid = tokens[0], marker = tokens[1] + ":" + tokens[2];
			Set<String> geneSymbol = new HashSet<String>(), entrez = new HashSet<String>();
			for (String geneSym: tokens[5].split("///")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				geneSymbol.add(geneSym);
			}
			for (String geneSym: tokens[6].split("///")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				geneSymbol.add(geneSym);
			}
			for (String geneSym: tokens[7].split("///")) {
				geneSym = geneSym.trim();
				if (geneSym.length() == 0 || geneSym.startsWith("-"))
					continue;
				entrez.add(geneSym);
			}
			for (String hugo: geneSymbol) {
				mrcaMarkerPairs.add(marker + "|" + hugo);
				mrcaRsIDPairs.add(rsid + "|" + hugo);
			}
			for (String hugo: entrez) {
				mrcaMarkerEntrezPairs.add(marker + "|" + hugo);
				mrcaRsIDEntrezPairs.add(rsid + "|" + hugo);
			}
		}
	}

	static void processAnnotation(String filename, TokenProcessor proc) {
		CsvReader reader = null;
		try {
			reader = new CsvReader(filename);
			long lineNo = 0;
			while(reader.readRecord()) {
				String[] tokens = reader.getValues();
				lineNo++;
				if (lineNo == 1) continue;
				proc.process(tokens);
				if (lineNo % 10000000 == 0) System.out.println(lineNo/1000000);
			};
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (reader != null)
				try {
					reader.close();
				} catch (Exception e) {}
		}
	}

	public static final void main(String[] args) {
		String inputFile = args[0];
		String outputFile = args[1];
		processAnnotation("/project/eqtl/validation/Lude_EQTL-hg19.txt", new LudeProcessor());
		processAnnotation("/project/eqtl/validation/MRCA_5percentFDR-hg19.txt", new MRCAProcessor());
		processAnnotation("/project/eqtl/validation/ncbi-lite.txt", new NCBIProcessor());
		CsvReader reader = null;
		PrintWriter writer = null;
		long lineNo = 0, valLude = 0, valMrca = 0, valNCBI = 0, valAny = 0;
		try {
			reader = new CsvReader(inputFile);
			writer = new PrintWriter(outputFile);
			while(reader.readRecord()) {
				String line = reader.getRawRecord();
				lineNo++;
				if (lineNo == 1) {
					writer.println(line+",ValidatedLude,ValidatedMRCA,ValidatedNCBI,ValidatedAny");
					continue;
				}
				String[] tokens = reader.getValues();
				String marker = tokens[2] + ":" + tokens[3], rsid = tokens[1];
				Set<String> geneSymbol = new HashSet<String>(), entrez = new HashSet<String>();
				for (String geneSym: tokens[15].split("\\|")) {
					geneSym = geneSym.trim();
					if (geneSym.length() == 0 || geneSym.startsWith("-"))
						continue;
					geneSymbol.add(geneSym);
				}
				for (String geneSym: tokens[16].split("\\|")) {
					geneSym = geneSym.trim();
					if (geneSym.length() == 0 || geneSym.startsWith("-"))
						continue;
					entrez.add(geneSym);
				}
				int
					lude = 0,
					mrca = 0,
					ncbi = 0;
				for (String geneSym : geneSymbol) {
					String pair = marker + "|" + geneSym;
					if (ludeMarkerPairs.contains(pair)) lude |= 1;
					if (mrcaMarkerPairs.contains(pair)) mrca |= 1;
					pair = rsid + "|" + geneSym;
					if (ludeRsIDPairs.contains(pair)) lude |= 1;
					if (mrcaRsIDPairs.contains(pair)) mrca |= 1;
					if (ncbiPairs.contains(pair)) ncbi |= 1;
				}
				for (String geneSym : entrez) {
					String pair = marker + "|" + geneSym;
					if (mrcaMarkerEntrezPairs.contains(pair)) mrca |= 1;
					pair = rsid + "|" + geneSym;
					if (mrcaRsIDEntrezPairs.contains(pair)) mrca |= 1;
					if (ncbiEntrezPairs.contains(pair)) ncbi |= 1;
				}
				valLude += lude; valMrca += mrca; valNCBI += ncbi;
				valAny = valAny + (lude + mrca + ncbi > 0 ? 1 : 0);
				writer.println(line + "," + lude + "," + mrca + "," + ncbi + "," + (lude + mrca + ncbi));
				if (lineNo % 10000000 == 0) System.out.println(lineNo/1000000);
			};
			System.out.println(lineNo + " lines were read.");
			System.out.println("Validated on Lude's = " + valLude);
			System.out.println("Validated on MRCA = " + valMrca);
			System.out.println("Validated on NCBI = " + valNCBI);
			System.out.println("Validated on at least one = " + valAny);
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
