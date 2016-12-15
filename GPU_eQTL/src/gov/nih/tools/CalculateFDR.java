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
import java.io.PrintWriter;

import static java.lang.Math.log10;
import static java.lang.Math.min;

/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class CalculateFDR {
	static double log10n = log10(39315185) + log10(283805);
	static final int one_b = 1000000000;
	static double[][] pval = new double[5][one_b];
	static long maxLineNo;

	public static final double get(long offset)
	{
		int
			offset1 = (int) (offset / one_b),
			offset2 = (int) (offset % one_b);
		return pval[offset1][offset2];
	}

	public static final void set(long offset, double v)
	{
		int
			offset1 = (int) (offset / one_b),
			offset2 = (int) (offset % one_b);
		pval[offset1][offset2] = v;
	}

	public static void main(String[] args)
	{
		LineNumberReader reader = null;
		long lineNo = 0;
		System.out.println("Usage: CalculateFDR numAllPossiblePairs filename.txt");
		int colIdx = -1;
		log10n = log10(Double.parseDouble(args[0]));
		String inputfile = args[1];
		int dotPos = inputfile.lastIndexOf('.');
		String
			prefix = inputfile.substring(0, dotPos),
			suffix = inputfile.substring(dotPos),
			outputfile = prefix + "-fdr" + suffix,
			curLine = null,
			header = null,
			tok[] = null;
		try {
			reader = new LineNumberReader(new FileReader(inputfile));
			do {
				curLine = reader.readLine();
				if (curLine == null)
					break;
				tok = curLine.split(",");
				if (header == null) {
					header = curLine;
					for (int i = 0; i < tok.length; i++)
						if ("log10P".equalsIgnoreCase(tok[i])) {
							colIdx = i;
							break;
						}
					System.out.println(tok[colIdx]);
					continue;
				}
				double curPval = Double.parseDouble(tok[colIdx]);
				set(lineNo, min(0, curPval + log10n - log10(lineNo)));
				lineNo++;
				if (lineNo % 10000000 == 0) System.out.println(lineNo/1000000);
			} while (true);
			System.out.println(lineNo + " lines were read.");
			reader.close(); reader = null;

			maxLineNo = lineNo;
			while(lineNo > 0) {
				lineNo--;
				set(lineNo, min(get(lineNo), get(lineNo + 1)));
			}

			reader = new LineNumberReader(new FileReader(inputfile));
			PrintWriter writer = new PrintWriter(outputfile);
			reader.readLine();
			writer.println(header + ",log10FDR");
			lineNo = 0;
			do {
				curLine = reader.readLine();
				if (curLine == null)
					break;
				writer.println(curLine + "," + get(lineNo));
				lineNo++;
			} while (true);
			reader.close(); reader = null;
			writer.flush(); writer.close();
			writer = null;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
