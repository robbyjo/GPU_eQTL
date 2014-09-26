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

import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.regex.Pattern;

class ResultFileReader {
	private String curLine = null;
	private File file;
	private LineNumberReader reader = null;
	private int colIdx = 0;

	public ResultFileReader(File fn, int idx) throws IOException {
		file = fn;
		colIdx = idx;
		reader = new LineNumberReader(new FileReader(file));
	}

	public void readNextLine() throws IOException
	{	curLine = reader.readLine(); }

	public String getCurrentLine()
	{	return curLine; }

	public double getCurrentValue()
	{
		if (curLine == null)
			return 1;
		String[] tokens = curLine.split(","); // $NON-NLS-1$
		return Double.parseDouble(tokens[colIdx]);
	}

	public void purge() {
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}

/**
 * Example script
 * 
 * @author Roby Joehanes
 */
public class ExternalSort {

	public static void main(String[] args)
	{
		System.out.println("External Sort (ascending), by Roby Joehanes");
		System.out.println("NOTE: Each individual files are assumed to be presorted in ascending order!");
		int colIdx = Integer.parseInt(args[0]);
		String
			inputFilePattern = args[1],
			outputFile = args[2];
		final Pattern fn_pattern = Pattern.compile(inputFilePattern.replaceAll("\\*", ".*").replaceAll("^.*" + File.separator, ""));
		File fileDir = new File(args[1]).getAbsoluteFile().getParentFile();
		File[] files = fileDir.listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return fn_pattern.matcher(name).matches();
			}
		});
		if (files == null)
			throw new RuntimeException("Not a valid directory!");
		int numfiles = files.length;
		if (numfiles == 0)
			throw new RuntimeException("Cannot find files according to the specification!");
		System.out.println("Sorting on column " + colIdx);
		System.out.println("Found " + numfiles + " files:");
		for (File file: files) {
			System.out.println(file.getAbsolutePath());
		}
		ResultFileReader[] readers = new ResultFileReader[numfiles];
		PrintWriter writer = null;
		double[] pval = new double[numfiles];
		int pmin_idx = 0;
		long lineNo = 0;
		try {
			writer = new PrintWriter(outputFile);
			for (int i = 0; i < numfiles; i++) {
				readers[i] = new ResultFileReader(files[i], colIdx);
				readers[i].readNextLine(); // Skip header
				if (i == 0)
					writer.println(readers[i].getCurrentLine());
				readers[i].readNextLine(); // Read the first line
			}

			for (int i = 0; i < numfiles; i++)
				pval[i] = readers[i].getCurrentValue();

			while (true) {
				pmin_idx = 0;
				for (int i = 1; i < numfiles; i++) {
					if (pval[i] < pval[pmin_idx]) pmin_idx = i;
				}
				if (pval[pmin_idx] == 1)
					break;
				lineNo++;
				writer.println(readers[pmin_idx].getCurrentLine());
				readers[pmin_idx].readNextLine();
				pval[pmin_idx] = readers[pmin_idx].getCurrentValue();
				if (lineNo % 1000000 == 0) System.out.println(lineNo);
			}
			writer.flush();
			writer.close(); writer = null;
			System.out.println(lineNo+" lines written");
			for (int i = 0; i < numfiles; i++) readers[i].purge();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
