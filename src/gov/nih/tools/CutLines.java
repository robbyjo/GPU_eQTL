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
package gov.nih.tools;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;

/**
 * @author Roby Joehanes
 */
public class CutLines {
	public static final void main(String[] args) {
		String inputFile = args[1];
		long cutoff = Long.parseLong(args[0]);
		if (cutoff <= 0) throw new RuntimeException();
		int
			dotPos = inputFile.lastIndexOf('.'),
			fileCount = 0;
		String
			prefix = dotPos >= 0 ? inputFile.substring(0, dotPos) : inputFile,
			suffix = dotPos >= 0 ? inputFile.substring(dotPos) : "",
			header = null;
		if (prefix.length() > 0) prefix += "-";
		LineNumberReader reader = null;
		PrintWriter writer = null;
		long lineNo = 0, lineCut = 0;
		try {
			reader = new LineNumberReader(new FileReader(inputFile));
			writer = null; //
			do {
				String line = reader.readLine();
				if (line == null)
					break;
				lineNo++;
				if (header == null) {
					header = line;
					continue;
				}
				if (writer == null || lineCut >= cutoff) {
					if (writer != null) {
						try { writer.close(); } catch (Exception e) {}
						fileCount++;
					}
					lineCut = 0;
					writer = new PrintWriter(prefix + fileCount + suffix);
					writer.println(header);
				}
				lineCut++;
				writer.println(line);
				if (lineNo % 10000000 == 0) System.out.println(lineNo/1000000);
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
