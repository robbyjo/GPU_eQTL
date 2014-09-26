package gov.nih.tools.old;

import java.io.FileReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.StringTokenizer;

/**
 * @author Roby Joehanes
 */
public class JoinLines {
	public static final void main(String[] args) {
		String inputfile = args[0];
		int dotPos = inputfile.lastIndexOf('.');
		String
			prefix = inputfile.substring(0, dotPos),
			suffix = inputfile.substring(dotPos),
			outputfile = prefix + "-joined" + suffix,
			curLine = null;
		System.out.println("Input file: " + inputfile);
		System.out.println("Output file: " + outputfile);
		LineNumberReader reader = null;
		PrintWriter writer = null;
		int lineNo = 0;
		try {
			reader = new LineNumberReader(new FileReader(inputfile));
			writer = new PrintWriter(outputfile);
			do {
				curLine = reader.readLine();
				if (curLine == null)
					break;
				lineNo++;
				if (lineNo % 10000 == 0) System.out.println(lineNo);
				StringTokenizer tok = new StringTokenizer(curLine);
				// Write first four tokens
				writer.print(tok.nextToken() + " " + tok.nextToken() + " " + tok.nextToken() + " " + tok.nextToken() + " ");
				do {
					writer.print(tok.nextToken() + tok.nextToken());
					if (!tok.hasMoreTokens()) break;
					writer.print(" ");
				} while (true);
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
