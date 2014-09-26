/*
 * Roby Joehanes
 * 
 * Copyright 2007 Roby Joehanes
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
package gov.nih.eqtl.io;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import gov.nih.utils.QStringUtils;

/**
 * Configuration saver / loader using INI-file format
 * INI-format is as follows:
 * 
 * "#" denotes the beginning of a comment.
 * The key/value pair is written as:
 * key = value
 * 
 * @author Roby Joehanes
 */
public class QIniLoader {
	private static final char
		cCommentSymbol = '#',
		cEqualSign = '=';

	private static final String sEqualSign = " = "; //$NON-NLS-1$

	public static final Map<String, String> load(Reader reader) throws IOException
	{
		HashMap<String, String> table = new HashMap<String, String>(); // The translation table
		LineNumberReader lineReader = new LineNumberReader(reader); // Convert to line-based reader
		String s;
		do
		{
			s = lineReader.readLine();
			if (s == null) // If the line is null, we're at the end of the file
				break;
			s = s.trim();
			int strLength = s.length();
			// If the line is empty or is a comment, ignore it. 
			if (strLength == 0 || s.charAt(0) == cCommentSymbol)
				continue;
			// Check the existence of an equal sign "=" 
			int equalPos = s.indexOf(cEqualSign);
			// If there is an equal sign, split the string into two:
			// The one before the equal sign as the key and
			// the one after the equal sign as the value
			if (equalPos == -1)
				break;
			String
				key = s.substring(0, equalPos).trim().toLowerCase(Locale.ENGLISH),  // Turkish Windows fix
				value = s.substring(equalPos + 1).trim();
			// Then put the key/value pair into the table
			table.put(key, QStringUtils.unescape(value));
		} while (true); // Repeat the whole stuff until the end of file
		// Return config table
		return table;
	}

	public static final void save(Map<String, String> cfg, Writer w) // throws IOException
	{
		PrintWriter writer = new PrintWriter(w);
		// Sort the available keys in the configuration
		Set<String> sortedKeys = new TreeSet<String>(cfg.keySet());
		for (String key: sortedKeys)
			writer.println(key.toLowerCase(Locale.ENGLISH) + sEqualSign + QStringUtils.escape(cfg.get(key)));  // Turkish Windows fix
	}
}
