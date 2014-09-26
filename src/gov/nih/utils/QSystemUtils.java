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
/*
 * Created on Nov 30, 2004
 */
package gov.nih.utils;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.lang.reflect.Field;
import java.net.URI;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.TreeSet;

import static gov.nih.utils.QStringUtils.parseVersionString;
import static gov.nih.utils.QStringUtils.sEq;
import static gov.nih.utils.QStringUtils.sLn;

/**
 * <P>A collection of system utilities.
 * 
 * <p>Some of the benchmarking code is adapted from Vladimir Roubtsov's code here:<br>
 * http://www.javaworld.com/javaworld/javatips/jw-javatip130.html
 * 
 * <P>The code is mainly for benchmarking. To do memory benchmark,
 * you can invoke one of the following:
 * <ul>
 * <li><tt>Object retval = memoryBenchmark(instance, "methodName",
 * new Class[] { ParamTypeClass1.class, ...}, new Object[]
 * {actualParam1, ...});</tt><br>Example:
 * <tt>String retval = (String) memoryBenchmark(myInstance, "myMethod",
 * new Class[] { int.class, String.class }, new Object[] { new Integer(1),
 * "duh" });</tt><br>
 * This one is for benchmarking the method <tt>String myMethod(int, String)</tt>
 * method of myInstance. So, instead of calling<br>
 * <tt>retval = myInstance.myMethod(1,"duh");</tt><br>
 * (which is not benchmarked), you can replace it with the benchmarked version.
 * For parameters with array types, you can specify like <tt>int[].class</tt>
 * in the paramTypes to indicate that the parameter type is an array of int.
 * If the method has no parameters, put null in the last two arguments.
 * </li>
 * <li><tt>Object retval = memoryBenchmark(ClassName.class, "staticMethodName",
 * new Class[] { ParamTypeClass1.class, ...}, new Object[] {actualParam1, ...});</tt><br>
 * This is similar to the first one except that this one is for static methods.
 * Example:<br>
 * <tt>String retval = (String) memoryBenchmark(MyClass.class, "myMethod",
 * new Class[] { int.class, String.class }, new Object[] { new Integer(1),
 * "duh" });</tt><br>
 * This is to benchmark the static method <tt>MyClass.myMethod(int, String)</tt>
 * </li>
 * </ul>
 * 
 * <P>For time benchmark, it's the same as memory benchmark. Just replace
 * <tt>memoryBenchmark</tt> with <tt>timeBenchmark</tt>.
 * 
 * <P>Note that in each call to <tt>memoryBenchmark</tt>, we invoke Java's
 * garbage collector (GC) aggressively so that this may be very slow for some
 * systems. This step is necessary to measure the true amount of used memory
 * without. If we don't invoke the GC, if the method runs long enough that
 * GC kicks in, it will definitely interfere with the memory measurement.
 * 
 * <P>When doing time benchmarks, we don't invoke the GC.
 * 
 * <P><strong>IMPORTANT: The benchmark result is printed through System.out!</strong>
 * <P><strong>IMPORTANT: The benchmark is not suitable for multi-threaded
 * benchmarking yet since if other threads are running, they are potentially
 * interfering with the benchmark results.</strong>
 * <P><strong>IMPORTANT: If you refactor this class, pay attention to getProgramPath
 * since it's package dependent!</strong>
 * 
 * @author Roby Joehanes
 */
public final class QSystemUtils
{
	/** Base memory **/
	public static final long mBaseMem = usedMemoryAfterGC();

	public static final String
		sNotEnoughMemoryErr = "Error.NotEnoughMemory"; //$NON-NLS-1$

	private static final int
		kNumTimesGCCall = 500,
		kNumTimesAggresiveGCCall = 4;
	// Delimiter used by jar URI
	private static final String
		sJarDelimiter = "!" + File.separator, //$NON-NLS-1$
		sFilePosDelimiter = "file:" + File.separator, //$NON-NLS-1$
		sSpaceDelimiter = "%20"; //$NON-NLS-1$

	private static final String
		sJava3DClassName = "javax.media.j3d.Canvas3D", //$NON-NLS-1$
		sJavaVersion = "java.version"; //$NON-NLS-1$

	private static SimpleDateFormat mDateFormatter; // If I initialized it here, it will give problems in Linux
	// I think this is a Java bug.

	public static final String
		sDotDotSlash = ".." + File.separator, //$NON-NLS-1$
		sLibDirectory = File.separator + "lib" + File.separator, //$NON-NLS-1$
		sBinDirectory = File.separator + "bin" + File.separator, //$NON-NLS-1$
		sx64 = "x64"; //$NON-NLS-1$

	public static final int
		kNumCPUCores = Runtime.getRuntime().availableProcessors(),
		kCurrentJavaVersion[] = QStringUtils.parseVersionString(System.getProperty(sJavaVersion));

	/**
	 * Original System.out
	 */
	public static final PrintStream mSystemOutput = System.out;

	/**
	 * Original System.err
	 */
	public static final PrintStream mSystemError = System.err;

	/**
	 * Redirected System.out / System.err. By default it has 256KB output buffer.
	 */
	public static final ByteArrayOutputStream mBufferedOutput = new ByteArrayOutputStream(256 * 1024);

	/**
	 * Run garbage collector aggressively (i.e. 4 times), even after
	 * the used memory stabilizes
	 */
	public static final void runGCAggressively()
	{
		// It helps to call Runtime.gc() using several method calls:
		for (int r = 0; r < kNumTimesAggresiveGCCall; ++r)
			runGC();
	}

	/**
	 * Run garbage collector until the used memory stabilizes
	 */
	public static final void runGC()
	{
		long
			usedMem1 = usedMemory(),
			usedMem2 = Long.MAX_VALUE;
		Runtime runtime = Runtime.getRuntime();
		for (int i = 0; (usedMem1 < usedMem2) && (i < kNumTimesGCCall); ++i)
		{
			runtime.runFinalization();
			runtime.gc();
			Thread.yield();

			usedMem2 = usedMem1;
			usedMem1 = usedMemory();
		}
	}

	/**
	 * Return used memory in bytes
	 * @return
	 */
	public static final long usedMemory()
	{
		Runtime runtime = Runtime.getRuntime();
		return runtime.totalMemory() - runtime.freeMemory();
	}

	/**
	 * Run garbage collector aggressively, then return the used memory in bytes.
	 * More accurate than usedMemory() since we force Java to release unused memory
	 * that's ready to be garbage collected.
	 * @return
	 */
	public static final long usedMemoryAfterGC()
	{
		runGCAggressively();
		Runtime runtime = Runtime.getRuntime();
		return runtime.totalMemory() - runtime.freeMemory();
	}

	/**
	 * Return memory statistics
	 * @param runGCFirst Runs garbage collector before getting the statistics or not.
	 * @return
	 */
	public static final QMemoryUsage getMemoryStatistics(boolean runGCFirst)
	{
		if (runGCFirst)
			QSystemUtils.runGCAggressively();
		MemoryMXBean memBean = ManagementFactory.getMemoryMXBean();
		MemoryUsage
			heapUsage = memBean.getHeapMemoryUsage(),
			nonHeapUsage = memBean.getNonHeapMemoryUsage();
		Runtime runtime = Runtime.getRuntime();
		long mem = runGCFirst ? usedMemoryAfterGC() : runtime.totalMemory() - runtime.freeMemory();
		return new QMemoryUsage(heapUsage.getUsed(), nonHeapUsage.getUsed(), mem - mBaseMem, heapUsage.getMax());
	}

	/**
	 * Get enumeration of current system. Example, Windows
	 * computers will return EPlatform.WINDOWS_x86
	 * @return
	 */
	public static final EPlatform getCurrentPlatform()
	{	return EPlatform.getCurrentPlatform(); }

	/**
	 * Get current date time. Example, if today is
	 * 23 June 2005, 18:26, it will output 200506231826
	 * @return
	 */
	public static final long getCurrentDateTime()
	{
		GregorianCalendar calendar = new GregorianCalendar();
		return calendar.get(Calendar.YEAR) * 10000000000L +
			calendar.get(Calendar.MONTH) * 100000000L +
			calendar.get(Calendar.DAY_OF_MONTH) * 1000000L +
			calendar.get(Calendar.HOUR_OF_DAY) * 10000L + 
			calendar.get(Calendar.MINUTE) * 100L +
			calendar.get(Calendar.SECOND);
	}

	/**
	 * Parse a path to an array of strings. Useful to
	 * parse classpath
	 * @param path
	 * @return
	 */
	public static final String[] parsePath(String path)
	{
		List<String> list = new ArrayList<String>();
		StringTokenizer tokenizer = new StringTokenizer(path, File.pathSeparator);
		while (tokenizer.hasMoreTokens())
			list.add(tokenizer.nextToken());
		return list.toArray(new String[list.size()]);
	}

	/**
	 * Get the directory where a particular class reside.<br>
	 * Example: If file MyClass.class is at /usr/local/bin/myprog/mypackage directory,
	 * it will return /usr/local/bin/myprog/mypackage
	 * @param c
	 * @return the path where that class reside
	 */
	public static final File getClassDirectoryFile(Class<?> c)
	{
		// Get the path that includes the class file
		String path = c.getResource(c.getSimpleName()+".class").getPath(); //$NON-NLS-1$
		// We're interested at the directory only. So, invoke getParentFile
		return new File(path).getParentFile();
	}

	/**
	 * Get the directory where a particular class reside.<br>
	 * Example: If file MyClass.class is at /usr/local/bin/myprog/mypackage directory,
	 * it will return /usr/local/bin/myprog/mypackage
	 * @param c
	 * @return the path where that class reside
	 */
	public static final String getClassDirectory(Class<?> c)
	{
		File f = getClassDirectoryFile(c);
		String fString = f.toString();
		int filePos = fString.indexOf(sFilePosDelimiter);
		if (filePos >= 0) 
			return fString.substring(filePos + sFilePosDelimiter.length());
			
		return f.getAbsolutePath();
	}

	/**
	 * Get the root directory where a particular class reside.<br>
	 * Example: If file MyClass.class is at /usr/local/bin/myprog/mypackage directory,
	 * it will return /usr/local/bin/myprog/
	 * @param c
	 * @return the path where that class reside
	 */
	public static final File getClassRootDirectoryFile(Class<?> c)
	{
		// Get the path that includes the class file
		String path = c.getResource(c.getSimpleName()+".class").getPath(); //$NON-NLS-1$
		// We're interested at the directory only. So, invoke getParentFile
		File dir = new File(path).getParentFile();
		Package p = c.getPackage();
		// we have package information, so we may need to append "../" as needed
		if (p == null) // this shouldn't happen
			return new File ("."); //$NON-NLS-1$
		StringTokenizer tokenizer = new StringTokenizer(p.getName(), "."); //$NON-NLS-1$
		int numTokens = tokenizer.countTokens();
		for (int i = 0; i < numTokens; i++)
			dir = dir.getParentFile();
		return dir;
	}

	/**
	 * Get the root directory where a particular class reside.<br>
	 * Example: If file MyClass.class is at /usr/local/bin/myprog/mypackage directory,
	 * it will return /usr/local/bin/myprog/
	 * @param c
	 * @return the path where that class reside
	 */
	public static final String getClassRootDirectory(Class<?> c)
	{
		File f = getClassRootDirectoryFile(c);
		String fString = f.getAbsolutePath();
		int filePos = fString.indexOf(sFilePosDelimiter);
		if (filePos >= 0)
		{
			int offset = EPlatform.isWindows() ? 0 : -1; // Windows messes up with directory naming
			return fString.substring(filePos + sFilePosDelimiter.length() + offset);
		}
		String path = f.getAbsolutePath();
		if (EPlatform.isWindows())
			path = unmangleSpace(path);
		return path;
	}

	/**
	 * Get current directory
	 * @return
	 */
	public static final String getCurrentDirectory()
	{
		return System.getProperty("user.dir"); //$NON-NLS-1$
	}

	/**
	 * Get user's home directory
	 * @return
	 */
	public static final String getHomeDirectory()
	{
		return System.getProperty("user.home"); //$NON-NLS-1$
	}

	/**
	 * Return the path of this program
	 * @return
	 */
	public static final String getProgramPath()
	{	return getProgramPath(QSystemUtils.class); }

	/**
	 * Return the path of a program that contains class c
	 * @param c
	 * @return
	 */
	public static final String getProgramPath(Class<?> c)
	{
		String filePath = getClassRootDirectory(c);
		//System.out.println("Orig filePath = " + filePath);

		if (!filePath.endsWith(File.separator)) // make sure it ends with a separator
			filePath += File.separator;
		//System.out.println("FilePath = " + filePath);

		int jarDelimPos = filePath.indexOf(sJarDelimiter);
		if (jarDelimPos >= 0)
		{
			// If in jar, assume in /lib directory.
			filePath = filePath.substring(0, filePath.lastIndexOf(File.separatorChar, jarDelimPos));
			filePath = new File(filePath).getParent();
			//System.out.println("is jar path, FilePath = " + filePath);
		}
		// If it is in /src or /bin directory, then
		// the root directory must be one level above it
		String ll = filePath.toLowerCase();
		if (ll.endsWith(sBinDirectory))
			filePath = new File(filePath).getParent();

		//System.out.println("FilePath = " + filePath);
		// I don't know why Windows need another set of space unmangling. -- RJ
		if (EPlatform.isWindows())
			filePath = unmangleSpace(filePath);
		//System.out.println("FilePath = " + filePath);

		return filePath;
	}

	/**
	 * Is the path inside or at a jar file?
	 * @param path
	 * @return
	 */
	public static final boolean isJarPath(String path)
	{	return path.indexOf(sJarDelimiter) >= 0; }

	/**
	 * Essentially a hack that changes spaces to %20 to avoid Windows bugs
	 * @param path
	 * @return
	 */
	public static final String unmangleSpace(String path)
	{
		int spaceDelimiter = path.indexOf(sSpaceDelimiter);
		return spaceDelimiter < 0 ? path :
			path.replaceAll(sSpaceDelimiter, " "); //$NON-NLS-1$
	}
	/**
	 * If the path is inside a jar file, this routine will split the
	 * path into two: The jar file's path and the path inside the jar
	 * file.
	 * @param path
	 * @return array of String with two elements: The first one is the
	 * jar file, the second one is the path inside. IF the path is NOT
	 * a jar path, it will return null 
	 */
	public static final String[] extractJarPath(String path)
	{
		int jarDelimiterPos = path.indexOf(sJarDelimiter);
		return jarDelimiterPos >= 0 ?
				new String[] { unmangleSpace(path.substring(0, jarDelimiterPos)),
				unmangleSpace(path.substring(jarDelimiterPos + sJarDelimiter.length()))}
		: null;
	}

	/**
	 * Convert backslash-delimited path into slash-delimited one
	 * @param pathhttp://java.sun.com/docs/books/tutorial/uiswing/layout/gridbag.html
	 * @return
	 */
	public static final String convertPath(String path)
	{
		URI
			curDir = new File(".").toURI(), //$NON-NLS-1$
			uri = new File(path).toURI();
		return curDir.relativize(uri).getPath();
	}
	
	/**
	 * Redirect System.out and System.err to default buffer
	 */
	public static final void redirectSystemOutErrDefault()
	{
		PrintStream ps = new PrintStream(mBufferedOutput);
		System.setOut(ps);
		System.setErr(ps);
	}

	/**
	 * Restore System.out back to their original state
	 */
	public static final void restoreSystemOut()
	{	System.setOut(mSystemOutput); }

	/**
	 * Restore System.err back to their original state
	 */
	public static final void restoreSystemErr()
	{	System.setErr(mSystemError); }

	/**
	 * Dump sytem properties
	 * @return
	 */
	public static final String dumpSystemProperties()
	{
		Properties prop = System.getProperties();
		StringBuilder buf = new StringBuilder();
		// FIXME This may break in Java 7. See Java Bug ID: 6804124 and 5045147 for details.
		// http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6804124
		// http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=5045147
		// I think this is okay since keyset in properties are strings.
		TreeSet<Object> sortedSet = new TreeSet<Object>();
		sortedSet.addAll(prop.keySet());
		for (Object key: sortedSet)
		{
			Object val = prop.get(key);
			buf.append(key + sEq + val + sLn);
		}
		return buf.toString();
	}

	/**
	 * Returns 1 if v1 > v2, 0 if v1 == v2, -1 if v1 < v2;
	 * @param v1
	 * @param v2
	 * @return
	 */
	public static final int compareVersion(String v1, String v2)
	{
		int[]
			thisVersion = parseVersionString(v1),
			newVersion = parseVersionString(v2);
		return compareVersion(thisVersion, newVersion);
	}

	/**
	 * Returns 1 if v1 > v2, 0 if v1 == v2, -1 if v1 < v2;
	 * @param v1
	 * @param v2
	 * @return
	 */
	public static final int compareVersion(int[] v1, int[] v2)
	{
		int verNo = v1.length;
		if (verNo > v2.length)
			verNo = v2.length;
		for (int i = 0; i < verNo; i++)
			if (v1[i] < v2[i])
				return -1;
			else if (v1[i] > v2[i])
				return 1;
		if (v1.length < v2.length)
			return -1;
		return 0;
	}

	/**
	 * Check Java 3D installation.
	 * @return true if the computer has Java 3D installed
	 */
	public static final boolean hasJava3D()
	{
		try
		{
			Class.forName(sJava3DClassName); //$NON-NLS-1$
			return true;
		}
		catch (Exception e)
		{	return false; }
	}

	/**
	 * Forcibly adding Java's library lookup directory.
	 * <P>From: http://forums.sun.com/thread.jspa?threadID=707176
	 * @param s
	 */
	public static final void forceAddLibDir(String s)
	{
		try
		{
			Field field = ClassLoader.class.getDeclaredField("usr_paths"); //$NON-NLS-1$
			field.setAccessible(true);
			String[] paths = (String[])field.get(null);
			for (int i = 0; i < paths.length; i++) {
				if (s.equals(paths[i])) {
					return;
				}
			}
			String[] tmp = new String[paths.length+1];
			System.arraycopy(paths,0,tmp,0,paths.length);
			tmp[paths.length] = s;
			field.set(null,tmp);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Format time in HH:MM:SS.xxx format
	 * @param time Time in milliseconds
	 * @return
	 */
	public static final String toTimeString(long time)
	{
		if (mDateFormatter == null)
			mDateFormatter = new SimpleDateFormat("HH:mm:ss.SSS"); //$NON-NLS-1$
		return mDateFormatter.format(new Date(time));
	}

	/**
	 * Returns true if the load is successful.
	 * @return
	 */
//	public static final boolean loadRLibrary(String dir)
//	{
//		try {
//			System.loadLibrary(sJRI);
//			return true;
//		} catch (UnsatisfiedLinkError e)
//		{}
//		if (!dir.endsWith(File.separator))
//			dir += File.separator;
//		dir += "bin" + File.separator; //$NON-NLS-1$
//		try {
//			forceAddLibDir(dir);
//			System.loadLibrary(sJRI);
//			return true;
//		} catch (UnsatisfiedLinkError e)
//		{	}
//		return false;
//	}

	private static final long numMSPerDay = 1000*60*60*24;

	/**
	 * Find the number of days between y1-m1-d1 and y2-m2-d2
	 * @param y1
	 * @param m1
	 * @param d1
	 * @param y2
	 * @param m2
	 * @param d2
	 * @return
	 */
	public static final int subtractDates(int y1, int m1, int d1, int y2, int m2, int d2)
	{
		return subtractDates(new GregorianCalendar(y1, m1, d1), new GregorianCalendar(y2, m2, d2));
	}

	/**
	 * Find the number of days between d1 and d2
	 * @param d1
	 * @param d2
	 * @return
	 */
	public static final int subtractDates(Calendar d1, Calendar d2)
	{
		if (d1.after(d2))
		{
			Calendar temp = d1; d1 = d2; d2 = temp;
		}
		long d3 = d2.getTimeInMillis() - d1.getTimeInMillis();
		return (int) (d3 / numMSPerDay);
	}

	/**
	 * Find the number of days between d1 and d2
	 * @param d1
	 * @param d2
	 * @return
	 */
	public static final int subtractDates(Date d1, Date d2)
	{
		if (d1.after(d2))
		{
			Date temp = d1; d1 = d2; d2 = temp;
		}
		long d3 = d2.getTime() - d1.getTime();
		return (int) (d3 / numMSPerDay);
	}

	public static final String[] exec(String[] cmd) throws IOException
	{
		Process p = Runtime.getRuntime().exec(cmd);
		BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
		BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		String s = null;
		StringBuilder
			stdout = new StringBuilder(),
			stderr = new StringBuilder();
		while ((s = stdInput.readLine()) != null) {
			stdout.append(s + sLn);
		}
		while ((s = stdError.readLine()) != null) {
			stderr.append(s + sLn);
		}
		return new String[] { stdout.toString(), stderr.toString() };
	}
}
