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
package gov.nih.utils;

import static gov.nih.utils.QStringUtils.parseVersionString;

import java.util.Locale;

/**
 * Enumeration to indicate system platform
 * 
 * @author Roby Joehanes
 */
public enum EPlatform
{
	WINDOWS_9X, WINDOWS_2000, WINDOWS_ME, WINDOWS_XP, WINDOWS_VISTA, WINDOWS_7, WINDOWS_CE, WINDOWS_OTHER,
	LINUX_x86, LINUX_x86_64, LINUX_PPC, LINUX_SPARC, LINUX_OTHER,
	SOLARIS_SPARC, SOLARIS_x86, SOLARIS_OTHER,
	MAC_OSX_PPC, MAC_OSX_x86, MAC_OSX_x86_64, MAC_CLASSIC,
	AIX, HPUX_PA_RISC, HPUX_IA64, HPUX_OTHER, DEC_ALPHA,
	OS2, BSD_x86, BSD_OTHER, OTHER;

	private static EPlatform mCurrentPlatform = null;
	private String mDescription;

	/**
	 * Automated system detection
	 * @return the enumerated platform
	 */
	public static EPlatform getCurrentPlatform()
	{
		if (mCurrentPlatform != null)
			return mCurrentPlatform;
		String osName = System.getProperty("os.name").toLowerCase(); //$NON-NLS-1$
		String arch = System.getProperty("os.arch").toLowerCase(); //$NON-NLS-1$
		mCurrentPlatform = OTHER;
		if (osName.startsWith("windows")) //$NON-NLS-1$
		{
			if (osName.indexOf("9") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_9X;
			else if (osName.indexOf("windows me") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_ME;
			else if (osName.indexOf("vista") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_VISTA;
			else if (osName.indexOf("2000") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_2000;
			else if (osName.indexOf("xp") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_XP;
			else if (osName.indexOf("windows 7") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_7;
			else if (arch.indexOf("arm") >= 0) //$NON-NLS-1$
				mCurrentPlatform = WINDOWS_CE;
			else
				mCurrentPlatform = WINDOWS_OTHER;
		}
		else if (osName.startsWith("linux")) //$NON-NLS-1$
		{
			if (arch.equals("x86") || (arch.startsWith("i") && arch.endsWith("86"))) //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				mCurrentPlatform = LINUX_x86;
			else if (arch.equals("x86_64")) //$NON-NLS-1$
				mCurrentPlatform = LINUX_x86_64;
			else if (arch.startsWith("sparc")) //$NON-NLS-1$
				mCurrentPlatform = LINUX_SPARC;
			else if (arch.startsWith("ppc") || arch.startsWith("power")) //$NON-NLS-1$ //$NON-NLS-2$
				mCurrentPlatform = LINUX_PPC;
			else
				mCurrentPlatform = LINUX_OTHER;
		}
		else if (osName.equals("mac os x")) //$NON-NLS-1$
		{
			if (arch.startsWith("ppc") || arch.startsWith("power")) //$NON-NLS-1$ //$NON-NLS-2$
				mCurrentPlatform = MAC_OSX_PPC;
			else if (arch.equals("x86") || (arch.startsWith("i") && arch.endsWith("86"))) //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				mCurrentPlatform = MAC_OSX_x86;
			else if (arch.equals("x86_64")) //$NON-NLS-1$
				mCurrentPlatform = MAC_OSX_x86_64;
		}
		else if (osName.equals("mac os")) //$NON-NLS-1$
			mCurrentPlatform = MAC_CLASSIC;
		else if (osName.equals("sun os") || osName.equals("solaris")) //$NON-NLS-1$ //$NON-NLS-2$
		{
			if (arch.equals("x86") || (arch.startsWith("i") && arch.endsWith("86"))) //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				mCurrentPlatform = SOLARIS_x86;
			else if (arch.startsWith("sparc")) //$NON-NLS-1$
				mCurrentPlatform = SOLARIS_SPARC;
			else
				mCurrentPlatform = SOLARIS_OTHER;
		}
		else if (osName.endsWith("bsd")) //$NON-NLS-1$
		{
			if (arch.equals("x86") || (arch.startsWith("i") && arch.endsWith("86"))) //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
				mCurrentPlatform = BSD_x86;
			else
				mCurrentPlatform = BSD_OTHER;
		}
		else if (osName.equals("aix")) //$NON-NLS-1$
			mCurrentPlatform = AIX;
		else if (osName.equals("hp-ux")) //$NON-NLS-1$
		{
			if (arch.equals("pa-risc") || arch.equals("pa_risc")) //$NON-NLS-1$ //$NON-NLS-2$
				mCurrentPlatform = HPUX_PA_RISC;
			else if (arch.startsWith("ia64")) //$NON-NLS-1$
				mCurrentPlatform = HPUX_IA64;
			else
				mCurrentPlatform = HPUX_OTHER;
		}
		else if (osName.indexOf("digital unix") >= 0) //$NON-NLS-1$
			mCurrentPlatform = DEC_ALPHA;
		else if (osName.startsWith("os/2")) //$NON-NLS-1$
			mCurrentPlatform = OS2;
		mCurrentPlatform.mDescription = osName + ":" + arch; //$NON-NLS-1$
		return mCurrentPlatform;
	}

	public static boolean isWindows()
	{
		return getCurrentPlatform().mDescription.startsWith("windows"); //$NON-NLS-1$
	}

	public static boolean isMac()
	{
		return getCurrentPlatform().mDescription.startsWith("mac"); //$NON-NLS-1$
	}

	/**
	 * I need this since Apple JDK seems to confuse between Mac and Unix
	 * @return
	 */
	public static boolean isMacOSX()
	{
		EPlatform p = getCurrentPlatform();
		return p == MAC_OSX_x86_64 || p == MAC_OSX_x86 || p == MAC_OSX_PPC;
	}

	public static boolean isLinux()
	{
		return getCurrentPlatform().mDescription.startsWith("linux"); //$NON-NLS-1$
	}

	public static boolean isUnix()
	{
		EPlatform curPlatform = getCurrentPlatform();
		return curPlatform == SOLARIS_x86 || curPlatform == SOLARIS_SPARC || curPlatform == SOLARIS_OTHER ||
		curPlatform == MAC_OSX_PPC || curPlatform == MAC_OSX_x86 || curPlatform == MAC_OSX_x86_64 || curPlatform == AIX ||
		curPlatform == BSD_x86 || curPlatform == BSD_OTHER || curPlatform == DEC_ALPHA || curPlatform == HPUX_IA64 ||
		curPlatform == HPUX_PA_RISC || curPlatform == HPUX_OTHER;
	}

	public static String getPlatformVersion()
	{
		return System.getProperty("os.version"); //$NON-NLS-1$
	}

	public static int[] getPlatformVersionParsed()
	{	return parseVersionString(getPlatformVersion()); }

	public static final boolean is64Bit()
	{
		String arch = System.getProperty("os.arch").toLowerCase(Locale.ENGLISH); //$NON-NLS-1$
		String sun_arch_string = System.getProperty("sun.arch.data.model"); //$NON-NLS-1$
		boolean arch64 = arch.equals("x86_64") || arch.equals("amd64"); //$NON-NLS-1$ //$NON-NLS-2$
		if (sun_arch_string == null)
			return arch64;
		return arch64 && Integer.parseInt(sun_arch_string) == 64;
	}

	@Override
	public String toString()
	{	return mDescription; }
}
