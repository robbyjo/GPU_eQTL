/*
 * Roby Joehanes
 * 
 * Copyright 2009 Roby Joehanes
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
package gov.nih.opencl;

import gov.nih.utils.QSystemUtils;

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.bridj.BridJ;

import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLPlatform;
import com.nativelibs4java.opencl.InfoName;
import com.nativelibs4java.opencl.JavaCL;

import static gov.nih.utils.QStringUtils.sEmpty;
/**
 * 
 * @author Roby Joehanes
 *
 */
public class QOpenCL
{
	public static final int kDefaultAlignment = 64;
	private static String sOpenCLLibrary = null;
	static {
		new JavaCL.OpenCLProbeLibrary();
		sOpenCLLibrary = BridJ.getNativeLibraryFile("OpenCL").toString();
	}

	public static final String getDetectedOpenCLLib() {
		return sOpenCLLibrary;
	}

	public static final void init() {
		JavaCL.listPlatforms();
	}

	public static final QOCLContextPool getGPUContextPool(boolean onlyAvailable, boolean only64bits)
	{	return new QOCLContextPool(getAllGPUContexts(onlyAvailable, only64bits)); }

	public static final QOCLContextPool getAllContextPool(boolean onlyAvailable)
	{	return new QOCLContextPool(getAllContexts(onlyAvailable)); }

	public static final QOCLContext[] getAllGPUContexts(boolean onlyAvailable, boolean only64bits)
	{
		List<QOCLContext> contexts = new ArrayList<QOCLContext>();
		try
		{
			for (CLPlatform platform: JavaCL.listPlatforms())
				for (CLDevice device: platform.listGPUDevices(onlyAvailable))
					if (device.isCompilerAvailable() && (!only64bits || device.isDoubleSupported()))
						contexts.add(new QOCLContext(device));
			int numContexts = contexts.size();
			return numContexts > 0 ? contexts.toArray(new QOCLContext[numContexts]) : null;
		}
		catch (UnsatisfiedLinkError ex)
		{}
		return null; 
	}

	public static final QOCLContext[] getAllContexts(boolean onlyAvailable)
	{
		List<QOCLContext> contexts = new ArrayList<QOCLContext>();
		try
		{
			for (CLPlatform platform: JavaCL.listPlatforms())
				for (CLDevice device: platform.listAllDevices(onlyAvailable))
					if (device.isCompilerAvailable())
						contexts.add(new QOCLContext(device));
			int numContexts = contexts.size();
			return numContexts > 0 ? contexts.toArray(new QOCLContext[numContexts]) : null;
		}
		catch (UnsatisfiedLinkError ex)
		{}
		return null; 
	}

	public static final CLDevice[] getAllGPUDevices(boolean onlyAvailable, boolean only64bits)
	{
		List<CLDevice> devices = new ArrayList<CLDevice>();
		try
		{
			for (CLPlatform platform: JavaCL.listPlatforms())
				for (CLDevice device: platform.listGPUDevices(onlyAvailable))
					if (device.isCompilerAvailable() && (!only64bits || device.isDoubleSupported()))
						devices.add(device);
			int numDevices = devices.size();
			return numDevices > 0 ? devices.toArray(new CLDevice[numDevices]) : null;
		}
		catch (UnsatisfiedLinkError ex)
		{}
		return null; 
	}

	/**
	 * Query methods that contains InfoName annotation. This is shamelessly
	 * taken from HardwareReport example of JavaCL
	 * @param c
	 * @return
	 */
	private static final Map<String, Method> listInfoMethods(Class<?> c)
	{
		Map<String, Method> mets = new TreeMap<String, Method>();
		for (Method met : c.getMethods())
		{
			InfoName name = met.getAnnotation(InfoName.class);
			if (name == null) continue;
			mets.put(name.value(), met);
		}
		return mets;
    }

	/**
	 * List platform information.<br>
	 * Usage:<br>
	 * <pre>import com.nativelibs4java.opencl.*;
	 * import qutils.opencl.QOpenCL;
	 * ...
	 * QOpenCL.initLibrary();
	 * CLPlatform[] platforms = JavaCL.listPlatforms();
	 * listPlatformInformation(platforms[0]);
	 * </pre>
	 * <br>
	 * <p>Most computer only have one platform, except for a Beowulf cluster.
	 * If you have multiple platforms, you can use a loop to query platform information.<br>
	 * This is shamelessly taken from <tt>HardwareReport</tt> example of JavaCL.
	 * @param c
	 * @return
	 */
	public static List<Map<String, Object>> listPlatformInformation(CLPlatform platform, boolean gpuOnly)
	{
		try {
			List<Map<String, Object>> ret = new ArrayList<Map<String, Object>>();
			Map<String, Method> platMets = listInfoMethods(CLPlatform.class);
			Map<String, Method> devMets = listInfoMethods(CLDevice.class);
			Map<String, Object> platInfos = new TreeMap<String, Object>();
			for (Map.Entry<String, Method> platMet : platMets.entrySet()) {
				try {
					platInfos.put(platMet.getKey(), platMet.getValue().invoke(platform));
				} catch (InvocationTargetException ex) {
					if (ex.getCause() instanceof UnsupportedOperationException)
						platInfos.put(platMet.getKey(), sEmpty);
					else
						throw ex;
				}
			}
			for (CLDevice device : gpuOnly ? platform.listGPUDevices(false) : platform.listAllDevices(false)) {
				Map<String, Object> devInfos = new TreeMap<String, Object>(platInfos);
				for (Map.Entry<String, Method> devMet : devMets.entrySet()) {
					try {
						devInfos.put(devMet.getKey(), devMet.getValue().invoke(device));
					} catch (InvocationTargetException ex) {
						if (ex.getCause() instanceof UnsupportedOperationException)
							devInfos.put(devMet.getKey(), sEmpty);
						else
							throw ex;
					}
				}
				ret.add(devInfos);
			}
			return ret;
		} catch (Throwable ex) {
			throw new RuntimeException(ex);
		}
	}

	public static final Map<CLPlatform, List<Map<String, Object>>> listAllPlatformInformation(boolean gpuOnly)
	{
		Map<CLPlatform, List<Map<String, Object>>> map = new TreeMap<CLPlatform, List<Map<String, Object>>>();
		try {
			for (CLPlatform platform: JavaCL.listPlatforms()) {
				map.put(platform, listPlatformInformation(platform, gpuOnly));
			}
			return map;
		} catch (Throwable ex) {
			throw new RuntimeException(ex);
		}
	}

	private static String toString(Object value, String lnTag) {
		if (value == null) {
			return "null";
		}
		Class<?> c = value.getClass();
		try {
			if (c.isArray()) {
				if (!c.getComponentType().isPrimitive()) {
					c = Object[].class;
				}
				value = Arrays.class.getMethod("toString", c).invoke(null, value);
			}
		} catch (Exception ex) {
			throw new RuntimeException(ex);
		}
		String s = value.toString();
		if (s.startsWith("[") && s.endsWith("]")) {
			s = s.substring(1, s.length() - 1);
			s = s.replaceAll(", ", lnTag);
		}
		return s;
	}

	private static String toHTMLTable(List<Map<String, Object>> list) {
		StringBuilder b = new StringBuilder();

		b.append("<table border='1'>\n");
		if (!list.isEmpty()) {
			Set<String> keys = list.get(0).keySet();

			b.append("\t<tr valign=\"top\">\n");
			b.append("\t\t<td></td>");
			for (Map<String, Object> device : list) {
				Object value = device.get("CL_DEVICE_NAME");
				b.append("<td><b>[" + toString(value, "<br>") + "]</b></td>\n");
			}
			b.append("\n");
			b.append("\t</tr>\n");

			for (String key : keys) {
				b.append("\t<tr valign=\"top\">\n");
				b.append("\t\t<td>" + key + "</td>");
				for (Map<String, Object> device : list) {
					Object value = device.get(key);
					b.append("<td>" + toString(value, "<br>") + "</td>");
				}
				b.append("\n");
				b.append("\t</tr>\n");
			}
		}

		b.append("</table>\n");
		return b.toString();
	}

	private static final String toText(List<Map<String, Object>> list) {
		StringBuilder b = new StringBuilder();
		int deviceNo = 1;
		for (Map<String, Object> device: list)
		{
			b.append("Device #"+ deviceNo + ": " + device.get("CL_DEVICE_NAME") + "\n\n");
			deviceNo++;
			for (String key: device.keySet())
				b.append(key + ": " + toString(device.get(key), ", ") + "\n");
			b.append("\n\n");
		}
		return b.toString();
	}

	public static final String dumpPlatformInformation(CLPlatform platform, boolean asHTML, boolean gpuOnly)
	{
		StringBuilder b = new StringBuilder();
		if (asHTML)
			b.append("<P>Platform: " + platform.toString() + "<br>\n\n");
		else
			b.append("Platform: " + platform.toString() + "\n\n");
		List<Map<String, Object>> list = listPlatformInformation(platform, gpuOnly);
		final String text = asHTML ? toHTMLTable(list) : toText(list);
		b.append(text);
		return b.toString();
	}

	public static final String dumpAllPlatformInformation(boolean asHTML, boolean gpuOnly)
	{
		StringBuilder b = new StringBuilder();
		if (asHTML)
			b.append("<html xmlns=\"http://www.w3.org/1999/xhtml\">\n<body>\n");
		for (CLPlatform platform : JavaCL.listPlatforms())
			b.append(dumpPlatformInformation(platform, asHTML, gpuOnly) + "\n\n");
		if (asHTML)
			b.append("</body></html>");
		return b.toString();
	}

	public static final int getNumNVidiaDevicesInLinux()
	{
		int numDevices = 0;
		String[] cmd = new String[] {
			"/bin/sh", //$NON-NLS-1$
			"-c", //$NON-NLS-1$
			"" //$NON-NLS-1$
		};
		try {
			cmd[2] = "/sbin/lspci | grep -i NVIDIA | grep \"3D controller\" | wc -l"; //$NON-NLS-1$
			String[] output = QSystemUtils.exec(cmd);
			if (output != null)
				numDevices += Integer.parseInt(output[0].trim());
		} catch (Exception e) {}
		try {
			cmd[2] = "/sbin/lspci | grep -i NVIDIA | grep \"VGA compatible controller\" | wc -l"; //$NON-NLS-1$
			String[] output = QSystemUtils.exec(cmd);
			if (output != null)
				numDevices += Integer.parseInt(output[0].trim());
		} catch (Exception e) {}
		return numDevices;
	}

	private static final boolean isReadWriteable(String filename)
	{
		File f = new File(filename);
		return f.canRead() && f.canWrite();
	}

	/**
	 * Returns true if the permission is okay.
	 */
	public static final boolean checkNVidiaDevicePermissionInLinux()
	{
		int numDevices = getNumNVidiaDevicesInLinux();
		if (!isReadWriteable("/dev/nvidiactl")) //$NON-NLS-1$
			return false;
		for (int i = 0; i < numDevices; i++)
			if (!isReadWriteable("/dev/nvidia" + i)); //$NON-NLS-1$
		return true;
	}
}
