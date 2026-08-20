/*
 * Copyright 2026 Roby Joehanes
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 */
package gov.nih.gpu.opencl;

import gov.nih.gpu.GpuBackend;
import gov.nih.gpu.GpuDevice;
import gov.nih.gpu.GpuException;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_device_id;
import org.jocl.cl_platform_id;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.jocl.CL.CL_DEVICE_NOT_FOUND;
import static org.jocl.CL.CL_DEVICE_TYPE_GPU;
import static org.jocl.CL.CL_SUCCESS;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetDeviceInfo;
import static org.jocl.CL.clGetPlatformIDs;
import static org.jocl.CL.clGetPlatformInfo;

/** Maintained JOCL 2.x binding to vendor OpenCL ICD drivers. */
public final class JoclGpuBackend implements GpuBackend {
	private static final int CL_PLATFORM_NOT_FOUND_KHR = -1001;

	static {
		CL.setExceptionsEnabled(false);
	}

	@Override
	public String getName() {
		return "opencl";
	}

	@Override
	public String getRuntimeDescription() {
		return "JOCL 2.0.6 using the system OpenCL ICD loader";
	}

	@Override
	public List<GpuDevice> discoverGpuDevices() {
		try {
			int[] platformCount = new int[1];
			int status = clGetPlatformIDs(0, null, platformCount);
			if (status == CL_PLATFORM_NOT_FOUND_KHR || platformCount[0] == 0) {
				return Collections.emptyList();
			}
			check(status, "clGetPlatformIDs");

			cl_platform_id[] platforms = new cl_platform_id[platformCount[0]];
			check(clGetPlatformIDs(platforms.length, platforms, null), "clGetPlatformIDs");
			List<GpuDevice> devices = new ArrayList<GpuDevice>();
			for (cl_platform_id platform : platforms) {
				int[] deviceCount = new int[1];
				status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, null, deviceCount);
				if (status == CL_DEVICE_NOT_FOUND || deviceCount[0] == 0) {
					continue;
				}
				check(status, "clGetDeviceIDs");
				cl_device_id[] ids = new cl_device_id[deviceCount[0]];
				check(clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, ids.length, ids, null), "clGetDeviceIDs");
				for (cl_device_id id : ids) {
					devices.add(new JoclGpuDevice(platform, id));
				}
			}
			return Collections.unmodifiableList(devices);
		} catch (UnsatisfiedLinkError e) {
			throw new GpuException("The system OpenCL loader could not be loaded. Install a vendor GPU driver with OpenCL support.", e);
		}
	}

	static void check(int status, String operation) {
		if (status != CL_SUCCESS) {
			throw new GpuException(operation + " failed: " + CL.stringFor_errorCode(status) + " (" + status + ")");
		}
	}

	static String getPlatformString(cl_platform_id platform, int parameter) {
		long[] size = new long[1];
		check(clGetPlatformInfo(platform, parameter, 0, null, size), "clGetPlatformInfo");
		byte[] value = new byte[(int) size[0]];
		check(clGetPlatformInfo(platform, parameter, value.length, Pointer.to(value), null), "clGetPlatformInfo");
		return decodeString(value);
	}

	static String getDeviceString(cl_device_id device, int parameter) {
		long[] size = new long[1];
		check(clGetDeviceInfo(device, parameter, 0, null, size), "clGetDeviceInfo");
		byte[] value = new byte[(int) size[0]];
		check(clGetDeviceInfo(device, parameter, value.length, Pointer.to(value), null), "clGetDeviceInfo");
		return decodeString(value);
	}

	static int getDeviceInt(cl_device_id device, int parameter) {
		int[] value = new int[1];
		check(clGetDeviceInfo(device, parameter, Sizeof.cl_int, Pointer.to(value), null), "clGetDeviceInfo");
		return value[0];
	}

	static long getDeviceLong(cl_device_id device, int parameter) {
		long[] value = new long[1];
		check(clGetDeviceInfo(device, parameter, Sizeof.cl_long, Pointer.to(value), null), "clGetDeviceInfo");
		return value[0];
	}

	static long getDeviceSizeT(cl_device_id device, int parameter) {
		long[] value = new long[1];
		check(clGetDeviceInfo(device, parameter, Sizeof.size_t, Pointer.to(value), null), "clGetDeviceInfo");
		return value[0];
	}

	static long[] getDeviceSizeTArray(cl_device_id device, int parameter, int length) {
		long[] value = new long[length];
		check(clGetDeviceInfo(device, parameter, (long) Sizeof.size_t * length, Pointer.to(value), null), "clGetDeviceInfo");
		return value;
	}

	private static String decodeString(byte[] value) {
		int length = value.length;
		while (length > 0 && value[length - 1] == 0) {
			length--;
		}
		return new String(value, 0, length, StandardCharsets.UTF_8);
	}
}
