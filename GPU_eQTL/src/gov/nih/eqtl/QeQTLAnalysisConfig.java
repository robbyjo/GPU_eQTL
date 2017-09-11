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
package gov.nih.eqtl;

import gov.nih.eqtl.io.QIniLoader;
import gov.nih.utils.QSystemUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Locale;
import java.util.Map;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QeQTLAnalysisConfig {
	private Map<String, String> mIni = null;
	
	public QeQTLAnalysisConfig(Map<String, String> ini)
	{	mIni = ini; }

	private static final String expandPath(String filename, String path)
	{
		if (filename == null)
			return null;
		File f = new File(filename);
		if (!f.isAbsolute())
			filename = path + filename;
		return filename;
	}

	public String getIniPath()
	{	return mIni.get("ini.path"); } //$NON-NLS-1$

	public String getFamilyFilename()
	{	return expandPath(mIni.get("family_file"), getIniPath()); } //$NON-NLS-1$

	public String getGenotypeFilename()
	{	return expandPath(mIni.get("genotype_file"), getIniPath()); } //$NON-NLS-1$

	public String getExpressionFilename()
	{	return expandPath(mIni.get("expression_file"), getIniPath()); } //$NON-NLS-1$

	public String getCovariateFilename()
	{	return expandPath(mIni.get("covariate_file"), getIniPath()); } //$NON-NLS-1$

	public String getPedigreeFilename()
	{	return expandPath(mIni.get("pedigree_file"), getIniPath()); } //$NON-NLS-1$

	public String getOutputFilename()
	{	return expandPath(mIni.get("output_file"), getIniPath()); } //$NON-NLS-1$

	public String[] getLibraryPaths()
	{
		String s = mIni.get("library_path"); //$NON-NLS-1$
		return s == null ? null : s.split(":"); //$NON-NLS-1$
	}

	public String getGenotypeFileFormat()
	{
		String s =  mIni.get("genotype_format");  //$NON-NLS-1$
		return s == null ? "tped" : s.toLowerCase(Locale.ENGLISH);  //$NON-NLS-1$
	}

	public boolean getGenotypeFileHeader()
	{
		String s =  mIni.get("genotype_file_header");  //$NON-NLS-1$
		return s == null ? true : Boolean.parseBoolean(s);
	}

	public boolean getSimplifyOutput()
	{
		String s =  mIni.get("simplify_output");  //$NON-NLS-1$
		return s == null ? false : Boolean.parseBoolean(s);
	}

	public boolean getRSqOutput()
	{
		String s =  mIni.get("rsq_only");  //$NON-NLS-1$
		return s == null ? false : Boolean.parseBoolean(s);
	}

	public boolean getPValOutput()
	{
		String s =  mIni.get("pval_only");  //$NON-NLS-1$
		return s == null ? false : Boolean.parseBoolean(s);
	}

	public String getGenotypeModel()
	{
		String s =  mIni.get("genotype_model");  //$NON-NLS-1$
		return s == null ? "additive" : s.toLowerCase(Locale.ENGLISH);  //$NON-NLS-1$
	}

	public String[] getFixedCovariates()
	{
		String s = mIni.get("covariate_fixed"); //$NON-NLS-1$
		return s == null ? null : s.split(" "); //$NON-NLS-1$
	}

	public String[] getRandomCovariates()
	{
		String s = mIni.get("covariate_random"); //$NON-NLS-1$
		return s == null ? null : s.split(" "); //$NON-NLS-1$
	}

	public String[] getFactorCovariates()
	{
		String s = mIni.get("covariate_factor"); //$NON-NLS-1$
		return s == null ? null : s.split(" "); //$NON-NLS-1$
	}

	public String getThresholdType()
	{
		String thresholdType = mIni.get("threshold"); //$NON-NLS-1$
		if (thresholdType == null)
			return "none"; //$NON-NLS-1$
		return thresholdType.split(" ")[0].trim(); //$NON-NLS-1$
	}

	public double getThresholdValue()
	{
		String thresholdType = mIni.get("threshold"); //$NON-NLS-1$
		if (thresholdType == null)
			return Double.NaN;
		return Double.parseDouble(thresholdType.split(" ")[1].trim()); //$NON-NLS-1$
	}

	public int getDFOffset()
	{
		String str = mIni.get("df_offset"); //$NON-NLS-1$
		if (str == null)
			return 0;
		return Integer.parseInt(str);
	}

	public void setDFOffset(int df)
	{
		mIni.put("block_size", String.valueOf(df)); //$NON-NLS-1$
	}

	public int getBlockSize()
	{
		String str = mIni.get("block_size"); //$NON-NLS-1$
		if (str == null)
			return 0;
		return Integer.parseInt(str);
	}

	public void setBlockSize(int s)
	{
		mIni.put("block_size", String.valueOf(s)); //$NON-NLS-1$
	}

	public int getNumThreads()
	{
		String str = mIni.get("num_threads"); //$NON-NLS-1$
		if (str == null)
			return QSystemUtils.kNumCPUCores;
		return Integer.parseInt(str);
	}

	public void setNumThreads(int s)
	{
		mIni.put("num_threads", String.valueOf(s)); //$NON-NLS-1$
	}

	public static final QeQTLAnalysisConfig loadConfig(String iniFilename) throws IOException
	{
		Map<String, String> ini = null;
		FileReader fr = null;
		if (iniFilename.contains("~"))
			iniFilename = iniFilename.replace("~",System.getProperty("user.home")); //$NON-NLS-1$ //$NON-NLS-2$
		File iniFile = new File(iniFilename).getAbsoluteFile();
		try {
			fr = new FileReader(iniFile);
			ini = QIniLoader.load(fr);
		} finally {
			try {
				if (fr != null) fr.close();
			} catch (Exception e) {}
		}
		String iniPath = iniFile.getParent();
		if (!iniPath.endsWith(File.separator))
			iniPath = iniPath + File.separator;
		ini.put("ini.path", iniPath); //$NON-NLS-1$
		return new QeQTLAnalysisConfig(ini);
	}
}
