/*
 * Roby Joehanes
 * 
 * Copyright 2011 Roby Joehanes
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
import gov.nih.gpu.GpuPrecision;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;
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
	{	mIni = new LinkedHashMap<String, String>(ini); }

	private static final String expandPath(String filename, String path)
	{
		if (filename == null)
			return null;
		File f = new File(filename);
		if (!f.isAbsolute())
			filename = (path == null ? "" : path) + filename;
		return filename;
	}

	public String get(String key)
	{ return mIni.get(key); }

	public void set(String key, String value)
	{
		if (value == null)
			mIni.remove(key);
		else
			mIni.put(key, value);
	}

	public Map<String, String> asMap()
	{ return new LinkedHashMap<String, String>(mIni); }

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

	public String getExpressionFileFormat()
	{
		String s = mIni.get("expression_format");
		return s == null ? "csv" : s.toLowerCase(Locale.ENGLISH);
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

	public boolean getValidateOnly()
	{
		String s = mIni.get("validate_only");
		return s != null && Boolean.parseBoolean(s);
	}

	public String getGenotypeModel()
	{
		String s =  mIni.get("genotype_model");  //$NON-NLS-1$
		return s == null ? "additive" : s.toLowerCase(Locale.ENGLISH);  //$NON-NLS-1$
	}

	public String getGenotypeField()
	{ return mIni.get("genotype_field"); }

	public String getGenotypeMissingPolicy()
	{
		String value = mIni.get("predictor_missing");
		if (value == null)
			value = mIni.get("genotype_missing");
		return value == null ? "mean" : value;
	}

	public QDataType getPredictorDataType()
	{ return QDataType.parse(mIni.get("predictor_type"), QDataType.GENOTYPE); }

	public QDataType getTraitDataType()
	{ return QDataType.parse(mIni.get("trait_type"), QDataType.EXPRESSION); }

	public QMissingValuePolicy getPredictorMissingPolicy()
	{ return QMissingValuePolicy.parse(getGenotypeMissingPolicy(), QMissingValuePolicy.MEAN); }

	public QMissingValuePolicy getTraitMissingPolicy()
	{ return QMissingValuePolicy.parse(mIni.get("trait_missing"), QMissingValuePolicy.PATTERN); }

	public int getPredictorFlankCount()
	{
		String value = mIni.get("predictor_flanks");
		if (value == null)
			return 1;
		int parsed = Integer.parseInt(value);
		if (parsed < 1)
			throw new IllegalArgumentException("predictor_flanks must be at least 1");
		return parsed;
	}

	public String getCovariateMissingPolicy()
	{
		String value = mIni.get("covariate_missing");
		return value == null ? "complete-samples" : value.trim().toLowerCase(Locale.ROOT).replace('_', '-');
	}

	public boolean getInspectMissingness()
	{ return getBoolean("inspect_missingness"); }

	public String getMissingnessQcOutputFilename()
	{ return expandPath(mIni.get("missingness_qc_output"), getIniPath()); }

	public String getMultiallelicPolicy()
	{ return mIni.get("multiallelic"); }

	public double getMinimumMaf()
	{ return getNonNegativeDouble("min_maf"); }

	public double getMinimumMac()
	{ return getNonNegativeDouble("min_mac"); }

	public String getVariantQcOutputFilename()
	{ return expandPath(mIni.get("variant_qc_output"), getIniPath()); }

	public String[] getFixedCovariates()
	{
		String s = mIni.get("covariate_fixed"); //$NON-NLS-1$
		return s == null ? null : s.split("\\s+"); //$NON-NLS-1$
	}

	public String[] getRandomCovariates()
	{
		String s = mIni.get("covariate_random"); //$NON-NLS-1$
		return s == null ? null : s.split("\\s+"); //$NON-NLS-1$
	}

	public String[] getFactorCovariates()
	{
		String s = mIni.get("covariate_factor"); //$NON-NLS-1$
		return s == null ? null : s.split("\\s+"); //$NON-NLS-1$
	}

	public String getThresholdType()
	{
		String thresholdType = mIni.get("threshold"); //$NON-NLS-1$
		if (thresholdType == null)
			return "none"; //$NON-NLS-1$
		return thresholdType.split("\\s+")[0].trim(); //$NON-NLS-1$
	}

	public double getThresholdValue()
	{
		String thresholdType = mIni.get("threshold"); //$NON-NLS-1$
		if (thresholdType == null)
			return Double.NaN;
		return Double.parseDouble(thresholdType.split("\\s+")[1].trim()); //$NON-NLS-1$
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
		mIni.put("df_offset", String.valueOf(df)); //$NON-NLS-1$
	}

	public int getBlockSize()
	{
		String str = mIni.get("block_size"); //$NON-NLS-1$
		if (str == null)
			return 0;
		int value = Integer.parseInt(str);
		if (value < 0)
			throw new IllegalArgumentException("block_size must not be negative");
		return value;
	}

	public void setBlockSize(int s)
	{
		mIni.put("block_size", String.valueOf(s)); //$NON-NLS-1$
	}

	public int getNumThreads()
	{
		String str = mIni.get("num_threads"); //$NON-NLS-1$
		if (str == null)
			return 0;
		int value = Integer.parseInt(str);
		if (value < 0)
			throw new IllegalArgumentException("num_threads must not be negative");
		return value;
	}

	public void setNumThreads(int s)
	{
		mIni.put("num_threads", String.valueOf(s)); //$NON-NLS-1$
	}

	public GpuPrecision getGpuPrecision()
	{ return GpuPrecision.parse(mIni.get("precision")); }

	public QResidualizationMode getResidualizationMode()
	{ return QResidualizationMode.parse(mIni.get("residualization")); }

	public int getGenotypeBlockRows()
	{ return getNonNegativeInt("genotype_block_rows"); }

	public int getExpressionBlockRows()
	{ return getNonNegativeInt("expression_block_rows"); }

	public String getGenotypeIdColumn()
	{ return mIni.get("genotype_id_column"); }

	public String getExpressionIdColumn()
	{ return mIni.get("expression_id_column"); }

	public String getCacheDirectory()
	{ return expandPath(mIni.get("cache_dir"), getIniPath()); }

	public boolean getRebuildCache()
	{ return getBoolean("rebuild_cache"); }

	public String getCheckpointDirectory()
	{ return expandPath(mIni.get("checkpoint_dir"), getIniPath()); }

	public boolean getResume()
	{ return getBoolean("resume"); }

	public boolean getKeepCheckpoints()
	{ return getBoolean("keep_checkpoints"); }

	public boolean getProfile()
	{ return getBoolean("profile") || mIni.get("profile_output") != null; }

	public String getProfileOutputFilename()
	{ return expandPath(mIni.get("profile_output"), getIniPath()); }

	private boolean getBoolean(String key)
	{
		String value = mIni.get(key);
		return value != null && Boolean.parseBoolean(value);
	}

	private int getNonNegativeInt(String key)
	{
		String value = mIni.get(key);
		if (value == null)
			return 0;
		int parsed = Integer.parseInt(value);
		if (parsed < 0)
			throw new IllegalArgumentException(key + " must not be negative");
		return parsed;
	}

	private double getNonNegativeDouble(String key)
	{
		String value = mIni.get(key);
		if (value == null)
			return 0;
		double parsed = Double.parseDouble(value);
		if (!Double.isFinite(parsed) || parsed < 0)
			throw new IllegalArgumentException(key + " must be a finite non-negative number");
		return parsed;
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
