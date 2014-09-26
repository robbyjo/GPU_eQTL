package gov.nih.exon.annotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import static gov.nih.utils.QFileUtils.*;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QAnnotation {
	static final int
		kDefaultBufferSize = 4 * 1024 * 1024;
	public static final String
		sListSeparatorRegex = "\\|", // $NON-NLS-1$
		sListSeparator = "|", // $NON-NLS-1$
		sMapSeparator = ":", // $NON-NLS-1$
		sNA = "N/A", // $NON-NLS-1$
		sTransClustID = "transcript_cluster_id", // $NON-NLS-1$
		sProbesetID = "probeset_id", // $NON-NLS-1$
		sGeneSymbol = "GeneSymbol", // $NON-NLS-1$
		sChromosome = "seqname", // $NON-NLS-1$
		sDescription = "GeneDescription", // $NON-NLS-1$
		sStrand = "strand", // $NON-NLS-1$
		sStart = "start", // $NON-NLS-1$
		sStop = "stop", // $NON-NLS-1$
		sEntrezID = "EntrezGeneID", // $NON-NLS-1$
		sChrPrefix = "chr"; // $NON-NLS-1$

	public int
		mTransClustID,
		mProbesetID,
		mEntrezID,
		mStart,
		mStop;
	public String
		mGeneSymbol,
		mChromosome, // Can be "X" or "Y" or "mitochondria"
		mDescription,
		mStrand;

	public QAnnotation() {}
	public QAnnotation(int probesetID, int transClustID, int entrezID, String geneSymbol, String chrom, String strand, int start, int stop, String desc) {
		mProbesetID = probesetID;
		mTransClustID = transClustID;
		mEntrezID = entrezID;
		mGeneSymbol = geneSymbol;
		mChromosome = chrom;
		mStrand = strand;
		mStart = start;
		mStop = stop;
		mDescription = desc;
	}

	protected static final String joinList(String[] strs) {
		StringBuilder buf = new StringBuilder();
		int n = strs.length;
		if (n == 0)
			return null;
		buf.append(strs[0]);
		for (int i = 1; i < n; i++)
			buf.append(sListSeparator + strs[i]);
		return buf.toString();
	}

	public Map<String, Object> exportToProperty() {
		Map<String, Object> map = new HashMap<String, Object>();
		if (mProbesetID != -1)
			map.put(sProbesetID, mProbesetID);
		if (mTransClustID != -1)
			map.put(sTransClustID, mTransClustID);
		if (mEntrezID != -1)
			map.put(sEntrezID, mEntrezID);
		if (mGeneSymbol != null)
			map.put(sGeneSymbol, mGeneSymbol);
		if (mDescription != null)
			map.put(sDescription, mDescription);
		if (mChromosome != null)
			map.put(sChromosome, mChromosome);
		if (mStrand != null)
			map.put(sStrand, mStrand);
		if (mStart != -1)
			map.put(sStart, mStart);
		if (mStop != -1)
			map.put(sStop, mStop);
		return map;
	}

	@Override
	public String toString()
	{	return mGeneSymbol + '@' + mChromosome; }

	public static final QAnnotation[] load(String filename) throws IOException
	{	return load(new BufferedReader(new FileReader(filename), kDefaultBufferSize), cCommaDelimiter); }

	public static final QAnnotation[] load(String filename, char[] delim) throws IOException
	{	return load(new BufferedReader(new FileReader(filename), kDefaultBufferSize), delim); }

	public static final QAnnotation[] load(Reader reader, char[] delim) throws IOException {
		String[][] annotStrs = null;
		annotStrs = readDelimitedFile(reader, delim, "#");
		int
			numAnnotations = annotStrs.length - 1,
			idxProbesetID,
			idxTransClustID,
			idxGeneSymbol,
			idxEntrezID,
			idxChroms,
			idxStrand,
			idxStart,
			idxStop,
			idxDescription;
		QAnnotation[] annots = new QAnnotation[numAnnotations];
		List<String> header = new ArrayList<String>();
		for (String hdr: annotStrs[0])
			header.add(hdr.toLowerCase(Locale.ENGLISH));

		// Figure out the index positions
		idxProbesetID = header.indexOf(sProbesetID.toLowerCase(Locale.ENGLISH));
		idxTransClustID = header.indexOf(sTransClustID.toLowerCase(Locale.ENGLISH));
		idxEntrezID = header.indexOf(sEntrezID.toLowerCase(Locale.ENGLISH));
		idxGeneSymbol = header.indexOf(sGeneSymbol.toLowerCase(Locale.ENGLISH));
		idxDescription = header.indexOf(sDescription.toLowerCase(Locale.ENGLISH));
		idxChroms = header.indexOf(sChromosome.toLowerCase(Locale.ENGLISH));
		idxStrand = header.indexOf(sStrand.toLowerCase(Locale.ENGLISH));
		idxStart = header.indexOf(sStart.toLowerCase(Locale.ENGLISH));
		idxStop = header.indexOf(sStop.toLowerCase(Locale.ENGLISH));
		int chrPrefixLength = sChrPrefix.length();

		if (idxProbesetID >= 0) {
			for (int i = 1; i <= numAnnotations; i++) {
				String[] s = annotStrs[i];
				annots[i-1] = new QAnnotation(Integer.parseInt(s[idxProbesetID]),
					idxTransClustID >= 0 ? Integer.parseInt(s[idxTransClustID]) : -1,
					idxEntrezID >= 0 ? Integer.parseInt(s[idxEntrezID]) : -1,
					idxGeneSymbol >= 0 ? s[idxGeneSymbol] : sNA,
					idxChroms >= 0 ? (s[idxChroms].startsWith(sChrPrefix) ? s[idxChroms].substring(chrPrefixLength) : s[idxChroms]) : sNA,
					idxStrand >= 0 ? s[idxStrand] : sNA,
					idxStart >= 0 ? Integer.parseInt(s[idxStart]) : -1,
					idxStop >= 0 ? Integer.parseInt(s[idxStop]) : -1,
					idxDescription >= 0 ? s[idxDescription] : sNA);
			}
		} else {
			for (int i = 1; i <= numAnnotations; i++) {
				String[] s = annotStrs[i];
				annots[i-1] = new QAnnotation(-1,
					idxTransClustID >= 0 ? Integer.parseInt(s[idxTransClustID]) : -1,
					idxEntrezID >= 0 ? Integer.parseInt(s[idxEntrezID]) : -1,
					idxGeneSymbol >= 0 ? s[idxGeneSymbol] : sNA,
					idxChroms >= 0 ? (s[idxChroms].startsWith(sChrPrefix) ? s[idxChroms].substring(chrPrefixLength) : s[idxChroms]) : sNA,
					idxStrand >= 0 ? s[idxStrand] : sNA,
					idxStart >= 0 ? Integer.parseInt(s[idxStart]) : -1,
					idxStop >= 0 ? Integer.parseInt(s[idxStop]) : -1,
					idxDescription >= 0 ? s[idxDescription] : sNA);
			}
		}
		return annots;
	}

	public static final Map<String, QAnnotation> indexByEntrezID(QAnnotation[] annots) {
		int n = annots.length;
		Map<String, QAnnotation> map = new HashMap<String, QAnnotation>(n);
		for (QAnnotation annot: annots)
			map.put(String.valueOf(annot.mEntrezID), annot);
		return map;
	}
}
