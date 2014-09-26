package gov.nih.exon.summarizer;

import gov.nih.parallel.IGenericParallelJob;
import gov.nih.parallel.QSynchronizedCounter;
import gov.nih.parallel.QThreadPool;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import qutils.QStringUtils;
import qutils.QSystemUtils;

import static gov.nih.exon.reader.QExonDataReader.readExonDataFile;
import static qutils.QFileUtils.writeText;
import static qutils.QStringUtils.sLn;
import static qutils.QStringUtils.sTab;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QSummarizer {
	public static final void summarize(Map<String, double[][]> geneData, ISummarizer summarizer, String delimiter, String saveFile)
		throws IOException {
		Set<String> geneIDSet = new TreeSet<String>();
		geneIDSet.addAll(geneData.keySet()); // Sorted
		String[] geneIDs = new String[geneIDSet.size()];
		geneIDSet.toArray(geneIDs);
		geneIDSet = null;
		// The parallel pattern is copied from my QeQTLAnalysisHotspot class
		int
			numIterations = geneIDs.length,
			numThreads = QSystemUtils.kNumCPUCores;
		double[][] summaries = new double[numIterations][];
		List<IGenericParallelJob> jobList = new ArrayList<IGenericParallelJob>();
		QSynchronizedCounter
			counter = new QSynchronizedCounter(0, numIterations),
			lock = new QSynchronizedCounter(0, numThreads);
		synchronized(lock) {
			for (int i = 0; i < numThreads; i++) {
				QSummarizerJob job = new QSummarizerJob(geneData, geneIDs, summarizer, summaries, lock, counter);
				jobList.add(job);
			}
			QThreadPool.mDefaultThreadPool.addAllJobs(jobList);
			while (lock.hasNext()) {
				try	{ lock.wait(); }
				catch (Exception e) {} // Interrupted exception is ignored
				//System.out.println(lock);
			}
		}

		StringBuilder buf = new StringBuilder();
		int geneNo = 0;
		for (String geneID: geneIDs) {
			buf.append(geneID); buf.append(delimiter);
			buf.append(QStringUtils.toString(summaries[geneNo++], delimiter));
			buf.append(sLn);
		}
		writeText(buf.toString(), saveFile);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			String fileName = args[0], ext, baseName;
			int pos = fileName.lastIndexOf('.');
			if (pos >= 0) {
				ext = fileName.substring(pos + 1);
				baseName = fileName.substring(0, pos);
			} else {
				baseName = fileName;
				ext = "";
			}
			long mem1 = QSystemUtils.usedMemoryAfterGC();
			Map<String, double[][]> geneData = readExonDataFile(fileName, sTab, true);
			long mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));
			ISummarizer[] summarizers = new ISummarizer[] {/*new QMeanSummarizer(),  QMedianSummarizer(),*/ new QPCASummarizer()};
			long time1 = System.currentTimeMillis();
			for (ISummarizer summarizer: summarizers) {
				String saveFile = baseName + '-' + summarizer.getPrefix() + '.' + ext;
				summarize(geneData, summarizer, sTab, saveFile);
			}
			long time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
