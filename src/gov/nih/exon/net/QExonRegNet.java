package gov.nih.exon.net;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.jgrapht.graph.SimpleGraph;

import qgraph.INode;
import qgraph.QDefaultEdge;
import qgraph.QDefaultNode;
import qgraph.QNetwork;
import qgraph.io.QGMLPlugin;
import qplugin.stats.fdr.QBenjaminiHochbergFDRPlugin;
import qplugin.stats.fdr.QFDRPlugin;
import qstats.glm.QGLM;
import qstats.glm.QPairwiseSimpleRegressionResult;
import gov.nih.table.IDataFilter;
import gov.nih.table.QDefaultColumnSelector;
import gov.nih.table.QTableData;
import gov.nih.utils.QFileUtils;
import gov.nih.utils.QSystemUtils;

import static qgraph.QGraphUtil.*;
import static gov.nih.utils.QDataUtils.translate;
import static gov.nih.utils.QFileUtils.readCSV;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableData;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;
import static gov.nih.utils.QStringUtils.cTabDelimiter;


/**
 * 
 * @author Roby Joehanes
 *
 */
public class QExonRegNet {
	static final int
		kAtlasID = 0,
		kUniqueID = 1,
		k2DBarCode = 2,
		kIncluded = 3,
		kSex = 4,
		kCaseType = 5,
		kPairNo = 6,
		kAB = 7,
		kMI = 8,
		kABI = 9,
		kCABG = 10,
		kPRCD = 11,
		kAge = 12,
		kSysBP = 13,
		kDiaBP = 14,
		kMedHyper = 15,
		kMedChol = 16,
		kSmoking = 17,
		kBMI = 18,
		kAlcohol = 19,
		kTotalChol = 20,
		kHDLChol = 21,
		kTriglyceride = 22,
		kMedDiabetes = 23,
		kGlucose = 24,
		includedCovars[] = new int[] { kSex, kAge, kSysBP, kDiaBP, kMedHyper, kMedChol, kSmoking, kBMI, kAlcohol, kTotalChol, kHDLChol, kTriglyceride };

	static double[][] removeUnusedCovar(double[][] covar, String[] covarNames, String[] colNames) {
		int
			numOrigCols = covarNames.length, // 450
			numCols = colNames.length, // 442
			numRows = covar[0].length; // 25
		double[][] result = new double[numRows][numCols];
		for (int i = 0; i < numOrigCols; i++) {
			String name = covarNames[i];
			for (int j = 0; j < numCols; j++) {
				if (name.equals(colNames[j])) {
					double[] curCovar = covar[i];
					for (int k = 0; k < numRows; k++)
						result[k][j] = curCovar[k];
					break;
				}
			}
		}
		return result;
	}

	static final double[][] getBatchCovar(String[] colNames) {
		int n = colNames.length;
		double[][] result = new double[5][n];

		for (int i = 0; i < n; i++) {
			int batchNo = Integer.parseInt(colNames[i].substring(1,3));
			if (batchNo > 1)
				result[batchNo - 2][i] = 1.0;
		}
		return result;
	}

	static final void mergeFiles() {
		int[][] net = new int[22011][];
		String dirName = "/Users/joehanesr/Desktop/regnet-gene/";
		try {
			for (String filename: new File(dirName).list()) {
				if (!filename.startsWith("rma.summary-"))
					continue;
				String[][] csv = QFileUtils.readCSV(dirName + filename);
				String sub = filename.substring(12, filename.lastIndexOf('.'));
				int
					numLines = csv.length,
					dashpt = sub.indexOf('-'),
					startPos = Integer.parseInt(sub.substring(0, dashpt)),
					endPos = Integer.parseInt(sub.substring(dashpt + 1)),
					delta = endPos - startPos;
				if (numLines != delta)
					System.out.println(numLines + " vs. " + delta);
				for (int i = 0; i < numLines; i++) {
					String[] curLine = csv[i];
					int
						numTokens = curLine.length,
						curGene[] = net[i + startPos] = new int[numTokens];
					for (int j = 0; j < numTokens; j++)
						curGene[j] = Integer.parseInt(curLine[j]);
				}
			}
			QFileUtils.writeText(simpleGraphToString(net), dirName + "regnet.txt");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	static final QNetwork<INode> annotateNet(String dataFilename, String netFilename, String annotFilename) throws IOException {
		//String annotFilename = "/Users/joehanesr/Desktop/annotation-gene.txt";
		String[][] annotations = readCSV(annotFilename);
		QTableData dataTable = readDelimitedFileAsTableData(dataFilename, cCommaDelimiter, "#", true, true); // num genes x num samples
		int[][] net = loadSimpleGraph(netFilename);
		System.out.println("Loaded.");
		return annotateNet(net, dataTable.getAllRowNames(), annotations, null, null, null, null);
	}

	static final QNetwork<INode> annotateNet(int[][] net, String[] transcriptIDs, String[][] annotations,
		double[] casePVals, double[] caseRSqs, double[] ctrlPVals, double[] ctrlRSqs) {
		String[] geneNames = translate(transcriptIDs, annotations, 0, 1);
		int numNodes = transcriptIDs.length;
		QDefaultNode[] nodes = new QDefaultNode[numNodes];
		SimpleGraph<INode, QDefaultEdge<INode>> graph = new SimpleGraph<INode, QDefaultEdge<INode>>((Class<? extends QDefaultEdge<INode>>) QDefaultEdge.class);

		boolean hasPVals = casePVals != null && caseRSqs != null && ctrlPVals != null && ctrlRSqs != null;
		int numPVals = casePVals.length;
		hasPVals = hasPVals && numPVals == caseRSqs.length && numPVals == ctrlPVals.length && numPVals == ctrlRSqs.length;
		if (net.length != numNodes)
			throw new RuntimeException();
		for (int i = 0; i < numNodes; i++) {
			String name = geneNames[i];
			nodes[i] = new QDefaultNode(Integer.parseInt(transcriptIDs[i]), name == null ? transcriptIDs[i] : name);
			graph.addVertex(nodes[i]);
		}
		for (int i = 0; i < numNodes; i++) {
			int[] parents = net[i];
			if (parents == null || parents.length == 0)
				continue;
			int numParents = parents.length;
			INode curNode = nodes[i];
			for (int j = 0; j < numParents; j++) {
				int paridx = parents[j];
				QDefaultEdge<INode> edge = graph.addEdge(nodes[paridx], curNode);
				if (hasPVals) {
					int idx = i > paridx ? paridx + i * (i - 1) / 2 : i + paridx * (paridx - 1) / 2;
					edge.setProperty("CasePVal", casePVals[idx]);
					edge.setProperty("CaseRSq", caseRSqs[idx]);
					edge.setProperty("CtrlPVal", ctrlPVals[idx]);
					edge.setProperty("CtrlRSq", ctrlRSqs[idx]);
				}
			}
		}
		return new QNetwork<INode>(graph);
	}

	static final int[][] doFDRCutoff(double[] pvals, QFDRPlugin plugin, double cutoff) {
		int
			numData = pvals.length,
			n = (int) Math.floor((1 + Math.sqrt(1 + 8 * numData)) / 2),
			result[][] = new int[n][];
		double[] fdrs = plugin.calculateFDR(pvals);
		List<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i < n; i++) {
			int idx = i * (i - 1) / 2;
			list.clear();
			for (int j = 0; j < i; j++)
				if (fdrs[idx + j] <= cutoff) list.add(j);
			int
				numSignif = list.size(),
				curResult[] = result[i] = new int[numSignif];
			for (int j = 0; j < numSignif; j++)
				curResult[j] = list.get(j);
		}
		return result;
	}

	public static void main(String[] args) {
		//annotateNet("/Users/joehanesr/Desktop/cc-resid-revised1-gene.txt",
		//		"/Users/joehanesr/Desktop/cc-resid-revised1-gene-diff.txt",
		//		"/Users/joehanesr/Desktop/annotation-gene.txt");
		try {
			String
				fileName = args[0],
				covarFileName = args[1],
				annotFile = args[2],
				//sStartIter = args[3],
				//sEndIter = args[4],
				sFToAddPVal = args[5],
				//sFToDropPVal = args[6],
				ext, baseName;
			int pos = fileName.lastIndexOf('.');
			if (pos >= 0) {
				ext = fileName.substring(pos + 1);
				baseName = fileName.substring(0, pos);
			} else {
				baseName = fileName;
				ext = "";
			}
			//int
				//startIter = Integer.parseInt(sStartIter);
			double
				fToAddPVal = Double.parseDouble(sFToAddPVal);

			IDataFilter filter = new IDataFilter() {
				@Override
				public boolean isRowSelected(int rowNo, double[] data, Object[] colNames) {
					return true;
				}

				@Override
				public boolean isColumnSelected(int colNo, Object colName, String[] colCategory) {
					return !"X02E03".equals(colName) && !"X02E04".equals(colName);
				}
				
				@Override
				public boolean isAllColumnsSelected() {
					return false;
				}
			};
			String strRep;
			long mem1, mem2, time1, time2; // Benchmark purposes

			mem1 = QSystemUtils.usedMemoryAfterGC();
			QTableData
				exonTable = readDelimitedFileAsTableData(fileName, cCommaDelimiter, "#", true, true), // num genes x num samples
				covarTable = readDelimitedFileAsTableData(covarFileName, cTabDelimiter, "#", true, true);
			String[][] annotations = readCSV(annotFile);
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));

			double[][]
				data = exonTable.createView(filter).extractData(),
				covar = removeUnusedCovar(covarTable.extractData(), covarTable.getRowNames(), exonTable.getColumnNames()),
				batchCovar = getBatchCovar(exonTable.getColumnNames());
			int
				numSubjects = covar[0].length,
				numGenes = data.length;
			if (numSubjects != data[0].length)
				throw new RuntimeException("Number of subjects in data and covariate table do not agree: " + data[0].length + " vs. " + numSubjects);
			System.out.println("Data loaded = " + numGenes + " x " + numSubjects);
			System.out.println("Covariate loaded = " + covar.length + " x " + numSubjects);
			System.out.println("Annotation loaded = " + annotations.length + " x " + annotations[0].length);
			assert(batchCovar != null);

			double[]
				mi = covar[9],
				abi = covar[10],
				cabg = covar[11],
				prcd = covar[12];
			BitSet cases = new BitSet(numSubjects);
			for (int i = 0; i < numSubjects; i++)
				if (mi[i] == 1 || abi[i] == 1 || cabg[i] == 1 || prcd[i] == 1)
					cases.set(i);
			BitSet ctrls = (BitSet) cases.clone();
			ctrls.flip(0, numSubjects);

			//*
			QPairwiseSimpleRegressionResult result1, result2;
			double[][]
				allcases = new QTableData(data).createView(new QDefaultColumnSelector(cases)).extractData(),
				allctrls = new QTableData(data).createView(new QDefaultColumnSelector(ctrls)).extractData();

			time1 = System.currentTimeMillis();
			result1 = QGLM.pairwiseSimpleRegression(allctrls, true);
			time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			//QFileUtils.writeText(QDataUtils.toStaggeredArray(result1.mPVals), baseName + "-ctrl-pvals." + ext);
			//QFileUtils.writeText(QDataUtils.toStaggeredArray(result1.mRSqs), baseName + "-ctrl-rsqs." + ext);
			int[][] result_ctrl = doFDRCutoff(result1.mPVals, new QBenjaminiHochbergFDRPlugin(), fToAddPVal);
			//strRep = simpleGraphToString(result_ctrl);
			//QFileUtils.writeText(strRep, baseName + "-ctrl." + ext);

			time1 = System.currentTimeMillis();
			result2 = QGLM.pairwiseSimpleRegression(allcases, true);
			time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			//QFileUtils.writeText(QDataUtils.toStaggeredArray(result2.mPVals), baseName + "-case-pvals." + ext);
			//QFileUtils.writeText(QDataUtils.toStaggeredArray(result2.mRSqs), baseName + "-case-rsqs." + ext);
			int[][] result_case = doFDRCutoff(result2.mPVals, new QBenjaminiHochbergFDRPlugin(), fToAddPVal);
			//strRep = simpleGraphToString(result_case);
			//QFileUtils.writeText(strRep, baseName + "-case." + ext);
			/*/

			int[][] result_ctrl = loadSimpleGraph("/Users/joehanesr/Desktop/cc-resid-revised1-gene-ctrl.txt");
			int[][] result_case = loadSimpleGraph("/Users/joehanesr/Desktop/cc-resid-revised1-gene-case.txt");
			System.out.println("Loaded.");
			/*/
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));
			time1 = System.currentTimeMillis();
			int[][] result_diff = calcDifferentialNetwork(result_case, result_ctrl);
			time2 = System.currentTimeMillis();
			System.out.println("Time = " + (time2 - time1));
			strRep = simpleGraphToString(result_diff);
			QFileUtils.writeText(strRep, baseName + "-diff." + ext);

			QNetwork<INode> network = annotateNet(result_diff, exonTable.getAllRowNames(), annotations,
				result2.mPVals, result2.mRSqs, result1.mPVals, result1.mRSqs);
			QGMLPlugin<INode> plugin = new QGMLPlugin<INode>();
			plugin.save(network, baseName+"-diff.gml");
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
