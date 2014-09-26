package gov.nih.exon.net;

import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.jgrapht.graph.SimpleGraph;

import qclustering.QClusterNode;
import qclustering.QClusterResult;
import qclustering.clusterer.flat.QFlatClusterResult;
import qclustering.io.qcf.QQCFParser;
import qgraph.INode;
import qgraph.QDefaultEdge;
import qgraph.QDefaultNode;
import qgraph.QNetwork;
import qgraph.io.QGMLPlugin;
import qplugin.network.ggm.QGGMUtils;

import gov.nih.table.QTableData;
import gov.nih.utils.QSystemUtils;
import gov.nih.exon.annotation.QAnnotation;

import static qclustering.QClusterNode.*;
import static qgraph.QNetwork.saveNetworks;
import static qstats.QStatsUtils.*;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableData;
import static gov.nih.utils.QFileUtils.writeText;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;
import static gov.nih.utils.QStringUtils.sLn;
import static gov.nih.utils.QStringUtils.sTab;

/**
 * @author Roby Joehanes
 *
 */
public class QExonGGM {
	static final int
		kFuzzyKDefaultNumClusters = 600,
		kFuzzyKDefaultNumCycles = 6,
		kIdealMaxNumChildrenPerCluster = 100,
		kMinClusterSize = 8,
		kMaxLevel = 2;
	static final String
		sGMLExtension = ".gml", // $NON-NLS-1$
		sNetworkDir = "network"; // $NON-NLS-1$
	static final int
		kMaskCL = 1,
		kMaskPAX = 2,
		kMaskPBMC = 4;

	static final void saveClusterResult(QClusterResult result, String filename) {
		try {
			QQCFParser.save(result, new FileOutputStream(filename));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Load cluster information from a file name
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static final QClusterNode[] loadClusters(String filename) throws IOException {
		QFlatClusterResult flatClusterResult = (QFlatClusterResult) QQCFParser.load(new FileInputStream(filename));
		return flatClusterResult.getClusters();
	}

	static final QNetwork<INode>[][] buildNetwork(QTableData data, QClusterNode[] allClusters) {
		double[][] dataArray = recenterData(data.extractData());
		int
			numRows = dataArray.length,
			numCols = dataArray[0].length,
			numParts = 3,
			onePart = numCols / numParts;
		double[][] partDataArray = new double[numRows][onePart];
		QNetwork[][] result = new QNetwork[numParts][];
		for (int partNo = 0; partNo < numParts; partNo++) {
			for (int rowNo = 0; rowNo < numRows; rowNo++)
				System.arraycopy(dataArray[rowNo], partNo * onePart, partDataArray[rowNo], 0, onePart);
			System.out.println("Subnetwork " + partNo);
			result[partNo] = QGGMUtils.buildNetwork(allClusters, partDataArray, null, null, 0.05, true, false);
			QSystemUtils.runGCAggressively();
		}
		return result;
	}

	private static final void compare(SimpleGraph<INode, QDefaultEdge<INode>> net1,
			SimpleGraph<INode, QDefaultEdge<INode>> net2, SimpleGraph<INode, QDefaultEdge<INode>> net3,
			int mask1, int mask2, int mask3, int[][] result) {
		for (QDefaultEdge<INode> edge: net1.edgeSet()) {
			INode
				src = edge.getSource(),
				dst = edge.getTarget();
			int
				srcID = src.getID(),
				dstID = dst.getID();
			if (srcID > dstID)
			{
				int temp = srcID;
				srcID = dstID;
				dstID = temp;
			}
			result[srcID][dstID] |= mask1;
			if (net2.getEdge(src, dst) != null)
				result[srcID][dstID] |= mask2;
			if (net3.getEdge(src, dst) != null)
				result[srcID][dstID] |= mask3;
		}
	}

	static final int[][] compareNetworks(QNetwork<INode>[][] networks, int numRows) {
		int
			numClusters = networks[0].length,
			result[][] = new int[numRows][numRows];
		for (int clusterNo = 0; clusterNo < numClusters; clusterNo++) {
			SimpleGraph<INode, QDefaultEdge<INode>>
				clNet = (SimpleGraph<INode, QDefaultEdge<INode>>) networks[0][clusterNo].getGraph(),
				paxNet = (SimpleGraph<INode, QDefaultEdge<INode>>) networks[1][clusterNo].getGraph(),
				pbmcNet = (SimpleGraph<INode, QDefaultEdge<INode>>) networks[2][clusterNo].getGraph();
			compare(clNet, paxNet, pbmcNet, kMaskCL, kMaskPAX, kMaskPBMC, result);
			compare(paxNet, clNet, pbmcNet, kMaskPAX, kMaskCL, kMaskPBMC, result);
			compare(pbmcNet, paxNet, clNet, kMaskPBMC, kMaskPAX, kMaskCL, result);
		}
		return result;
	}

	static final void saveComparison(int[][] result, String filename, String[] rowNames) throws IOException {
		int numRows = result.length;
		StringBuilder buf = new StringBuilder();
		String[] msg = new String[] {
			"",
			"The following edges are found in CL, but not in other cell types",
			"The following edges are found in PAX, but not in other cell types",
			"The following edges are found in CL and PAX, but not in PBMC",
			"The following edges are found in PBMC, but not in other cell types",
			"The following edges are found in CL and PBMC, but not in PAX",
			"The following edges are found in PBMC and PAX, but not in CL",
			"The following edges are found all cell types"
		};
		for (int k = 1; k < 8; k++) {
			buf.append(msg[k]);
			buf.append(sLn);
			for (int i = 1; i < numRows; i++) {
				for (int j = 0; j < i; j++) {
					if (result[j][i] == k) {
						buf.append(rowNames[j] + "--" + rowNames[i]);
						buf.append(sLn);
					}
				}
			}
			buf.append(sLn);
		}

		for (int i = 1; i < numRows; i++) {
			for (int j = 0; j < i; j++) {
				if (result[j][i] >= 8 || result[j][i] < 0) {
					System.out.println("Should never happened!");
				}
			}
		}
		writeText(buf.toString(), filename);
	}

	static final void annotateClusterNodes(QClusterNode[] clusters, QAnnotation[] annotations) {
		Map<String, QAnnotation> annotMapIndexedByEntrezID = QAnnotation.indexByEntrezID(annotations);
		Set<String> seenIDsWithoutAnnot = new HashSet<String>();
		for (QClusterNode cluster: clusters) {
			for (QClusterNode child: cluster.getChildren()) {
				Map<String, Object> prop = child.getProperties();
				if (prop == null || prop.size() == 0) {
					QAnnotation annot = annotMapIndexedByEntrezID.get(child.getName());
					if (annot == null) {
						String id = child.getName();
						seenIDsWithoutAnnot.add(id);
						System.out.println("No annotation found for gene " + id);
					} else {
						child.setProperties(annot.exportToProperty());
						child.setName(annot.mGeneSymbol);
					}
				}
			}
		}
	}

	static final void reannotateNetwork(String networkDir, String basename) {
		String dir = new File(basename).getParent();
		if (!dir.endsWith(File.separator))
			dir += File.separator;
		dir += networkDir + File.separator;
		QGMLPlugin<INode> plugin = new QGMLPlugin<INode>();
		for (File file: new File(dir).listFiles()) {
			String filename = file.getAbsolutePath();
			if (file.isDirectory() || filename.endsWith(".DS_Store")) // $NON-NLS-1$
				continue;
			try {
				QNetwork<INode> net = plugin.load(filename);
				for (INode vertex: net.getGraph().vertexSet()) {
					QDefaultNode node = (QDefaultNode) vertex;
					Object sym = node.getProperty(QAnnotation.sGeneSymbol);
					if (sym != null)
						node.setName(sym);
				}
				plugin.save(net, filename);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	static final String statNetwork(String networkDir, String basename) {
		String
			dir = new File(basename).getParent(),
			filename = new File(basename).getName();
		if (!dir.endsWith(File.separator))
			dir += File.separator;
		dir += networkDir + File.separator;
		File[] files = new File(dir).listFiles(new FileFilter() {
			@Override
			public boolean accept(File pathname) {
				return pathname.isFile() && pathname.getAbsolutePath().endsWith(sGMLExtension);
			}
		});

		int
			numSubnets = 3,
			numClusters = files.length / numSubnets,
			count[] = new int[8];
		assert (numClusters * numSubnets == files.length);
		files = null;
		filename = dir + filename;
		QGMLPlugin<INode> plugin = new QGMLPlugin<INode>();
		StringBuilder buf = new StringBuilder();
		buf.append("Cluster#" + sTab + "CLedges" + sTab + "PAXedges" + sTab + "PBMCedges" + sTab);
		buf.append("CLonly" + sTab + "PAXonly" + sTab + "CL_PAX" + sTab + "PBMConly" + sTab + "PBMC_CL" + sTab + "PBMC_PAX" + sTab + "All" + sLn);
		for (int clusterNo = 0; clusterNo < numClusters; clusterNo++) {
			String ext = clusterNo + sGMLExtension;
			Arrays.fill(count, 0);
			buf.append(clusterNo + sTab);
			try {
				SimpleGraph<INode, QDefaultEdge<INode>>
					clNet = (SimpleGraph<INode, QDefaultEdge<INode>>) plugin.load(filename + "-0-" + ext).getGraph(), // $NON-NLS-1$
					paxNet = (SimpleGraph<INode, QDefaultEdge<INode>>) plugin.load(filename + "-1-" + ext).getGraph(), // $NON-NLS-1$
					pbmcNet = (SimpleGraph<INode, QDefaultEdge<INode>>) plugin.load(filename + "-2-" + ext).getGraph(); // $NON-NLS-1$
				buf.append(clNet.edgeSet().size() + sTab +
					paxNet.edgeSet().size() + sTab +
					pbmcNet.edgeSet().size());
				/*
				 * 1 = CL only
				 * 2 = PAX only
				 * 3 = CL & PAX
				 * 4 = PBMC only
				 * 5 = PBMC & CL
				 * 6 = PBMC & PAX
				 * 7 = ALL
				 */
				for (QDefaultEdge<INode> edge: clNet.edgeSet()) {
					INode
						src = edge.getSource(),
						dst = edge.getTarget();
					boolean
						b1 = paxNet.getEdge(src, dst) != null,
						b2 = pbmcNet.getEdge(src, dst) != null;
					if (b1 & b2) // ALL
						count[7]++;
					else if (b1) // CL & PAX
						count[3]++;
					else if (b2) // CL & PBMC
						count[5]++;
					else
						count[1]++; // CL only
				}
				for (QDefaultEdge<INode> edge: paxNet.edgeSet()) {
					INode
						src = edge.getSource(),
						dst = edge.getTarget();
					if (clNet.getEdge(src, dst) == null) {
						if (pbmcNet.getEdge(src, dst) != null) // PBMC & PAX
							count[6]++;
						else // PAX only
							count[2]++;
					}
				}
				for (QDefaultEdge<INode> edge: pbmcNet.edgeSet()) {
					INode
						src = edge.getSource(),
						dst = edge.getTarget();
					if (clNet.getEdge(src, dst) == null && paxNet.getEdge(src, dst) == null) // PBMC only
						count[4]++;
				}
				for (int k = 1; k < 8; k++)
					buf.append(sTab + count[k]);
			} catch (Exception e) {
				e.printStackTrace();
			}
			buf.append(sLn);
		}
		return buf.toString();
	}

	public static void main(String[] args) {
		try {
			String
				fileName = args[0],
				ext, baseName;
			int pos = fileName.lastIndexOf('.');
			if (pos >= 0) {
				ext = fileName.substring(pos + 1);
				baseName = fileName.substring(0, pos);
			} else {
				baseName = fileName;
				ext = "";
			}
			//reannotateNetwork(sNetworkDir, baseName);
			//System.exit(0);
			String clusterFilename = baseName + "-cluster600l2." + ext;
			QTableData data = readDelimitedFileAsTableData(fileName, cCommaDelimiter, "#", true, true); // num genes x num samples
			QAnnotation[] annotations = QAnnotation.load(args[1]);
			QClusterNode[] clusters;
			/*
			clusters = qclustering.clusterer.flat.fuzzy.QFuzzyKPlugin.hierarchicalFuzzyCluster(data, kFuzzyKDefaultNumClusters, kFuzzyKDefaultNumCycles, baseName, ext,
					kIdealMaxNumChildrenPerCluster, kMinClusterSize, 0, kMaxLevel);
			QFlatClusterResult flatClusterResult = new QFlatClusterResult(clusters);
			//QFlatClusterResult flatClusterResult = cluster(recenterData(data), kFuzzyKDefaultNumClusters, kFuzzyKDefaultNumCycles, EFuzzyUnassignedRow.GROUP_INTO_ONE_CLUSTER);
			//clusters = flatClusterResult.getClusters();
			saveClusterResult(flatClusterResult, clusterFilename);
			/*/
			clusters = loadClusters(clusterFilename);
			//*/
			annotateClusterNodes(clusters, annotations);
			printClusterStats(clusters);
			String stat = statNetwork(sNetworkDir, baseName);
			writeText(stat, baseName + "-comparison-comp-600l2." + ext);
			System.out.println(stat);
			System.exit(0);
			QNetwork<INode>[][] networks = buildNetwork(data, clusters);
			saveNetworks(networks, sNetworkDir, baseName, sGMLExtension);
			int[][] compResult = compareNetworks(networks, data.getNumberOfRows());
			saveComparison(compResult, baseName + "-comparison600l2." + ext, data.getRowNames());
			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
