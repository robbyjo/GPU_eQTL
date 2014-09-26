package gov.nih.exon.net;


import java.util.BitSet;

import qplugin.stats.fdr.QBenjaminiHochbergFDRPlugin;
import qstats.glm.QGLM;
import qstats.glm.QGLMResult;
import qstats.glm.QGLMResultBasic;
import gov.nih.parallel.ISimpleParallelJob;
import gov.nih.table.QTableData;
import gov.nih.utils.QFileUtils;
import gov.nih.utils.QStringUtils;
import gov.nih.utils.QSystemUtils;

import static qstats.QStatsUtils.calcPOfF;
import static gov.nih.parallel.QParallelUtils.*;
import static gov.nih.utils.QFileUtils.readDelimitedFileAsTableData;
import static gov.nih.utils.QStringUtils.cCommaDelimiter;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QGeneRegNetCC implements ISimpleParallelJob {
	static final double[]
		caseness = new double[] {
		0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
		0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
		0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1,
		0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0,
		0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1,
		0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0,
		1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0,
		1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0,
		1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0,
		1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
		1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0,
		1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1,
		1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0,
		0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
		1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0,
		1, 0, 1, 1, 0, 0, 1 },
		sbp = new double[] {
		115, 120, 148, 141, 154, 124, 125, 146, 113, 130, 115, 113, 132, 117, 123, 141,  96, 107,
		164, 146, 132, 153,  92,  96, 109, 124, 137, 126, 125, 131, 157, 114, 137, 102, 124, 164,
		112, 111, 113, 115, 106, 128, 159, 147, 122, 148, 118, 120, 119, 148, 158, 115, 108, 140,
		113, 111, 115, 129, 139, 104, 104, 119, 164, 143, 122, 108, 122, 150, 140, 128, 132, 146,
		127, 136, 192, 137, 139, 160,  96, 113, 147, 144, 143, 126, 159, 144, 128, 146, 106, 130,
		127, 143, 122, 106, 105,  90, 115, 117, 106, 116, 146, 129, 165, 128, 139, 117, 123, 132,
		161, 136, 110, 149, 143, 105, 118,  95, 131, 127, 117, 130, 133, 139, 127, 123, 149, 155,
		137, 138, 131, 115, 145, 130, 170, 120, 139, 158, 141, 148, 125, 119, 143, 126, 127, 144,
		112, 143, 144, 136, 123, 116, 157, 162, 129, 127, 135, 147, 127, 129, 132, 181, 130, 138,
		138, 115, 106, 146, 131,  99, 133, 116, 111, 112, 108, 115, 139, 124, 118, 133, 116, 146,
		116, 111, 150, 115, 135, 141, 130, 123, 113, 178, 127, 122, 137, 145, 100, 115, 147, 138,
		121, 106, 152,  97, 127, 125, 157, 125, 116, 171, 130, 137, 131,  94, 141, 145, 145, 111,
		113, 130, 144, 112, 141, 191, 152, 117, 111, 135, 120, 118, 109, 116, 129, 145, 105, 110,
		117,  99, 141, 133, 103, 123, 156, 137, 152, 121, 139, 127, 110, 123, 116, 142, 163, 119,
		130, 113, 136, 138, 116, 159, 120, 116, 135, 126, 145, 123,  97, 128, 104, 136, 111, 121,
		125, 136, 124, 133, 130, 133, 125, 109, 112, 112, 139, 140, 136, 139, 137, 107, 118, 129,
		 94, 168, 125,  96, 129, 117, 135, 114, 128, 138, 115,  98, 108, 114, 179, 120, 143, 111,
		122, 124, 165, 117, 104, 141, 170, 118, 143, 115, 156, 122, 104, 130, 117, 140, 133, 127,
		154, 152, 110, 112, 122, 128, 136, 144, 130, 125, 136, 130, 129, 124, 137, 121, 125, 155,
		126, 105, 155, 155, 128, 107, 118, 124, 143, 158, 127, 144, 142, 135, 155, 129, 118, 145,
		131, 146, 121, 119, 115, 133, 121, 103, 123, 101, 114, 120, 130,  95, 136, 126, 105, 134,
		 89, 115, 165,  98, 146, 141, 141, 113, 119, 130, 116, 150, 132, 107, 125, 148, 199, 152,
		115, 125, 126, 133,  94, 164, 109, 142, 118, 123, 131, 158, 144, 147, 107, 122, 136, 132,
		126, 158, 124, 129, 157, 129, 134, 125, 145, 137, 131, 144, 143, 119, 140, 132, 131, 116,
		134, 154, 153, 141, 199, 110, 166, 123, 134, 127 },
		dbp = new double[] {
		63, 52, 87, 70, 78, 46, 65, 85, 68, 83, 74, 68, 75, 68, 68, 81, 59, 67, 79, 70, 84, 68, 41,
		54, 69, 75, 70, 65, 62, 85, 68, 76, 82, 73, 69, 63, 61, 69, 65, 77, 75, 77, 70, 68, 77, 70,
		60, 77, 76, 87, 87, 79, 60, 74, 71, 78, 73, 82, 71, 47, 68, 81, 68, 51, 78, 72, 70, 75, 82,
		60, 68, 73, 74, 79, 73, 94, 77, 86, 35, 63, 80, 53, 64, 67, 79, 71, 71, 74, 59, 69, 75, 59,
		56, 62, 68, 54, 76, 75, 58, 67, 57, 86, 53, 61, 67, 68, 70, 66, 85, 87, 67, 87, 69, 47, 68,
		62, 78, 64, 62, 73, 76, 74, 39, 67, 60, 79, 73, 61, 57, 58, 89, 63, 57, 64, 57, 61, 72, 67,
		72, 70, 82, 57, 67, 69, 47, 86, 60, 73, 61, 63, 91, 92, 85, 82, 76, 70, 66, 67, 68, 60, 82,
		90, 83, 84, 66, 84, 73, 56, 86, 74, 76, 82, 74, 70, 77, 77, 64, 54, 79, 72, 64, 60, 81, 58,
		77, 78, 62, 84, 66, 90, 50, 66, 79, 82, 52, 68, 74, 73, 53, 67, 78, 49, 73, 69, 59, 68, 75,
		74, 74, 82, 71, 60, 57, 68, 82, 75, 73, 65, 77, 71, 74, 102, 89, 67, 74, 90, 82, 47, 70, 72,
		71, 88, 49, 73, 73, 49, 73, 68, 68, 79, 80, 54, 73, 60, 77, 74, 65, 59, 77, 67, 75, 73, 59,
		78, 73, 83, 69, 85, 71, 81, 72, 65, 80, 78, 52, 85, 62, 62, 69, 52, 77, 85, 82, 90, 73, 74,
		62, 77, 70, 65, 68, 77, 72, 78, 80, 50, 48, 79, 59, 86, 81, 58, 68, 86, 88, 75, 78, 59, 69,
		52, 64, 75, 66, 70, 74, 68, 82, 80, 74, 70, 70, 73, 92, 62, 72, 64, 60, 53, 68, 72, 65, 74,
		90, 57, 80, 88, 56, 78, 68, 73, 66, 70, 78, 60, 65, 47, 64, 85, 73, 79, 70, 82, 67, 65, 84,
		67, 60, 57, 71, 72, 74, 82, 53, 88, 72, 68, 75, 61, 77, 79, 50, 78, 61, 73, 57, 72, 58, 64,
		77, 58, 79, 68, 70, 58, 65, 60, 52, 81, 62, 61, 91, 64, 75, 79, 70, 62, 80, 68, 71, 82, 57,
		59, 44, 68, 82, 67, 68, 79, 61, 74, 53, 73, 76, 80, 79, 80, 68, 94, 87, 77, 77, 70, 67, 64,
		57, 75, 79, 60, 91, 64, 73, 73, 78, 85, 72, 68, 75, 70, 75, 71, 58, 69, 79, 82, 76, 61, 61,
		75, 71, 69, 85, 75},
		caseness3Disease = new double[] {
		0, 1, 1, 0, 1, 0, 1, 0, 0, -1, 1, 0, 1, 0, 0, 1, 1, 0, 1, -1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, -1, -1, 1, 1, 1, 0, 0,
		1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, -1, -1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0,
		0, 1, 0, 1, -1, -1, 1, 1, 1, 1, 0, 0, 1, 1, 0, -1, -1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0,
		1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, -1, -1, 0, 1, 0, 1, 0, 1,
		1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, -1, 1, 0, 1, 0, 1, 0, 1, -1, -1, 1,
		-1, 1, -1, 0, 1, 0, 1, -1, -1, -1, -1, 0, -1, 1, 0, 0, 1, 0, 1, -1, -1, 0, 1, 1, 0, -1, -1, 1, 0, 0, 1, 1, 0, 0, 1, -1, -1, -1, -1, 1, 0,
		0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, -1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, -1, 0, 0, 0, -1, -1, 1, 0, -1, 1, 1, 1, 1, 0, 1, 0, 1,
		0, 1, 0, -1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, -1, 1, 1, 1, 0, 1, 1, 0, 1, 0, -1, 0, 0, 1, 1, 0, 0, 1, -1, -1, 0, 1, -1, -1, 1, 0,
		0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, -1, -1, 0, 1, -1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, -1, 1, -1, 0, 1, 0, 1, 0, 1, 0, -1, 1, -1,
		1, 0, 1, 0, 1, -1, 0, 0, 0, -1, -1, -1, -1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, -1, -1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, -1, -1, 1, 0,
		0, 1, 1, 0, -1, -1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, -1, 1, -1, 0, 0, 1};

	// A bitset indicating which genes are annotated (=1) 
	static final long[]
		geneAnnotationBitset = new long[] {
		0xFFFFFFFFFFFFFFFFL, 0xFFF7FFFFFFFFFFFFL, 0xFFFFFFFFFFFBFFFFL, 0xFFFFFFFFFFFFFEFFL, 0xDFFFFFFFFFFF7FFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFBL, 0xFFFFFFFFFFFFFFFDL, 0xDFFFFDFFFFFFFFFFL, 0xFFFFEFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFAFFFFFFFFFFDFFL, 0xFFFFFFF7FFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xF7FFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFEFFFL, 0xF7FFFFFFFF7EFFFFL, 0xFFBFFFFFFFFFFFFFL, 0xFFFFFFFEFFDFFFFFL,
		0xFFF77FFFFFF7FFFFL, 0xFFFFBFFFFFFFFFFFL, 0xFFFFFFFFFFE7FFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFDFFFFFFFL, 0xFFFFF7EFFFFFFFFFL,
		0xFEFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFBFFFL, 0xFFFF7FFFFFFFFFFFL, 0xFFFFFEFFFFFBFFFFL,
		0xFFDFBFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFEFL, 0xFFFFFFFFFFFFFFFFL, 0xFDFDFFFFFFFFFFFFL,
		0xFFFFFFFDFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFDFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFF7FFFFFFFFFEFFL, 0xFFF7FFFFFFFFFFFFL,
		0xFFFFFEFFFFFBFFFFL, 0xFFFFFFFFFFFFFF7FL, 0xFFFFFFFFFFFFFF7FL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFF7FFFFFL, 0xFFEFFFDF7FFFFFEFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFF7EFFFFFF7FL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFBFFFFL,
		0xFFFF7FFFFFF7FFFFL, 0xFFDFFFFFFFFFFFEFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFEFFFL, 0xFFFFFF6FFFFFFFFFL, 0xFF7FFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xDFFFFFFFFFFF7FFFL,
		0xDFDFFFFFFFFDFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFBFFFFFFFF7FFFFFL, 0xBFFFFFFEFFFFFFFFL, 0xFFFFFFFFFFFFFBFFL, 0xFFFFFFEFFFFFFFFFL,
		0xFFFFFFFFFDFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFDL, 0xFF7FFFFFFFFFFFFFL, 0xFFFFFFFFF7DFFFF7L,
		0xFF7FFFFDFFFFFFFFL, 0xFFFFEFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFDL, 0xFFFFFFFFFFFFFFFFL, 0xEFFFFFFFFFFFFFFFL, 0xFFFFFFFFBFFFFFFFL, 0xFFFFFFFFEFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFDFFFFL, 0xFFFFFFFFFFF7FFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFDFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFF7FFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFDFFFFFFFL, 0xFFFFFFFFFFFFFEFFL, 0xFFFFFFFFFFFDEFFFL,
		0xEFFFFFFFFFFFFFFFL, 0xFFFFFBFFFFFBDFFFL, 0xF7FFFFFFFFFFEFFFL, 0xFFFFFFBFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFBFFFFFFFFFFFL,
		0xFFFFFFFDFBFFFFFFL, 0xFFF7FFFFFFFFDFFFL, 0xFFFBEFFFFFFFFFFFL, 0xA7FFFFFFFDFFFFFFL, 0xFFFFFFFBFFFFFFFFL, 0xDFFFFFFFFFFFFFFFL,
		0x7BFF7FFEFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFEFFFFFFFFBFFL,
		0xFFFFFFFFFFEFF7FFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFF7FFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFDFFFFFFFF7FFL,
		0xFFFFBFFFF7FFEFFDL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFDFFFFFFL, 0xFFFFFFFFFFFDFFBFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFDFFFFFFFFFFL, 0xFFFFFFFFEFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFDFFFFFFFFL, 0xFFFFF7FFFFBFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFEFFEFFFFFL, 0xFFFFFFFBFFFFFFFDL, 0xFFFFFFBFFFFFFFBFL, 0xFEFFFFFFFFFFFFFFL, 0xFBFFFFFFFFFFFFFFL, 0xFBDEFFFFFFF7FFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFF7FFDFFFFFFL, 0xFFFFFFD7FFFFFFFFL, 0xF7FFFFFFFFFFFFFFL, 0xFFFDFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFDFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFDFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFDFFFFFBFL, 0xFFFFFFFDBFFFFFFBL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFEFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFDFFFFFFFFFFL, 0xFFFFFFDFFFFFFFFFL, 0xFFFFFFFFFFF7FFFFL,
		0xFEFDFDFFFFFFFFFFL, 0x7FFFFFFFFFFFFFFFL, 0xFDFFEDFFFFFFFFFFL, 0xFFFBFFFFFFFFFFFFL, 0xFFBFFFFFFFFFFFFFL, 0xFFFDFFFFFFBFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFDFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFEFF7FFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFEFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFBDFFL, 0xFFFDDFFFDFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFEFFEFDFL, 0xFDFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFF7FFL,
		0xFFFFFFFF7FFFFFFFL, 0xFF7FFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFF7FFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFF7FFFFFFFFFFL,
		0xFFFBFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFF7FFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xDFFFFFFFFFFFFFFFL,
		0xF7FFFFFFFFFBFFFFL, 0xFEFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFEL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xDFFFEFFFFFFFFFFFL, 0xFFFFF7FFFFFFFFBFL, 0xFFFFFF77FFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFDFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFF7FL, 0xFFFFFFFFFFEFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFDFFFFFFFFFFBFL, 0xFFFFFFBFFFFFFFFFL,
		0xFFFFFBFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFBFFFFFFFFL, 0xFFFFFFFFFFFFFFF7L, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFBFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFDFDFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFF7L, 0x7FFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL,
		0xFFFFFFFFFFFF7FEFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFDFFFFFFFFFFFFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFBFFFDFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFEFL, 0xFFFFFFFFF7FFFFFFL, 0xFFF3FFFFFFFFFFEFL,
		0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFEFFFFFFFFFL, 0xFFFFFFFFFFFFF7EFL,
		0xFFFFFFFFBFFFFFFFL, 0xFFFFFEFBFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFFFFF7FFFFFFL, 0xFF7FFFFFFFFFFFFFL, 0xFFFFFFFFFFDFFFFFL,
		0xBFFFFFFFFFFFFFFFL, 0xFFFFFBFFFFFFFFFFL, 0xFFFFFFFFFFFFFFFFL, 0xFFFFFF8000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L, 0x0000000000000000L,
		0x0000000000000000L, 0x0000000000000000L
		};
	private double[][]
		pvals, data;
	private double[][] result;

	private double designMat1[][], designMat2[][], y[], yy[][] = new double[1][], covars[][];

	public QGeneRegNetCC(double[][] dt, double[][] rs, double[][] ps, double[][] cov) {
		data = dt;
		result = rs;
		pvals = ps;
		covars = cov;
		int
			len = covars[0].length,
			numCovars = covars.length;
		designMat1 = new double[len][2*numCovars + 2];
		designMat2 = new double[len][2];
		yy[0] = y = new double[len];
		for (int i = 0; i < len; i++) {
			designMat1[i][0] = designMat2[i][0] = 1;
			for (int j = 0; j < numCovars; j++)
				designMat1[i][j + 2] = covars[j][i];
		}
	}

	@Override
	public boolean execute(int iterNo) {
		int
			numGenes = data.length,
			len = covars[0].length,
			numCovars = covars.length;
		result[iterNo][iterNo] = 1;
		pvals[iterNo][iterNo] = 1;
		QGLMResult glm1 = null, glm2 = null;
		double[] x = data[iterNo];
		for (int i = 0; i < len; i++) {
			designMat1[i][1] = designMat2[i][1] = x[i];
			double xi = x[i];
			for (int j = 0; j < numCovars; j++)
				designMat1[i][j + 3] = covars[j][i] * xi;
			//designMat1[i][3] = caseness[i] * x[i];
		}
		for (int geneNo = 0; geneNo < numGenes; geneNo++) {
			if (geneNo == iterNo) continue;
			x = data[geneNo];
			for (int i = 0; i < len; i++)
				y[i] = x[i];
			if (glm1 == null || glm2 == null) {
				glm1 = QGLM.multipleRegression(designMat1, y);
				glm2 = QGLM.multipleRegression(designMat2, y);
			} else {
				glm1 = glm1.recalculate(yy);
				glm2 = glm2.recalculate(yy);
			}
			glm1.calcMoreStats();
			glm2.calcMoreStats();
			QGLMResultBasic glm_diff = glm1.calcDiffStatistic(glm2);
			result[geneNo][iterNo] = (float) glm_diff.mRSq;
			pvals[geneNo][iterNo] = 1 - calcPOfF(glm_diff.mF, glm_diff.mDFReg, glm_diff.mDFErr); 
		}
		if (iterNo % 100 == 0) System.out.println(iterNo);
		return true;
	}

	@Override
	public QGeneRegNetCC clone() {
		return new QGeneRegNetCC(data, result, pvals, covars);
	}

	static final int[] countGraduation(double[][] result, double[] graduation) {
		int len = graduation.length, n = result.length;
		int[] counts = new int[graduation.length];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) continue;
				double val = result[i][j];
				for (int k = 0; k < len; k++)
					if (val >= graduation[k]) counts[k]++;
			}
		}
		return counts;
	}

	static final double[][] filterDataByColumn(double[][] data, BitSet bitset) {
		int
			newCols = bitset.cardinality(),
			oldCols = data[0].length,
			nrows = data.length;
		double[][] newData = new double[nrows][newCols];
		for (int rowNo = 0; rowNo < nrows; rowNo++) {
			double[]
				dataRow = data[rowNo],
				newDataRow = newData[rowNo];
			int newColNo = 0;
			for (int colNo = 0; colNo < oldCols; colNo++) {
				if (bitset.get(colNo))
					newDataRow[newColNo++] = dataRow[colNo];
			}
		}
		return newData;
	}

	static final double[][] filterData(double[][] data, BitSet rowBitset, BitSet colBitset) {
		int
			newCols = colBitset.cardinality(),
			newRows = rowBitset.cardinality(),
			newRowNo = 0,
			oldCols = data[0].length,
			nrows = data.length;
		double[][] newData = new double[newRows][newCols];
		for (int rowNo = 0; rowNo < nrows; rowNo++) {
			if (rowBitset.get(rowNo)) {
				double[]
					dataRow = data[rowNo],
					newDataRow = newData[newRowNo];
				int newColNo = 0;
				for (int colNo = 0; colNo < oldCols; colNo++) {
					if (colBitset.get(colNo))
						newDataRow[newColNo++] = dataRow[colNo];
				}
				newRowNo++;
			}
		}
		return newData;
	}

	static final BitSet getAnnotationBitset(long[] arr) {
		int n = arr.length;
		BitSet bs = new BitSet(n * 64);
		bs.set(0, n * 64, true);
		for (int i = 0; i < n; i++) {
			long bits = arr[i];
			for (int j = 0; j < 64; j++) {
				if ((bits & 0x8000000000000000L) == 0L)
					bs.clear(i*64+j);
				bits = bits << 1;
			}
		}
		return bs;
	}

	public static void main(String[] args) {
		//annotateNet("/Users/joehanesr/Desktop/cc-resid-revised1-gene.txt",
		//		"/Users/joehanesr/Desktop/cc-resid-revised1-gene-diff.txt",
		//		"/Users/joehanesr/Desktop/annotation-gene.txt");
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

			long mem1, mem2, time1, time2; // Benchmark purposes

			mem1 = QSystemUtils.usedMemoryAfterGC();
			QTableData dataTable = readDelimitedFileAsTableData(fileName, cCommaDelimiter, "#", true, true); // num genes x num samples
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));

			BitSet
				annotBitset = getAnnotationBitset(geneAnnotationBitset),
				bitset = new BitSet(caseness.length);
			bitset.set(0, caseness.length, true);
			for (int i = 0; i < caseness3Disease.length; i++)
				bitset.set(i, caseness3Disease[i] >= 0);
			//System.out.println(annotBitset.cardinality());
			double[][] data = dataTable.extractData();
			if (data[0].length != bitset.cardinality())
				data = filterDataByColumn(data, bitset);
			int
				numSubjects = data[0].length,
				numGenes = data.length;
			System.out.println("Data loaded = " + numGenes + " x " + numSubjects);
			double[][]
				rsqs = new double[numGenes][numGenes],
				pvals = new double[numGenes][numGenes];
			double[] rsq_graduation = new double[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};
			int[] rsq_counts;
			mem2 = QSystemUtils.usedMemoryAfterGC();
			System.out.println("Memory used = " + (mem2 - mem1));

			System.out.println("Case-control");
			time1 = System.currentTimeMillis();
			executeSimpleParallelJob(new QGeneRegNetCC(data, rsqs, pvals, filterDataByColumn(new double[][] {caseness}, bitset)), numGenes);
			time2 = System.currentTimeMillis();
			System.out.println("Regression time = " + (time2 - time1));
			time1 = System.currentTimeMillis();
			QFileUtils.writeText(rsqs, baseName + "-rsq-cc3." + ext);
			QFileUtils.writeText(pvals, baseName + "-pvals-cc3." + ext);
			time2 = System.currentTimeMillis();
			System.out.println("Saving time = " + (time2 - time1));
			rsq_counts = countGraduation(rsqs, rsq_graduation);
			System.out.println(QStringUtils.toString(rsq_graduation));
			System.out.println(QStringUtils.toString(rsq_counts));
			rsqs = null;
			int numAnnotatedGenes = annotBitset.cardinality();
			double[][] fdr = new double[numAnnotatedGenes][numAnnotatedGenes];
			pvals = filterData(pvals, annotBitset, annotBitset);
			QBenjaminiHochbergFDRPlugin.calculateFDR(pvals, fdr);
			QFileUtils.writeText(fdr, baseName + "-fdrnamedonly-cc3." + ext);

			System.exit(0);
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
