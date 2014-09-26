package gov.nih.exon.coulter;

import java.io.IOException;

import gov.nih.utils.QStringUtils;

import static gov.nih.utils.QFileUtils.writeText;
import static gov.nih.utils.QStringUtils.sTab;

/**
 * 
 * @author Roby Joehanes
 *
 */
public class QCoulterDataUtils {
	static final void saveResult(String filename, double[][] result, String[] rowNames) throws IOException {
		String str = QStringUtils.toDelimitedString(result, rowNames, null, sTab);
		writeText(str, filename);
	}

	static final void saveResult(String filename, int[][] result, String[] rowNames) throws IOException {
		String str = QStringUtils.toDelimitedString(result, rowNames, null, sTab);
		writeText(str, filename);
	}
}
