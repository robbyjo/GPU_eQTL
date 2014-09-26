package gov.nih.tools.old;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import com.csvreader.CsvReader;

/**
 * @author Roby Joehanes
 */
public class Validation {
	static ExecutorService threadPool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() + 1);
	static final long unknownPos = Long.MIN_VALUE;
	static String[] SNPChrs, TxChrs;
	static long[] SNPPos, TxStart, TxStop;
	static boolean[] isCis;
	static ResultParams params1, params2;
	static int numLinesOnFile1 = -1, counter = 0;
	static int[] distances;
	static PrintWriter writer;

	static final ResultParams parseParams(String[] params) {
		return new ResultParams(Integer.parseInt(params[0]),
			Integer.parseInt(params[1]),
			Integer.parseInt(params[2]),
			Integer.parseInt(params[3]),
			Integer.parseInt(params[4]),
			params[5].length() == 0 ? -1 : Integer.parseInt(params[5]),
			params[6]);
	}

	public static final int countNumLines(String filename) throws IOException {
	    InputStream is = new BufferedInputStream(new FileInputStream(filename));
	    try {
	        byte[] c = new byte[1024];
	        int count = 0;
	        int readChars = 0;
	        boolean empty = true;
	        while ((readChars = is.read(c)) != -1) {
	            empty = false;
	            for (int i = 0; i < readChars; ++i) {
	                if (c[i] == '\n') {
	                    ++count;
	                }
	            }
	        }
	        return (count == 0 && !empty) ? 1 : count;
	    } finally {
	        is.close();
	    }
	}

	static final int[] parseDistances(String[] dist) {
		int[] d = new int[dist.length];
		for (int i = 0; i < dist.length; i++)
			d[i] = Integer.parseInt(dist[i]);
		return d;
	}

	public static final void main(String[] args) {
		// 1000000,100000 Cis_1M 15,1,2,19,20,21,MRCA_5percentFDR-hg19.txt  18,2,3,11,13,14,eqtl-gene-1000g-result-p1e-4-na33-irsq3-maf01-fdr.txt
		distances = parseDistances(args[0].split(","));
		String cisToken = args[1];
		params1 = parseParams(args[2].split(","));
		params2 = parseParams(args[3].split(","));

		System.out.println("Processing file: " + params1.filename);
		try {
			numLinesOnFile1 = countNumLines(params1.filename);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		SNPChrs = new String[numLinesOnFile1];
		TxChrs = new String[numLinesOnFile1];
		SNPPos = new long[numLinesOnFile1];
		TxStart = new long[numLinesOnFile1];
		TxStop = params1.idxTxStop < 0 ? null : new long[numLinesOnFile1];
		isCis = new boolean[numLinesOnFile1];
		CsvReader reader = null;
		long lineNo = 0;
		try {
			reader = new CsvReader(params1.filename);
			reader.setTrimWhitespace(true);
			while(reader.readRecord()) {
				lineNo++;
				String[] tokens = reader.getValues();
				if (lineNo == 1) {
					System.out.println("Reading from columns: " + tokens[params1.idxIsCis] +
						"," + tokens[params1.idxSNPChr] + "," + tokens[params1.idxSNPPos] +
						"," + tokens[params1.idxTxChr] + "," + tokens[params1.idxTxStart] +
						(params1.idxTxStop > 0 ? "," + tokens[params1.idxTxStop] : ""));
					continue; // Skip first line
				}
				int idx = ((int) lineNo) -1;
				isCis[idx] = tokens[params1.idxIsCis].equals(cisToken);
				SNPChrs[idx] = tokens[params1.idxSNPChr];
				SNPPos[idx] = Long.parseLong(tokens[params1.idxSNPPos]);
				TxChrs[idx] = tokens[params1.idxTxChr];
				try {
					TxStart[idx] = Long.parseLong(tokens[params1.idxTxStart]);
				} catch (Exception ex) { TxStart[idx] = unknownPos; }
				if (params1.idxTxStop > 0)
					try {
						TxStop[idx] = Long.parseLong(tokens[params1.idxTxStop]);
					} catch (Exception ex) { TxStop[idx] = unknownPos; }
			}
			System.out.println(lineNo + " lines were read from file " + params1.filename);
			reader.close(); reader = null;
		} catch (Exception e) {
			e.printStackTrace();
			if (reader != null)
				try {
					reader.close();
				} catch (Exception e1) {}
		}

		System.out.println("Processing file: " + params2.filename);
		String outputFilename = params2.filename + "-" + params1.filename + ".out.txt";
		reader = null;
		writer = null;
		lineNo = counter = 0;
		try {
			reader = new CsvReader(params2.filename);
			reader.setTrimWhitespace(true);
			writer = new PrintWriter(outputFilename);
			while(reader.readRecord()) {
				String line = reader.getRawRecord();
				String[] tokens = reader.getValues();
				if (lineNo == 0) {
					lineNo++;
					writer.println(line + ",IsValidated");
					System.out.println("Reading from columns: " + tokens[params2.idxIsCis] +
						"," + tokens[params2.idxSNPChr] + "," + tokens[params2.idxSNPPos] +
						"," + tokens[params2.idxTxChr] + "," + tokens[params2.idxTxStart] +
						(params2.idxTxStop > 0 ? tokens[params2.idxTxStop] : ""));
					continue; // Skip first line
				}
				lineNo++;
				int curIsCis = Integer.parseInt(tokens[params2.idxIsCis]);
				String curSNPChr = tokens[params2.idxSNPChr];
				long curSNPPos = Long.parseLong(tokens[params2.idxSNPPos]);
				String curTxChr = tokens[params2.idxTxChr];
				long curTxPos = Long.parseLong(tokens[params2.idxTxStart]);
				if (lineNo % 100000 == 0) System.out.println(lineNo/1000);
				if (lineNo - counter > 1000000) Thread.sleep(2000);
				threadPool.execute(new ValidationThread(curIsCis, curSNPChr, curSNPPos, curTxChr, curTxPos, line));
			}
			System.out.println("Waiting for all threads to finish.");
			threadPool.shutdown();
			try {
				threadPool.awaitTermination(3, TimeUnit.DAYS);
			} catch (InterruptedException exx) {}
			System.out.println(lineNo + " lines were read from file " + params2.filename);
			reader.close(); reader = null;
			writer.close(); writer = null;
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		} finally {
			if (reader != null)
				try {
					reader.close();
				} catch (Exception e) {}
			if (writer != null)
				try {
					writer.close();
				} catch (Exception e) {}
		}

		System.exit(0);
	}
}

class ResultParams {
	int idxSNPChr, idxSNPPos, idxTxChr, idxTxStart, idxTxStop, idxIsCis;
	String filename;
	ResultParams() {}
	ResultParams(int idxIsCis, int idxSNPChr, int idxSNPPos, int idxTxChr, int idxTxStart, int idxTxStop, String filename) {
		this.idxIsCis = idxIsCis;
		this.idxSNPChr = idxSNPChr;
		this.idxSNPPos = idxSNPPos;
		this.idxTxChr = idxTxChr;
		this.idxTxStart = idxTxStart;
		this.idxTxStop = idxTxStop;
		this.filename = filename;
	}
}

class ValidationThread implements Runnable {
	int curIsCis;
	String curSNPChr, curTxChr;
	long curSNPPos, curTxPos;
	String line;

	ValidationThread(int curIsCis, String curSNPChr, long curSNPPos, String curTxChr, long curTxPos, String line) {
		this.curIsCis = curIsCis;
		this.curSNPChr = curSNPChr;
		this.curSNPPos = curSNPPos;
		this.curTxChr = curTxChr;
		this.curTxPos = curTxPos;
		this.line = line;
	}

	public void run() {
		boolean validated = false;
		for (int i = 0; i < Validation.numLinesOnFile1; i++) {
			if (((curIsCis == 1) == Validation.isCis[i]) &&
				Math.abs(curSNPPos - Validation.SNPPos[i]) <= Validation.distances[curIsCis] &&
				Math.abs(curTxPos - Validation.TxStart[i]) <= Validation.distances[curIsCis] &&
				curSNPChr.equals(Validation.SNPChrs[i]) &&
				curTxChr.equals(Validation.TxChrs[i])) {
				validated = true;
				break;
			}
		}
		synchronized(Validation.writer) {
			Validation.writer.println(line + (validated ? ",1" : ",0"));
			Validation.writer.flush();
			Validation.counter++;
		}
	}
}