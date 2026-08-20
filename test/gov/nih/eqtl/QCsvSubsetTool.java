/* Copyright 2011 Roby Joehanes; GNU GPL version 3.0. */
package gov.nih.eqtl;

import static gov.nih.utils.QStringUtils.cCommonDelimiter;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;

import com.csvreader.CsvReader;
import com.csvreader.CsvWriter;

import gov.nih.eqtl.io.QDelimitedMatrixSource;

/** Test-only utility for creating an untracked row subset without changing sample columns. */
public final class QCsvSubsetTool {
    private QCsvSubsetTool() { }

    public static void main(String[] args) throws IOException {
        if (args.length != 3)
            throw new IllegalArgumentException("Usage: input.csv output.csv number-of-data-rows");
        Path input = Path.of(args[0]);
        Path output = Path.of(args[1]);
        int maximumRows = Integer.parseInt(args[2]);
        Files.createDirectories(output.toAbsolutePath().getParent());
        CsvReader reader = null;
        CsvWriter writer = null;
        try {
            reader = new CsvReader(QDelimitedMatrixSource.openReader(input));
            reader.setDelimiter(cCommonDelimiter);
            reader.setTrimWhitespace(true);
            writer = new CsvWriter(output.toString(), ',', Charset.defaultCharset());
            int records = 0;
            while (records <= maximumRows && reader.readRecord()) {
                writer.writeRecord(reader.getValues(), true);
                records++;
            }
            if (records <= 1)
                throw new IOException("Input contains no data rows: " + input);
        } finally {
            if (reader != null)
                reader.close();
            if (writer != null)
                writer.close();
        }
    }
}
