/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.io.IOException;
import java.nio.file.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
class QPooledCisIndexTest {
    @Test void inclusiveWindowsPreserveSourceOrderAndChromosomes(@TempDir Path root)throws Exception {
        Path file=root.resolve("positions.tsv");Files.writeString(file,"ID\tCHROM\tPOS\nt.1\tchr1\t105\n");
        var index=new QPooledScoreRunner.CisIndex(file,new String[]{"1:110:A:G","2:105:A:G","chr1:100:A:C","1:111:A:T","1:105:A:C"},5);
        assertArrayEquals(new int[]{0,2,4},index.rows("t.1"));
        assertThrows(IOException.class,()->index.rows("t"));
    }
    @Test void windowsDoNotOverflowAtCoordinateLimits(@TempDir Path root)throws Exception {
        Path file=root.resolve("positions.tsv");Files.writeString(file,"ID\tCHROM\tPOS\nt\t1\t"+Long.MAX_VALUE+"\n");
        var index=new QPooledScoreRunner.CisIndex(file,new String[]{"1:1:A:G","1:"+Long.MAX_VALUE+":A:C"},Long.MAX_VALUE);
        assertArrayEquals(new int[]{0,1},index.rows("t"));
    }
    @Test void malformedOrAmbiguousCoordinatesFail(@TempDir Path root)throws Exception {
        Path file=root.resolve("positions.tsv");
        for(String text:new String[]{"ID,CHROM,POS\nt,1,1\n","ID\tCHROM\tPOS\nt\t1\t0\n","ID\tCHROM\tPOS\nt\tchr\t1\n","ID\tCHROM\tPOS\nt\t1\t1\nt\t2\t1\n"}){
            Files.writeString(file,text);assertThrows(IOException.class,()->new QPooledScoreRunner.CisIndex(file,new String[]{"1:1:A:G"},0));
        }
        Files.writeString(file,"ID\tCHROM\tPOS\nt\t1\t1\n");
        assertThrows(IOException.class,()->new QPooledScoreRunner.CisIndex(file,new String[]{"rs1"},0));
        assertThrows(IllegalArgumentException.class,()->new QPooledScoreRunner.CisIndex(file,new String[]{"1:1:A:G"},-1));
    }
}
