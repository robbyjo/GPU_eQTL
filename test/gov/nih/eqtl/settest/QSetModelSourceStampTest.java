/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl.settest;

import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QSetModelSourceStampTest {
    private void replace(Path path,String text) throws Exception {
        var stamp=Files.getLastModifiedTime(path);Path temporary=Files.createTempFile(path.getParent(),"replacement",".txt");
        Files.writeString(temporary,text);Files.setLastModifiedTime(temporary,stamp);
        Files.move(temporary,path,StandardCopyOption.REPLACE_EXISTING);
    }
    @Test void unchangedBytesRemainPortableAcrossFileIdentityChanges(@TempDir Path root) throws Exception {
        Path path=root.resolve("pedigree.tsv");String text="ID\tSIRE\tDAM\na\t0\t0\n";Files.writeString(path,text);
        var signer=new QSetModelSourceStamp();String before=signer.capture(path);
        replace(path,text);
        assertEquals(before,signer.capture(path));assertEquals(before,new QSetModelSourceStamp().capture(path));
        assertTrue(before.contains(":sha256="));
    }
    @Test void sameSizeAndTimeReplacementStillDetectsDifferentContent(@TempDir Path root) throws Exception {
        Path path=root.resolve("pedigree.tsv");Files.writeString(path,"ID\tSIRE\tDAM\na\t0\t0\n");
        var signer=new QSetModelSourceStamp();String before=signer.capture(path);
        replace(path,"ID\tSIRE\tDAM\nb\t0\t0\n");
        assertNotEquals(before,signer.capture(path));assertNotEquals(before,new QSetModelSourceStamp().capture(path));
    }
    @Test void pathsAndTimestampsStillBindCheckpointIdentity(@TempDir Path root) throws Exception {
        Path a=root.resolve("a"),b=root.resolve("b");Files.writeString(a,"same");Files.copy(a,b,StandardCopyOption.COPY_ATTRIBUTES);
        var signer=new QSetModelSourceStamp();String before=signer.capture(a);
        assertNotEquals(before,signer.capture(b));
        Files.setLastModifiedTime(a,java.nio.file.attribute.FileTime.fromMillis(Files.getLastModifiedTime(a).toMillis()+2000));
        assertNotEquals(before,signer.capture(a));
    }
}
