/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.eqtl;
import static org.junit.jupiter.api.Assertions.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class QPooledNullProcessTest {
    @Test void chromosomeProcessesSerializeOneNullFit(@TempDir Path root)throws Exception {
        List<Process> processes=new ArrayList<>();
        try {
            for(String name:List.of("a","b")) {
                List<String> args=new ArrayList<>(List.of(Path.of(System.getProperty("java.home"),"bin","java").toString(),"--enable-native-access=ALL-UNNAMED"));
                for(String property:List.of("java.io.tmpdir","org.bytedeco.javacpp.cachedir","eqtl.cpu.threads"))if(System.getProperty(property)!=null)args.add("-D"+property+"="+System.getProperty(property));
                args.addAll(List.of("-cp",System.getProperty("java.class.path"),Child.class.getName(),root.toString(),name));
                processes.add(new ProcessBuilder(args).redirectErrorStream(true).redirectOutput(root.resolve(name+".log").toFile()).start());
            }
            long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(45);
            while(!(Files.exists(root.resolve("a.ready"))&&Files.exists(root.resolve("b.ready")))) {
                assertTrue(System.nanoTime()<deadline,"Child startup timed out");
                for(Process child:processes)assertTrue(child.isAlive(),"Child failed before null fit");Thread.sleep(10);
            }
            Files.writeString(root.resolve("go"),"go");
            for(Process child:processes){assertTrue(child.waitFor(45,TimeUnit.SECONDS));assertEquals(0,child.exitValue());}
            Set<String> results=Set.of(Files.readString(root.resolve("a.result")),Files.readString(root.resolve("b.result")));
            assertEquals(Set.of("computed:1","cache:0"),results);
        }finally{for(Process child:processes)if(child.isAlive()){child.destroyForcibly();child.waitFor();}}
    }
    public static final class Child {
        public static void main(String[] args)throws Exception {
            Path root=Path.of(args[0]);String name=args[1];var d=QPooledPedigreeModelTest.data();var model=QPooledPedigreeModelTest.model(d);
            var geometry=model.geometry(java.util.stream.IntStream.range(0,96).toArray());Files.writeString(root.resolve(name+".ready"),"ready");
            long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(45);while(!Files.exists(root.resolve("go"))){if(System.nanoTime()>deadline)throw new IllegalStateException("Gate timeout");Thread.sleep(10);}
            var fit=QPooledNullCache.fit(geometry,d.y(),"same-trait",root.resolve("shared-nulls"),false);
            Files.writeString(root.resolve(name+".result"),(fit.reused?"cache":"computed")+":"+model.fitExecutions());
        }
    }
}
