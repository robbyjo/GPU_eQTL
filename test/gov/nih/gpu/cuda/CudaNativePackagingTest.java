/* Copyright 2026 Roby Joehanes; GNU GPL version 3. */
package gov.nih.gpu.cuda;

import static org.junit.jupiter.api.Assertions.*;
import java.io.InputStream;
import org.junit.jupiter.api.Test;

/** No GPU/runtime required: cross-platform JNI payloads must survive dependency resolution. */
class CudaNativePackagingTest {
    @Test void linuxAndWindowsJniBindingsArePresentRegardlessOfBuildHost() throws Exception {
        for (String name : new String[]{"JCudaDriver", "JCudaRuntime", "JNvrtc",
                "JNvPTXCompiler", "JCublas", "JCublas2"}) {
            check("/lib/lib" + name + "-12.6.0-linux-x86_64.so",
                    new byte[]{0x7f, 'E', 'L', 'F'});
            check("/lib/" + name + "-12.6.0-windows-x86_64.dll",
                    new byte[]{'M', 'Z'});
        }
    }

    private static void check(String resource, byte[] magic) throws Exception {
        try (InputStream input = CudaNativePackagingTest.class.getResourceAsStream(resource)) {
            assertNotNull(input, "Missing CUDA JNI resource: " + resource);
            assertArrayEquals(magic, input.readNBytes(magic.length), resource);
        }
    }
}
