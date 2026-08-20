# GPU eQTL

GPU-accelerated Java software for eQTL and related QTL analyses. This codebase originated in 2011–2013 and is being modernized incrementally while preserving its double-precision statistical behavior.

## Current status

The application targets Java 17 and now has selectable GPU backends behind a vendor-neutral API:

| Setting | Implementation | Intended hardware |
| --- | --- | --- |
| `auto` (default) | CUDA-first discovery plus OpenCL fallback | Mixed machines; avoids running the same NVIDIA card twice |
| `cuda` | JCuda 12.6.0 and full-double-precision cuBLAS DGEMM | NVIDIA GPUs |
| `opencl` | JOCL 2.0.6 and the production OpenCL C kernel | NVIDIA, AMD, or Intel GPUs whose driver exposes OpenCL with FP64 |

Automatic discovery uses CUDA for a usable NVIDIA device and also includes distinct FP64 OpenCL devices from other vendors. If CUDA cannot initialize, NVIDIA OpenCL remains available as a fallback. A native HIP/ROCm backend is not included yet; AMD cards currently use JOCL/OpenCL.

The real-valued eQTL calculation remains double precision and has been validated through both CUDA and OpenCL on an NVIDIA GeForce RTX 2080 against the same CPU reference. Most client Intel Iris Xe GPUs do not provide the native FP64 capability required by this analysis, so they may be reported by `--printgpuinfo` but are excluded from execution. The legacy categorical-SNP path remains disabled; the CUDA backend currently implements `eqtlReal` only.

The run configuration is still the legacy INI format. Matrix streaming, input-ID alignment, categorical covariate encoding, SNP interactions, and forward selection remain on the modernization roadmap in `AGENTS.md` and `SESSION_HISTORY.md`.

## Build and test

Install a Java 17 or newer JDK. Maven itself does not have to be installed because the repository includes a checksum-pinned wrapper.

On Windows:

```powershell
.\mvnw.cmd clean package
```

On Linux or macOS:

```sh
./mvnw clean package
```

The runnable jar is `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar`. It includes the Java dependencies and the JCuda/JOCL JNI bindings; it does not bundle GPU drivers or vendor runtimes. The CUDA backend needs a compatible NVIDIA driver plus CUDA runtime/cuBLAS installation. The OpenCL backend needs the vendor's OpenCL ICD.

Hardware-independent tests always run. CUDA and OpenCL numerical tests skip cleanly when their respective FP64 runtime/device is unavailable.

## Eclipse

Use an Eclipse installation with Maven Integration for Eclipse (m2e) and a Java 17 or newer JDK. Import this directory with **File → Import → Maven → Existing Maven Projects**.

If the project was already open before the Maven metadata was added, select the project and run **Maven → Update Project…** (`Alt+F5`), enable **Force Update of Snapshots/Releases**, and then run **Project → Clean…**. Eclipse should show **Maven Dependencies** in the project tree and use the `JavaSE-17` execution environment.

The custom multi-delimiter JavaCSV jar remains attached directly with its source archive; JCuda, JOCL, JDistlib, Commons Compress, and JUnit are supplied by Maven.

## Inspect and select GPUs

Automatic selection is the default:

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo
```

Force CUDA or OpenCL in PowerShell (quote the `-D` argument):

```powershell
java '-Deqtl.gpu.backend=cuda' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo
java '-Deqtl.gpu.backend=opencl' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo
```

In POSIX shells, the same property does not need PowerShell quoting:

```sh
java -Deqtl.gpu.backend=cuda -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo
```

Recent Java runtimes may print a warning when a JNI binding loads. Add `--enable-native-access=ALL-UNNAMED` before `-jar` to suppress it. Java 17 remains the compilation target.

A manual, transfer-inclusive comparison is available after the build. Its results are indicative rather than a substitute for benchmarking a representative analysis:

```powershell
java --enable-native-access=ALL-UNNAMED -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.gpu.GpuBackendBenchmark
```

## Run an analysis

With automatic backend selection:

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar path\to\analysis.ini
```

To force a backend, place the quoted `-Deqtl.gpu.backend=...` option before `-jar`, as shown above. The existing INI keys and input formats are unchanged. The former `library_path` key is no longer used; runtime discovery is handled by the CUDA/OpenCL backends.

See `AGENTS.md` for contributor constraints and the ordered roadmap. See `SESSION_HISTORY.md` for timestamped changes, verification, and known limitations.