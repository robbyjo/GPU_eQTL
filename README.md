# GPU eQTL

GPU-accelerated Java software for eQTL and related QTL analyses. This codebase originated in 2011–2013 and is being modernized incrementally while preserving its double-precision statistical behavior.

## Current status

The application now targets Java 17 and uses JOCL 2.0.6 through a backend-neutral GPU layer. The current backend runs OpenCL C kernels against the OpenCL runtime installed with the GPU driver. It has been validated on an NVIDIA RTX 2080; AMD and Intel devices require vendor drivers that expose an OpenCL GPU with compiler and FP64 support.

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

The runnable jar is `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar`. Tests that do not require hardware always run; the numerical production-kernel test skips cleanly if no suitable FP64 OpenCL GPU is available.

## Eclipse

Use an Eclipse installation with Maven Integration for Eclipse (m2e) and a Java 17 or newer JDK. Import this directory with **File → Import → Maven → Existing Maven Projects**.

If the project was already open before the Maven metadata was added, select the project and run **Maven → Update Project…** (`Alt+F5`), enable **Force Update of Snapshots/Releases**, and then run **Project → Clean…**. Eclipse should show **Maven Dependencies** in the project tree and use the `JavaSE-17` execution environment.

The custom multi-delimiter JavaCSV jar remains attached directly with its source archive; JOCL, JDistlib, Commons Compress, and JUnit are supplied by Maven.

## Inspect GPUs

```sh
java -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo
```

If no device is reported, install or update the GPU vendor's driver/OpenCL runtime. JOCL supplies Java/native bindings but does not install the vendor's GPU driver.

## Run an analysis

```sh
java -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar path/to/analysis.ini
```

The existing INI keys and input formats are unchanged in this first modernization stage. The former `library_path` key is no longer used; GPU runtime discovery is handled by the operating system's OpenCL loader.

See `AGENTS.md` for contributor constraints and the ordered roadmap. See `SESSION_HISTORY.md` for timestamped changes, verification, and known limitations.
