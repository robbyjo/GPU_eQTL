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

The legacy INI format remains supported, and every analysis setting now also has a command-line form. Headered CSV matrices are validated before computation: blank or duplicate row/sample IDs are fatal, genotype and expression samples are explicitly reordered, and a covariate table may bridge different ID namespaces. Text covariates are automatically one-hot encoded using the lexicographically first level as the reference; numeric-looking factors can be forced with `--factor-covariates`. Rank-deficient covariate models stop before analysis.

Headered CSV genotype and expression matrices can also be processed in bounded-RAM blocks. The first run creates indexed FP64 cache files after validation, reordering, covariate residualization, and standardization. Later runs with the same source metadata, sample order, and covariate projection reuse those prepared rows instead of reparsing or reprocessing the CSV. Each cache row has an integrity checksum. SNP interactions, forward selection, and indexed VCF/BCF input remain on the roadmap.

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

Force CUDA or OpenCL with the application argument:

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --gpu-backend cuda --printgpuinfo
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --gpu-backend opencl --printgpuinfo
```

The equivalent system property remains supported in PowerShell (quote the `-D` argument):

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

The equivalent migration form is `--config`; command-line options override values read from the INI file:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --gpu-backend cuda `
  --threads 4
```

A run no longer needs an INI file:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype genotype.csv `
  --expression expression.csv `
  --covariates covariates.csv `
  --output results.csv `
  --fixed-covariates Age,Batch,PC1,PC2 `
  --threshold pval 1e-4 `
  --block-size 10240 `
  --threads 4
```

If genotype and expression headers use different identifiers, the program searches for covariate columns whose unique values exactly match each header. Ambiguous cases must be resolved explicitly:

```powershell
  --genotype-id-column NWDID --expression-id-column TORID
```

Use `--validate-only` to scan all row/sample IDs, resolve alignment, encode covariates, check model rank, and report degrees of freedom without running the association analysis.

### Bounded-RAM CSV mode

Set either block-row option to enable the streamed path; an omitted block-row value falls back to `--block-size`. Values are rounded up to the 16-row GPU tile boundary.

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --genotype-block-rows 10240 `
  --expression-block-rows 10240
```

The same settings can be added to an INI file as `genotype_block_rows` and `expression_block_rows`. ID bridge columns are `genotype_id_column` and `expression_id_column`. Plain CSV, `.csv.gz`, and `.csv.bz2` sources are supported by the metadata/block reader. CSV headers and a first-column row identifier are required for validation and streaming; legacy headerless genotype CSV files continue through the old full-memory loader.

By default, caches are placed in `.gpu-eqtl-cache` beside the output file. Choose another location, preferably fast local SSD storage, with:

```powershell
  --cache-dir D:\eqtl-cache
```

Use `--rebuild-cache` after an intentional forced rebuild. Cache signatures include the input path, size, modification time, row/sample metadata, column permutation, and covariate projection. A changed signature creates a different cache rather than silently reusing incompatible prepared values. Old cache versions are not automatically deleted. For the complete WHI chromosome example, allow roughly 6–7 GB of additional uncompressed cache storage.

### Checkpoint and resume

Bounded-RAM analyses write one atomic result part per genotype block. Only a completely written part is considered finished, and final output is assembled in genotype-block order. On success the temporary checkpoint directory is removed. If execution stops, rerun the identical command with `--resume`:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --genotype-block-rows 10240 `
  --expression-block-rows 10240 `
  --resume
```

The default directory is `<output-file>.checkpoint`. `--checkpoint-dir DIR` selects a different location, and `--keep-checkpoints` retains completed parts after successful assembly. Resume is rejected when input/cache signatures, thresholds, degrees of freedom, block sizes, or output modes do not match. The corresponding INI keys are `cache_dir`, `rebuild_cache`, `checkpoint_dir`, `resume`, and `keep_checkpoints`.

Run `--help` for the complete argument list. See `TODO.md` for the remaining modernization work. The former `library_path` key is no longer used; runtime discovery is handled by the CUDA/OpenCL backends.

See `AGENTS.md` for contributor constraints and the ordered roadmap. See `SESSION_HISTORY.md` for timestamped changes, verification, and known limitations.
