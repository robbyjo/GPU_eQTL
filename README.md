# GPU eQTL

GPU-accelerated Java software for eQTL and related QTL analyses. This codebase originated in 2011–2013 and is being modernized incrementally. FP64 remains the default scientific mode; FP32 is available as an explicit performance/compatibility choice.

## Current status

The application targets Java 17 and now has selectable GPU backends behind a vendor-neutral API:

| Setting | Implementation | Intended hardware |
| --- | --- | --- |
| `auto` (default) | CUDA-first discovery plus OpenCL fallback | Mixed machines; avoids running the same NVIDIA card twice |
| `cuda` | JCuda 12.6.0 and cuBLAS DGEMM/SGEMM | NVIDIA GPUs |
| `opencl` | JOCL 2.0.6 and the production OpenCL C kernel | NVIDIA, AMD, or Intel GPUs exposed by an OpenCL driver |

Automatic discovery uses CUDA for a usable NVIDIA device and also includes distinct OpenCL devices from other vendors. It filters those devices for the requested precision. If CUDA cannot initialize, NVIDIA OpenCL remains available as a fallback. A native HIP/ROCm backend is not included yet; AMD cards currently use JOCL/OpenCL.

The real-valued eQTL calculation defaults to FP64 and is validated through both CUDA and OpenCL against the same CPU reference. `precision = fp32` or `--precision fp32` converts only the prepared GPU input tiles to single precision and uses SGEMM/a `float` OpenCL kernel; covariate adjustment, standardization, cache values, effects, test statistics, and p-values remain FP64. This permits OpenCL execution on Intel Iris Xe devices that lack native FP64, provided their installed driver exposes a usable OpenCL GPU. The legacy categorical-SNP path remains disabled; the CUDA backend currently implements `eqtlReal` only.

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

Hardware-independent tests always run. CUDA and OpenCL FP64/FP32 numerical tests skip cleanly when their respective runtime/device is unavailable.

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

## Precision, automatic sizing, and multiple GPUs

FP64 is the default. Select FP32 explicitly in an INI file or on the command line:

```ini
precision = fp32
```

```powershell
  --precision fp32
```

FP32 changes results slightly and should not be substituted silently for an established FP64 pipeline. On the representative 128-SNP by 128-trait WHI subset (2,005 samples; 16,384 associations), CUDA FP32 versus FP64 had maximum absolute differences of `5.59e-9` in R-squared, `1.03e-8` in effect, `9.70e-7` in t, and `2.51e-6` in log10 p. Both modes selected the same 11 associations at p <= `1e-4`. A larger 4,096-SNP by 8,192-trait study (33,554,432 associations) found one FP32-only result at that reporting boundary; among the 48,948 common reported results, maximum absolute differences were `9.73e-7` in R-squared, `5.71e-7` in effect, `1.11e-4` in t, and `9.90e-4` in log10 p. These are bounded accuracy studies, not guarantees for every matrix or threshold; revalidate important borderline associations in FP64.

When `block_size`/`--block-size` is omitted or zero, the application reads total VRAM and maximum-allocation limits from every selected GPU. It chooses one block size supported by the least-capable device, accounting for both input tiles, the quadratic result tile, FP32 versus FP64 element size, VRAM headroom, the per-allocation limit, and Java's array-size limit. It also caps the target result allocation at 1 GiB because consuming all available VRAM with one quadratic tile is usually slower and much harder on host memory. When the genotype matrix is smaller than the memory limit, the recommendation retains up to four genotype jobs per selected GPU so one oversized tile does not eliminate packing/post-processing overlap. The choice and its limiting reason are printed before computation; an explicit block size remains an override.

When `num_threads`/`--threads` is omitted or zero, bounded-RAM mode uses up to four workers per GPU to overlap cache reads, packing, GPU execution, and result processing; full-memory mode uses up to two. The bounded-RAM recommendation estimates prepared inputs, packed precision-specific inputs, and the result array per worker, retains JVM-heap headroom, is capped by required genotype jobs, and leaves one CPU core free. Explicit values remain supported, with warnings for idle GPUs, CPU oversubscription, or a streamed configuration estimated to consume more than 75% of the available JVM heap.

All distinct usable GPUs are opened, with one exclusive context per device. Jobs reserve contexts from a shared pool, so a multi-GPU system uses multiple devices concurrently. Automatic discovery suppresses only the duplicate OpenCL representation of an NVIDIA device already selected through CUDA. On mixed-capacity systems, automatic block sizing follows the smallest selected GPU. Use `--gpu-backend cuda` or `opencl` when you intentionally want only one backend family.

The old `lambda` INI setting no longer controls memory or result storage and is ignored with a compatibility notice. It should be removed from new configurations.

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
  --precision fp64
```

Omitting `--block-size` and `--threads` in this example enables the device-aware recommendations described above.

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

Prepared rows are encoded and decoded in checksummed bulk records. This preserves the existing cache format and per-row CRC validation while avoiding a file operation and checksum update for every individual double. Caches created by the earlier version remain reusable.

### Performance profiling

Add `--profile` to print phase totals for metadata/alignment, cache creation and reads, host packing, GPU-context wait, GPU buffer setup, upload, computation, download, CPU statistics/result writing, and final assembly. Worker-phase totals can overlap and therefore need not sum to wall time. Add `--profile-output FILE` to write the same measurements as CSV; specifying a profile output also enables profiling.

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --profile `
  --profile-output run-profile.csv
```

The INI equivalents are `profile = true` and `profile_output = run-profile.csv`. GPU profiling adds synchronization needed to separate upload, compute, and download time, so use the reported breakdown diagnostically and compare complete wall times for throughput decisions.

### Checkpoint and resume

Bounded-RAM analyses write one atomic result part per genotype block. Only a completely written part is considered finished, and final output is assembled in genotype-block order. On success the temporary checkpoint directory is removed. If execution stops, rerun the identical command with `--resume`:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --genotype-block-rows 10240 `
  --expression-block-rows 10240 `
  --resume
```

The default directory is `<output-file>.checkpoint`. `--checkpoint-dir DIR` selects a different location, and `--keep-checkpoints` retains completed parts after successful assembly. Resume is rejected when input/cache signatures, precision, thresholds, degrees of freedom, block sizes, or output modes do not match. The corresponding INI keys are `cache_dir`, `rebuild_cache`, `checkpoint_dir`, `resume`, and `keep_checkpoints`.

Run `--help` for the complete argument list. See `TODO.md` for the remaining modernization work. The former `library_path` key is no longer used; runtime discovery is handled by the CUDA/OpenCL backends.

See `AGENTS.md` for contributor constraints and the ordered roadmap. See `SESSION_HISTORY.md` for timestamped changes, verification, and known limitations.
