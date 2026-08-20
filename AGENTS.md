# GPU eQTL contributor guide

This file applies to the entire `NIH-Project` subtree. Read `SESSION_HISTORY.md` before making changes and append a timestamped entry after any material work.

## Project purpose

This is a Java application for GPU-accelerated eQTL and related QTL analyses. The scientific contract is more important than raw speed: FP64 is the default and must preserve double-precision results, sample ordering, degrees of freedom, thresholds, and output identifiers. FP32 is an explicit approximate mode and must always be identified as such.

## Supported build

- Use Java 17 or newer. Java 17 is the compilation target.
- Use the checked-in Maven wrapper: `./mvnw test` on POSIX systems or `.\mvnw.cmd test` on Windows.
- In Eclipse, use Maven Integration for Eclipse (m2e), import as an Existing Maven Project, and run Maven → Update Project after pom.xml changes.
- Build the runnable jar with `./mvnw clean package` or `.\mvnw.cmd clean package`.
- The packaged application is `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar`.
- Print detected GPU information with `java -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo`.
- Do not commit `target/`, `bin/`, downloaded dependencies, GPU drivers, or vendor SDKs.

The Maven build currently compiles the eQTL application and its supporting NIH utility packages. `src/gov/nih/exon` is deliberately outside the build because it depends on the absent QGeneric/qstats/qplugin projects. Do not silently add it to the build; either restore those dependencies or extract and test that code first.

The tracked `lib/javacsv-src.zip` contains a project-specific JavaCSV fork that accepts multiple delimiters. Maven unpacks it as a generated source. Do not replace it with the Maven Central JavaCSV artifact without adding equivalent parsing tests and preserving that behavior.

## GPU architecture and rules

- Analysis code depends only on interfaces in `src/gov/nih/gpu`.
- `auto` is the default backend. It prefers a usable CUDA device over the duplicate NVIDIA OpenCL device, while retaining distinct OpenCL devices from other vendors and filtering them for the requested precision. Users can force `cuda` or `opencl` with the `eqtl.gpu.backend` system property.
- JOCL-specific code belongs in `src/gov/nih/gpu/opencl` and must not leak JOCL handles into `gov.nih.eqtl`.
- The `opencl` backend uses JOCL 2.0.6 and the machine's OpenCL ICD. A usable device must be available and provide an OpenCL compiler; FP64 mode additionally requires native FP64 support.
- CUDA-specific code belongs in `src/gov/nih/gpu/cuda` and must not leak JCuda/cuBLAS handles into `gov.nih.eqtl`. The `cuda` backend uses JCuda 12.6.0 and cuBLAS DGEMM/SGEMM, requires a compatible system CUDA runtime, and currently supports `eqtlReal` only.
- NVIDIA GPUs can use CUDA or OpenCL. AMD and Intel GPUs currently use OpenCL when the installed driver exposes a usable GPU for the requested precision. A future native HIP/ROCm, Vulkan-compute, or Level Zero implementation must implement `GpuBackend`, `GpuDevice`, and `GpuContext` instead of branching throughout the analysis code.
- `GpuTuning` owns automatic block and thread recommendations. Block sizing must account for both input allocations, the quadratic output allocation, numeric width, total VRAM, maximum allocation, and the least-capable selected GPU. Do not restore the obsolete `lambda` host-RAM saturation heuristic.
- Keep each `GpuContext` exclusive to one worker. Reserve it immediately before submission, release it in `finally`, and close the pool after workers finish.
- Reuse device allocations where safe. For partial final tiles, submit only rounded active work dimensions; never read or interpret padded cells as results.
- Kernel compilation and GPU API failures must stop the analysis. Never continue after merely logging such an error.

## Scientific verification

- Hardware-independent unit tests must run without a physical GPU. Hardware tests must skip cleanly when no suitable runtime/device is installed.
- GPU numerical tests must compare against a CPU calculation with an explicit tolerance. Keep at least one test exercising the real production kernel.
- Before changing matrix layout, padding, work sizes, covariate residualization, or statistical formulas, add a small deterministic reference dataset and compare identifiers, pair counts, R-squared values, p-values, and degrees of freedom.
- FP32 must remain opt-in. Prepared caches and CPU statistical work stay FP64; checkpoint signatures must include precision. Any expansion of FP32 behavior requires an explicit accuracy study and user approval.
- Treat missing, duplicated, reordered, or mismatched sample IDs as fatal unless an explicit alignment policy is selected and reported.

## Data and compatibility

- The positional legacy INI interface remains supported. `QeQTLCommandLine` also supports `--config` plus overrides or an argument-only run; keep both paths covered by compatibility tests.
- Headered CSV input goes through `QDelimitedMatrixSource`, which scans metadata and identifiers before exposing reordered row blocks. Blank/duplicate identifiers and mismatched field counts are fatal.
- `QCovariateTable` may bridge different genotype/expression ID columns, automatically encodes text covariates with a deterministic reference level, and rejects rank-deficient models before computation.
- `genotype_block_rows` or `expression_block_rows` enables bounded-RAM CSV analysis. `QBinaryMatrixCache` stores aligned, residualized, standardized FP64 rows with an index and per-row checksum. Its signature must cover source metadata, sample ordering, and the covariate projection; never reuse a mismatched cache.
- `QAnalysisCheckpoint` writes one atomic part per genotype block. Resume must validate the complete analysis signature, and final output must be assembled in genotype-block order. Never treat a `.partial` file as complete.
- An omitted/zero `block_size` and `num_threads` selects the tested GPU-aware recommendations. Explicit values remain overrides. Preserve one exclusive context per discovered GPU so multi-GPU execution continues through `GpuContextPool`.
- Prepared caches remove repeated CSV parsing and covariate preprocessing, but matrix cache blocks may still be reread for the cross-product schedule. Profile OS cache behavior, disk reads, and batching before claiming an overall speed improvement. Plain, gzip, and bzip2 CSV streams are supported.
- `--validate-only` performs the metadata, ID-alignment, covariate, rank, and degrees-of-freedom checks without association computation.
- Preserve the repository's GNU GPL version 3 headers and original author attribution.
- Avoid unrelated reformatting in legacy source files; mixed historical line endings can otherwise obscure reviews.

## Modernization roadmap

Use this dependency-aware order unless the user asks otherwise:

Completed foundations: deterministic fixtures and ID validation/reordering; categorical-covariate encoding and rank checks; CLI plus legacy INI compatibility; bounded CSV blocks; reusable indexed prepared caches; atomic checkpoint/resume with deterministic streamed output assembly; opt-in FP32; and device-aware block/thread recommendations with multi-GPU context pooling.

The remaining ordered work is maintained in `TODO.md`. Profile the cache-backed schedule before kernel changes; add indexed VCF/BCF before forward selection; and treat categorical-SNP repair separately from categorical covariates.

## Session records

Append to `SESSION_HISTORY.md`; never rewrite prior entries. Each entry must include:

- An ISO-8601 timestamp with timezone.
- The baseline commit and goal.
- Decisions and files changed.
- Exact verification commands and test outcomes, including GPU/driver details for hardware tests.
- Known limitations, compatibility notes, and the next recommended step.
