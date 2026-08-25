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
- Print detected compute-backend information with `java -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` (`--printgpuinfo` remains an alias).
- Do not commit `target/`, `bin/`, downloaded dependencies, GPU drivers, or vendor SDKs.

The Maven build currently compiles the eQTL application and its supporting NIH utility packages. `src/gov/nih/exon` is deliberately outside the build because it depends on the absent QGeneric/qstats/qplugin projects. Do not silently add it to the build; either restore those dependencies or extract and test that code first.

The tracked `lib/javacsv-src.zip` contains a project-specific JavaCSV fork that accepts multiple delimiters. Maven unpacks it as a generated source. Do not replace it with the Maven Central JavaCSV artifact without adding equivalent parsing tests and preserving that behavior.

## Compute architecture and rules

- Analysis code depends only on interfaces in `src/gov/nih/gpu`.
- `auto` is the default backend. It prefers eligible CUDA/OpenCL devices, suppresses a duplicate NVIDIA OpenCL representation when CUDA is usable, and falls back to CPU only when no eligible GPU remains. It never mixes CPU into an otherwise usable GPU pool. Users can force `cuda`, `opencl`, or `cpu` with `--backend` or the compatibility `eqtl.gpu.backend` property.
- JOCL-specific code belongs in `src/gov/nih/gpu/opencl` and must not leak JOCL handles into `gov.nih.eqtl`.
- The `opencl` backend uses JOCL 2.0.6 and the machine's OpenCL ICD. A usable device must be available and provide an OpenCL compiler; FP64 mode additionally requires native FP64 support.
- CUDA-specific code belongs in `src/gov/nih/gpu/cuda` and must not leak JCuda/cuBLAS handles into `gov.nih.eqtl`. The `cuda` backend uses JCuda 12.6.0 and cuBLAS DGEMM/SGEMM, requires a compatible system CUDA runtime, and currently supports `eqtlReal` only.
- CPU-specific code belongs in `src/gov/nih/gpu/cpu`. The `cpu` backend can use oneMKL 2026.1 on x64 Windows/Linux evaluation builds, bundled OpenBLAS 0.3.34 on supported desktop platforms, or a pure-Java FP64/FP32 fallback. `eqtl.cpu.blas=auto|mkl|openblas|java` selects the engine; explicitly selected native-engine failures are fatal. `auto` prefers an available oneMKL runtime, then OpenBLAS, then Java. One exclusive CPU context owns the process-global BLAS thread setting (`eqtl.cpu.threads`), so do not create nested BLAS worker pools.
- The project is GPLv3 with the narrow oneMKL linking permission in `LICENSE_EXCEPTION` for portions owned by Roby Joehanes. Intel oneMKL remains under the Intel Simplified Software License. Every jar must contain `META-INF/LICENSE`, `META-INF/LICENSE_EXCEPTION`, and all files from `THIRD_PARTY_LICENSES`; never strip Intel oneMKL/OpenMP or incorporated third-party notices from an MKL artifact.
- NVIDIA GPUs can use CUDA or OpenCL. AMD and Intel GPUs currently use OpenCL when the installed driver exposes a usable GPU for the requested precision. A future native HIP/ROCm, Vulkan-compute, or Level Zero implementation must implement `GpuBackend`, `GpuDevice`, and `GpuContext` instead of branching throughout the analysis code.
- `GpuTuning` owns automatic block and thread recommendations. GPU block sizing must account for both input allocations, the quadratic output allocation, numeric width, total VRAM, maximum allocation, the least-capable selected GPU, and enough genotype jobs to pipeline work when the matrix is smaller than the memory limit. CPU block sizing instead uses available JVM heap and a bounded host result target. Streamed worker tuning must retain JVM-heap headroom for prepared inputs, packed inputs, and result arrays. Do not restore the obsolete `lambda` host-RAM saturation heuristic.
- Keep each `GpuContext` exclusive to one worker. Reserve it immediately before submission, release it in `finally`, and close the pool after workers finish.
- Fixed-effect QR/rank decisions stay on the CPU. GPU residualization implements `Y - (Y Q) Q^T` without materializing an `N x N` projector. Cache Q once per exclusive context, cap concurrent cache preparation for JVM-heap headroom, release projection buffers before association allocation, and keep `--residualization cpu` as the reproducibility path.
- Reuse device allocations where safe. For partial final tiles, submit only rounded active work dimensions; never read or interpret padded cells as results.
- Kernel compilation and GPU API failures must stop the analysis. Never continue after merely logging such an error.

## Scientific verification

- Hardware-independent unit tests must run without a physical GPU. Hardware tests must skip cleanly when no suitable runtime/device is installed.
- Accelerated numerical tests must compare against a scalar CPU calculation with an explicit tolerance. Keep tests exercising the real production CUDA/OpenCL operations, native oneMKL/OpenBLAS, and the portable Java fallback.
- Before changing matrix layout, padding, work sizes, covariate residualization, or statistical formulas, add a small deterministic reference dataset and compare identifiers, pair counts, R-squared values, p-values, and degrees of freedom.
- FP32 must remain opt-in. GPU residualization follows the requested precision, while QR/rank decisions, prepared-cache encoding, standardization, and CPU statistical work stay FP64. Cache/checkpoint signatures must distinguish approximate preprocessing. Any expansion of FP32 behavior requires an explicit accuracy study and user approval.
- Treat missing, duplicated, reordered, or mismatched sample IDs as fatal unless an explicit alignment policy is selected and reported.
- Literal sample-ID prefix removal must be opt-in and audited. `strict` alignment remains the default; `covariate-subset` may exclude only matrix-header extras and must still find every unique covariate sample in both matrices. Reject blank normalized IDs and prefix-induced collisions.
- Missingness must be audited before replacement or deletion. Keep predictor, trait, and selected-covariate policies explicit in the QC output. Pattern-wise trait deletion must recompute the covariate QR/rank, threshold, effective N, and residual degrees of freedom for every exact sample mask.
- `local-pattern` genotype filling is an unphased nearest-flanking-dosage proxy, not haplotype/reference-panel imputation. It is valid only for a declared genotype predictor, must respect chromosome boundaries when `CHROM:POS` metadata are present, must warn when trusting row order, and must retain row-mean fallback behavior in tests. Mean remains the default predictor policy.

## Data and compatibility

- The positional legacy INI interface remains supported. `QeQTLCommandLine` also supports `--config` plus overrides or an argument-only run; keep both paths covered by compatibility tests.
- Headered CSV input goes through `QDelimitedMatrixSource`, which scans metadata and identifiers before exposing reordered row blocks. Blank/duplicate identifiers and mismatched field counts are fatal.
- VCF/BCF frequency, HWE, missingness, singleton/doubleton, and monomorphic QC must be finalized on the aligned analysis columns, not matrix-only header samples. The HTSJDK reader/decoder remains single-threaded; bounded workers may evaluate distinct already-decoded variants, but their decisions, duplicate checks, QC rows, and retained identifiers must be consumed in source order. Do not recompute global QC on later source passes. Preserve the additive biallelic variant semantics and include the selected sample permutation in prepared-cache signatures.
- Resumable variant-QC state is scientific state. Its signature must cover every input, resolved-field, filtering, and aligned-sample-order choice that can change EAF/MAF/MAC/HWE or inclusion. Commit only complete source-ordered parts atomically, validate an incomplete checkpoint against the decoded source prefix, and never reuse an incompatible checkpoint by guessing.
- Indexed region state must also cover the index identity and normalized region/set definitions. Merge only query intervals, not set membership; retain deterministic source-contig order, reject source/index mutation during a run, and never fall back to a full sequential scan when a user requested a region.
- Aligned variant QC persists dosage sum, squared-dosage sum, called/missing counts, and diploid hard-call counts. Exact trait-pattern analysis may reuse a checksummed aligned raw FP64 cache and pattern-specific sufficient-statistics/prepared caches. Cohort HWE and aligned MAF/MAC remain the default QC contract; a pattern-specific frequency filter must be explicit and reported because it changes the variant set across traits.
- `QCovariateTable` may bridge different genotype/expression ID columns, automatically encodes text covariates with a deterministic reference level, and rejects rank-deficient models before computation.
- `genotype_block_rows` or `expression_block_rows` enables bounded-RAM CSV analysis. `QBinaryMatrixCache` stores aligned, residualized, standardized FP64 rows with an index and per-row checksum. Bulk row I/O must preserve the version-1 byte format and checksum validation so existing caches remain compatible. Its signature must cover source metadata, sample ordering, covariate projection, and non-legacy preprocessing precision/backend; never reuse a mismatched cache.
- `QAnalysisCheckpoint` writes one atomic part per genotype block. Resume must validate the complete analysis signature, and final output must be assembled in genotype-block order. Never treat a `.partial` file as complete.
- An omitted/zero `block_size` and `num_threads` selects backend-aware recommendations. Explicit values remain overrides. Preserve one exclusive context per discovered GPU so multi-GPU execution continues through `GpuContextPool`, and one CPU context so OpenBLAS parallelism is not nested.
- Prepared caches remove repeated CSV parsing and covariate preprocessing, but matrix cache blocks may still be reread for the cross-product schedule. Profile OS cache behavior, disk reads, and batching before claiming an overall speed improvement. Plain, gzip, and bzip2 CSV streams are supported.
- `trait_cache = auto|memory|disk` changes only residency of already aligned/residualized/standardized FP64 trait rows. Memory mode must retain bounded row arrays, share them read-only across workers, reserve heap for active worker buffers, and preserve disk-cache/checkpoint signatures and exact output bytes.
- `--validate-only` performs the metadata, ID-alignment, covariate, rank, and degrees-of-freedom checks without association computation.
- `--profile` records host-observed phase timings; backend timing must remain opt-in because separating GPU upload, compute, and download may require synchronization. Concurrent worker totals can overlap and must not be represented as additive wall-clock percentages.
- Preserve the repository's GNU GPL version 3 headers, the oneMKL exception notice where applicable, and original author attribution.
- Avoid unrelated reformatting in legacy source files; mixed historical line endings can otherwise obscure reviews.

## Modernization roadmap

Use this dependency-aware order unless the user asks otherwise:

Completed foundations: deterministic fixtures and ID validation/reordering; categorical-covariate encoding and rank checks; CLI plus legacy INI compatibility; bounded CSV blocks; reusable indexed prepared caches; atomic checkpoint/resume with deterministic streamed output assembly; opt-in FP32; backend-aware block/thread recommendations with multi-GPU context pooling; phase profiling; reusable CUDA/OpenCL fixed-effect block residualization; and an OpenBLAS/pure-Java CPU analysis backend.

The remaining ordered work is maintained in `TODO.md`. Profile the cache-backed schedule before kernel changes; add indexed VCF/BCF before forward selection; and treat categorical-SNP repair separately from categorical covariates.

## Session records

Append to `SESSION_HISTORY.md`; never rewrite prior entries. Each entry must include:

- An ISO-8601 timestamp with timezone.
- The baseline commit and goal.
- Decisions and files changed.
- Exact verification commands and test outcomes, including GPU/driver details for hardware tests.
- Known limitations, compatibility notes, and the next recommended step.
