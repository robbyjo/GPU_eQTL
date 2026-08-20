# GPU eQTL contributor guide

This file applies to the entire `NIH-Project` subtree. Read `SESSION_HISTORY.md` before making changes and append a timestamped entry after any material work.

## Project purpose

This is a Java application for GPU-accelerated eQTL and related QTL analyses. The scientific contract is more important than raw speed: changes must preserve double-precision results, sample ordering, degrees of freedom, thresholds, and output identifiers.

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
- JOCL-specific code belongs in `src/gov/nih/gpu/opencl` and must not leak JOCL handles into `gov.nih.eqtl`.
- The current `opencl` backend uses JOCL 2.0.6 and the machine's OpenCL ICD. A usable device must be available, provide an OpenCL compiler, and support FP64.
- NVIDIA, AMD, and Intel GPUs can use this backend when their installed vendor driver exposes OpenCL. A future native CUDA, HIP/ROCm, Vulkan-compute, or Level Zero implementation should implement `GpuBackend`, `GpuDevice`, and `GpuContext` instead of branching throughout the analysis code.
- Keep each `GpuContext` exclusive to one worker. Reserve it immediately before submission, release it in `finally`, and close the pool after workers finish.
- Reuse device allocations where safe. For partial final tiles, submit only rounded active work dimensions; never read or interpret padded cells as results.
- Kernel compilation and GPU API failures must stop the analysis. Never continue after merely logging such an error.

## Scientific verification

- Hardware-independent unit tests must run without a physical GPU. Hardware tests must skip cleanly when no suitable runtime/device is installed.
- GPU numerical tests must compare against a CPU calculation with an explicit tolerance. Keep at least one test exercising the real production kernel.
- Before changing matrix layout, padding, work sizes, covariate residualization, or statistical formulas, add a small deterministic reference dataset and compare identifiers, pair counts, R-squared values, p-values, and degrees of freedom.
- Do not lower the computation from `double` to `float` without an explicit accuracy study and user approval.
- Treat missing, duplicated, reordered, or mismatched sample IDs as fatal unless an explicit alignment policy is selected and reported.

## Data and compatibility

- The current run interface is a legacy INI file. Keep it working until the command-line replacement has compatibility tests and a documented deprecation path.
- The current loaders materialize the genotype and expression matrices in RAM. New loaders should expose metadata separately from block iteration so validation can happen before computation.
- Preserve the repository's GNU GPL version 3 headers and original author attribution.
- Avoid unrelated reformatting in legacy source files; mixed historical line endings can otherwise obscure reviews.

## Modernization roadmap

Use this dependency-aware order unless the user asks otherwise:

1. Add deterministic end-to-end reference fixtures, sample/gene/SNP ID validation, duplicate detection, and explicit reordering rules.
2. Add categorical-covariate parsing and automatic reference-level one-hot encoding, with rank/collinearity checks. Repair and test the currently disabled categorical-SNP path separately.
3. Replace the INI-only entry point with full command-line arguments, while supporting `--config` for the old files during migration.
4. Introduce block readers for SNP and expression matrices, then stream blocks through the existing GPU scheduler without changing result ordering.
5. Profile transfers, packing, kernel occupancy, partial tiles, multi-GPU load balancing, and output writing before making further performance changes.
6. Add SNP interaction analysis with explicit model/degree-of-freedom tests and bounded block-pair scheduling.
7. Add forward selection on top of an indexed/re-readable genotype source; do not require retaining every SNP in RAM.

## Session records

Append to `SESSION_HISTORY.md`; never rewrite prior entries. Each entry must include:

- An ISO-8601 timestamp with timezone.
- The baseline commit and goal.
- Decisions and files changed.
- Exact verification commands and test outcomes, including GPU/driver details for hardware tests.
- Known limitations, compatibility notes, and the next recommended step.
