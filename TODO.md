# GPU eQTL modernization to-do list

This list contains work remaining after the correctness, CLI, bounded-RAM cache, checkpoint, optional FP32, and device-aware tuning milestones. Preserve the scientific and verification rules in `AGENTS.md`.

## Next priority: profile and optimize the cache-backed schedule

- Measure metadata scanning, first-run cache creation, cache reads, host packing, GPU transfer, kernel/cuBLAS time, CPU post-processing, and final output assembly on representative sample and block sizes.
- Measure whether repeated prepared-block reads are satisfied by the operating-system file cache. If not, batch genotype jobs so one expression cache block can feed several GPU submissions, or add a bounded shared expression-block cache.
- Benchmark the new device-aware block and worker recommendations against representative end-to-end workloads; tune the 1-GiB output-tile target only with evidence. Include partial final tiles and real multi-GPU machines.
- Benchmark FP32 versus FP64 throughput on NVIDIA, Intel, and AMD hardware and expand the accuracy study beyond the representative WHI subset, especially around reporting thresholds.
- Reduce legacy residualization thread/log overhead for small blocks and buffer output without changing FP64 results.
- Add cache inspection and safe stale-cache pruning commands; never delete caches merely because they are not used by the current run.

## Input and output formats

- Add indexed VCF/VCF.gz and BCF genotype input, with explicit `DS`/`GT` policy, multiallelic handling, missing-value rules, region selection, and exact header-sample alignment.
- Extend validation to family and pedigree identifiers if those inputs remain supported.
- Add optional compressed result output and decide whether the legacy full-memory scheduler should also assemble deterministic block-ordered output.

## Statistical analyses

- Add SNP interaction analysis with an explicit model definition, degrees-of-freedom tests, deterministic CPU references, and bounded block-pair scheduling.
- Add forward selection on top of an indexed/re-readable genotype source, including stopping rules, conditional-model degrees of freedom, and reproducible tie handling.
- Repair and separately verify the disabled categorical-SNP analysis path. Automatic categorical covariates do not make that path correct.

## GPU backends and release work

- Validate on AMD FP64 hardware before considering a native HIP/ROCm backend; retain JOCL/OpenCL as the current AMD path.
- Validate the JOCL FP32 path on an Intel Iris Xe and reassess Level Zero only if it materially improves a verified FP32 workflow or supports scientifically usable FP64 hardware.
- Consolidate shaded-jar metadata warnings, document supported driver/runtime combinations, and prepare reproducible release artifacts and checksums.
