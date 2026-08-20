# GPU eQTL modernization to-do list

This list contains work remaining after the correctness, CLI, bounded-RAM cache, and checkpoint milestones. Preserve the scientific and verification rules in `AGENTS.md`.

## Next priority: profile and optimize the cache-backed schedule

- Measure metadata scanning, first-run cache creation, cache reads, host packing, GPU transfer, kernel/cuBLAS time, CPU post-processing, and final output assembly on representative sample and block sizes.
- Measure whether repeated prepared-block reads are satisfied by the operating-system file cache. If not, batch genotype jobs so one expression cache block can feed several GPU submissions, or add a bounded shared expression-block cache.
- Compare block sizes and worker counts without exceeding GPU allocation limits. Include partial final tiles and multiple GPUs.
- Reduce legacy residualization thread/log overhead for small blocks and buffer output without changing double-precision results.
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
- Reassess Intel Level Zero only on Intel hardware with scientifically usable FP64.
- Consolidate shaded-jar metadata warnings, document supported driver/runtime combinations, and prepare reproducible release artifacts and checksums.
