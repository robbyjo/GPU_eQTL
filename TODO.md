# GPU eQTL modernization to-do list

This list contains work remaining after the correctness, CLI, bounded-RAM cache, checkpoint, optional FP32, device-aware tuning, and first cache/scheduler profiling milestones. Preserve the scientific and verification rules in `AGENTS.md`.

## Next priority: indexed VCF/VCF.gz and BCF input

- Add indexed VCF/VCF.gz and BCF genotype input, with explicit `DS`/`GT` policy, multiallelic handling, missing-value rules, region selection, and exact header-sample alignment.
- Expose metadata and block iteration through the same source abstraction used by bounded CSV input so validation, preparation caches, checkpoint/resume, and future forward selection do not depend on the file format.
- Add deterministic fixtures for allele orientation, multiallelic records, missing dosage/genotype values, compressed/indexed access, sample reordering, and malformed headers before representative-data testing.

## Performance follow-up

- Use `--profile` on the complete WHI chromosome run to determine whether rereading the full prepared expression cache once per genotype block remains material after bulk-row I/O. Only then consider expression-block sharing or a different loop/checkpoint schedule.
- Validate automatic block/worker tuning on real multi-GPU machines and on larger 5,100–5,700-sample cohorts. Revisit the 1-GiB output target and four-worker pipeline only with those measurements.
- Profile GPU fixed-effect residualization at 5,100–5,700 samples, higher covariate ranks, and on real multi-GPU/Intel/AMD systems. The 2,005-sample NVIDIA study was cache-I/O-bound even though projection compute itself was fast; do not claim a universal wall-time gain.
- Benchmark FP32 versus FP64 on Intel and AMD hardware. The larger NVIDIA WHI study found one reporting-threshold classification difference, so retain FP64 verification for borderline FP32 findings.
- Add cache inspection and safe stale-cache pruning commands; never delete caches merely because they are not used by the current run.

## Other input and output work

- Extend validation to family and pedigree identifiers if those inputs remain supported.
- Add optional compressed result output and decide whether the legacy full-memory scheduler should also assemble deterministic block-ordered output.

## Statistical analyses

- Add SNP interaction analysis with an explicit model definition, degrees-of-freedom tests, deterministic CPU references, and bounded block-pair scheduling.
- Add forward selection on top of an indexed/re-readable genotype source, including stopping rules, conditional-model degrees of freedom, and reproducible tie handling.
- Repair and separately verify the disabled categorical-SNP analysis path. Automatic categorical covariates do not make that path correct.

## Further goal: cohort-aware joint analysis

- Add a cohort-partitioned joint-analysis mode for FHS, WHI, CARDIA, JHS, and similar studies. Each cohort must accept its own fixed-effect design, including cohort-specific technical covariates, batches, exams, ancestry covariates, and intercept; do not require covariate columns or coefficients to be shared across cohorts.
- Project both genotype and expression through the same cohort-specific nuisance model before concatenating samples. Treat literal stacking of transformed sample columns and blockwise accumulation of score/information statistics as equivalent implementations, and add a deterministic test proving that equivalence.
- Support covariance-aware pre-whitening where observations are not independent: an aligned pedigree/GRM block for FHS-like familial cohorts and a subject-level correlation model for CARDIA-like repeated exams. Do not treat subtraction of phenotype BLUPs alone as independent residualization. For repeated measures, also support a documented one-phenotype-per-subject aggregation policy when a stable average eQTL effect is the intended estimand; never absorb the time-invariant genotype effect with subject fixed effects.
- Keep ordinary fixed-effect projection for independent cohorts such as WHI and JHS. Model a small number of batches or cohorts as categorical fixed effects rather than automatically estimating random-intercept variances.
- Define the cross-cohort estimand and scaling explicitly. Preserve allele dosage for a per-allele effect; require harmonized effect alleles and comparable expression units or documented cohort-specific precision scaling. Do not silently substitute equal-cohort, sample-size, or standardized-genotype weighting for the information implied by the joint model.
- Stream cohort contributions by SNP/expression tile so the program writes only combined results rather than complete per-cohort QTL files. Track cohort availability/effective sample information and optionally accumulate a common-effect heterogeneity statistic for every tile; run cohort-specific or leave-one-cohort-out diagnostics only for retained findings.
- Add deterministic CPU fixtures covering different covariate sets, unequal cohort sizes and residual variances, allele-frequency differences, missing cohort/variant combinations, related families, repeated subjects, rank and degrees-of-freedom accounting, stacked-versus-streamed equality, and failure on ID or allele mismatches before adding GPU acceleration.

## GPU backends and release work

- Validate on AMD FP64 hardware before considering a native HIP/ROCm backend; retain JOCL/OpenCL as the current AMD path.
- Validate the JOCL FP32 path on an Intel Iris Xe and reassess Level Zero only if it materially improves a verified FP32 workflow or supports scientifically usable FP64 hardware.
- Consolidate shaded-jar metadata warnings, document supported driver/runtime combinations, and prepare reproducible release artifacts and checksums.
