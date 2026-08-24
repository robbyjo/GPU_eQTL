# GPU eQTL modernization to-do list

This list contains work remaining after the correctness, CLI, bounded-RAM cache, checkpoint, optional FP32, device-aware tuning, and first cache/scheduler profiling milestones. Preserve the scientific and verification rules in `AGENTS.md`.

## Missing-data production follow-up

- Cross-check the new common missingness audit on representative CSV, VCF.gz, and BCF cohorts, including covariate-driven sample removal, and compare exact row/sample/pattern counts with independent R or bcftools summaries.
- For exact trait-pattern deletion, define and audit pattern-specific handling when a rare predictor becomes monomorphic after that pattern's additional sample removal. Global VCF QC now uses the final aligned cohort, but an individual trait mask can be smaller still.
- Add a deterministic end-to-end CPU association reference for exact trait-pattern deletion and compare identifiers, effective N, residual DF, effects, and p-values with the GPU path. The first implementation is exact but deliberately conservative: it rebuilds prepared predictor rows once per pattern, groups output by pattern, and disables checkpoint resume.
- Add resumable pattern scheduling and reusable pattern-specific predictor caches without weakening cache signatures. Investigate batching compatible masks only if it preserves each trait's exact sample set, projection, rank, threshold, and output N/DF.
- Benchmark `local-pattern` genotype imputation against row-mean and standard reference-panel imputation. Define maximum distance, minimum comparable flanks/donors, phased-input behavior, and richer per-imputation QC before treating the proxy as a production default. Mean remains the default.
- Consider an explicit global complete-sample policy across predictor, trait, and selected covariates. Keep it distinct from `exclude-row` (feature removal) and exact trait-pattern deletion (per-pattern sample removal).

## VCF/VCF.gz and BCF milestone delivered; production follow-up remains

- Delivered sequential plain/gzip VCF and BCF 2.1/2.2 input through the common metadata/block source, with exact header-sample order, canonical variant IDs, explicit `DS`/`GT`, missing and multiallelic policies, MAF/MAC filters, monomorphic/singleton/doubleton classification, EAF/MAF/MAC, exact biallelic HWE, and an atomic variant-QC report.
- Keep the deterministic VCF.gz/BCF fixture and optional independently generated BCF test. Before a production run, validate a representative cohort VCF.gz and BCF against `bcftools query` counts, allele statistics, chosen field, sample order, and several hand-checked variants.
- Add tabix/CSI region selection and true random-access variant blocks. The current first pass scans the whole source for QC and each preparation pass reads sequentially; a prepared cache is indexed, but the original VCF/BCF source is not yet region-queryable in this application.
- Add an explicit multiallelic split policy only after tests define per-ALT DS/GT projection, other-ALT calls, canonical IDs, allele counts, and HWE behavior. Current choices are `exclude` and `error`.

## Performance follow-up

- The prepared-trait residency layer now supports `auto`, `memory`, and `disk`: memory reads the once-residualized FP64 trait cache once and shares immutable rows across genotype workers; auto reserves heap for worker buffers and falls back to disk. Profile all three modes on the complete WHI chromosome and 4,746-sample batch1234 inputs before changing the default heuristic or loop/checkpoint schedule.
- Validate automatic block/worker tuning on real multi-GPU machines and on larger 5,100–5,700-sample cohorts. Revisit the 1-GiB output target and four-worker pipeline only with those measurements.
- Profile GPU fixed-effect residualization at 5,100–5,700 samples, higher covariate ranks, and on real multi-GPU/Intel/AMD systems. The 2,005-sample NVIDIA study was cache-I/O-bound even though projection compute itself was fast; do not claim a universal wall-time gain.
- Benchmark FP32 versus FP64 on Intel and AMD hardware. The larger NVIDIA WHI study found one reporting-threshold classification difference, so retain FP64 verification for borderline FP32 findings.
- Add cache inspection and safe stale-cache pruning commands; never delete caches merely because they are not used by the current run.

## Other input and output work

- Extend validation to family and pedigree identifiers if those inputs remain supported.
- Add optional compressed result output and decide whether the legacy full-memory scheduler should also assemble deterministic block-ordered output.
- Generalize legacy association column labels (`Rs_ID`, `ProbesetID`) through an opt-in output schema for methylation, proteomics, and other declared matrix roles without breaking existing parsers.

## Statistical analyses

- Add gene/region set ingestion and deterministic CPU references for burden, SKAT, and SKAT-O before GPU acceleration. Define allele orientation, MAF/MAC masks, missingness, variant weights, covariate projection, set overlap, null-model reuse, p-value method, and failure behavior for monomorphic or rank-deficient sets.
- Accelerate the set-test linear algebra in bounded variant/set blocks, reuse the existing multi-GPU context pool, and compare every GPU statistic to the CPU reference. SKAT-O must define its rho grid and adjusted omnibus p-value; do not report the minimum unadjusted component p-value as SKAT-O.
- Stream only retained set/trait results and support restartable deterministic scheduling so burden/SKAT analyses do not require all variant-by-trait results on disk. Defer relatedness-aware set tests until the cohort covariance/null-model design below is implemented.
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
