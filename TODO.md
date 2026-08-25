# GPU eQTL modernization to-do list

This list contains work remaining after the correctness fixtures, CLI, sample alignment, categorical covariates, bounded-RAM caches, checkpointing, optional FP32, backend-aware tuning, missingness policies, VCF/BCF input, aligned-sample variant QC, association progress, reusable VCF/BCF preprocessing, initial profiling, and OpenBLAS/pure-Java CPU-backend milestones. Preserve the scientific and verification rules in `AGENTS.md`.

Recommended dependency order after closing the production VCF/VCF.gz milestone:

`scalable trait-pattern scheduling -> set-test foundation -> burden -> GPU burden -> SKAT -> SKAT-O`

Bounded CPU-output profiling and release validation should proceed alongside that main sequence. Cohort covariance models must be established before relatedness-aware set tests.

## 0. Completed: production VCF/VCF.gz and QC milestone

- Delivered: resumable aligned-sample QC with atomically committed ordered parts, signatures covering input metadata/policies/sample order, source-prefix or indexed-boundary validation after interruption, and direct reuse of completed scans. Production-test interruption/restart on a representative chromosome and record checkpoint storage/reread overhead.
- Delivered: approximately 15-second association progress based on actual variant–trait comparisons, including resumed-work accounting, percentage, elapsed time, throughput, and ETA across full-memory, streamed, multi-backend, and trait-pattern paths.
- Delivered: `--preprocess-only` aligned VCF/BCF QC and checksummed raw-dosage cache creation without backend initialization or association; matching later analyses discover and reuse the cache. Production chr22 validation reused a completed 15,979,340-record QC checkpoint in 71.81 seconds, opened the 440,406-row genotype cache in 0.96 seconds and the 114,406-row expression cache in 0.25 seconds, and retained 37,824,081,886 bytes (35.23 GiB) across three matrix-cache files.
- Delivered for VCF/VCF.gz: low-allocation worker parsing of raw FORMAT/GT/DS/FT columns, header-index sample selection, four bounded in-flight records per worker, HTSJDK-reference byte-equivalence tests, and pipeline diagnostics; this replaced a slower first attempt that materialized full HTSJDK genotype objects concurrently. A 4,746-sample/1,500-variant bounded benchmark plateaued at eight workers, so automatic mode is capped at 8 while explicit overrides remain available. Unrestricted unsorted-header VCF checkpoint resumes avoid HTSJDK's eager indexed decoder by validating the completed raw prefix before returning to parallel QC; indexed region scans retain direct seeks. The corrected production resume reached approximately 17,600 records/s versus 1,000–1,100 records/s on the eager indexed path and approximately 1,300 records/s in the older implementation.
- Closure decision: the biallelic additive VCF/VCF.gz workflow is production-hardened for the current use case. A short external bcftools/R QC audit remains desirable release evidence but does not block subsequent work. Production BCF validation is explicitly deferred because no representative BCF input is currently available; BCF reader-thread genotype expansion remains a known performance limitation.
- Delivered baseline for `--trait-cache memory` on the 4,746-sample batch1234 chr22 analysis: 50,385,088,836 variant–trait comparisons completed in 1,995.79 seconds (33.26 minutes; 25.25 million comparisons/s wall-clock), with 12,096 GPU-compute calls totaling 1,544.93 worker-seconds and CPU result/filter/write calls totaling 1,269.69 worker-seconds. Profile `auto` and `disk` only if their operational tradeoff is needed; do not repeat the expensive full analysis solely to complete a mode sweep.

## 1. Validate and package the completed CPU analysis backend

- Delivered in the current worktree: `--backend auto|cuda|opencl|cpu`; GPU-first automatic CPU fallback; one FP64/FP32 CPU context; OpenBLAS and pure-Java engines; association and residualization operations; separate BLAS/application worker controls; host-heap block sizing; CPU pipeline sizing; native/fallback numerical tests; and one universal shaded jar containing Windows x64, Linux x64/ARM64, and macOS x64/ARM64 OpenBLAS natives.
- Run representative full analyses through forced `openblas` and `java` modes. Compare identifiers, partial tiles, thresholds, effects, p-values, checkpoint assembly, and exact output order with CUDA/OpenCL FP64 results; record throughput and memory, but do not imply CPU is generally faster than a suitable GPU.
- Smoke-test the shaded jar on every bundled platform and publish checksums. If the universal jar is inconveniently large, add reproducible classifier-specific release profiles while retaining one simple extract-and-run choice.
- Delivered: selectable oneMKL 2026.1 FP64/FP32 association and residualization on x64 Windows/Linux, automatic oneMKL -> OpenBLAS -> Java fallback, numerical reference tests, isolated platform profiles, a GPLv3 Section 7 oneMKL linking exception, and exact packaged Intel oneMKL/OpenMP plus incorporated-code notices. Benchmark it against forced OpenBLAS on representative cohort tiles and full runs before choosing it automatically for any release; the default release remains free of Intel-licensed oneMKL native binaries.

## 2. Deferred indexed/BCF and region extensions

- Delivered: tabix-indexed BGZF VCF and Tribble-indexed VCF/BCF interval queries; repeatable inline and TSV region/set definitions; explicit one-based/BED coordinates; deterministic contig aliases/order; merged query intervals with overlapping set membership; empty-set reporting; source/index/region checkpoint signatures; mutation checks; and direct indexed seek to an interrupted QC boundary.
- Production BCF validation is deferred until a representative BCF file is available. Keep synthetic/hardware-independent BCF tests active in the meantime; this does not block VCF.gz association or the initial rare-variant set-test foundation.
- Add a pure-Java standard CSI decoder (or a carefully packaged cross-platform htslib bridge) for ordinary bcftools `.csi` files. HTSJDK 5 currently handles tabix/Tribble indexes but not standard variant CSI; unsupported CSI must remain a clear pre-analysis failure, never a silent sequential fallback.
- Validate representative indexed VCF.gz and BCF queries, region unions, set memberships, and interruption/restart against bcftools on production copies before depending on them for set tests or forward selection.
- Add an explicit multiallelic split policy only after tests define per-ALT DS/GT projection, other-ALT calls, canonical IDs, allele counts, missingness, and HWE behavior. Current choices remain `exclude` and `error`.

## 3. Rare-variant burden, SKAT, and SKAT-O analyses

### 3a. Shared set-test foundation

- Define gene/region set input, effect-allele orientation, allele harmonization, MAF/MAC masks, missingness policies, variant weights, set overlap, covariate projection, null-model reuse, and failure behavior for monomorphic, empty, or rank-deficient sets.
- Begin with continuous traits and deterministic FP64 CPU fixtures. Compare identifiers, retained variants, weights, effective N, degrees of freedom, statistics, and p-values against an independent implementation before GPU acceleration.
- Stream variant/set and trait blocks and write only requested or retained set/trait results. Add restartable deterministic scheduling so set tests do not require complete variant-by-trait QTL files on disk.

### 3b. Burden tests first

- Implement weighted and unweighted burden tests on the CPU first. Burden testing should validate the complete region-ingestion, filtering, weighting, residualization, missingness, and output pipeline before the more complex kernel tests.
- Add GPU/optimized-CPU acceleration only after exact CPU-reference agreement. Reuse prepared covariate/null-model state and bounded variant/set tiles.

### 3c. SKAT

- Implement a deterministic CPU SKAT reference, including the weighted kernel, covariance/eigenvalue construction, and a defined quadratic-form p-value method with explicit convergence/fallback behavior.
- Accelerate the large repeated set-by-trait score and covariance operations through CPU BLAS and the existing multi-GPU context pool. Keep the relatively small eigenvalue and p-value work on the CPU unless profiling proves otherwise.

### 3d. SKAT-O

- Define and test the rho grid, burden/SKAT mixture construction, correlation between component tests, and the adjusted omnibus p-value. Never report the minimum unadjusted component p-value as SKAT-O.
- Compare every optimized/GPU statistic and final adjusted p-value with the deterministic CPU reference over common, rare, singleton, empty, collinear, overlapping, and missing-data fixtures.
- Defer relatedness-aware burden/SKAT/SKAT-O until the cohort covariance/null-model design below is implemented and verified.

## 4. Missing-data production hardening

- Delivered: exact subset mean/EAF/MAF/MAC/squared-dosage statistics; pattern-specific mean genotype filling; per-pattern monomorphic skipping; optional `--frequency-scope pattern`; compact pattern QC summaries; a checksummed aligned raw dosage cache; and persistent statistics/prepared caches keyed by source, alignment, policy, and sample mask.
- Delivered: a deterministic end-to-end FP64 CPU association reference for exact trait-pattern deletion covering identifiers, effective N, residual degrees of freedom, R-squared, effects, t statistics, and p-values.
- Delivered: durable ordered pattern-result/QC groups, nested genotype-block resume, signature validation, atomic final assembly, completed-group reuse, retained-checkpoint support, and a failure/restart test proving byte-identical output.
- Delivered: fail-fast preflight reporting pattern count, observed-N range, repeated predictor-preparation amplification, upper-bound prepared numeric storage, and structurally non-positive-DF patterns. `--max-trait-patterns` defaults to 256 and requires an explicit decision for larger schedules.
- Delivered: checksummed, atomically written exact matrix missingness scans keyed by source path/size/mtime, parser/cache tag, dimensions, selected sample identities, and exact selected order. Matching reruns regenerate the human-readable audit without rereading the full matrix; corrupt or incompatible scans rebuild.
- Delivered first scalable genotype-block-outer scheduler: one globally prepared padded trait cache; per-mask FP64 covariate QR/rank/DF/threshold state; GPU/CPU sufficient-statistics products for exact predictor projection and variance; exact per-mask mean/zero predictor imputation; pattern-monomorphic skipping; one compute-context-exclusive worker per selected device; deterministic genotype-block/trait-block/variant/original-estimable-trait output; signed block resume; and audited `error|skip` handling for unestimable masks. Deterministic tests match the established pattern-outer VCF implementation for complete, mean-imputed, and zero-imputed predictors and prove byte-identical restart.
- Delivered bounded production validation on batch1234: a tabix-indexed 10-kb chr22 slice retained 146 variants; the exact 31-column preflight processed 16,862 masks with 15 CPU workers in 28.3 seconds and found 735 unestimable masks covering 2,997 traits. The 6,656,382,998-byte global cache retained 175,108 traits. CUDA FP64 completed 25,565,768 comparisons in 27 seconds (approximately 0.94 million/s) and wrote 58,038 retained results; the first run took 241.4 seconds including metadata and one-time cache work. An identical resume reused variant QC, the signed 3.95-GB expression scan, raw predictor cache, global trait cache, and all association comparisons, finishing in 78.3 seconds (27.5 seconds of which was deterministic rank/model preflight).
- Do not launch the full incomplete-trait chr22 analysis yet. The bounded profile required 338 backend calls, 34.0 GiB of uploads, 3.7 GiB of downloads, and 8.76 GPU-compute seconds; end-to-end association throughput was far below the complete-trait baseline. Add an on-device pattern-statistics finalization/reduction (or an equivalent compact exact formulation) so `X'g` intermediates and repeated design masks are not transferred for every genotype block. Reuse device allocations/design state, profile larger 256–2,048-variant blocks, and retain FP64 pattern-outer/reference agreement before authorizing the expensive full run.
- Add durable per-pattern aggregate variant-QC counts to genotype-outer (no-call, pattern-monomorphic, and optional pattern MAF/MAC) before enabling `--frequency-scope pattern` there. FP32 genotype-outer also remains intentionally disabled pending an accuracy study; FP64 is the current scientific path.
- Benchmark `local-pattern` genotype imputation against row-mean and standard reference-panel imputation. Define maximum distance, minimum comparable flanks/donors, phased-input behavior, and richer per-imputation QC before treating the proxy as a production default. Mean remains the default.
- Consider an explicit global complete-sample policy across predictor, trait, and selected covariates. Keep it distinct from `exclude-row` feature removal and exact trait-pattern deletion.

## 5. Cohort-aware joint analysis and covariance models

- Add a cohort-partitioned joint-analysis mode for FHS, WHI, CARDIA, JHS, and similar studies. Each cohort must accept its own fixed-effect design, including cohort-specific technical covariates, batches, exams, ancestry covariates, and intercept.
- Project both genotype and expression through the same cohort-specific nuisance model before combining samples. Treat literal stacking of transformed sample columns and blockwise accumulation of score/information statistics as equivalent implementations, and add a deterministic test proving that equivalence.
- Support covariance-aware pre-whitening where observations are not independent: an aligned pedigree/GRM block for FHS-like familial cohorts and a subject-level correlation model for CARDIA-like repeated exams. Do not treat subtraction of phenotype BLUPs alone as independent residualization.
- For repeated measures, support a documented one-phenotype-per-subject aggregation policy when a stable average QTL effect is the intended estimand. Never absorb a time-invariant genotype effect with subject fixed effects.
- Keep ordinary fixed-effect projection for independent cohorts such as WHI and JHS. Model a small number of batches or cohorts as categorical fixed effects rather than automatically estimating random-intercept variances.
- Define the cross-cohort estimand and scaling explicitly. Preserve allele dosage for a per-allele effect; require harmonized effect alleles and comparable expression units or documented cohort-specific precision scaling. Do not silently substitute equal-cohort, sample-size, or standardized-genotype weighting.
- Stream cohort contributions by variant/trait tile and write only combined results. Track cohort availability/effective sample information and optionally accumulate a common-effect heterogeneity statistic; run cohort-specific or leave-one-cohort-out diagnostics only for retained findings.
- Add deterministic CPU fixtures covering different covariate sets, unequal cohort sizes and residual variances, allele-frequency differences, missing cohort/variant combinations, related families, repeated subjects, rank/DF accounting, stacked-versus-streamed equality, and failure on ID or allele mismatches before GPU acceleration.

## 6. Other statistical analyses

- Add SNP interaction analysis with an explicit model definition, degrees-of-freedom tests, deterministic CPU references, pruning controls, and bounded block-pair scheduling.
- Add forward selection on top of indexed/re-readable genotype sources, including stopping rules, conditional-model degrees of freedom, cache reuse, and reproducible tie handling.
- Repair and separately verify the disabled categorical-SNP analysis path. Automatic categorical covariates do not make that path correct.

## 7. Performance and operational follow-up

- Validate automatic GPU block/worker tuning on real multi-GPU machines and on larger 5,100–5,700-sample cohorts. Revisit the 1-GiB output target and four-worker pipeline only with representative measurements.
- Profile GPU fixed-effect residualization at 5,100–5,700 samples, higher covariate ranks, and on real multi-GPU/Intel/AMD systems. The 2,005-sample NVIDIA study was cache-I/O-bound even though projection compute was fast; do not claim a universal wall-time gain.
- Benchmark FP32 versus FP64 on Intel and AMD hardware. The larger NVIDIA WHI study found one reporting-threshold classification difference, so retain FP64 verification for borderline FP32 findings.
- Add cache inspection and safe stale-cache pruning commands; never delete caches merely because they are not used by the current run.

## 8. Other input and output work

- Extend validation to family and pedigree identifiers if those inputs remain supported.
- Add optional compressed result output and decide whether the legacy full-memory scheduler should also assemble deterministic block-ordered output.
- Generalize legacy association labels such as `Rs_ID` and `ProbesetID` through an opt-in output schema for methylation, proteomics, and other declared matrix roles without breaking existing parsers.

## 9. GPU backend validation and release work

- Validate JOCL/OpenCL on AMD FP64 hardware before considering a native HIP/ROCm backend. Add HIP/ROCm only if representative end-to-end measurements justify its Java bridge and packaging burden.
- Validate the JOCL FP32 path on an Intel Iris Xe. Reassess Level Zero only if it materially improves a verified FP32 workflow or supports scientifically usable FP64 hardware.
- Consolidate shaded-jar metadata warnings, document supported Java/driver/runtime combinations, and produce reproducible per-platform artifacts, checksums, licenses, and release notes.
