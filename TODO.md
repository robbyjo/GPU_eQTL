# GPU eQTL modernization to-do list

This list contains work remaining after the correctness fixtures, CLI, sample alignment, categorical covariates, bounded-RAM caches, checkpointing, optional FP32, backend-aware tuning, missingness policies, VCF/BCF input, aligned-sample variant QC, association progress, reusable VCF/BCF preprocessing, initial profiling, and OpenBLAS/pure-Java CPU-backend milestones. Preserve the scientific and verification rules in `AGENTS.md`.

Recommended dependency order after closing the production VCF/VCF.gz milestone:

`set-test foundation -> burden -> GPU burden -> SKAT -> SKAT-O`

Bounded CPU-output profiling and release validation should proceed alongside that main sequence. Cohort covariance models must be established before relatedness-aware set tests.

## 0. Completed: production VCF/VCF.gz and QC milestone

- Delivered: resumable aligned-sample QC with atomically committed ordered parts, signatures covering input metadata/policies/sample order, source-prefix or indexed-boundary validation after interruption, and direct reuse of completed scans. Production-test interruption/restart on a representative chromosome and record checkpoint storage/reread overhead.
- Delivered: approximately 15-second association progress based on actual variant–trait comparisons, including resumed-work accounting, percentage, elapsed time, throughput, and ETA across full-memory, streamed, multi-backend, and trait-pattern paths.
- Delivered: `--preprocess-only` aligned VCF/BCF QC and checksummed raw-dosage cache creation without backend initialization or association; matching later analyses discover and reuse the cache. Production chr22 validation reused a completed 15,979,340-record QC checkpoint in 71.81 seconds, opened the 440,406-row genotype cache in 0.96 seconds and the 114,406-row expression cache in 0.25 seconds, and retained 37,824,081,886 bytes (35.23 GiB) across three matrix-cache files.
- Delivered for VCF/VCF.gz: low-allocation worker parsing of raw FORMAT/GT/DS/FT columns, header-index sample selection, four bounded in-flight records per worker, HTSJDK-reference byte-equivalence tests, and pipeline diagnostics; this replaced a slower first attempt that materialized full HTSJDK genotype objects concurrently. A 4,746-sample/1,500-variant bounded benchmark plateaued at eight workers, so automatic mode is capped at 8 while explicit overrides remain available. Unrestricted unsorted-header VCF checkpoint resumes avoid HTSJDK's eager indexed decoder by validating the completed raw prefix before returning to parallel QC; indexed region scans retain direct seeks. The corrected production resume reached approximately 17,600 records/s versus 1,000–1,100 records/s on the eager indexed path and approximately 1,300 records/s in the older implementation.
- Closure decision: the biallelic additive VCF/VCF.gz workflow is production-hardened for the current use case. Production indexed VCF/BCF validation against bcftools 1.24 is now recorded below; BCF reader-thread genotype expansion remains a known performance limitation.
- Delivered baseline for `--trait-cache memory` on the 4,746-sample batch1234 chr22 analysis: 50,385,088,836 variant–trait comparisons completed in 1,995.79 seconds (33.26 minutes; 25.25 million comparisons/s wall-clock), with 12,096 GPU-compute calls totaling 1,544.93 worker-seconds and CPU result/filter/write calls totaling 1,269.69 worker-seconds. Profile `auto` and `disk` only if their operational tradeoff is needed; do not repeat the expensive full analysis solely to complete a mode sweep.

## 1. Validate and package the completed CPU analysis backend

- Delivered in the current worktree: `--backend auto|cuda|opencl|cpu`; GPU-first automatic CPU fallback; one FP64/FP32 CPU context; OpenBLAS and pure-Java engines; association and residualization operations; separate BLAS/application worker controls; host-heap block sizing; CPU pipeline sizing; native/fallback numerical tests; and one universal shaded jar containing Windows x64, Linux x64/ARM64, and macOS x64/ARM64 OpenBLAS natives.
- Delivered: forced OpenBLAS, portable Java, and CUDA FP64 completed the same 513-variant by 129-trait, 4,746-sample analysis with partial predictor/trait tiles and a p=0.05 threshold. All retained the same 5,115 ordered rows and exact IDs/N/DF; maximum OpenBLAS/Java log10-p differences from CUDA were `3.55e-14`/`2.70e-13`. Every backend assembled two checkpoint parts and reproduced byte-identical backend-specific output on resume. CUDA compute was 0.168 seconds versus 0.429 OpenBLAS and 0.640 Java on the RTX 2080 fixture.
- Smoke-test the shaded jar on every bundled platform and publish checksums. If the universal jar is inconveniently large, add reproducible classifier-specific release profiles while retaining one simple extract-and-run choice.
- Delivered: selectable oneMKL 2026.1 FP64/FP32 association and residualization on x64 Windows/Linux, automatic oneMKL -> OpenBLAS -> Java fallback, numerical reference tests, isolated platform profiles, a GPLv3 Section 7 oneMKL linking exception, and exact packaged Intel oneMKL/OpenMP plus incorporated-code notices. A production-derived 8,193-by-1,025 tile schedule over 4,746 samples compared two timed runs per engine: mean compute was 3.111 seconds for oneMKL and 3.506 for OpenBLAS, while mean analysis wall was 5.062 and 4.922 seconds respectively. Outputs matched in IDs/order/N/DF with FP64-scale differences. The mixed end-to-end result does not justify a release-default change; the default release remains free of Intel-licensed oneMKL native binaries.

## 2. Deferred indexed/BCF and region extensions

- Delivered: tabix-indexed BGZF VCF and Tribble-indexed VCF/BCF interval queries; repeatable inline and TSV region/set definitions; explicit one-based/BED coordinates; deterministic contig aliases/order; merged query intervals with overlapping set membership; empty-set reporting; source/index/region checkpoint signatures; mutation checks; and direct indexed seek to an interrupted QC boundary.
- Delivered production indexed-region validation against bcftools 1.24 on a 4,059-record chr22 slice converted to BCF. Overlapping intervals yielded 3,321 canonical records; tabix VCF and Tribble-indexed BCF matched bcftools IDs/order and memberships exactly, matched AC/AN/NS exactly after the 4,746-sample alignment, agreed in displayed MAF/HWE to `5e-7`, and produced byte-identical application QC. An injected BCF interruption left two durable parts, sought directly to the indexed boundary on restart, and reassembled the same QC bytes.
- Add a pure-Java standard CSI decoder (or a carefully packaged cross-platform htslib bridge) for ordinary bcftools `.csi` files. HTSJDK 5 currently handles tabix/Tribble indexes but not standard variant CSI; unsupported CSI must remain a clear pre-analysis failure, never a silent sequential fallback.
- Standard CSI support was considered after that validation. Keep it deferred until a pure-Java decoder or carefully packaged htslib bridge can preserve the same deterministic contract; `.csi` region requests now fail with an explicit supported-index message and never substitute a full scan.
- Add an explicit multiallelic split policy only after tests define per-ALT DS/GT projection, other-ALT calls, canonical IDs, allele counts, missingness, and HWE behavior. Current choices remain `exclude` and `error`.

## 3. Rare-variant burden, SKAT, and SKAT-O analyses

### 3a. Shared set-test foundation

- Delivered initial CPU contract/reference: strict tab-separated variant-to-set definitions; exact REF/ALT and explicit effect-allele orientation without implicit swaps/strand flips; positive explicit/default-one weights; deterministic overlapping memberships; aligned-cohort min/max MAF/MAC masks; mean/zero/error genotype missingness; reusable full-rank continuous-trait FP64 null models; and explicit error/skip behavior for absent, empty, and post-projection monomorphic sets. Definitions have a deterministic content signature.
- Delivered deterministic eight-sample FP64 fixtures covering retained IDs/order, overlapping membership, REF-effect orientation, non-unit weights, mean-filled dosage, MAF exclusion, effective N, residual DF, R-squared, effect, t, log10 p, empty/monomorphic audit states, allele mismatch, duplicate membership, missing-dosage failure, and covariate-rank failure. R-squared/effect/t values were independently calculated with NumPy QR before being anchored to the Java reference.
- Delivered: aligned CSV and VCF/BCF production adapters, native one-based sliding windows with configurable size/stride, exact ALT-effect adaptation of automatic/custom memberships, stable set-audit/result schemas, and CLI/INI modes for burden, SKAT, and SKAT-O.
- Delivered: bounded set/trait tiling, retained-result filtering, and signature-bound atomic checkpoint/restart without materializing complete variant-by-trait QTL files.
- Delivered: end-to-end automatic-window CSV/VCF/BCF byte equality for burden, SKAT, and SKAT-O plus byte-identical CSV checkpoint resume. Sliding runs reuse a version-1 checksummed aligned raw cache and a derived signature/file-identity-bound row-offset sidecar; requested rows retain source order and CRC validation, adjacent rows are read sequentially in bulk, and trusted buffers transfer without a duplicate dosage copy.
- Delivered full chr22 CPU burden profiles on 428,629 variants, 4,746 samples, two traits, and explicit 256-set tiles. Warm-cache 10-kb/5-kb, 50-kb/10-kb, and 100-kb/50-kb grids emitted 7,453/3,808/765 windows in 182.14/229.10/148.35 seconds and observed 3.86/6.30/7.91 GB JVM heap at tile boundaries. Heap-aware omitted/zero `--set-block-size` now accounts for aligned samples, method workspaces, and actual membership density. On 100-kb/50-kb it selected 70, completed in 172.15 seconds at a 5.41-GB observed peak, reproduced the same result-row multiset and byte-identical audit, and preserved the 16.29-GB raw cache.

### 3b. Burden tests first

- Delivered production weighted/unweighted FP64 burden execution. The batched optimized-CPU set-by-trait path reuses prepared null state and is tested directly against the scalar reference over partial set/trait tiles.
- Direct GPU burden products remain optional future work after profiling demonstrates value over the batched CPU path.

### 3c. SKAT

- Delivered deterministic FP64 CPU SKAT with weighted projected kernels, batched score/covariance products, CPU eigenvalues, exact single/equal-eigenvalue chi-square cases, bounded Imhof inversion, and an explicitly reported Satterthwaite fallback.
- Direct multi-GPU score/covariance products remain optional future work; keep eigenvalue and p-value work on the CPU.

### 3d. SKAT-O

- Delivered a configurable rho grid (default `0,0.25,0.5,0.75,1`), weighted burden/SKAT mixture construction, shared-Gaussian correlated parametric-null replicates, fixed seed/signature, and adjusted omnibus p-value. The minimum component p-value is audit-only.
- Delivered additional statistical hardening fixtures for singleton exact scaled chi-square behavior, nearly collinear variants, extreme finite weights, an explicitly forced bounded-Imhof fallback, and SKAT-O Monte Carlo resolution/determinism. Retain these as gates before enabling any direct GPU set-test product.
- Defer relatedness-aware burden/SKAT/SKAT-O until the cohort covariance/null-model design below is implemented and verified.

## 4. Missing-data production hardening

- Delivered: exact subset mean/EAF/MAF/MAC/squared-dosage statistics; pattern-specific mean genotype filling; per-pattern monomorphic skipping; optional `--frequency-scope pattern`; compact pattern QC summaries; a checksummed aligned raw dosage cache; and persistent statistics/prepared caches keyed by source, alignment, policy, and sample mask.
- Delivered: a deterministic end-to-end FP64 CPU association reference for exact trait-pattern deletion covering identifiers, effective N, residual degrees of freedom, R-squared, effects, t statistics, and p-values.
- Delivered: durable ordered pattern-result/QC groups, nested genotype-block resume, signature validation, atomic final assembly, completed-group reuse, retained-checkpoint support, and a failure/restart test proving byte-identical output.
- Delivered: fail-fast preflight reporting pattern count, observed-N range, repeated predictor-preparation amplification, upper-bound prepared numeric storage, and structurally non-positive-DF patterns. `--max-trait-patterns` defaults to 256 and requires an explicit decision for larger schedules.
- Delivered: checksummed, atomically written exact matrix missingness scans keyed by source path/size/mtime, parser/cache tag, dimensions, selected sample identities, and exact selected order. Matching reruns regenerate the human-readable audit without rereading the full matrix; corrupt or incompatible scans rebuild.
- Delivered first scalable genotype-block-outer scheduler: one globally prepared padded trait cache; per-mask FP64 covariate QR/rank/DF/threshold state; GPU/CPU sufficient-statistics products for exact predictor projection and variance; exact per-mask mean/zero predictor imputation; pattern-monomorphic skipping; one compute-context-exclusive worker per selected device; deterministic genotype-block/trait-block/variant/original-estimable-trait output; signed block resume; and audited `error|skip` handling for unestimable masks. Deterministic tests match the established pattern-outer VCF implementation for complete, mean-imputed, and zero-imputed predictors and prove byte-identical restart.
- Delivered bounded production validation on batch1234: a tabix-indexed 10-kb chr22 slice retained 146 variants; the exact 31-column preflight processed 16,862 masks with 15 CPU workers in 28.3 seconds and found 735 unestimable masks covering 2,997 traits. The 6,656,382,998-byte global cache retained 175,108 traits. CUDA FP64 completed 25,565,768 comparisons in 27 seconds (approximately 0.94 million/s) and wrote 58,038 retained results; the first run took 241.4 seconds including metadata and one-time cache work. An identical resume reused variant QC, the signed 3.95-GB expression scan, raw predictor cache, global trait cache, and all association comparisons, finishing in 78.3 seconds (27.5 seconds of which was deterministic rank/model preflight).
- Delivered: backend-neutral compact pattern-statistics plans plus CUDA/OpenCL on-device finalization. Each exclusive context uploads the packed genotype aggregates once per genotype block, reuses device allocations across pattern batches, retains `X'g` intermediates on-device, and returns only replacement, residual sum of squares, and filled sum of squares per pattern/variant. CPU retains a portable reference fallback. CUDA uses a fail-fast NVRTC-compiled finalizer; OpenCL compiles the equivalent FP64 kernel through the selected ICD.
- Delivered bounded CSV block sweep on the first 2,048 rows of `chr22-filtered-batch1234-dose-meanimp.csv` against all 16,862 production trait masks and 175,108 estimable traits. Blocks 256/512/1,024/2,048 retained the same 876,057 identifiers and N/DF; 512/1,024/2,048 were numerically identical by identifier, while 256 differed only by FP64 roundoff (maximum `7.67e-13` T and `1.23e-13` log10 p) with no threshold-set changes. Analysis-wall times were 303.76/230.43/234.16/240.42 seconds; 512 rows is the current measured choice on the RTX 2080 fixture, not a universal default. Download was 6.11 GiB for every run, including association matrices; the compact pattern component was approximately 0.76 GiB across eight 256-row blocks instead of roughly 25 GiB of legacy `X'g` intermediates. A full 440,406-variant incomplete-trait run remains an approximately half-day operation at the measured rate, so launch it only when that cost is acceptable; further work should reduce repeated expanded-design upload/GEMM work rather than re-expand host output.
- Delivered durable per-pattern aggregate genotype-outer variant QC. Compact CUDA/OpenCL/CPU finalizers now return called counts and dosage sums with filled sufficient statistics; successful genotype blocks atomically commit matching checksummed QC parts, and deterministic assembly reports input/included/no-call/monomorphic/missing-genotype plus candidate MAF/MAC counts per pattern. Fresh and resumed aggregates match the established pattern-outer oracle. Candidate pattern thresholds are deliberately audit-only until filtering is enabled; `--frequency-scope pattern` and FP32 genotype-outer remain intentionally disabled.
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
