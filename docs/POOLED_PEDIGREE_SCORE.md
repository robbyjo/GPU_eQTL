# Pooled pedigree-block score analysis

This opt-in single-variant mode fits one Gaussian REML null per trait, with an
additive pedigree covariance in one named cohort and a separate residual
variance in every cohort. It uses exact observed-sample masks, FP64 score
products and chi-square(1) score inference. The existing legacy mixed runner
and related-set spectral fitter remain separate.

```sh
java -Xmx32g -jar gpu-eqtl.jar \
  --genotype pooled.bcf --expression pooled-jointtmmfpkm-ratio-winsor.csv.gz \
  --covariates pooled-covariates.csv \
  --genotype-id-column SampleID --expression-id-column SampleID \
  --fixed-covariates Cohort,Age,Sex,WBC,LY_PER,MO_PER,EO_PER,BA_PER,FHS_LevySet_X1,FHS_LevySet_X2,FHS_LevySet_X4,WHI_Batch_x2,CARDIA_visit_Y25,PC1,PCE1 \
  --pedigree fhs-pedigree.tsv --pedigree-id-column PedigreeID \
  --mixed-model score --mixed-null-fitter pedigree-block \
  --residual-group-column Cohort --pedigree-group FHS \
  --trait-missing pattern --predictor-missing mean \
  --unestimable-trait-patterns skip \
  --trait-positions transcript-tss.tsv --mixed-cis-window 1000000 \
  --mixed-null-cache-dir shared-pooled-nulls \
  --reml-max-iterations 1000 --backend cuda --precision fp64 \
  --keep-checkpoints --output pooled-score.csv
```

The PC list is abbreviated in this example; supply every prespecified PC.
`Cohort` is a categorical fixed effect. Nested batch indicators are zero outside
their cohort: FHS indicators encode the nonreference `Levy_Set` levels, WHI's
indicator encodes `Batch=x2`, and CARDIA's encodes `Collection_visit=Y25`.
Use the actual observed levels and inspect rank rather than inferring a
numeric trend from batch labels. These labels are example configuration, not
hardcoded engine behavior.

The covariance is `blockdiag(a A_FHS + e_FHS I, e_WHI I, e_CARDIA I)`.
`A` is twice kinship, not a family random-intercept kernel. The pedigree only
needs individuals in the selected pedigree group and their required ancestors.
All analysis sample IDs remain unique. Other cohorts are modeled as independent;
their duplicates/close relatives must be resolved before selecting this model.

Pedigree relationships and covariance blocks are constructed once. Per-mask
submatrices are diagonalized without a pooled dense projection; one exact mask
geometry is held in memory. Large-family rotations and repeated weighted design
calculations use the resolved native CPU algebra when available. No eigenvalues
are truncated to approximate the model; roundoff-sized negatives may be clamped.
Each trait receives its own variance fit. Score information includes the full
fixed-effect correction, and missing genotype means use the trait's observed
participants. Zero and one are observed trait values. Missing parser sentinels
and NaN exclude participants. A mask losing a nuisance factor level retains the
remaining design column space with its actual rank. Unestimable traits fail
unless explicitly skipped; the paired `.pooled.tsv` audit records their status,
observed cohort counts, rank, variance components and convergence diagnostic.

The raw aligned genotype cache is indexed. With cis options, variant identifiers
must have a `CHROM:POS:` prefix. Trait positions are an exact-ID TSV headed
`ID`, `CHROM`, `POS`, with 1-based coordinates; callers must provide the intended
strand-aware TSS. No transcript-version matching is inferred. Without the cis
options, all supplied variants are tested. Atomic result/audit trait parts and
input/model/backend signatures protect compatible restarts. This mode has a
distinct checkpoint identity from existing mixed results.

Association products use the orthonormal columns of the weighted fixed design.
Score information is the weighted predictor norm minus the squared projection
norm, without multiplying by the potentially ill-conditioned inverse information
matrix. When subtraction would remove at least 99% of the raw norm, the equivalent
explicit residual norm and residual score are evaluated and reprojected once to
remove roundoff along the fixed design. This numerical branch changes neither
the score law nor the existing `information <= 1e-12 * raw norm` collinearity rule.
Nonfinite products and materially invalid information still fail explicitly.

## Identifiability and convergence

The response is projected onto the fixed-design residual space and normalized
by its residual sample SD. Component bounds are 1e-10 to 1e10 in these normalized
variance units. At most two deterministic starts share the total
`--reml-max-iterations` budget. Acceptance requires the maximum absolute projected
log-variance score to be at most 1e-6. Estimates and likelihood are returned in
original response units. This optimizer and its bounds are separate from the
legacy fitter. Numerical nonconvergence is fatal even when unestimable traits
are explicitly skipped.

Each likelihood evaluation whitens the fixed design and response, then uses
thin Householder QR for the residual and information determinant. This avoids
squaring the weighted design's condition number when a cohort residual variance
reaches its lower bound. The calculation uses native LAPACK or the portable Java
reference without constructing an `N x N` projector. A restart satisfying the
score criterion is accepted even if an earlier nonstationary fit has a slightly
smaller objective due to floating-point roundoff. A line search that cannot move
any parameter proceeds to the remaining start instead of spending its budget on
identical iterates.

The line search computes the REML objective change directly from covariance
ratios, using `expm1`/`log1p` and a small Cholesky factor of the relative
fixed-effect information. This avoids subtracting independently rounded full
likelihoods, which can reject useful steps near a variance boundary and stall
above the score criterion. With `A = orth(sqrt(W) X)`, `T = W_new/W - I` and
the current whitened residual `e`, the residual-sum-of-squares change is
`e' T e - (A' T e)' (I + A' T A)^-1 (A' T e)`. The information determinant
change is `log det(I + A' T A)`. No pooled `N x N` matrix is formed.
The Armijo rule and roundoff allowance, total iteration budget, component bounds
and strict `1e-6` convergence criterion are unchanged; an unresolved fit remains
fatal. Absolute likelihoods and final score statistics retain their existing
weighted-QR and orthogonal-score calculations.

If a selected pedigree covariance is exactly `c I`, its additive and residual
components are not separately identifiable. The fitter uses their identifiable
total with the exact summed original bounds and labels it `total:FHS`.
`unidentified-additive:FHS=0` is bookkeeping, not evidence of absent genetic
variance. Masks without FHS use the corresponding exact reduction. Other
covariance dependencies after fixed-effect projection fail explicitly.

## Shared nulls and restart integrity

Use distinct output/checkpoint paths for chromosome jobs and the same explicit
`--mixed-null-cache-dir` to share compatible nulls. By default this directory is
the output's cache directory plus `pooled-nulls`. Each immutable null has a
SHA-256 payload checksum and a key covering the trait ID, response bits, ordered
sample/pedigree/group identities, exact mask, fixed-design geometry, rotations,
bounds, algorithm, iteration budget and resolved CPU algebra. Predictor data are
absent because the null contains no predictor. Stable OS locks prevent duplicate
cross-process fitting/writing. Corruption fails visibly. `--rebuild-cache`
recomputes without overwriting a prior immutable fit. Geometry must still be
constructed to validate a shared null, so this cache does not eliminate all
per-mask eigendecomposition work.

Result and audit parts carry separate model-bound SHA-256 records. Resume rejects
missing digest indexes, missing/modified committed parts and incompatible
inputs/model/backend, including trait-position changes. Uncommitted parts are
retained outside the checkpoint in its sibling `.uncommitted` directory before
recomputation; retained audits use that directory's `audit` subdirectory, so
cleanup preserves all interrupted evidence. One OS lock owns the checkpoint.
Both final outputs must assemble before default cleanup; `--keep-checkpoints`
retains committed parts. Output order is trait then original predictor order.
This implementation uses pooled result signature version 6 and the new
`pooled-null-v5-relative-reml` algorithm key. Pooled results/checkpoints from
v1.9.2.4 and earlier require fresh output paths. Earlier immutable null files
remain preserved under their old keys and are not reused by the changed fitter;
the binary null format is unchanged. Preserve all earlier outputs and diagnostics.
Existing legacy/spectral signatures and aligned raw-cache formats are unchanged.
The released validation-only checks for exact IDs, observed-mask rank/DF and
covariance identifiability remain available without fitting nulls.

The VCF/BCF reader's aligned-sample `hwe_p` diagnostic is distinct from an external
full-cohort `P_HWE` annotation. Join that independently prepared annotation using
the exact allele-aware variant key and retain its reference sample definition.
This fitter does not recalculate HWE per transcript or introduce an HWE filter.

## Validation and limits

### Optional genome-wide trait batching

For complete predictors, add `--mixed-trait-batch-size 32` (maximum 256) to a
genome-wide `--mixed-null-fitter pedigree-block` score run. The default is 1,
which preserves the existing runner and result-checkpoint identity. Batching
requires `--predictor-missing error` and cannot be combined with cis-only
`--trait-positions`. Each exact trait sample mask keeps its own geometry and
per-trait REML null; a mask change flushes the current batch.

The batched runner caches the original genotype centering and pedigree rotation
under `centered-rotated-v1` in the cache directory. Its identity includes the raw
input/sample order and exact geometry/algebra. Per-key process locks protect
creation, and existing corrupted rows fail their checksums rather than being
silently rebuilt. Budget disk space for an additional FP64 genotype matrix for
each mask. High-entropy missingness can make this unsuitable.

Score columns and weighted squared norms use bounded FP64 matrix products.
Nearly collinear/boundary cases use the original explicit residual and
reprojection calculation. Effect/SE/statistic range checks still run for every
pair; the chi-square tail is skipped only safely below the output threshold.
Selected pairs still pass the original log-p cutoff. The null fit, covariance
bounds, `1e-6` convergence criterion and chi-square(1) inference are unchanged.

Use fresh result/checkpoint paths when enabling batching or changing its size.
Compatible immutable nulls can be reused. Results and audits remain separate
checksummed atomic trait parts; an interrupted batch can leave some completed
traits, and resume reconstructs missing paired audits while retaining partial
evidence. Batches are submitted sequentially through an exclusive backend
context; this option does not add a concurrent multi-GPU scheduler. Explicit
heap/device checks reject oversized tiles or batches with instructions to
reduce them. Heap admission retries once after requesting garbage collection,
so unreachable prior tiles are not mistaken for live allocations; the same
headroom requirement still applies. Up to 256 checkpoint writers are opened
iteratively, avoiding a call stack proportional to the trait batch size.
Every writer must close successfully before new parts receive digest commits.
Their combined character buffers are capped at 16 MiB, with encoder/writer
headroom included in heap admission; GPU admission counts device matrices only.
Performance depends on input parsing, cache reads, trait batch
size, null fitting and available hardware; no full-cohort timing is promised.

### Reference checks

Synthetic tests compare complete/incomplete masks, likelihoods, gradients and
scores against independent base-R dense covariance Cholesky and finite
differences. They cover response scaling, inbred proportional-identity bounds,
unequal component likelihood changes against dense-R references for complete
and incomplete masks in both directions,
covariance confounding, native/scalar large-family rotations, inclusive cis
boundaries, null reuse/corruption, raw-row CRC validation and result/audit
resume/corruption. Actual CUDA, NVIDIA OpenCL and native CPU product checks
include partial tiles. These checks do not establish full-cohort throughput or
production-server performance.

Supported combinations are single-SNP Gaussian score tests, one observation per
pedigree individual, one pedigree group and FP64; no GRM/additional random kernels
or supplied-covariance cohort mode. Predictor missing policies are error/mean/zero;
frequency filters use aligned participants, with pattern frequency filtering
rejected. Trait-pattern scheduling is handled by this runner. Null fits and
mask geometry remain CPU work. Many distinct masks or one large connected
pedigree can remain expensive. PC reconstruction, cohort-wide HWE annotation and
scientific calibration are analysis-workflow responsibilities.

The independent reference generator is `scripts/generate-pooled-reference.R`.
It uses base R dense covariance Cholesky, rather than the production block
eigendecomposition, and commits synthetic fixtures only.

The standalone `test/gov/nih/eqtl/PooledSyntheticBenchmark.java` is an optional
engineering benchmark, excluded from ordinary JUnit execution. Compile it with
`javac -cp target/gpu-eqtl-1.9.2.7-all.jar -d target/benchmark-classes test/gov/nih/eqtl/PooledSyntheticBenchmark.java`
and run `gov.nih.eqtl.PooledSyntheticBenchmark FRESH_OUTPUT_DIRECTORY` with
`target/benchmark-classes` plus the shaded JAR on the Java classpath, native CPU
algebra selected, and an available CUDA runtime. It uses 8,955 synthetic samples,
55 design columns, a 512-person synthetic pedigree component, complete/incomplete
masks, a shared-null hit and a scalar/CUDA product comparison. It reports its
observed timings as JSON; those timings are not an end-to-end cohort estimate.
