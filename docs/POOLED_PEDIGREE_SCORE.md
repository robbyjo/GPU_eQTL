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

## Identifiability and convergence

The response is projected onto the fixed-design residual space and normalized
by its residual sample SD. Component bounds are 1e-10 to 1e10 in these normalized
variance units. At most two deterministic starts share the total
`--reml-max-iterations` budget. Acceptance requires the maximum absolute projected
log-variance score to be at most 1e-6. Estimates and likelihood are returned in
original response units. This optimizer and its bounds are separate from the
legacy fitter. Numerical nonconvergence is fatal even when unestimable traits
are explicitly skipped.

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
missing/modified committed parts and incompatible inputs/model/backend, including
trait-position changes. Uncommitted parts are retained in a sibling
`.uncommitted` directory before recomputation. One OS lock owns the checkpoint.
Both final outputs must assemble before default cleanup; `--keep-checkpoints`
retains committed parts. Output order is trait then original predictor order.
The algorithm namespace differs from the initial draft and from other mixed
modes. Existing legacy signatures remain unchanged.

The VCF/BCF reader's aligned-sample `hwe_p` diagnostic is distinct from an external
full-cohort `P_HWE` annotation. Join that independently prepared annotation using
the exact allele-aware variant key and retain its reference sample definition.
This fitter does not recalculate HWE per transcript or introduce an HWE filter.

## Validation and limits

Synthetic tests compare complete/incomplete masks, likelihoods, gradients and
scores against independent base-R dense covariance Cholesky and finite
differences. They cover response scaling, inbred proportional-identity bounds,
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
`javac -cp target/gpu-eqtl-1.9.1.6-all.jar -d target/benchmark-classes test/gov/nih/eqtl/PooledSyntheticBenchmark.java`
and run `gov.nih.eqtl.PooledSyntheticBenchmark FRESH_OUTPUT_DIRECTORY` with
`target/benchmark-classes` plus the shaded JAR on the Java classpath, native CPU
algebra selected, and an available CUDA runtime. It uses 8,955 synthetic samples,
55 design columns, a 512-person synthetic pedigree component, complete/incomplete
masks, a shared-null hit and a scalar/CUDA product comparison. It reports its
observed timings as JSON; those timings are not an end-to-end cohort estimate.
