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
  --trait-positions transcript-tss.tsv --mixed-cis-window 1000000 \
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

The participant pedigree graph and its connected covariance blocks are built
once. Each observed mask subsets and diagonalizes those blocks without a pooled
dense projection. Storage is quadratic in each connected pedigree block;
one large connected pedigree can still require substantial memory. The most
recent mask geometry is cached, with no unbounded cache of distinct masks.
Each trait receives its own variance fit. Score information includes the full
fixed-effect correction, and missing genotype means use the trait's observed
participants. Zero and one are observed trait values. Missing parser sentinels
and NaN exclude participants. A mask losing a nuisance factor level retains the
remaining design column space with its actual rank. Unestimable traits fail
unless explicitly skipped; the paired `.pooled.tsv` audit records their status,
observed cohort counts, rank, variance components and convergence diagnostic.

Each component is bounded between `1e-10` and `1e10` times the trait's
fixed-effect residual mean square. If the observed pedigree kernel is exactly
identity, only its sum with that cohort's residual variance is estimable; the
sum retains bounds `2e-10` to `2e10`, and additive variance is reported as zero.
Other zero or dependent covariance components after fixed-effect projection
are explicitly unestimable (normalized Gram eigenvalue tolerance `1e-10`).
Convergence requires maximum projected log-variance score at most `1e-6`.
Score products use normalized response units; effects, standard errors,
variance components and likelihoods are converted back to the supplied units.
Unrepresentable FP64 values and numerical convergence failures stop the run.

The raw aligned genotype cache is indexed. With cis options, variant identifiers
must have a `CHROM:POS:` prefix. Trait positions are an exact-ID TSV headed
`ID`, `CHROM`, `POS`, with 1-based coordinates; callers must provide the intended
strand-aware TSS. No transcript-version matching is inferred. Without the cis
options, all supplied variants are tested. Atomic result/audit trait parts and
input/model/backend signatures protect compatible restarts. This mode has a
distinct checkpoint identity from existing mixed results. Release v1.9.2 uses
pooled signature version 2; draft PR #20 version-1 result directories must not
be reused. Aligned raw genotype caches and existing legacy/spectral signatures
are unchanged. `--validate-only` checks exact cis IDs, observed-mask rank/DF and
covariance identifiability without fitting nulls or writing association results.

## Validation and practical limits

Tests compare complete and incomplete masks with independent base-R dense
REML/GLS references, response scaling from `1e-100` to `1e100`, identity
reduction, projected covariance confounding, partial FP64 product tiles, inclusive
cis selection, indexed-cache integrity and byte-identical result/audit resume.
Native CPU, NVIDIA CUDA and OpenCL products are compared with scalar references.
Release-wide and packaged validation is recorded in the v1.9.2 release receipt.

This opt-in model has synthetic validation; no full-cohort speedup or production
analysis is claimed. Null fits are not persisted for reuse across chromosome
jobs. Many distinct masks require repeated eigendecompositions and fits; measure
time and memory on a representative subset before scheduling a large analysis.
Only one cohort can carry pedigree covariance; other cohorts must be independent.
FP32, GRMs, extra random kernels, set tests and per-predictor REML are unsupported
in this path. Existing mixed and spectral set-test modes remain available.

The independent reference generator is `scripts/generate-pooled-reference.R`.
It uses base R dense covariance Cholesky, rather than the production block
eigendecomposition, and commits synthetic fixtures only.
