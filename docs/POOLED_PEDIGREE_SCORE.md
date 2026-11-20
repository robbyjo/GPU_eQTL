# Pooled pedigree-block score analysis (development draft)

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

Per-mask pedigree components are diagonalized without a pooled dense projection.
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

## Validation and work remaining

The initial targeted tests compare complete and incomplete masks with independent
base-R dense REML/GLS references, response scaling, an identity reduction,
partial product tiles, direct cis selection and byte-identical result/audit
resume. Native CPU and NVIDIA OpenCL products passed locally; the CUDA test
skipped when its runtime was unavailable. These are synthetic checks.

This is a development draft, not a release or a demonstrated full-cohort speedup.
Before production use:

- Precompute the participant pedigree graph once, avoiding repeated all-pair
  relationship discovery for distinct transcript masks.
- Review covariance identifiability and component-bound handling for degenerate
  masks, and broaden independent numerical/gradient and malformed-input tests.
- Add reusable per-trait null-fit persistence across chromosome jobs and measure
  the tradeoffs for many distinct missingness masks.
- Complete full regression, packaged CLI/CUDA checks, representative timing and
  memory measurements, source/checkpoint review and user documentation.

The independent reference generator is `scripts/generate-pooled-reference.R`.
It uses base R dense covariance Cholesky, rather than the production block
eigendecomposition, and commits synthetic fixtures only.
