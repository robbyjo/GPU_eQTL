# Related samples, exact trait deletion, and FP64 set products

Burden, SKAT, and SKAT-O support fitted pedigree/GRM adjustment, exact observed-sample trait masks, and CUDA/OpenCL/CPU FP64 matrix products. These are continuous Gaussian-trait tests. A ratio of **0 or 1 is observed**; only the matrix parser's missing sentinel or NaN removes a sample.

```sh
java -Xmx12g -jar target/gpu-eqtl-1.9.1-all.jar \
  --genotype dosage.csv --expression isoform-ratios.csv \
  --covariates covariates.csv --fixed-covariates Age,Sex,PC1,PC2 \
  --variant-sets genes.tsv --analysis skat-o \
  --relationship-matrix grm.tsv --mixed-model score \
  --trait-missing pattern --predictor-missing mean \
  --backend cuda --precision fp64 --skat-o-method integration \
  --threshold pval 1e-8 --keep-checkpoints --output related-skato.csv
```

Use `--pedigree pedigree.tsv` instead of `--relationship-matrix` for the existing ID/SIRE/DAM pedigree format. Both default to `--mixed-model score`, and either can accompany `--random-covariates` grouping terms. Specifying both a pedigree and GRM is rejected. `--mixed-model reml` remains a single-predictor refit mode and is rejected for set tests. Omit relationship/random terms for fixed-effect set tests. Existing set definitions, sliding windows, indexed VCF/BCF regions, filtering, explicit filling/exclusion policies, and legacy INI input remain available.

## GRM contract and fitted-null calculation

The dense UTF-8 TSV starts with `ID`, followed by unique sample column labels. Every row has its sample ID followed by one finite numeric entry per matrix column. Row/column order can differ; matching is by exact ID. All columns must have exactly one corresponding row. Extra GRM individuals are permitted. Every aligned genotype ID must be present, after any explicitly selected predictor-ID prefix removal. No automatic name normalization, pedigree-ID mapping, scaling to kinship, or allele inference is applied to a GRM.

```text
ID	s1	s2	s3
s1	1	0.5	0
s2	0.5	1	0
s3	0	0	1
```

The matrix must be symmetric and positive semidefinite within the documented numerical checks in `QMixedModel.readRelationship` (symmetry tolerance `1e-12` times the largest diagonal magnitude; eigenvalue tolerance `1e-10` times that magnitude). A residual identity component is fitted alongside it. The pedigree path supplies the numerator relationship matrix, twice kinship. Matrix scale is absorbed by its fitted variance component, so record the supplied matrix's definition in study provenance.

For each trait, REML fits `V = sum_k theta_k K_k` using the existing pinned JLinAlg CPU fitter. With fixed design `X`, the fitted null supplies

```text
P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1
U = W G' P y
C = W G' P G W
```

`G` contains filled declared-effect dosages and `W` the positive definition weights. Both scores and their covariance include fitted relatedness. Genotype columns are first projected off the ordinary fixed-effect span to reduce cancellation; this leaves `U` and `C` mathematically unchanged because `P X = 0`.

Related burden uses `(1'U)^2 / (1'C1)` with chi-square(1) inference and reports the score effect `(1'U)/(1'C1)`. Fixed-effect burden retains its regression t test, effect, R-squared, and `N-rank(X)-1` DF. SKAT uses `U'U` with the nonzero eigenvalues of `C`. Fixed-effect kernel scores use the null residual variance estimated with `N-rank(X)` DF. Related kernel scores already contain the fitted variance through `P` and must **not** be divided by an additional residual variance. Kernel output `DF` records null residual DF, not the degrees of freedom of a single chi-square distribution. Fitted-null inference remains asymptotic.

## Exact transcript sample masks

`--trait-missing pattern` is the default. Alignment and selected-covariate missingness are resolved first. Each transcript then selects its ordered observed samples. The engine subsets **every covariance basis on both axes before fitting**, reconstructs the design QR/rank/DF, and fits a new trait null. Subsetting a full-cohort fitted `P` would give a different model and is never used.

Called genotype means and monomorphic/no-call decisions are recomputed in that mask. Mean filling uses only the mask's called samples; zero filling refers to zero copies of the declared effect allele. Set MAF/MAC filtering retains the aligned-cohort default. `--frequency-scope pattern` explicitly applies set MAF/MAC limits within the observed mask, including for CSV set analyses; this changes set membership across transcripts. VCF/BCF global QC and its separate `--min-maf`/`--min-mac` semantics remain upstream. No-call or monomorphic mask-specific members are audited separately.

Preflight checks all exact masks before result checkpoint creation. `--max-trait-patterns` defaults to 256. Rank-deficient or non-positive-DF masks fail, or are audited and omitted with `--unestimable-trait-patterns skip`. A constant residual response, REML convergence failure, invalid covariance, or unresolved probability calculation is fatal even with set-degeneracy skipping. `--validate-only` checks alignment, matrix/pedigree validity and mask-specific rank/DF without fitting a null or initializing a backend.

One fitted null is resident at a time, reused only for that transcript's set tiles. Another transcript with the same mask gets its own REML fit. There is no persistent fitted-null cache to confuse different trait/model states.

## Actual backend work

The engine uses `GpuContext`/`GpuContextPool` and the production FP64 `eqtlReal` product operation. CUDA uses cuBLAS DGEMM; OpenCL uses its compiled FP64 kernel; CPU uses its selected native BLAS or Java engine.

| Operation | Execution |
| --- | --- |
| Weighted burden aggregation | Selected backend |
| Fixed-effect genotype projection products | Selected backend; `--residualization cpu` forces scalar projection |
| Fitted `P` times genotype rows | Selected backend |
| Genotype–response scores and genotype score covariance | Selected backend |
| QR/rank, mask construction, dosage filling/weighting, REML | CPU FP64 |
| Eigenvalues, quadratic-form tails, rho-grid integration | CPU FP64 |

Matrix products use at most 128-row/column tiles, reduced by `GpuTuning` for the least-capable device's allocation/total-memory limits or the CPU heap budget. They round only active kernel dimensions and extract only active cells. A context is reserved immediately before submission and released in `finally`. All eligible selected contexts participate through the pool, with sequential submissions; this is not a concurrent multi-GPU speedup claim. Each run prints device/runtime identity and **executed call counts by backend and operation**. Resume can legitimately execute zero products. A forced GPU error stops the run rather than substituting CPU results. Use `java -Deqtl.cpu.blas=java -jar ... --backend cpu ...` and the separate scalar reference for verification.

## SKAT-O joint-tail integration

The default `--skat-o-method integration` retains the configured rho grid, including exact rho 1. For each component, `Q_rho = (1-rho) U'U + rho (1'U)^2`. The minimum component p-value sets a rejection threshold for **each** component by inverting the same quadratic-form tail used to evaluate that component. The reported p-value is the probability of their correlated union, not an unadjusted minimum.

Let `U = C^(1/2) Z`, `Z = R D`, with `R^2 ~ chi-square(rank(C))` independent of uniform sphere direction `D`. The union event conditional on direction is

```text
R^2 >= min_rho critical_rho / Q_rho(C^(1/2) D).
```

The algorithm integrates this radial survival probability directly. It never counts rare exceedances or subtracts an almost-one joint CDF from one. Rank one is exact. Rank two uses adaptive angular quadrature. Equal nonzero eigenvalues use one-dimensional conditional chi-square integration. Other spectra use 16 independently shifted Halton direction rules, with four standard errors across replicates and two successive refinement checks. These are numerical error **estimates**, not certified probability bounds. The code also checks the theoretical minimum-p/union bounds. FP64 tails below `1e-280` are explicitly outside the supported omnibus range.

- `--skat-o-relative-tolerance`: default `0.005`, allowed `1e-6` to `0.05`.
- `--skat-o-max-evaluations`: default `1048576`; a work limit, not a p-value resolution.
- `--skat-o-seed`: default `20260827`; controls randomized quadrature or optional simulation.
- `P_value_method` records the method, estimated error, evaluations, and seed.
- Failure to meet convergence/work/union checks stops the current part. No Bonferroni, moment-matching, minimum-p, or simulation fallback is substituted.

`--skat-o-method simulation` remains available with `--skat-o-simulations` (default 10,000), preserving correlated Gaussian calibration and the `(extreme+1)/(simulations+1)` resolution limit. Version-4 score-space simulation is deterministic for its seed, but its random draws/order differ from the former sample-space implementation; old seeded output is not promised byte-identical.

## Results, restart, and limits

The first 15 result columns retain their names and meanings. `Sample_mask` and `Null_model` are appended; the latter reports fitted variances, likelihood, iterations, convergence and singularity diagnostics for related tests. Audits now contain one row **per set and trait**, including effective N, exact-mask hash, and mask-specific no-call/monomorphic counts. Output order is input trait order, then definition set order. The SHA-256 mask identifies ordered aligned column positions; aligned sample IDs/order and source identity are signed separately.

Each version-4 checkpoint part atomically stores both result and audit records for one trait/set tile. Signatures cover ordered source/alignment/design/definition state, trait content, resolved relationship covariance bases, fitter/dependency version, numerical settings, thresholds and backend/device/runtime identity. Result and audit assembly reproduces the same bytes after interrupted-tile recovery. Changed masks, trait values, GRM state, or settings reject reuse. Source metadata changes before part commit or assembly are fatal. Existing raw caches remain compatible; **use a new set-result checkpoint directory** for pre-v4 analyses. Preserve historical results rather than overwriting their provenance.

The dense REML workspace and GRM PSD check have heap guards. Related set-tile recommendations reserve more heap for null fitting. A very large individual set can still exceed memory or make its eigenvalue decomposition expensive. The trait-outer schedule rereads genotype tiles for every transcript, including transcripts sharing a mask. Large, high-entropy isoform analyses have no new end-to-end performance guarantee. Broad spectra can still exceed the existing positive-gamma work bound, or fail omnibus quadrature convergence. Sparse REML, binary/GLMM outcomes, supplied cohort GLS with set tests, and FP32 set tests remain unsupported. Ordinary single-variant mixed-model exact deletion is still separate future work.

## Independent validation

Fixtures and replayable generators are under `test/resources/related-set-reference` and `scripts/generate-*-reference.R`. They require no R installation at test runtime.

- Base R `optim`, Cholesky and direct GLS matrix products independently fit pedigree/GRM nulls on 48 or 43 samples. Checks cover all score/covariance entries, burden/SKAT statistics, a repeated mask with a different response, subset mean filling, and observed ratio endpoints. Maximum scaled score/covariance discrepancy on this fixture: `7.16e-7` (gate `3e-6`).
- Rank-two omnibus tails use independent conditional-normal interval integration with R `stats::integrate` (QUADPACK), covering positive/negative correlation and adjusted p-values down to `2.22e-23`. Maximum relative discrepancy: `1.28e-8`.
- Equal-spectrum ranks 3, 8 and 20 use independent conditional chi-square integration through minimum component p `1e-12` (relative gate `2e-6`).
- General rank-three component quantiles use **CompQuadForm 1.4.4 Farebrother** and independent **cubature 2.1.4.1 hcubature** for the joint event, down to minimum p `1e-10`. Maximum observed relative discrepancy: `2.77e-4`, below the default `0.005` integration tolerance. These checks validate the implemented joint event; they do not claim numerical equivalence to every R SKAT approximation/version.
- Real CUDA/OpenCL and native CPU operations are compared with scalar FP64 on partial tiles of 1, 17 and 129 variants, both fixed and related models, including fitted projections, burden aggregation, scores, covariance and SKAT p-values. Hardware checks skip only when no eligible runtime/device exists.

These are bounded correctness fixtures, not a universal small-tail accuracy certification. A converged numerical integral does not establish biological model adequacy or finite-sample calibration for every rare-variant design.
