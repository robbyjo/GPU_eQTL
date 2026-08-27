# Rare-variant set-test foundation

This document freezes the scientific contracts implemented by the deterministic FP64 references and the production burden, SKAT, and SKAT-O modes under `gov.nih.eqtl.settest`.

## Explicit variant-set input

`QVariantSetTable` reads a UTF-8 tab-separated file with a header. Required columns are:

| Column | Contract |
| --- | --- |
| `SET_ID` | Nonblank set/gene/region identifier. First appearance defines deterministic set order. |
| `VARIANT_ID` | Exact canonical identifier expected from the aligned variant source. Within-set order follows the file. |
| `REF` | Exact source reference allele, compared case-insensitively after uppercasing. |
| `ALT` | Exact single alternate allele. Multiallelic comma-separated values are rejected. |
| `EFFECT_ALLELE` | Must equal the declared REF or ALT. ALT uses the source dosage; REF uses `2 - ALT dosage`. |
| `WEIGHT` | Optional finite positive weight. Blank or omitted means `1`. |

One variant may occur in several sets; this is the supported overlapping-set representation. Repeating the same variant within one set is fatal. The normalized ordered content has a SHA-256 signature. Set IDs and variant IDs may not contain tabs, commas, or semicolons.

Allele harmonization is intentionally conservative. Declared REF/ALT must match the source REF/ALT exactly. The reference does not infer strand complements, swap REF/ALT, or guess from allele frequency. A mismatch is fatal before statistics are produced.

For indexed VCF/BCF region runs, omitting `--variant-sets` adapts the aligned variant-QC `region_sets` memberships. Each retained biallelic row becomes an exact REF/ALT definition with ALT as the effect allele and weight one. Declared region order, overlapping memberships, and empty declared regions are retained. An explicit TSV takes precedence when supplied.

## Frequency, missingness, and weights

MAF and MAC are calculated from finite additive ALT dosages in the final aligned analysis samples, before replacement. `QSetTestPolicy` has inclusive minimum and maximum MAF/MAC bounds. Its unfiltered factory uses MAF `[0, 0.5]` and MAC `[0, infinity]`; a rare-variant maximum is never silently assumed.

Supported predictor-missingness policies are:

- `error`: any missing dosage in an otherwise retained variant is fatal.
- `mean`: fill with the called-sample mean dosage of the declared effect allele.
- `zero`: fill with zero copies of the declared effect allele.

Frequency-filtered variants do not enter a burden. Absent definitions and empty or post-covariate monomorphic sets each have explicit `error` or `skip` behavior. Audits retain requested, absent, frequency-excluded, and included counts plus included variant IDs and final status.

## Reusable continuous-trait null model

`QSetTestNullModel` accepts an already aligned covariate design and continuous trait rows. Fixed-effect QR/rank decisions remain FP64 on the CPU. A rank-deficient design, non-positive residual degrees of freedom, missing/non-finite trait or covariate value, or zero-variance residual trait is fatal.

The null object prepares each trait once as `Y - (Y Q) Q^T`, records its residual standard deviation, and exposes the same state to every set. For one tested burden predictor, residual degrees of freedom are `N - rank(X) - 1`.

## Scalar burden reference

For each retained variant and aligned sample, `QBurdenReference` adds

`weight * declared-effect-allele dosage`

in deterministic file order. It projects the resulting burden through the same covariate Q, standardizes it in FP64, and uses the production `QeQTLStatistics` conversion for R-squared, per-weighted-dosage effect, signed t statistic, and log10 two-sided p-value. Unweighted burden is the exact special case in which every weight is one.

Output order is set first-appearance order followed by input trait order. Every result records set ID, trait ID, retained variant count, effective N, residual DF, R-squared, effect, t, and log10 p.

The production burden path prepares every resident set burden once and batches its FP64 set-by-trait products. Tests compare this path directly with the scalar reference.

## SKAT and SKAT-O

SKAT residualizes every declared-effect dosage row against the fixed-effect Q matrix without genotype standardization. Positive definition weights multiply kernel columns. The statistic is the weighted score quadratic form divided by the FP64 null residual variance; covariance eigenvalues come from the weighted projected-genotype Gram matrix. The null variance uses `N - rank(X)` degrees of freedom.

Quadratic-form survival probabilities use deterministic Imhof numerical inversion with fixed tolerance, block, and iteration limits. A single positive eigenvalue and equal positive eigenvalues use their exact scaled chi-square distributions. Failure of the bounded Imhof convergence test is reported as `satterthwaite-fallback` in the output; it is never silent.

SKAT-O defaults to rho `0,0.25,0.5,0.75,1`. Each component combines `(1-rho)` times the SKAT kernel with `rho` times the weighted burden kernel. The same seeded Gaussian vector drives all rho components in each parametric-null replicate, preserving their correlation. The adjusted p-value is `(extreme + 1)/(simulations + 1)` and the output includes the minimum component p-value only as an audit field; it is not reported as the omnibus result. Defaults are 10,000 simulations and seed 20260827, both configurable.

## Production scheduling and output

Use `--analysis burden|skat|skat-o` with `--variant-sets FILE`, or with indexed VCF/BCF region definitions. `--set-block-size` (default 256) bounds resident set definitions/dosages and `--expression-block-rows` bounds resident traits. The scheduler writes one atomic checkpoint part for each deterministic set-tile/trait-tile pair. Resume validates source identities, aligned sample orders, covariate design, definitions, masks, thresholds, tile sizes, method, and SKAT-O settings before reusing a part.

The stable CSV result columns are `Set_ID,Trait_ID,Method,Variants,N,DF,Statistic,R_squared,Effect,T,P_value,log10P,P_value_method,Minimum_component_P,Rho_component_P_values`. Inapplicable fields are `NaN` or blank. The companion TSV audit records requested, absent, frequency-excluded, and included variants, status, retained IDs, definition signature, and method.

## Current limits

- Set tests require FP64 and continuous traits. Exact trait-pattern deletion is not yet supported in these modes.
- Optimized CPU products are production-supported; direct multi-GPU set-test products remain a profiling-led extension and must match these references before becoming selectable.
- Relatedness-aware set tests remain deferred until the cohort covariance/null-model design is implemented.
