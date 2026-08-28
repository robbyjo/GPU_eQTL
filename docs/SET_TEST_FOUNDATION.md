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

## Automatic sliding windows

`--window-size BP` generates chromosome-local sliding-window sets directly from retained variant identifiers; no region TSV is needed. `--window-stride BP` controls the distance between starts and defaults to the window size. Both values must be positive and the stride cannot exceed the size. The grid is one-based and anchored at coordinate 1: window starts are `1 + k * stride`, endpoints are inclusive, and only windows containing at least one retained variant are emitted. Set IDs are `CHROM:START-END`, ordered by source-contig first appearance and then increasing start. Variants within a window retain source order.

Automatic windows work with CSV and VCF/BCF sources whose canonical row identifiers encode `CHROM:POS:REF:ALT`; additional colon-separated identifier fields are allowed. Membership uses POS. Every automatic member uses exact ALT effect orientation and weight one. Use an explicit `--variant-sets` TSV for custom weights or memberships, or indexed `--region`/`--regions-file` definitions for custom genomic intervals. Automatic windows are mutually exclusive with those custom-definition options.

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

Use `--analysis burden|skat|skat-o` with `--window-size BP`, `--variant-sets FILE`, or indexed VCF/BCF custom region definitions. Every expression-matrix row is a tested phenotype; selected fixed covariates are adjustment variables and are not candidate phenotypes. An omitted or zero `--set-block-size` selects a heap-aware resident set tile from the aligned sample count, method workspaces, and actual set-membership density; a positive value is an explicit override. `--expression-block-rows` bounds resident traits. The resolved set size is part of the checkpoint signature. The scheduler writes one atomic checkpoint part for each deterministic set-tile/trait-tile pair. Resume validates source identities, aligned sample orders, covariate design, definitions, masks, thresholds, tile sizes, method, and SKAT-O settings before reusing a part.

Automatic sliding-window runs build or reuse the version-1 checksummed aligned FP64 raw cache before set construction. A separately atomic and checksummed row-offset sidecar is bound to the raw-cache signature, byte length, modification time, and row count. Corrupt or incompatible sidecars are rebuilt; selected raw rows still undergo their version-1 CRC check. Resident set tiles retrieve only requested rows in source order, coalesce adjacent rows into sequential reads, bulk-decode each big-endian FP64 row, and transfer the trusted row buffer into the internal set-test value without a second dosage copy. Existing version-1 raw cache bytes remain compatible.

Full chr22 CPU burden profiling used 428,629 variants, 4,746 aligned samples, two expression traits, `Age`, FP64, OpenBLAS, an 8-GiB heap, and set/trait tile sizes 256/1. The completed warm-cache measurements were:

| Window / stride | Nonempty windows | Set blocks | Wall seconds | Selected FP64 bytes | Observed JVM heap at tile boundaries |
| --- | ---: | ---: | ---: | ---: | ---: |
| 10 kb / 5 kb | 7,453 | 30 | 182.14 | 16,345,186,032 | 3,861,495,624 |
| 50 kb / 10 kb | 3,808 | 15 | 229.10 | 16,451,762,208 | 6,297,342,968 |
| 100 kb / 50 kb | 765 | 3 | 148.35 | 16,315,419,120 | 7,905,198,712 |

The first CSV run additionally spent substantial time building the 16.29-GB aligned raw cache; preserve and reuse it for repeated grids. Heap depends on the number of unique variants spanned by one set tile, not just the number of windows. On the 100-kb/50-kb grid, automatic sizing selected 70 sets from a 3,956,932,400-byte worst-tile estimate under a 4,026,531,840-byte budget. It completed in 172.15 seconds with a 5,406,173,208-byte observed heap peak, versus 7,905,198,712 bytes with the former 256-set default; result rows were identical as an unordered set and the audit was byte-identical. Explicit sizing remains available for reproducibility and experiments. These measurements justify the indexed/bulk source path; they do not imply that one biological grid is preferable to another.

The stable CSV result columns are `Set_ID,Trait_ID,Method,Variants,N,DF,Statistic,R_squared,Effect,T,P_value,log10P,P_value_method,Minimum_component_P,Rho_component_P_values`. Inapplicable fields are `NaN` or blank. The companion TSV audit records requested, absent, frequency-excluded, and included variants, status, retained IDs, definition signature, and method.

## Current limits

- Set tests require FP64 and continuous traits. Exact trait-pattern deletion is not yet supported in these modes.
- Optimized CPU products are production-supported; direct multi-GPU set-test products remain a profiling-led extension and must match these references before becoming selectable.
- Relatedness-aware set tests remain deferred until the cohort covariance/null-model design is implemented.
