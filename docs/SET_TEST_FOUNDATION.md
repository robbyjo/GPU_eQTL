# Rare-variant set-test foundation

This document freezes the scientific contracts implemented by the initial deterministic FP64 CPU reference. These classes are a developer/reference API under `gov.nih.eqtl.settest`; they are not yet a supported command-line analysis mode.

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

## Deliberately remaining work

- Connect the contract to aligned CSV and VCF/BCF readers and to existing indexed region/set memberships.
- Define a stable user-facing set audit/result schema and CLI/INI options.
- Stream bounded variant/set and trait blocks and add signature-bound checkpoint/restart.
- Add production weighted/unweighted burden execution only after the end-to-end input/output path matches this scalar reference.
- Do not begin GPU burden, SKAT, or SKAT-O acceleration until those boundaries are deterministic and tested.
