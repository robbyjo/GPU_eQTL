# REML, pedigree, JLinAlg and SKAT follow-up

This implements the follow-up to the [initial audit](AUDIT_2026-09-09.md),
against baseline `218b7b17d718c9c8bff3e3611d05d1e1c277775b`.
[MIXED_MODELS.md](MIXED_MODELS.md) documents the current interface and limits.

## Delivered

- Selective JLinAlg 0.3.0 adoption for Gaussian REML and pedigree construction,
  with modern JDistlib 0.10.2 for new inference. Existing fixed-effect QR and
  the GPU association engine stay in place.
- `--mixed-model score`/`genesis`: one REML null per trait, full projected
  score variance, chi-square(1) inference. This implements GENESIS-style null
  reuse, not its optional approximate diagonal score-SE method.
- `--mixed-model reml`: refit each predictor, with residual, Satterthwaite or
  Kenward–Roger Wald inference. Both modes support combined pedigree and
  random-intercept grouping components.
- Strict pedigree/sample checks, convergence/variance diagnostics, bounded
  predictor tiles, backend-neutral FP64 projection and atomic trait checkpoints.
- SKAT tails now use exact special cases, rank-two angular quadrature and a
  cached positive gamma series with a remainder/work limit. Unresolved spectra
  fail explicitly. SKAT-O inverts the same tail law. Set checkpoints advance
  to signature version 3.

## Correctness

Final default suite: **155 tests, 150 passed, five optional skips, zero errors
or failures**. Four skips are optional oneMKL checks; one requires external
bcftools. The separate oneMKL profile passes all **33 selected tests**.

Six independent R references cover grouping, pedigree and combined covariance,
under null-score and full-refit inference. Base R Cholesky/REML optimization
generates expected values; the grouping refit also agrees with `nlme`. Maximum
checked beta/SE/statistic/DF/log-p/likelihood difference: **3.48e-8**.
Tests also cover response units, predictor scale/offset, convergence failures,
pedigree errors/order, degeneracy, reordered samples and checkpoint integrity.

The 30-case independent CompQuadForm Davies grid has **zero fallbacks** and
maximum absolute p error **1.22e-8**, within R's requested 1e-7 accuracy.
The old fallback's maximum error was 0.0526. Rank-two tails down to about
**2.54e-219** agree with independent 80-digit mpmath quadrature. Monotonicity,
inverse-tail consistency, scale invariance and unresolved-spectrum failure
are covered.

The observed discrepancy was local to GPU eQTL's approximation. Historical R
SKAT versions also had documented tail/fallback issues, discussed by
[Wu, Guan and Pankow (2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761292/).
These results do not establish an R SKAT bug as the cause here.

Actual CUDA, NVIDIA OpenCL and CPU projections match scalar references,
including 1/17/33-column partial tiles. Hardware: i9-14900K, 32 logical
processors, RTX 4090, driver 610.88/API 13.3; JCuda 12.6 with staged CUDA 12
libraries (cuBLAS 12.9.2.10, runtime/NVRTC 12.9.79/12.9.86). The host has two
4090s; this is not a two-device scaling study.

## Timings

Same deterministic CPU harness; medians in milliseconds. Production fixture:
16 variants, eight traits, 512 samples; SKAT-O uses 99 simulations and rho
0, 0.5, 1. These are microbenchmarks, not production throughput.

| Operation | Original baseline | Initial audit | Current |
| --- | ---: | ---: | ---: |
| Fixed preparation | 2.1487 | 0.2815 | 0.2978 |
| Burden production | 2.6944 | 0.8560 | 0.8222 |
| SKAT production | 21.4314 | 2.7902 | 1.0091 |
| SKAT-O production | 1186.8624 | 68.0480 | 53.4428 |

On the 48-sample, 32-predictor combined-covariance fixture, including the null
fit, score mode takes **5.2844 ms** versus **100.8443 ms** for per-predictor
REML. These modes use different inference procedures; this is not a benchmark
against GENESIS or evidence of large-cohort scaling.

## Reproduction and remaining limits

See [verification.json](mixed-model-2026-09-09/verification.json) for commands
and artifact hash. The standard clean package build and cold dependency fetch
succeed. Deliberate dependency corruption fails SHA-256 verification before
extraction; the verified artifact was restored. Packaged score/REML commands
reproduce the combined-covariance reference. GPL, oneMKL and upstream
JLinAlg/JDistlib/posterior notices are present in the jar. Upstream bytecode is
compatible with the Java 17 compilation target.

Limits: dense REML state/fitting; sequential mixed-score tile submission;
no fitted mixed-model set tests or exact trait-pattern deletion; no measured
private-cohort/multi-GPU scaling. Some spectra exceed the positive-series work
limit and fail. The truncation bound does not certify all rounding error.
Satterthwaite/Kenward–Roger option propagation and finite inference are tested;
the six independent R association references use residual-DF Wald inference.
Next step: profile a representative large pedigree/cohort before adopting
sparse or native factorizations.
