# GPU eQTL 1.9.1.1

This patch release incorporates [PR #1](https://github.com/robbyjo/GPU_eQTL/pull/1), correcting failures on difficult SKAT eigenvalue spectra and stalled REML null fits. Java 17 or newer is required; FP64 remains the scientific default.

## Numerical fixes

- When the checked positive-gamma series cannot resolve a spectrum within its work limit, tilted characteristic-function inversion evaluates the same weighted chi-square distribution. It directly integrates small tails, checks truncation and estimated numerical error, and retains an explicit evaluation limit. SKAT-O quantiles invert this same distribution. This is not a moment-matched or saddlepoint approximation; unresolved calculations still fail.
- Recoverable REML null line-search failures retry in orthonormal fixed-effect coordinates and residual-SD units. The original model, variance bounds, convergence tolerance and total iteration budget are preserved. A likelihood plateau permits a refinement step only when the projected score decreases, and success still requires the original score tolerance. The fitted projection uses the orthonormal design to reduce cancellation. Ordinary per-predictor refits retain their existing convergence behavior.
- Recovered set nulls identify `fit-method=qr-scaled-retry`. The pinned library's coefficient-space result adapter represents unavailable original-coordinate score diagnostics as zero; those adapter values are not an additional convergence check. Actual convergence is checked before adaptation.
- The JLinAlg download URL follows the upstream rename to OmicsMatic. Pinned library bytes, checksum and dependency versions are unchanged.

## Upgrade compatibility

**Use fresh set-test and mixed-model result checkpoint directories.** Set execution signatures advance to version 5 and ordinary mixed-model signatures to version 3. Retain older results and checkpoints for provenance. Aligned raw dosage caches remain compatible. Output columns, sample alignment rules, missingness policies, and opt-in FP32 behavior are unchanged.

## Verification

- The merged source passed the full Maven wrapper suite on Windows x64: **190 tests, 185 passed, five optional skips, no failures or errors**. Skips cover four separate oneMKL-profile checks and the external bcftools fixture.
- Real CUDA and OpenCL calculations ran on an NVIDIA GeForce RTX 2080, driver 621.14, alongside OpenBLAS and portable Java references. The build uses JDK 25 and targets Java 17. The clean release build and exact packaged-JAR checks are recorded in the published manifest and [verification receipt](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/release-1.9.1.1/verification.json).
- Synthetic regressions include independent Davies/Imhof references, analytic upper tails near `7e-218`, direct lower-tail checks near `2.66e-37`, scale invariance and same-distribution quantile inversion. REML checks cover balanced closed-form estimates/direct GLS, response offsets and small variances, covariate scaling, active bounds, and exhausted iteration budgets. Related-set tests cover CPU/GPU agreement, exact sample masks, and deterministic checkpoint recovery.

General quadrature and floating-point error estimates are not certified interval bounds. Difficult spectra and fits can still reach explicit work or convergence limits; dense REML/relationship/eigen memory limits remain. This release check does not establish full-cohort performance, universal tail calibration, concurrent multi-GPU scaling, other-vendor behavior, or fresh Linux/macOS runtime validation. The merged PR separately records bounded Linux/H200 failure replays performed before merge.

## Download

Download `gpu-eqtl-1.9.1.1-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.1.1 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.1). The universal JAR bundles OpenBLAS for Windows x64, Linux x64/ARM64 and macOS x64/ARM64; GPU drivers and CUDA/OpenCL runtimes remain external requirements.

See the [changelog](../CHANGELOG.md), [related-set guide](RELATED_SET_TESTS.md), [mixed-model guide](MIXED_MODELS.md), and [set-test numerical details](SET_TEST_FOUNDATION.md).
