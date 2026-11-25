# Changelog

For release summaries and upgrade instructions, see the website [Updates page](https://robbyjo.github.io/GPU_eQTL/updates.html).

## [1.9.3](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.3) — 2026-11-25

- Accelerate window SKAT eigendecomposition with native CPU LAPACK values-only mode, retaining eigenvectors for SKAT-O and the portable FP64 reference. Preserve full score covariance, PSD/rank checks and the selected pedigree spectral REML null.
- Hand expensive positive-gamma series to the existing same-law, error-checked tilted contour inversion after a preferred 2,048-term budget. Keep exact special cases, relative accuracy and explicit failure limits.
- Overlap bounded CPU p-value workers (`--threads`) with subsequent context-owned matrix products, preserving source-ordered atomic result/audit checkpoints. Use fresh set result/checkpoint paths after upgrading; compatible raw caches and null identities are unchanged. See [resource budgeting](docs/SET_TEST_PERFORMANCE.md).

## [1.9.2.8](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.8) — 2026-11-24

- Fix pooled REML stalls when a variance component reaches a bound. Project constraint scores before BFGS directions and secant updates, resetting curvature when the free set changes. Keep the model, variance bounds, FP64 precision, total iteration budget and strict `1e-6` projected-score convergence criterion. Numerical nonconvergence remains fatal.
- Use pooled result signature version 7 and shared-null key `pooled-null-v6-projected-bfgs`. Upgrading from v1.9.2.7 or earlier requires fresh pooled result/checkpoint paths and refitted nulls; retain old immutable state. Legacy/set identities and raw/rotated genotype-cache formats are unchanged.

## [1.9.2.7](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.7) — 2026-11-23

- Add opt-in `--mixed-trait-batch-size 2..256` for genome-wide pedigree-block score scans with complete predictors. Traits sharing an exact observed mask share one centered/rotated genotype read and batched FP64 score/norm products. Keep per-trait REML fits, variance bounds, strict convergence, explicit residual safeguards and chi-square(1) inference. Default size 1 keeps the previous runner and checkpoint identity; batched results require fresh checkpoints and bind the batch size. Immutable compatible nulls remain reusable.

## [1.9.2.6](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.6) — 2026-11-23

- Include both Linux x64 and Windows x64 JCuda/cuBLAS JNI bindings in every build, independent of the build host. Add hardware-independent payload tests and release-JAR checks so a Windows-built universal JAR cannot silently omit its Linux CUDA bindings. The compatible NVIDIA driver and system CUDA/cuBLAS/NVRTC runtime remain separate requirements; numerical code and checkpoint identities are unchanged.

## [1.9.2.5](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.5) — 2026-11-23

- Prevent pooled REML line-search stalls near variance boundaries by computing the likelihood change directly from relative covariance and fixed-effect information, instead of subtracting two nearly equal full objectives. Preserve the model, bounds, iteration budget, Armijo allowance and strict `1e-6` projected-score criterion.
- Pooled result signature version 6 and a new shared-null algorithm key require fresh pooled results/checkpoints and refitted nulls when upgrading from v1.9.2.4 or earlier. Preserve earlier immutable state; legacy/set identities and aligned raw-cache formats are unchanged.

## [1.9.2.4](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.4) — 2026-11-22

- Stabilize pooled association information at cohort variance boundaries. GPU and scalar products now use the orthonormal weighted design; close subtractions use a reprojected residual norm. This fixes invalid negative information and inaccurate estimates after successful null convergence, preserving the model, fitted null, score law and collinearity threshold.
- Pooled result signature version 5 requires fresh result/checkpoint paths from v1.9.2.3 and earlier. The v4 shared-null key and immutable format are unchanged: the stable score basis is reconstructed from the verified null. Legacy/set signatures and raw-cache formats are unchanged.

## [1.9.2.3](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.3) — 2026-11-22

- Stabilize pooled pedigree REML at cohort variance boundaries with weighted Householder QR. Accept score-converged restarts despite objective roundoff and stop zero-length line-search steps. Preserve the statistical model, component bounds, total iteration budget and `1e-6` score criterion.
- Pooled result signature version 4 and a new shared-null key require fresh pooled result/checkpoint paths when upgrading from v1.9.2.2 or earlier. Preserve earlier outputs; legacy/set signatures and aligned raw-cache formats are unchanged.

## [1.9.2.2](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.2) — 2026-11-22

- Read explicit BCF header `IDX` values for FILTER/INFO/FORMAT and contigs, including sparse and reordered dictionaries. Reject conflicting or undeclared indices instead of crashing or decoding another field. Native bcftools synthetic fixtures cover GT/DS, sample order, missing calls, non-PASS sites and indexed traversal.
- BCF-derived QC, raw/prepared caches and result checkpoints have a new reader identity. Preserve previous BCF results and use fresh result/checkpoint paths when upgrading; CSV and VCF identities are unchanged.

## [1.9.2.1](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2.1) — 2026-11-20

- Reuse immutable, checksummed pooled null fits across chromosome processes with `--mixed-null-cache-dir`; per-key process locks avoid duplicate fits. Native weighted-design and large-family products preserve FP64 scalar/dense-R references.
- Verify pooled result/audit parts with model-bound SHA-256 records and checkpoint ownership. Preserve interrupted parts outside cleanup directories and reject lost digest indexes. Proportional inbred-identity reductions retain summed bounds, and deterministic starts share the requested iteration budget.
- Pooled signature version 3 requires fresh result/checkpoint directories when upgrading from v1.9.2. Legacy/spectral signatures and aligned raw-cache formats remain unchanged. See the [release notes](docs/RELEASE_NOTES_1.9.2.1.md) and [pooled guide](docs/POOLED_PEDIGREE_SCORE.md).

## [1.9.2](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.2) — 2026-11-20

- Add opt-in pooled pedigree-block single-variant score analysis: one pedigree cohort, fitted residual variances by cohort, exact observed masks, FP64 score products and optional indexed cis windows.
- Review fixes preserve inference across trait units, reject projected covariance confounding, retain identity variance-sum bounds, reuse the participant pedigree graph, validate exact cis metadata and protect indexed cache reads. See the [release notes](docs/RELEASE_NOTES_1.9.2.md) and [pooled guide](docs/POOLED_PEDIGREE_SCORE.md).
- Pooled draft version-1 checkpoints require fresh result directories. Existing legacy/spectral signatures and aligned raw-cache formats remain unchanged.

## [1.9.1.6](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.6) — 2026-11-20

- Resolve extreme SKAT tail accuracy failures for spectra with many weak eigenvalues using a spectrum-aware integration-by-parts remainder bound. Retain every positive eigenvalue, the same weighted chi-square law and SKAT-O inversion, the `2e-11` relative target, work limits and explicit unresolved-case failures. Accepted version-6 checkpoints and aligned raw dosage caches remain compatible under the existing input/model/backend checks.
- Add independent high-precision synthetic tail references down to approximately `2.12e-306`, a closed-form convolution check, scale-invariance tests and extreme-tail quantile inversion. See the [release notes](docs/RELEASE_NOTES_1.9.1.6.md) and [numerical method](docs/SET_TEST_FOUNDATION.md).

## [1.9.1.5](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.5) — 2026-11-19

- Stabilize spectral REML profile comparisons and score accumulation across trait units; retain the strict convergence check and use a new spectral checkpoint identity.

## [1.9.1.4](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.4) — 2026-11-19

- Add explicit `--set-null-fitter spectral` for continuous pedigree-adjusted set tests. Reuse one bounded exact-mask error-contrast decomposition while fitting each response's original two-component REML objective and preserving genotype score covariance. Include portable/native FP64 references, exact identity handling, actual score diagnostics and signature-incompatible spectral result checkpoints. The default legacy fitter and SKAT tail laws are unchanged. See [spectral set-null documentation](docs/SPECTRAL_SET_NULL.md).

## [1.9.1.3](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.3) — 2026-11-17

- Make set-test model-input fingerprints portable across shared-filesystem mounts by signing file contents instead of node-local device/inode keys. Retain local mutation checks and backend/runtime binding. Use fresh set-result checkpoints for the version-6 signature; raw dosage caches remain compatible.

## [1.9.1.2](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.2) — 2026-11-17

- Retain pedigree-adjusted set-test traits whose observed pedigree covariance is exactly identity by fitting the identifiable total genetic-plus-residual variance. The combined component is explicit in null diagnostics, and its bounds preserve the original covariance family.
- Bind the new set-null handling to the fitted-model checkpoint signature. Use fresh related-set checkpoints; raw dosage caches remain compatible.

## [1.9.1.1](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.1) — 2026-11-16

- Resolve difficult SKAT eigenvalue spectra with a tilted characteristic-function integral when the checked positive-gamma series reaches its work or resolution limit. This evaluates the same distribution, retains explicit error/work checks, and supplies the same tail to SKAT-O quantile inversion.
- Retry recoverable REML null line-search failures in QR/scaled coordinates, preserving the model, variance bounds, convergence tolerances and total iteration budget. Refine likelihood plateaus only with a decreasing projected score and require the original score tolerance. Use the orthonormal design in the fitted projection to reduce cancellation.
- Add independent synthetic Davies/Imhof, closed-form tiny-tail, balanced REML and direct GLS regressions. Version-5 set and version-3 mixed result signatures require fresh checkpoint directories; aligned raw dosage caches remain compatible.
- Repair the renamed upstream JLinAlg artifact URL without changing the pinned library bytes or checksum.

## [1.9.1](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1) — 2026-11-12

- Added pedigree- and genetic-relationship-adjusted burden, SKAT, and SKAT-O tests, fitting a Gaussian REML null per trait and incorporating its covariance in both genotype scores and score variance.
- Added FP64 CUDA/OpenCL execution for burden aggregation, genotype projection, scores, and score-covariance products, with operation reporting and scalar CPU reference checks.
- Replaced default SKAT-O simulation with joint rho-grid numerical integration, supporting stringent significance thresholds with explicit convergence failures. Correlated-null simulation remains selectable.
- Added exact observed-sample deletion for missing transcript ratios, preserving valid zero/one measurements. Null models, relationship submatrices, genotype filling, and sample-dependent quantities follow each trait's exact sample mask.
- Added independent R/CompQuadForm references, small-tail and hardware comparisons, checkpoint recovery tests, and the [related-set guide](docs/RELATED_SET_TESTS.md).
- Set outputs retain their first 15 columns and append `Sample_mask` and `Null_model`; results are trait-major and audits are per set/trait. Version-4 set checkpoints require a fresh checkpoint directory when upgrading. Raw dosage caches remain compatible.

See the [release notes](docs/RELEASE_NOTES_1.9.1.md) for validation scope and limitations.

## Earlier releases

- [1.9.0.3](docs/RELEASE_NOTES_1.9.0.3.md): unified OpenBLAS GEMM/thread-control binding and clarified requested-thread diagnostics.
- [1.9.0.2](docs/RELEASE_NOTES_1.9.0.2.md): scalar verification of small exact-pattern correlation excursions.
- [1.9.0.1](docs/RELEASE_NOTES_1.9.0.1.md): see the archived release notes.
- [1.9.0](docs/RELEASE_NOTES_1.9.0.md): see the archived release notes.
