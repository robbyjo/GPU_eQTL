# Changelog

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
