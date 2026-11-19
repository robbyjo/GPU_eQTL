# GPU eQTL 1.9.1.5

This patch fixes spectral REML fits that could fail convergence because rounding excluded an accurately stationary profile candidate. It includes [PR #13](https://github.com/robbyjo/GPU_eQTL/pull/13).

Profile candidates are now compared after removing a trait-unit-dependent constant using the pilot residual variance. Compensated summation improves the profiled scale, objective and log-variance scores. The covariance family, variance bounds, original-scale reported likelihood, 65,536-evaluation work limit and strict 1e-6 projected-score threshold are unchanged. The fix does not alter genotype score covariance, FP64 defaults, dependencies or SKAT/SKAT-O tail laws.

Continue to select this fitter with `--set-null-fitter spectral` for continuous pedigree burden, SKAT or SKAT-O analyses. The default remains `legacy`. See the [spectral method and CLI guide](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/SPECTRAL_SET_NULL.md).

**Use a fresh spectral result/checkpoint directory when upgrading from v1.9.1.4.** The spectral fitter signature advances to `spectral-error-contrast-reml-v3-scaled-compensated-profile-resolved-engine`. Retain previous outputs and checkpoints for provenance. Legacy fitter identities and aligned raw dosage-cache formats are unchanged; all other input, setting and backend/runtime compatibility checks still apply.

The new synthetic 512-sample regression tests trait scales from 1e-100 to 1e100, checking convergence, variance scaling and the original-scale likelihood transformation. It reproduced a projected-score failure against the released v1.9.1.4 JAR and passed with the corrected source. The full reviewed suite passed 213 tests: 208 passed, five optional skips, zero failures or errors. Required hosted Linux/Java 17 CI also passed. Local validation used Windows with an NVIDIA GeForce RTX 2080, driver 621.14, and the Java/OpenBLAS/CUDA/OpenCL backends. Final exact-JAR and publication evidence is recorded in `docs/release-1.9.1.5/verification.json` after publication.

This release task used synthetic validation; it did not rerun private cohort analyses, remote diagnostics or production jobs. Dense-memory and explicit numerical convergence/work limits remain. Fresh multi-node, macOS, other-vendor and multi-GPU validation is not claimed.

Download the universal JAR with `SHA256SUMS` and `RELEASE-MANIFEST.txt`, and verify the checksum before use. Java 17 or newer is required. OpenBLAS natives are bundled for Windows x64, Linux x64/ARM64 and macOS x64/ARM64; CUDA/OpenCL require compatible installed runtimes.
