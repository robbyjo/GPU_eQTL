# GPU eQTL 1.9.1.4

This release adds an opt-in spectral REML null fitter for continuous pedigree-adjusted burden, SKAT and SKAT-O tests, from [PR #10](https://github.com/robbyjo/GPU_eQTL/pull/10). Select it with:

```text
--analysis skat --pedigree pedigree.tsv --mixed-model score --set-null-fitter spectral
```

Use `burden` or `skat-o` for the other supported set methods. The equivalent INI setting is `set_null_fitter = spectral`; command-line arguments override it. The default remains `legacy`. See the [complete command and method documentation](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/SPECTRAL_SET_NULL.md).

The spectral implementation reuses one exact ordered-mask decomposition and fits each trait's variances separately. It preserves the two-component pedigree-plus-residual covariance family, variance bounds, full genotype score covariance, FP64 arithmetic and SKAT/SKAT-O tail laws. Geometry and REML run on the CPU; the selected CPU, CUDA or OpenCL backend performs the set matrix products. Positive eigenvalues are retained without low-rank truncation. Unresolved fits fail explicitly.

Spectral mode requires a pedigree and score mode, with no additional random kernels or GRM. Ordinary single-predictor models, per-predictor REML refits and FP32 are unsupported for this selector. Output reports `fit-method=spectral-reml` and the actual maximum projected score. The spectral search uses its own documented evaluation limit; `--reml-max-iterations` continues to govern legacy fitting.

Use separate result/checkpoint directories for spectral and legacy analyses. Spectral checkpoints bind the resolved native OpenBLAS or portable Java geometry engine, and a mid-analysis engine change is fatal. Pre-release development spectral checkpoints are also incompatible with the released fitter identity. The existing legacy signature and aligned raw dosage-cache formats are unchanged; other input, setting and backend/runtime compatibility checks still apply.

Reviewed source validation passed 212 tests: 207 passed, five optional skips, no failures or errors, plus required hosted Linux/Java 17 CI. Tests cover dense REML/GLS references, native/portable geometry, exact identity, ordered masks, per-response fits, CLI/INI selection, invalid combinations, checkpoint identity and resume. Local hardware validation used a Windows NVIDIA GeForce RTX 2080 with driver 621.14. Final exact-JAR, cross-backend CLI and public-asset evidence is recorded in `docs/release-1.9.1.4/verification.json` after publication.

This is synthetic correctness validation, not a full-cohort throughput claim. Dense memory and numerical work limits remain. Private deployment pilots and production analyses were not rerun by this release task; fresh multi-node, other-vendor, macOS and multi-GPU validation is not claimed.

Download the universal JAR with `SHA256SUMS` and `RELEASE-MANIFEST.txt`, and verify the checksum before use. Java 17 or newer is required. OpenBLAS natives are bundled for Windows x64, Linux x64/ARM64 and macOS x64/ARM64; CUDA/OpenCL require compatible installed runtimes.
