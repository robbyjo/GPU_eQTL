# Spectral pedigree nulls for continuous set tests

`--set-null-fitter spectral` selects an FP64 spectral REML null for continuous burden, SKAT and SKAT-O analyses with one pedigree covariance and an independent residual component. The default remains `legacy`. Other random-effect combinations and ordinary single-predictor mixed analyses retain their existing fitter.

```powershell
java -Xmx12g --enable-native-access=ALL-UNNAMED -jar target/gpu-eqtl-1.9.2-all.jar `
  --genotype genotype.csv --expression traits.csv --covariates covariates.csv `
  --fixed-covariates Age,Sex --genotype-id-column ID --expression-id-column ID `
  --analysis skat --variant-sets sets.tsv --pedigree pedigree.tsv --mixed-model score `
  --set-null-fitter spectral --trait-missing pattern --output spectral-skat.csv `
  --checkpoint-dir spectral-skat.checkpoint --keep-checkpoints
```

Choose `--analysis burden` or `--analysis skat-o` for the other supported set tests. Repeat the same command with `--resume` only after a compatible spectral checkpoint exists. Omit `--set-null-fitter` or specify `--set-null-fitter legacy` to use the existing fitter. Do not reuse a legacy result/checkpoint directory for spectral output.

For a positional INI or `--config analysis.ini`, add `set_null_fitter = spectral` alongside `analysis = skat`, `pedigree_file = pedigree.tsv` and `mixed_model = score`. Command-line `--set-null-fitter legacy` or `spectral` overrides that INI setting. Unknown values, ordinary single-predictor analyses, GRM-only models, additional random kernels, per-predictor REML mode and FP32 spectral requests are rejected before analysis.

For every exact ordered observed-sample mask, the implementation builds an orthonormal error-contrast basis B for the selected fixed effects, diagonalizes B'KB, and profiles the original restricted likelihood over the genetic/residual variance ratio. The total scale is optimized analytically subject to the original individual-component variance bounds. All positive covariance eigenvalues are retained. This is the same two-component covariance family, without low-rank truncation, shared heritability, or expression-only residualization.

One geometry is cached per model instance and exact ordered mask. Each response receives its own fitted variances. Replacing the cached mask releases the previous geometry; there is no unbounded cross-trait cache. A mask whose pedigree covariance is exactly identity retains the existing one-total-variance reduction and its summed original bounds. General dependent covariance bases are not silently merged.

The score response is P*y. Genotypes are transformed by a factor F satisfying F*F'=P, so their score covariance is (G*F)*(G*F)'. The existing CPU/CUDA/OpenCL product interface executes these products; the SKAT and SKAT-O tail laws remain unchanged. LAPACK geometry uses bundled OpenBLAS with the configured CPU thread count and a process-wide lock. `eqtl.cpu.blas=java` selects the portable reference. The full OpenBLAS thread setting is restored after each operation scope; the separate association BLAS engine's settings are preserved.

Profile candidates are compared after removing the response-unit constant using the pilot residual variance. Compensated accumulation is used for the profiled scale, objective and log-variance scores, preserving stationary-root selection when traits have very different units. The final reported restricted likelihood remains on its original scale, and the 1e-6 projected-score check and original variance bounds are retained. This numerical-profile version has a distinct spectral checkpoint identity; use a fresh spectral result/checkpoint directory after upgrading from v1.9.1.4.

The profile search has an explicit 65,536-evaluation work bound and requires a maximum projected log-variance score no larger than 1e-6. These are spectral-search evaluations, not iterations of the legacy Fisher-scoring optimizer; `--reml-max-iterations` continues to govern the legacy fitter. Invalid covariance, failed decomposition, unsupported models, nonfinite values and unresolved fits remain fatal. Rank, N, residual DF and the existing `--unestimable-trait-patterns` policy are unchanged.

Outputs identify `fit-method=spectral-reml` and report the actual maximum projected score. The JLinAlg coefficient-space adapter carries genuine original-coordinate coefficients, covariance and residuals, but represents unavailable adapter log-variance scores as zero; the independently checked spectral score is therefore reported separately. The spectral model signature includes the fitter/factor algorithm and the resolved native OpenBLAS or portable Java geometry engine, including when `eqtl.cpu.blas=auto` is requested. A change of geometry engine during fitting is fatal. Legacy and spectral result checkpoints are incompatible and must use separate directories. Compatible aligned raw dosage caches may be reused. A resumed spectral run preserves byte-identical committed output parts.

Changing the numerical optimizer does not mean previous converged fits were invalid. If a study requires all output from one implementation, run every trait in a fresh result/checkpoint namespace and keep the previous outputs for comparison. Benchmark the complete pipeline before interpreting fitter speedups as cohort completion-time improvements.

The tests compare native and portable geometry, independent dense REML/GLS results, fixed-effect coefficients/covariance, identity reduction, exact mask ordering, per-response fits, weighted SKAT windows, numerical output, and signed checkpoint/resume incompatibility. Tests and committed fixtures are synthetic. Private cohort inputs and row-level results are excluded from the repository.
