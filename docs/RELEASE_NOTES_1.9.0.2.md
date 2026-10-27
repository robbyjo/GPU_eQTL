# GPU eQTL 1.9.0.2

This patch addresses a genotype-outer exact-pattern scan abort at correlation `-1.000000061191563` and improves association progress estimates. Java 17 or newer is required.

## Correlation verification

The reported overshoot of approximately `6.12e-8` exceeds the existing FP64 `1e-8` allowance. Independently accumulated products and subtraction of covariate projection energy can lose accuracy when residual predictor variance is small.

For finite excursions greater than `1e-8` and at most `1e-6` beyond either end of `[-1, 1]`, genotype-outer exact-pattern analysis now verifies the association with scalar CPU QR. It uses the same observed samples, covariate design, and pattern-specific predictor replacements, calculates residual norms directly, and reprojects the prepared trait. The verified correlation and scales determine threshold inclusion, effect, t statistic, and p value. This preserves finite inference for a nearly perfect fit instead of automatically declaring a perfect fit.

The `1e-6` window is a verification limit, not an expanded clamping tolerance. Correlation and relative predictor/trait scales must agree within that limit. Disagreement, larger excursions, and non-finite statistics still stop the block. Errors identify predictor, trait, pattern, N, and DF. Successful blocks report their scalar verification count. Ordinary FP64 validation and opt-in FP32 behavior are unchanged.

## Progress

Genotype-outer scans show genotype row count times original expression row count as an approximate upper bound from the start, with approximate percentage and ETA. The discovered active comparison count remains separate, and completion reports the exact eligible count. Reused checkpoint comparisons remain excluded from current-run throughput.

## Resume and download

Download `gpu-eqtl-1.9.0.2-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.0.2 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.0.2).

Keep the failed run's log, caches, and checkpoints. Replace the JAR and use the same inputs, options, and checkpoint directory with `--resume`. Previously accepted correlations retain their existing calculation, so compatible completed 1.9.0.1 fixed-effect blocks remain reusable; incomplete blocks are retried. Cache/checkpoint formats, result schema/order, sample alignment, missingness policies, and degrees of freedom are unchanged. JLinAlg remains pinned to 0.3.6 and JDistlib to 0.10.2; mixed-model checkpoints from before the 1.9.0.1 dependency update remain incompatible.

## Validation and limits

The full Maven suite passes 166 tests: 161 passed, five optional skips, zero failures or errors. Nine targeted regressions also pass against the packaged JAR. Coverage includes both signs, finite near-perfect inference, independent scalar projection, exact masks, mean/zero replacement, thresholds, failure atomicity, and deterministic resume. CUDA/OpenCL and Java/OpenBLAS production operations remain compared with scalar references. Required license/notice files are retained.

Validation host: Windows x64, JDK 25 targeting Java 17, NVIDIA RTX 2080 driver 617.42, CUDA runtime 12.6 with Toolkit 12.9 libraries on the process PATH. The production dataset has not been rerun; scalar disagreement still requires investigation. Linux/macOS native loading, a separate Java 17 runtime, and production-scale behavior remain unvalidated. Release build identity, checksums, and verification details are recorded in the release manifest and SESSION_HISTORY.md.
