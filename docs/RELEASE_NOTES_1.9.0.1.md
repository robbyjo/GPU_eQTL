# GPU eQTL 1.9.0.1

This patch fixes an association scan abort near a perfect correlation. Java 17 or newer is required.

## Fix

The FP64 correlation validator now clamps finite excursions up to `1e-8` beyond either end of `[-1, 1]`. The previous `1e-10` allowance rejected the reported genotype-outer value `1.0000000001232048`. Independently accumulated cross-products and residual variances can disagree by this small amount through rounding and subtraction cancellation. Clamping happens before R-squared and t/p inference. In-range correlations are unchanged; larger excursions and NaN/infinity still stop the run. FP32 retains its existing `1e-4` allowance.

Regression tests cover the exact reported value and its negative, the FP64 boundary, unchanged valid correlations, rejection of larger/non-finite errors, FP32 compatibility, and perfect-correlation inference (`R-squared=1`, signed infinite t, `log10P=-Infinity`).

## Build and download

Download `gpu-eqtl-1.9.0.1-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.0.1 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.0.1). GPU runtime requirements are unchanged; see the README for CUDA setup.

The original JLinAlg v0.3.0 download is no longer available upstream. The build now extracts its unchanged `org.jlinalg` classes from the published GPU eQTL v1.9.0 universal archive, pinned by SHA-256 `5a457224db20a87737e797c9cc315d492b29930ab3496dcea811fc59f4b3d64b`. JDistlib remains pinned at 0.10.2. No numerical-library version is upgraded.

## Validation and limits

The full Maven suite passes: 159 tests, 154 passed, five optional skips, zero failures/errors. Validation ran on Windows x64, JDK 25 targeting Java 17, NVIDIA RTX 2080 / driver 617.42, with CUDA 12.9 libraries on the process PATH. CUDA/OpenCL hardware tests and portable Java/OpenBLAS comparisons passed. Packaged checks are recorded in SESSION_HISTORY.md and the release manifest. The supplied production error excerpt identifies the failing check; the restricted production analysis has not been rerun by an agent. Retain its logs and incomplete checkpoint state for the audit trail, and use the patched JAR for the user-executed continuation. This patch does not change input/cache/checkpoint schemas or signatures. Linux/macOS native loading and production-scale behavior require validation in those environments.
