# GPU eQTL 1.9.1.2

This patch release incorporates [PR #2](https://github.com/robbyjo/GPU_eQTL/pull/2), retaining pedigree-adjusted set-test traits whose observed-sample pedigree covariance is exactly identity. Java 17 or newer is required; FP64 remains the scientific default.

## Exact-identity pedigree nulls

After missing transcript measurements are removed, a pedigree mask can contain only unrelated, noninbred individuals. Its additive-genetic and residual covariance matrices are then both identity. Their separate variances cannot be identified, but their sum can be estimated.

The set-test null now fits one `total-genetic-plus-residual` component for this exact two-component case. Summing the original component bounds preserves the same covariance family. The trait and its observed samples are retained; results and audits explicitly report `identity-pedigree=combined-total-variance`.

Recognition requires correctly sized covariance matrices and exact identity entries. Near-identity matrices, additional covariance terms and ordinary single-predictor fits retain their existing behavior. Convergence tolerances and iteration limits remain enforced; genuine fit failures still stop execution.

## Upgrade compatibility

**Use fresh related-set result checkpoint directories.** The fitted-model signature changes to identify this rule. Aligned raw dosage caches remain compatible. The set output schema, sample alignment, missingness, allele and threshold policies are unchanged. Dependency versions and packaged native libraries are unchanged.

## Verification and limits

- The exact PR source passed the full Maven wrapper suite on Windows x64: **195 tests, 190 passed, five optional skips, no failures or errors**. Skips cover four separate oneMKL-profile checks and the external bcftools fixture.
- New synthetic tests compare total variance, coefficients, likelihood, projection and scores against closed-form REML/OLS/direct GLS. They also cover observed principal submatrices, invalid dimensions, near-identity rejection, ordinary-fit behavior, convergence failures, trait retention and byte-identical checkpoint recovery.
- Real CUDA/OpenCL checks ran on NVIDIA GeForce RTX 2080, driver 621.14, alongside OpenBLAS and portable Java references. The build uses JDK 25 with Java 17 as its class target. Clean release build and exact packaged-JAR results are recorded in the published manifest and [verification receipt](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/release-1.9.1.2/verification.json).

This fix does not perform general covariance-rank reduction or approximate merging of covariance bases. Dense model/eigendecomposition limits and explicit numerical work/convergence checks remain. This release review uses synthetic inputs; the PR separately records prior bounded Linux/A100 comparisons. No new private-cohort rerun, fresh Linux/macOS runtime validation, other-vendor check or concurrent multi-GPU scaling study is claimed here.

## Download

Download `gpu-eqtl-1.9.1.2-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.1.2 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1.2). The universal JAR bundles OpenBLAS for Windows x64, Linux x64/ARM64 and macOS x64/ARM64. GPU drivers and CUDA/OpenCL runtimes remain external requirements.

See the [changelog](../CHANGELOG.md) and [related-set guide](RELATED_SET_TESTS.md) for further details.
