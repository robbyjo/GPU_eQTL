# GPU eQTL 1.9.2.7

This release includes [PR #41](https://github.com/robbyjo/GPU_eQTL/pull/41), adding opt-in trait batching for genome-wide pooled pedigree scores. Traits with the same exact observed sample mask share centered/rotated predictor reads and wider FP64 score and weighted-norm products.

Add `--mixed-trait-batch-size 32` to a genome-wide `--mixed-null-fitter pedigree-block` score run with `--predictor-missing error`. Supported batch sizes are 1 through 256. The default remains 1 and uses the existing runner. Batching cannot be combined with cis-only trait positions, predictor imputation or FP32. See the [method, restrictions and storage costs](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/POOLED_PEDIGREE_SCORE.md#optional-genome-wide-trait-batching).

Each trait retains its own REML null, variance bounds, strict convergence criterion and chi-square(1) inference. Nearly collinear and covariance-boundary cases use the original explicit residual calculation. Result identifiers, ordering, counts, degrees of freedom and output thresholds retain their checks.

**When enabling batching or changing its size, use fresh result/checkpoint paths.** Batched identities bind the configured size. Compatible immutable nulls remain reusable. Runs that omit the new option keep the v1.9.2.6 default checkpoint identity; existing input/model/backend compatibility checks still apply. Earlier pooled and BCF upgrade requirements remain unchanged. Preserve earlier outputs and state.

The rotation cache stores another FP64 predictor matrix for each exact mask. Budget disk space accordingly; many masks can make batching unsuitable. Null fitting, parsing and cache bandwidth can still dominate. Batches run sequentially through an exclusive context and do not add concurrent multi-GPU scheduling.

The reviewed local suite passed 270 tests: 265 passed and five optional skips, with zero failures/errors, including actual RTX 2080 CUDA/OpenCL, native OpenBLAS and portable Java references. Tests include maximum-size batch products, independent within-cohort GLS boundary references, masks, partial tiles, thresholds, atomic resume/recovery and corruption rejection. A review correction bounds the combined checkpoint buffers and includes them in heap admission; maximum-batch checkpoint recovery also passes in a 128 MiB JVM. The release verification record distinguishes exact-JAR checks from the PR author's earlier Linux CPU/A100 prototype timings. No full-cohort throughput or production deployment is claimed.
