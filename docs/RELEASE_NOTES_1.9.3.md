# GPU eQTL 1.9.3

This release includes [PR #47](https://github.com/robbyjo/GPU_eQTL/pull/47), accelerating set-test numerical work while preserving FP64 models, full score covariance, masks, allele/weight declarations, sample counts, degrees of freedom and thresholds.

- CPU LAPACK computes window eigenvalues for SKAT and eigenvectors for SKAT-O; the portable reference and existing PSD/rank checks remain available.
- Expensive positive gamma series move to the existing checked, same-law contour inversion after a preferred 2,048 terms. The original 131,072-term series remains recovery if early inversion fails. Accuracy/work limits, exact special cases and explicit failures are preserved; no moment matching is introduced.
- Bounded CPU result workers overlap p-values with later matrix products. Native algebra and exclusive device contexts stay on the producer; compact immutable metadata and results are consumed in source order before atomic checkpoint commits. Review added cooperative cancellation to real simulation/quadrature paths so failure cleanup stops outstanding work.

**Set upgrade requirement:** use fresh set result/checkpoint paths from v1.9.2.8 or earlier, including burden, SKAT and SKAT-O. Numerical execution signatures change. Compatible raw dosage caches and fitted-null identities remain unchanged. Changing only `--threads` result workers permits compatible resumes. Pooled single-variant identities are unchanged from v1.9.2.8; earlier pooled and BCF migrations still apply when skipping versions.

Budget the producer's native BLAS team plus CPU result workers across all concurrent JVMs. For an eight-core allocation, an example is `-Deqtl.cpu.threads=4` before `-jar` with `--threads 4` among analysis options. Worker limits are admission bounds, not complete memory estimates. SKAT-O still needs quadratic workspace for portable rho-dependent eigenproblems. See [resource guidance and reproducible component benchmarks](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/SET_TEST_PERFORMANCE.md).

The reviewed local suite passed 283 tests: 278 passed and five optional skips, with zero failures/errors. Native/portable eigenspectra and SKAT-O scales, independent extreme-tail references, full-series recovery, real simulation cancellation, ordering/atomic failure and CPU/GPU resumes pass. Local execution used JDK25 targeting Java17, RTX2080 CUDA/OpenCL, native OpenBLAS and portable Java. The release verification record identifies checks run against the exact official JAR.

Component timings do not establish full-cohort speedup. No production deployment, controlled cohort analysis, fresh Linux/other-vendor GPU or concurrent multi-GPU throughput is claimed.
