# GPU eQTL 1.9.2.6

This release includes [PR #38](https://github.com/robbyjo/GPU_eQTL/pull/38), fixing Linux CUDA JNI libraries missing from Windows-built JARs. Every build now includes the six JCuda/cuBLAS JNI bindings for each of Linux x64 and Windows x64. Hardware-independent tests check the payloads and binary signatures; both release scripts check all twelve entries in the final JAR.

Install the compatible NVIDIA driver and system CUDA/cuBLAS runtime separately. Genotype-outer exact-pattern CUDA analysis also needs NVRTC. Bundled JNI bindings do not replace these runtimes. See the [installation guide](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/USER_GUIDE.md#download-and-first-run).

Numerical code, statistical models, precision defaults and cache/checkpoint identities are unchanged from v1.9.2.5. Existing exact compatibility checks still apply. Upgrading from v1.9.2.4 or earlier still requires fresh pooled results/checkpoints and refitted nulls, as described in the [v1.9.2.5 release notes](https://github.com/robbyjo/GPU_eQTL/blob/master/docs/RELEASE_NOTES_1.9.2.5.md). Preserve earlier state.

The reviewed local suite passed 261 tests: 256 passed and five optional skips, with zero failures/errors, including actual RTX 2080 CUDA/OpenCL operations, native OpenBLAS and portable Java references. The PR author also validated an earlier Windows-built candidate on Linux/Java 17 with an A100 and CUDA 12.6; that evidence does not describe the final release artifact. The release verification record identifies checks run against the official JAR.

This packaging repair does not change pedigree scan performance. No production deployment, controlled cohort analysis, full-cohort throughput, or fresh other-vendor/multi-GPU validation is claimed.
