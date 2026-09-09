# GPU eQTL 1.9.0

This release publishes the current modernization and numerical fixes as **1.9.0**. Version 2.0 is reserved until broader testing is complete. Java 17 or newer is required; FP64 remains the default and FP32 remains explicitly approximate.

## Download and CUDA setup

Download `gpu-eqtl-1.9.0-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.0 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.0). The JAR includes Java dependencies and portable OpenBLAS natives. NVIDIA drivers and CUDA libraries are installed separately.

For CUDA, install the runtime libraries from [NVIDIA CUDA Toolkit 12.6](https://developer.nvidia.com/cuda-12-6-0-download-archive), including cuBLAS/cuBLASLt and NVRTC/builtins. Follow the [README installation instructions](https://github.com/robbyjo/GPU_eQTL/blob/v1.9.0/README.md#install-the-cuda-126-runtime-libraries) for driver requirements and Windows/Linux library paths. A newer NVIDIA driver is compatible; Toolkit 13 libraries alone do not replace the CUDA 12 dependencies.

```sh
java -jar gpu-eqtl-1.9.0-universal.jar --printbackendinfo
```

## Changes

- GENESIS-style mixed-model scores reuse a fitted REML null per trait. This mode is the default when pedigree or random effects are supplied; `--mixed-model reml` retains full per-association REML refitting with selectable Wald degrees of freedom.
- Selective JLinAlg 0.3.0 and JDistlib 0.10.2 adoption supports fitted Gaussian mixed models and pedigree construction while retaining the existing GPU association engine.
- SKAT tails use exact special cases and error-checked mixture evaluation. Unresolved spectra stop explicitly rather than silently substituting the old moment approximation. This fixes a local approximation issue; it does not establish an R SKAT bug.
- Corrected GLS whitening/residual handling and optimized compound-symmetry preparation, rare-variant kernel reuse, and small CPU matrix operations.
- Fixed-effect scans retain streamed/cache-backed input, exact trait-pattern deletion, deterministic checkpoint/resume, explicit alignment/QC, CUDA/OpenCL acceleration, and a CPU fallback. Burden, SKAT, and SKAT-O retain FP64 set-test contracts.
- Build metadata, application banner, release scripts, and current command examples now identify version 1.9.0.

## Validation and limits

The release is verified on Windows x64 with two RTX 4090 GPUs and NVIDIA driver 610.88, using CUDA 12.6 runtime 12.6.37, cuBLAS 12.6.0.22, and NVRTC 12.6.20. The universal artifact includes OpenBLAS natives for other desktop platforms, but Linux/macOS native loading is not validated by this Windows release run. Java 25 builds Java 17-targeted bytecode; a separate Java 17 runtime smoke remains unperformed on this host.

The full Maven suite and packaged backend/license checks are the publication gates. Independent R and multiprecision comparisons and measured synthetic benchmarks are documented in the [mixed-model/SKAT verification report](https://github.com/robbyjo/GPU_eQTL/blob/v1.9.0/docs/REML_SKAT_FOLLOWUP_2026-09-09.md). Synthetic benchmarks do not establish full-cohort or multi-GPU speedups.

Fitted mixed models currently use dense CPU REML state. Related-sample mixed-model set tests, exact trait-pattern REML, GLMMs, random slopes, categorical-SNP repair, interactions, and forward selection remain unsupported. Very wide SKAT spectra may exceed the explicit convergence limit. Existing incompatible cohort/set checkpoints must be rebuilt; current set checkpoints use signature v3. Broader platform and representative large-cohort testing remains necessary before version 2.0.
