# GPU eQTL 1.9.0.3

This patch fixes CPU OpenBLAS thread control and clarifies its diagnostics. Java 17 or newer is required.

## OpenBLAS thread control

Earlier releases loaded separate native OpenBLAS instances for matrix multiplication and thread control. On Linux, `-Deqtl.cpu.threads=N` could change the nolapack instance while the full-library GEMM instance retained its initial setting. In an isolated reproduction with `OPENBLAS_NUM_THREADS=1`, requesting two threads left GEMM at one despite the application reporting two.

Version 1.9.0.3 binds all FP64/FP32 ordinary and right-transposed GEMM operations and their thread control to the same nolapack instance. The backend description labels the count as **requested**: JavaCPP's helper echoes the request rather than querying the effective native count, and OpenBLAS may impose a smaller limit.

The bundled Linux x64 OpenBLAS 0.3.34 binary retains `MAX_THREADS=64`. This release does not enable a 96-thread BLAS team. Application `--threads` remains separate from the native BLAS pool. Exact-pattern CPU analysis still uses one exclusive context/outer worker and includes serial Java phases between BLAS calls. Budget concurrent JVMs together; higher CPU utilization does not by itself establish higher throughput.

## Download and compatibility

Download `gpu-eqtl-1.9.0.3-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.0.3 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.0.3).

The universal JAR includes OpenBLAS natives for Windows x64, Linux x64/ARM64, and macOS x64/ARM64. CUDA/OpenCL requirements and the separate oneMKL profiles are unchanged. Dependencies remain JLinAlg 0.3.6, JDistlib 0.10.2, OpenBLAS 0.3.34 and JavaCPP 1.5.14.

Cache/checkpoint formats, fixed-effect signature rules, identifiers, sample ordering, missingness policies, degrees of freedom, thresholds and result schemas are unchanged. Compatible completed fixed-effect blocks remain reusable with identical analysis settings and `--resume`; incomplete parts are retried. Fitted mixed-model checkpoints also bind the backend description, so the changed OpenBLAS description makes older OpenBLAS mixed-model checkpoints incompatible. Keep those original checkpoints and start a fresh mixed-model checkpoint directory instead of forcing reuse.

Using the correctly selected native thread count can change FP64 accumulation order, so newly computed results may differ at floating-point roundoff scale. Keep the original logs/checkpoints, and change the executable only when starting a new invocation. Publishing this release does not update or restart existing runs.

## Validation and limits

A fresh-JVM regression rejects the duplicate full/nolapack binding and checks a deterministic FP64 product. Production FP64/FP32 GEMM, residualization and compact-pattern operations remain compared with scalar references. Clean release verification runs the full Maven suite, checks packaged licenses and prints backend information; source identity and artifact SHA-256 are recorded in the release manifest.

The exact patched engine class was also tested in isolated Linux x64 OpenJDK 21 processes while three production jobs continued. A native JNI getter queried the actual OpenBLAS instance: the patch honored two-thread and one-thread settings with one library loaded, while the original two-library implementation reproduced the failure. Ordinary/right-transposed FP64/FP32 products and output padding matched scalar references within `1e-12`/`2e-5` tolerances. Linux synthetic tests passed 21 checks with one optional oneMKL-profile skip, including six exact-pattern analysis/checkpoint/resume cases explicitly using OpenBLAS. These tests ran at low priority, restricted to two logical CPUs, with private temporary/native-cache paths; the installed JAR and existing jobs were unchanged.

The release build host is Windows x64, JDK 25 targeting Java 17, NVIDIA GeForce RTX 2080 driver 617.42. Linux validation covers isolated numerical and synthetic analysis tests, not production-cohort throughput or 96-thread scaling. macOS/Linux ARM64 native loading and a separate Java 17 runtime remain unvalidated. No performance gain is claimed for the binding fix.
