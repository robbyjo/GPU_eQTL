# GPU eQTL 1.9.1

This release extends rare-variant burden, SKAT, and SKAT-O analyses with relatedness adjustment, FP64 GPU matrix calculations, stringent-tail SKAT-O integration, and exact deletion of missing transcript measurements. Java 17 or newer is required.

## Relatedness and missing transcript measurements

Use `--pedigree` or `--relationship-matrix` to fit a Gaussian REML null for each trait. Both genotype scores `U = G'Py` and their covariance `C = G'PG` use the fitted covariance. Related burden uses chi-square(1) inference; ordinary fixed-effect burden retains its t statistic, effect, R-squared, and degrees-of-freedom contract.

With `--trait-missing pattern`, each transcript is tested on its observed samples. Ratios of zero and one remain valid. Relationship matrices are subset before fitting; covariate rank, degrees of freedom, genotype filling, and degeneracy decisions are recomputed for the mask. Distinct traits sharing a mask still receive distinct fitted nulls. Aligned-cohort frequency filters remain the default; pattern-specific frequencies require explicit selection.

## FP64 GPU execution

CUDA and OpenCL now execute substantial burden aggregation, genotype projection (including fitted `P*G`), score, and score-covariance products through the existing backend interfaces. Logs identify actual backend operation counts. QR/rank decisions, REML, eigendecomposition, and p-value calculations remain on the CPU. The scalar FP64 reference and `--residualization cpu` fixed-projection path remain available. Set tests require FP64.

## SKAT-O significance

The default `--skat-o-method integration` integrates the correlated union across the rho grid using the same component distribution for tail evaluation and inversion. Exact rank-one, adaptive rank-two, and conditional equal-spectrum cases are supplemented by independently shifted directional quadrature with analytic radial survival integration. It has no simulation-count p-value floor. Tolerance, evaluation budget, and seed are configurable; unresolved integrals fail explicitly. `--skat-o-method simulation` retains correlated-null simulation for comparison.

General directional-quadrature error estimates are estimates, not certified bounds. Difficult high-rank or broad-spectrum sets can exhaust the work budget, and tails below `1e-280` are unsupported. An unsuccessful convergence check must not be interpreted as a valid significant result.

## Compatibility

Legacy INI and command-line workflows remain supported. Set result columns 1–15 are preserved, with `Sample_mask` and `Null_model` appended; output order is now trait-major, and audit records are per set/trait. Version-4 set checkpoints bind sample masks, fitted-model state, integration settings, and backend identity. **Start a fresh set-test checkpoint directory when upgrading**; do not force reuse of older set result parts. Aligned raw dosage cache formats remain compatible.

Dependencies and native packaging are unchanged. The OpenBLAS thread-control fix from 1.9.0.3 is retained, including the bundled Linux x64 library's 64-thread build limit. This release does not replace installed JARs or alter running analyses.

## Validation and limits

- Independent base-R REML/direct GLS references cover pedigree/GRM models, full and partial masks, and distinct responses sharing a mask. Maximum scaled score/covariance discrepancy was `7.16e-7` against a `3e-6` tolerance.
- Independent R conditional-normal quadrature covers adjusted SKAT-O p-values down to `2.22e-23`, with worst relative discrepancy `1.28e-8`. Equal-spectrum references cover ranks 3, 8, and 20. General rank-three CompQuadForm 1.4.4 Farebrother/cubature references cover component minima through `1e-10`, with worst relative discrepancy `2.77e-4` under the default 0.005 tolerance.
- Real CUDA/OpenCL matrix operations are compared with scalar FP64 results, including partial 1/17/129-variant tiles. Tests also cover missing-sample policies, mismatched cache state, convergence failures, and deterministic checkpoint recovery. The implementation suite passed 177 tests with five optional skips (four oneMKL-profile checks and one external bcftools fixture).
- Release verification runs the full Maven suite from a clean commit, checks packaged notices, and runs backend discovery. The published manifest records the source commit and SHA-256. Hardware exercised: Windows x64, JDK 25 targeting Java 17, NVIDIA GeForce RTX 2080, driver 621.14, CUDA runtime 12.6, and OpenCL 621.14; OpenBLAS and portable Java CPU references also ran.

This is synthetic correctness validation, not production-cohort throughput or all-spectrum calibration. Dense REML/relationship matrices and eigendecomposition retain their memory limits; genotype blocks are reread per trait. Concurrent multi-GPU performance, other GPU vendors, Linux/macOS execution of this release, and a separate Java 17 runtime were not validated. Binary/GLMM outcomes, supplied cohort-GLS set models, ordinary single-variant mixed-model exact deletion, and FP32 sets remain unsupported.

## Download and documentation

Download `gpu-eqtl-1.9.1-universal.jar`, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` from the [v1.9.1 release](https://github.com/robbyjo/GPU_eQTL/releases/tag/v1.9.1). The universal JAR bundles OpenBLAS natives for Windows x64, Linux x64/ARM64, and macOS x64/ARM64; GPU drivers and CUDA/OpenCL runtimes remain external requirements. Bundled platforms do not imply that native loading was tested on every host.

See the [changelog](../CHANGELOG.md), [related-set guide](RELATED_SET_TESTS.md), [set definitions](SET_TEST_FOUNDATION.md), and [command-line reference](COMMAND_LINE_REFERENCE.md).
