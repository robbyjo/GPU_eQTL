# Session history

Append-only record of material modernization work. Times use ISO-8601 with the local timezone.

## 2026-08-20T10:54:03.9481519-04:00 — JavaCL replacement and build foundation

### Baseline and review

- Baseline commit: `7cf681e1f4e76b74d773075911a210fcf4824420` (`7cf681e1`, 2026-08-10, merge by robbyjo).
- Reviewed the legacy Eclipse source layout, bundled jars, application startup, INI configuration, matrix loaders, multi-GPU scheduling, OpenCL kernels, real-valued and categorical SNP jobs, and output flow.
- The repository had no Maven/Gradle build, no test suite, and no user/contributor documentation. JavaCL RC1/RC3, BridJ, Commons Compress 1.3, and JDistlib 0.3.6b were supplied as local jars.
- Confirmed the major follow-on gaps: INI-only invocation, whole-matrix RAM loading, no sample-ID alignment check, no automatic categorical-covariate encoding, disabled categorical-SNP execution, no SNP interactions, and no forward selection.

### Decisions

- Replaced JavaCL with JOCL 2.0.6 behind a small vendor-neutral GPU API (`GpuBackend`, `GpuDevice`, and `GpuContext`). JOCL runs the existing OpenCL C FP64 kernel through the system OpenCL ICD, giving one initial path for NVIDIA, AMD, and Intel drivers while leaving an explicit seam for future CUDA, HIP/ROCm, Vulkan, or Level Zero backends.
- Kept double precision mandatory. Device discovery now requires an available device, a compiler, and FP64 support, and reports platform, vendor, driver, API, memory, and work-group metadata.
- Made kernel compilation and worker execution fail fast. Worker exceptions now propagate to a distinct nonzero analysis exit code. Contexts and buffers have explicit lifecycles, context reservation is protected by `try/finally`, and device buffers are reused across matrix blocks.
- Reduced wasted execution on the last partial matrix block by launching only the rounded active SNP/expression dimensions. Padded host buffers remain for the unchanged kernel layout.
- Added a Java 17 Maven build and a checksum-pinned Maven 3.9.16 wrapper. Dependencies now come from Maven Central where compatible: JOCL 2.0.6, Commons Compress 1.28.0, JDistlib 0.4.5, and JUnit Jupiter 5.13.4.
- Retained the project's multi-delimiter JavaCSV fork by generating sources from `lib/javacsv-src.zip`; the Central artifact is not behaviorally equivalent.
- Excluded `src/gov/nih/exon` from this build because its QGeneric/qstats/qplugin dependencies are not present.

### Changes

- Added `pom.xml`, `mvnw`, `mvnw.cmd`, and `.mvn/wrapper/maven-wrapper.properties`.
- Added the backend-neutral implementation under `src/gov/nih/gpu` and the JOCL implementation under `src/gov/nih/gpu/opencl`.
- Migrated `QeQTLAnalysis`, `QeQTLSNPJobReal`, and `QeQTLSNPJobCat` to the new API.
- Updated JDistlib package imports in the eQTL core and standalone tools, removed two missing external utility references, and brought all 83 non-exon/generated Java sources into the build.
- Removed the obsolete `src/gov/nih/opencl` JavaCL adapters and both tracked `lib/javacl-*-shaded.jar` files.
- Updated the Eclipse classpath to Java 17 plus Maven dependencies, and ignored Maven build output.
- Added lifecycle, device-filtering, and production-kernel numerical tests under `test/gov/nih/gpu`.

### Verification

- `.\mvnw.cmd clean package` — successful Java 17-targeted compilation of 83 sources, five-test execution, and runnable shaded jar creation.
- `.\mvnw.cmd test` — 5 tests passed, 0 failed, 0 errored, 0 skipped.
- The numerical integration test compiled and ran the production `eqtlReal` OpenCL kernel and matched a CPU dot product within `1e-12`.
- `java -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — successful, exit code 0.
- Hardware used for this session: NVIDIA GeForce RTX 2080; NVIDIA CUDA OpenCL platform; driver 610.88; OpenCL 3.0 CUDA; FP64 available; maximum allocation 2,147,401,728 bytes; maximum work-group size 1,024.
- A current Java runtime emits a native-access warning while JOCL loads its JNI library. It does not occur as a failure in the Java 17 target and may be suppressed on newer Java runtimes with `--enable-native-access=ALL-UNNAMED`.

### Known limitations and next work

- No representative end-to-end eQTL input/output fixture exists yet. The real kernel is numerically tested, but loading, covariate residualization, thresholding, and output should be locked down with a small reference dataset before deeper changes.
- Only the OpenCL/JOCL backend is implemented. Vendor-specific CUDA/HIP/Vulkan/Level Zero backends can now be added without changing analysis scheduling, but should be justified by profiling and tested against the same CPU reference.
- The application still accepts one INI path rather than full command-line arguments. The obsolete INI `library_path` setting is no longer used for GPU loading.
- Genotype and expression matrices are still read fully into RAM.
- Sample IDs are not yet checked/aligned across genotype, expression, pedigree, and covariate inputs.
- Categorical covariates are not automatically one-hot encoded. The categorical-SNP branch is still deliberately blocked by the legacy `This branch is broken` guard.
- The legacy JavaCSV source emits two removal warnings, JOCL's OpenCL 1.2-compatible queue constructor is deprecated by newer OpenCL APIs, and the shaded jar reports overlapping `META-INF` license/manifest/module descriptors; none blocks the verified build, but packaging notices should be consolidated before a formal release.
- Recommended next session: add a deterministic end-to-end fixture and ID matching/validation first, then categorical covariates and the command-line interface. Those correctness boundaries should precede chunked matrix loading, interactions, and forward selection.

## 2026-08-20T11:32:49.0715718-04:00 — Eclipse Maven dependency integration

### Baseline and goal

- Baseline commit: `7cf681e1f4e76b74d773075911a210fcf4824420`.
- Goal: make the existing Eclipse project resolve the Maven-managed JOCL, JDistlib, Commons Compress, and JUnit libraries instead of reporting unresolved imports.
- Root cause: `.classpath` referenced the m2e Maven container, but `.project` had neither the Maven builder nor Maven nature. Eclipse compiler settings also targeted Java 21 while Maven targets Java 17, the test source folder was absent, and the custom JavaCSV fork needed an explicit Eclipse library attachment.

### Decisions and changes

- Added the `org.eclipse.m2e.core.maven2Builder` builder and `org.eclipse.m2e.core.maven2Nature` nature to `.project`.
- Expanded `.classpath` with Maven-derived attributes, the `test` source folder, Java 17, Maven Dependencies, `target/classes` output, and the tracked `lib/javacsv.jar` with `lib/javacsv-src.zip` attached as source.
- Aligned `.settings/org.eclipse.jdt.core.prefs` to Java 17 with release compilation enabled and added `.settings/org.eclipse.m2e.core.prefs`.
- Added Eclipse import/update instructions to `README.md` and the m2e requirement to `AGENTS.md`.
- Verification exposed legacy `jdistlib.*` imports in `SVALite.java`; restored them to `net.sourceforge.jdistlib.*` to match JDistlib 0.4.5.

### Verification

- `.\mvnw.cmd clean test dependency:build-classpath '-Dmdep.outputFile=target/eclipse-maven-classpath.txt'` — successful; 83 sources compiled and 5 tests passed with 0 failures, errors, or skips.
- The resolved Maven classpath contains JOCL 2.0.6, JDistlib 0.4.5, Commons Compress 1.28.0, JUnit Jupiter 5.13.4, and their transitive dependencies.
- `lib/javacsv.jar` contains both `com/csvreader/CsvReader.class` and `CsvWriter.class`; its tracked source archive is attached in Eclipse.
- Parsed `.project` and `.classpath` as XML after editing.

### Compatibility and next step

- Eclipse must include Maven Integration for Eclipse (m2e) and a Java 17 or newer JDK. An already-open workspace must run **Maven → Update Project…** (`Alt+F5`) once, followed by **Project → Clean…**, so Eclipse reloads the new nature and dependency container.
- `src/gov/nih/exon` remains excluded because its external QGeneric/qstats/qplugin dependencies are still absent.
- Next: refresh the project in Eclipse and confirm the Problems view is clean. If Maven Dependencies still does not appear, verify that m2e is installed and `JavaSE-17` is configured under Installed JREs.
## 2026-08-20T12:20:05.8903634-04:00 — CUDA/cuBLAS backend and automatic multi-vendor selection

### Baseline and goal

- Baseline commit: `7cf681e1f4e76b74d773075911a210fcf4824420`.
- Goal: determine whether native vendor GPU libraries are preferable to JOCL and add selectable, automatically detected backends for the available NVIDIA/Intel-oriented hardware without changing the double-precision real-valued eQTL interface.

### Decisions and changes

- Added JCuda and JCublas 12.6.0 Maven dependencies. The shaded jar includes their Java classes and platform JNI bindings; NVIDIA's driver, CUDA runtime, and cuBLAS remain system prerequisites.
- Added `CudaGpuBackend`, `CudaGpuDevice`, and `CudaGpuContext` under `src/gov/nih/gpu/cuda`. The context maps the existing row-major `eqtlReal` matrix multiplication to full-double-precision cuBLAS DGEMM as `C^T = B^T A^T`, reuses device buffers, clears unused output capacity, and preserves the existing normalization and output stride.
- Kept CUDA/cuBLAS types behind `GpuBackend`, `GpuDevice`, and `GpuContext`; no vendor handles leaked into `gov.nih.eqtl`.
- Added `AutoGpuBackend` and made `auto` the default. It discovers CUDA first, suppresses only a matching usable NVIDIA OpenCL duplicate, and retains distinct OpenCL devices so mixed-vendor systems can use them together. Explicit `cuda`, `opencl`, and `jocl` selections remain available through `-Deqtl.gpu.backend=...`.
- Kept JOCL as the AMD and Intel path. A native HIP/ROCm backend was not added because there is no AMD validation hardware in scope and the official rocBLAS/hipBLAS interfaces are native C/C++, requiring a separately tested Java bridge.
- Added automatic-selection/fallback tests, expanded the production numerical integration test to exercise CUDA, OpenCL, and auto over multiple traits/SNPs with partial tile capacity, and added a manual transfer-inclusive backend benchmark.
- Updated `pom.xml`, `README.md`, and `AGENTS.md` with backend requirements, selection commands, Intel FP64 limitations, and current compatibility.

### Verification

- `.\mvnw.cmd test` — successful; 10 tests passed, 0 failed, 0 errored, 0 skipped. The hardware suite ran the production `eqtlReal` operation through auto, CUDA/cuBLAS, and JOCL/OpenCL, and every checked matrix cell matched the CPU calculation within `1e-11`.
- `.\mvnw.cmd clean package` — successful Java 17-targeted build; 10 tests passed and `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar` was created with the JCuda, JCublas, and JOCL JNI artifacts included.
- `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — exit code 0; `auto` selected the CUDA representation and omitted the duplicate NVIDIA OpenCL representation.
- `java '-Deqtl.gpu.backend=cuda' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — exit code 0; explicit CUDA selection succeeded.
- `java '-Deqtl.gpu.backend=opencl' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — exit code 0; explicit JOCL/OpenCL selection succeeded.
- `java --enable-native-access=ALL-UNNAMED -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.gpu.GpuBackendBenchmark` — exit code 0. Transfer-inclusive median CUDA/OpenCL times were 10.635/11.753 ms for 128 rows and a 2048x2048 tile, 6.617/6.662 ms for 512 rows and a 1024x1024 tile, and 6.227/6.488 ms for 2048 rows and a 512x512 tile. These microbenchmarks show only a modest, shape-dependent CUDA advantage and are not a promise for complete-analysis runtime.
- Hardware: NVIDIA GeForce RTX 2080, compute capability 7.5, CUDA driver API 13.3, CUDA Runtime 12.6, CUDA maximum allocation reported as 8,589,606,912 bytes, and FP64 available. The same card's OpenCL path used NVIDIA driver 610.88, OpenCL 3.0 CUDA, maximum allocation 2,147,401,728 bytes, and FP64 available.

### Known limitations, compatibility, and next step

- The CUDA backend currently implements the active real-valued `eqtlReal` operation only. The legacy categorical-SNP branch remains deliberately blocked and would require a separately verified CUDA operation or explicit OpenCL selection after that path is repaired.
- CUDA is not universally faster than OpenCL. Automatic discovery prefers CUDA for duplicate NVIDIA devices, but users can force OpenCL and should benchmark representative sample counts/tile sizes when throughput matters.
- Most client Intel Iris Xe devices lack native FP64 and therefore cannot satisfy this project's scientific execution contract. They may appear in diagnostics but are filtered from analysis. No reduced-precision fallback was introduced.
- AMD GPUs currently use JOCL/OpenCL. Native HIP/ROCm or rocBLAS support remains a future backend that needs an AMD test machine and the same CPU-reference tests.
- Current Java releases may emit a restricted-native-access warning unless `--enable-native-access=ALL-UNNAMED` is supplied; Java 17 remains supported.
- No representative end-to-end input/output fixture exists yet. The next recommended work remains deterministic reference fixtures plus sample/gene/SNP ID validation before categorical covariates, CLI replacement, or streaming loaders.

## 2026-08-20T12:49:00.7068131-04:00 — Standalone Git repository migration

### Baseline and goal

- Combined-repository baseline: `7af8b27df2ab434ac33ec94168d7b7abd4276052` on `master`, with the working-tree root at `D:\git` and the old combined remote `ssh://robbyjo@desktop/c:/git`.
- Project-only split baseline: `7202e53b45d9c2332787ad1a6f408d98438a9d19` (`Allow CUDA backend as well`).
- Goal: separate NIH-Project from the multi-project repository and its corrupt unrelated `Analysis-Pipeline` index entry, preserve NIH-Project history and current files, and establish independent local and server repositories without deleting the old combined copy.

### Decisions and changes

- Confirmed that the combined repository had only `master` and no tags, then extracted the `NIH-Project` subtree into a project-rooted history. The split examined 2,697 combined commits and retained 31 project commits, beginning with `bfd4f88414467e0972eb1792ae7a401b66bdd894` (`Initial commit`).
- Created the standalone bare remote at `D:\git\NIH-Project.git` on DESKTOP (`100.96.31.61`) and seeded its `master` branch through SSH.
- Set the canonical remote URL to `robbyjo@desktop:D:/git/NIH-Project.git`. The SCP-style path is required because the `ssh://.../d:/...` spelling is interpreted with an invalid leading slash by the Windows SSH Git service.
- Created the canonical standalone local checkout at `D:\projects\NIH-Project`. Its repository root is now the project directory itself, so Git operations no longer inspect the combined `D:\git` index.
- Overlaid the previous `D:\git\NIH-Project` working files while excluding `.git`, `target`, and `bin`. After Git normalization, the overlay had no content differences from the split tip; the apparent modifications were LF/CRLF stat noise only.
- Left the original `D:\git` working tree and the old `C:\git` combined remote unchanged as recovery copies. No old repository, branch, or project directory was deleted.
- Source-tree change for this migration: appended this record to `SESSION_HISTORY.md`.

### Verification

- `git subtree split --prefix=NIH-Project --branch standalone-master master` in a temporary clone — successful; produced `7202e53b45d9c2332787ad1a6f408d98438a9d19` with 31 commits and project-relative root paths.
- `ssh -o BatchMode=yes robbyjo@100.96.31.61 "git init --bare --initial-branch=master D:/git/NIH-Project.git"` — successful.
- `git push 'robbyjo@desktop:D:/git/NIH-Project.git' standalone-master:master` — successful; created remote `master`.
- `git clone 'robbyjo@desktop:D:/git/NIH-Project.git' 'D:\projects\NIH-Project'` — successful; the initial standalone checkout tracked `origin/master` and contained 31 commits.
- `git diff --exit-code --ignore-space-at-eol` after overlay — exit code 0. Re-indexing produced no staged diff, confirming that current tracked content matched the split tip.
- `.\mvnw.cmd clean package` in `D:\projects\NIH-Project` — successful; 87 production sources and 5 test sources compiled, 10 tests passed with 0 failures/errors/skips, and the shaded runnable jar was created.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` in the standalone checkout — exit code 0. Auto selected the NVIDIA GeForce RTX 2080 CUDA backend; driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64 available.

### Known limitations, compatibility, and next step

- Existing clones of the old combined repository are not automatically redirected. New clones should use `robbyjo@desktop:D:/git/NIH-Project.git`.
- The old combined repository still contains its unrelated corrupt `Analysis-Pipeline` cache entry. Use `D:\projects\NIH-Project` for future NIH-Project development rather than `D:\git\NIH-Project`.
- The migration preserved `master`; there were no other branches or tags to migrate.
- Next: open `D:\projects\NIH-Project` as the Codex/Eclipse workspace and continue the correctness-first roadmap with deterministic end-to-end fixtures and ID validation.

## 2026-08-20T14:57:31.0343291-04:00 — ID-safe CSV analysis, categorical covariates, CLI, and bounded-RAM blocks

### Baseline and goal

- Baseline commit: `d5b2d49428ddb7dfc3507111f91a4eaa43438544` (`Document standalone repository migration`) in the canonical checkout at `D:\projects\NIH-Project`.
- Goal: implement roadmap items 1–3—deterministic scientific fixtures and ID validation, automatic categorical-covariate encoding, and a full command-line interface—plus partial genotype/expression loading for CSV while preserving the legacy INI invocation and double-precision production results.
- The supplied WHI files under `D:\Research\topmed\sqtl\sqtl-whi` were used read-only for schema and validation checks. No participant data, identifiers, or WHI-derived output was copied into the repository.

### Decisions and changes

- Added `QDelimitedMatrixSource`, a metadata-first, re-readable row-block source for headered plain/gzip/bzip2 delimited matrices. It validates field counts and rejects blank or duplicate sample and row identifiers before computation. Column permutations are applied while parsing each numeric block.
- Added `QCovariateTable` and `QSampleAlignment`. Genotype and expression headers may use different ID namespaces; matching covariate columns are inferred only when their unique value sets exactly equal the matrix headers, or selected explicitly with `genotype_id_column` / `expression_id_column`. Missing, duplicate, unmatched, or ambiguous IDs are fatal. Covariate-row order is the canonical sample order and reordering counts are reported.
- Replaced the headered-CSV covariate path with mixed-type parsing. Non-numeric selected covariates are automatically categorical, the lexicographically first level is the reference, and the remaining levels are one-hot encoded. `covariate_factor` / `--factor-covariates` forces numeric-looking categories. Interactions between encoded covariates are supported, missing values are fatal, and a rank-deficient design matrix stops before analysis.
- Added `QeQTLCommandLine`. A positional INI file remains valid; `--config` permits command-line overrides; and argument-only runs default to CSV. Added file, covariate, ID bridge, threshold, model, block, thread, output, backend, debug, validation, and streamed-row options. `--gpu-backend auto|cuda|opencl` is applied before lazy GPU runtime creation. Fixed `QeQTLAnalysisConfig.setDFOffset`, which incorrectly wrote `block_size`.
- Added `--validate-only`, which performs complete metadata/identifier scans, alignment, covariate encoding, rank checking, and degree-of-freedom calculation without initializing a GPU or writing analysis results.
- Added `QeQTLPreprocessor`, `QeQTLStatistics`, and `QeQTLStreamedJobReal`. Supplying `genotype_block_rows` or `expression_block_rows` enables a bounded-RAM real-valued CSV path. Each block uses the same covariate residualization, sample standardization, FP64 kernel/cuBLAS operation, threshold, effect, t, p-value, and output identifier formulas as the full-memory path. Source block capacities are rounded to the 16-row GPU tile boundary, GPU contexts remain worker-exclusive, queued genotype blocks are bounded by worker count, and output failures terminate the analysis.
- Added the synthetic `test/resources/eqtl-reference` fixture with intentionally permuted and differently named genotype/expression samples plus a text categorical covariate. Added tests for duplicate identifiers, exact bridging/permutations, automatic one-hot encoding, full rank, INI/CLI compatibility, and fixed reference identifiers, pair count, degrees of freedom, R-squared, effects, t statistics, and log10 p-values.
- Updated `README.md` and `AGENTS.md` with the current CLI, validation rules, categorical behavior, compression support, bounded-RAM settings, compatibility path, and performance caveat.

### Verification

- `.\mvnw.cmd clean package` — successful Java 17 build of 94 production sources; 17 tests passed with 0 failures, errors, or skips; created `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar`.
- The hardware integration tests exercised `eqtlReal` through auto, CUDA/cuBLAS, and JOCL/OpenCL; all production-kernel cells matched the CPU calculation within `1e-11`.
- Argument-only full-memory fixture run: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\eqtl-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --output target\reference-full.csv --threshold none 0 --block-size 16 --threads 1` — successful; six association rows.
- Bounded-RAM fixture run with the same arguments plus `--genotype-block-rows 16 --expression-block-rows 16` and output `target\reference-stream.csv` — successful. Full and streamed output files were byte-for-byte identical, SHA-256 `70FF175A10C5B006BD5458A6947EC8B6960A4D54A6BD239FAFAE8E809F9939AB`.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --gpu-backend cuda --printgpuinfo` — successful; the new application argument selected CUDA before discovery.
- `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar ... --validate-only` on the synthetic fixture — successful without GPU initialization and without creating the requested output.
- `java --enable-native-access=ALL-UNNAMED -jar D:\projects\NIH-Project\target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --config sqtl-chr1a.ini --validate-only` from the WHI data directory — successful in 49.718 seconds. Reported 250,000 SNPs, 136,840 expression traits, and 2,005 aligned samples; inferred genotype `NWDID` and expression `TORID`; both matrices already had canonical covariate order; covariate rank was 28; regression/error/offset degrees of freedom were 28/1976/7.
- Hardware: NVIDIA GeForce RTX 2080, driver 610.88, compute capability 7.5, CUDA driver API 13.3, CUDA Runtime 12.6, plus the same card's NVIDIA OpenCL FP64 path.
- `git diff --check` — clean apart from Git's informational LF-to-CRLF warnings on existing Windows-oriented files.

### Known limitations, compatibility, and next step

- A positional headered-CSV INI run remains compatible but now performs a complete metadata/ID scan before loading or streaming. Headerless legacy genotype CSV and TPED continue through the old full-memory path and do not receive the new ID validation; `--validate-only` requires headered CSV matrices.
- The bounded-RAM scheduler rereads and residualizes the expression CSV once per genotype block. This strictly bounds matrix RAM but can increase disk and CPU work substantially, especially for gzip/bzip2 input. It should be profiled on representative data before being selected solely for speed.
- With multiple workers, result chunks retain the legacy scheduler's worker-interleaved ordering. The byte-identical full/stream verification used one worker; identifiers and statistics are invariant, but a deterministic global sort is not yet imposed.
- Modern alignment currently covers genotype, expression, and covariate tables. Family/pedigree identifier validation remains in the legacy path. The categorical-SNP execution branch remains disabled; automatic categorical covariates do not repair that separate model.
- The next recommended work is representative profiling of CSV scans, repeated expression reads/residualization, packing, GPU transfer/occupancy, and output writing. Based on those measurements, add an indexed reusable binary cache or indexed VCF/BCF genotype source before interactions or forward selection.

### 2026-08-20T14:59:50.4033764-04:00 verification addendum

- Added a final safety check requiring rounded streamed row capacities not to exceed `block_size`, preventing an oversized host/device result allocation.
- `.\mvnw.cmd test` — successful after that check; 17 tests passed with 0 failures, errors, or skips, including all three available GPU integration backends.
- `.\mvnw.cmd -q -DskipTests package` — successful; refreshed the shaded runnable jar with the final verified sources.

## 2026-08-20T17:00:10.8833586-04:00 — Prepared-matrix cache, representative verification, and checkpoint/resume

### Baseline and goal

- Baseline commit: `9dfed723295cd3f4d212ce91c3ca1a20f62e533e` in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: complete the next three priorities requested by the user: replace repeated CSV parsing/preprocessing with a reusable indexed bounded-RAM representation, verify it on a representative WHI subset, and add safe checkpoint/resume. Move all later modernization work to an explicit to-do list.
- The supplied WHI files under `D:\Research\topmed\sqtl\sqtl-whi` were used read-only. Temporary subsets, caches, checkpoints, and result files were created only under ignored `target\whi-subset`; no WHI participant data or derived result was added to Git.

### Decisions and files changed

- Added `src/gov/nih/eqtl/io/QBinaryMatrixCache.java`. The first bounded-RAM run writes aligned, covariate-residualized, standardized double-precision rows to an indexed binary cache; later matching runs read arbitrary blocks without reparsing or re-residualizing the CSV. Each row includes its identifier, original standard deviation, FP64 values, and CRC32 integrity checksum. Cache construction uses a temporary file followed by an atomic replacement.
- Cache signatures cover the format version, matrix kind, normalized source path, source size and modification time, scanned row/sample metadata, alignment permutation, and every double bit of the covariate projection. Incompatible or corrupt cache headers are not reused. Cache filenames include the first 20 hexadecimal signature characters; stale caches are retained for explicit future pruning rather than deleted automatically.
- Added `src/gov/nih/eqtl/QAnalysisCheckpoint.java`. Bounded-RAM output is written as one independently buffered `.partial` file per genotype block and promoted to `.part` only after successful close. Resume validates an analysis signature covering both prepared caches, block sizes, threshold, degrees of freedom, and output mode, skips complete genotype blocks, and deterministically assembles the final output in genotype-block order through an atomic output replacement.
- Checkpoints default to `<output>.checkpoint` and are cleaned after a successful assembly. `--keep-checkpoints` retains them; `--resume` is required when a checkpoint directory exists; and a mismatched manifest is fatal. Added `--cache-dir`, `--rebuild-cache`, `--checkpoint-dir`, `--resume`, and `--keep-checkpoints`, with matching INI keys.
- Updated `QeQTLAnalysis` and `QeQTLStreamedJobReal` to schedule cached genotype blocks and indexed expression blocks through the existing FP64 GPU context pool. The GPU operation and statistical formulas were not changed. The full-memory analysis path remains available and unchanged.
- Added hardware-independent cache/checkpoint tests in `test/gov/nih/eqtl/io/QBinaryMatrixCacheTest.java` and `test/gov/nih/eqtl/QAnalysisCheckpointTest.java`, expanded command-line parsing tests, and added the generic test-only `QCsvSubsetTool` used to create reproducible row-limited validation inputs without embedding WHI identifiers or data.
- Updated `README.md` and `AGENTS.md`, and added `TODO.md`. Remaining work is now explicitly grouped into cache-backed scheduling/profile optimization; indexed VCF/VCF.gz/BCF and output formats; SNP interactions, forward selection, and categorical-SNP repair; and AMD/Intel backend validation plus release packaging.

### Verification

- `.\mvnw.cmd test` — successful; 20 tests passed with 0 failures, errors, or skips. The hardware integration tests exercised the production real-valued operation through auto, CUDA/cuBLAS, and JOCL/OpenCL and matched the CPU calculation within the existing `1e-11` tolerance.
- `.\mvnw.cmd -q clean package` — successful; produced `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` with all 20 tests passing.
- Synthetic full-memory reference run: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\eqtl-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --output target\reference-full.csv --threshold none 0 --block-size 16 --threads 1` — successful; six association rows.
- Synthetic cache/checkpoint run used the same command with output `target\reference-stream.csv` plus `--genotype-block-rows 16 --expression-block-rows 16`; it was byte-for-byte identical to the full-memory output, SHA-256 `70FF175A10C5B006BD5458A6947EC8B6960A4D54A6BD239FAFAE8E809F9939AB`. A second run logged reuse of both caches. A retained, fully complete checkpoint resumed without submitting GPU block work and produced the same output.
- Representative input creation used `java -cp target\test-classes;target\classes gov.nih.eqtl.QCsvSubsetTool <source> <destination> 128` for the genotype and expression matrices. The ignored fixtures contained 128 SNPs, 128 expression traits, and all 2,005 sample columns; the original covariate file was used directly.
- WHI full-memory reference command: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --config D:\Research\topmed\sqtl\sqtl-whi\sqtl-chr1a.ini --genotype target\whi-subset\genotype-128.csv --expression target\whi-subset\expression-128.csv --covariates D:\Research\topmed\sqtl\sqtl-whi\sqtl-whi-mastermat-withpcs-2005.csv --genotype-id-column NWDID --expression-id-column TORID --output target\whi-subset\full.csv --threshold none 0 --block-size 64 --threads 4` — successful. It reported zero sample reorders, covariate rank 28, error degrees of freedom 1976, and offset degrees of freedom 7.
- The WHI cache/checkpoint command used the same arguments with output `target\whi-subset\cached.csv` plus `--genotype-block-rows 32 --expression-block-rows 32 --cache-dir target\whi-subset\cache --checkpoint-dir target\whi-subset\checkpoint --keep-checkpoints`. It created four genotype checkpoint parts. The full and cached results each contained 16,384 associations and their rows had zero differences after sorting.
- Rerunning that WHI command with `--resume` found all four completed blocks, submitted no GPU block work, and reproduced SHA-256 `F2C74FCA2DCA1EC9F297431B3E35FB2CC5336EBE8694BB054C00CD47236A076D`. After deleting exactly ignored temporary parts `block-00000002.part` and `block-00000003.part` to simulate interruption, `--resume` reported `2 / 4` completed, recomputed only genotype row offsets 64 and 96, restored four parts, and produced an exactly identical final output.
- Hardware for GPU verification: NVIDIA GeForce RTX 2080, compute capability 7.5, NVIDIA driver 610.88, CUDA driver API 13.3, CUDA Runtime 12.6, and an FP64-capable NVIDIA OpenCL path.
- `git diff --check` — no whitespace errors; only Git's informational LF-to-CRLF working-copy warnings were emitted for Windows-oriented files.

### Known limitations, compatibility, and next step

- A run still performs the complete source metadata/ID scan before cache selection so alignment and duplicate-ID validation cannot be bypassed. The signature deliberately uses path, size, and modification time rather than hashing every byte of multi-gigabyte source files; preserving size and timestamp while changing content could evade source-change detection.
- The prepared expression cache can still be reread once per genotype block. This removes repeated text decompression/parsing and covariate preprocessing but does not guarantee a speedup when prepared blocks exceed the operating-system file cache. Profile representative block sizes, worker counts, storage, host packing, transfers, GPU time, and output assembly before changing the schedule.
- Prepared caches require additional disk space. The complete WHI example is expected to require roughly 6–7 GB. Old signature-specific caches are not automatically pruned.
- Checkpoint parts are atomic but do not yet carry independent content hashes. They protect against ordinary process interruption because `.partial` files are never accepted, while storage corruption outside that protocol is not independently detected.
- Cache/checkpoint mode currently applies to the headered CSV bounded-RAM real-valued path. Full-memory runs, legacy headerless input, TPED, and the disabled categorical-SNP path do not gain resume support.
- Next recommended work is the first section of `TODO.md`: profile and optimize the cache-backed schedule without changing results. Indexed VCF/VCF.gz/BCF support should precede forward selection; SNP interactions remain a separately verified statistical milestone.

## 2026-08-20T17:35:46.5647026-04:00 — Optional FP32, VRAM-aware sizing, automatic workers, and multi-GPU preservation

### Baseline and goal

- Baseline commit: `eb85f6ac534617233b17211735b535a66f72f598` (`RAM optimization, adding checkpoint / resume`) in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: add an explicitly selectable FP32 calculation, replace the obsolete `lambda`/RAM-saturation concept with GPU-memory-aware block sizing when `block_size` is omitted, choose a useful worker count when `num_threads` is omitted, and preserve use of every distinct compatible GPU.
- User approval for an optional float calculation was explicit. FP64 remains the default; the FP32 accuracy study and its limitations are recorded below.

### Decisions and files changed

- Added `GpuPrecision` and the `precision = fp64|fp32` / `--precision fp64|fp32` setting. FP64 remains the default. FP32 converts only prepared input tiles to `float`; source parsing, sample alignment, covariate projection, residualization, standardization, cache storage, standard deviations, and downstream statistics remain `double`.
- Extended `GpuContext` with precision-aware compilation and an FP32 execution method. `CudaGpuContext` now uses cuBLAS SGEMM for FP32 and DGEMM for FP64; `JoclGpuContext` uses the same production kernel with `DATATYPE float` or `double`. Both reuse precision-sized device buffers and fail when the compiled and requested precisions disagree.
- FP32 mode allows available OpenCL GPUs without native FP64, including a compatible Intel Iris Xe driver/device. FP64 discovery still rejects non-FP64 devices. CUDA-first automatic discovery suppresses the duplicate NVIDIA OpenCL representation even when the NVIDIA device is usable only for FP32, while retaining distinct Intel/AMD OpenCL devices.
- Added total global-memory reporting to `GpuDevice` and `--printgpuinfo`. `GpuTuning` calculates an automatic block size from the least-capable selected device, numeric width, padded sample count, two input buffers, quadratic output buffer, total VRAM, maximum single-allocation size, Java array limit, 12.5%/minimum-256-MiB VRAM headroom, and a 1-GiB target output allocation. Choices are rounded to the 16-row tile or 512-row large-block granularity. For the 2,005-sample, 250,000-SNP, 136,840-trait WHI dimensions and an 8-GiB RTX 2080, the tested recommendations are 11,264 rows for FP64 and 16,384 for FP32.
- An omitted or zero `num_threads` now remains distinguishable from an explicit value. `GpuTuning` recommends up to one worker per GPU for bounded-RAM jobs and up to two per GPU for full-memory jobs, capped by schedulable jobs and all but one CPU core. Streamed recommendations use the number of genotype-block jobs rather than the larger cross-product iteration count. Explicit values are retained with warnings for CPU oversubscription or unused GPUs.
- Preserved multi-GPU execution: initialization opens one context for every distinct precision-compatible device and `GpuContextPool` hands each exclusively to a worker. Mixed devices use the block size supported by the smallest device. No vendor handle was introduced into analysis code.
- Removed the unused `kLambda` argument/constant and the legacy commented RAM-saturation calculation. Existing INI files containing `lambda` remain accepted but receive a notice that the setting is ignored; new configurations should omit it.
- Added a finite/range check for GPU correlations. Tiny precision-dependent excursions above absolute one are clamped only within an explicit tolerance (`1e-4` for FP32, `1e-10` for FP64); larger or non-finite results stop the analysis.
- Precision is included in the checkpoint analysis signature, so FP32 cannot resume FP64 parts or vice versa. Prepared matrix caches remain reusable FP64 data for either mode.
- Added `GpuTuningTest`, expanded CUDA-first duplicate filtering, CLI/config, cache-signature, and production-kernel integration tests, and updated `README.md`, `AGENTS.md`, and `TODO.md`.

### Verification

- `.\mvnw.cmd clean package` — successful Java 17 build of 98 production sources and 13 test sources; 26 tests passed with 0 failures, errors, or skips; created `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`.
- `\.\mvnw.cmd -Dtest=GpuKernelIntegrationTest test` — successful; six production-kernel tests exercised FP64 and FP32 through automatic, CUDA/cuBLAS, and JOCL/OpenCL paths. FP64 retained the existing `1e-11` tolerance; FP32 used an explicit `max(2e-5, |expected| * 5e-6)` CPU-reference tolerance.
- Hardware-independent tuning tests verified 11,264-row FP64 and 16,384-row FP32 recommendations for an 8-GiB device at WHI dimensions, selection of 8,192 rows when a second device had a 512-MiB maximum allocation, and worker recommendations that use multiple GPUs without selecting every CPU core.
- Synthetic FP64 command omitted block/thread settings: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\eqtl-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --threshold none 0 --output target\precision\fp64.csv` — successful; selected block size 16 and one worker automatically.
- The same synthetic command with `--precision fp32` produced six matching identifier pairs. Maximum FP32-versus-FP64 absolute differences were `9.01e-8` in R-squared, `4.00e-8` in effect, `2.87e-7` in t, and `1.47e-7` in log10 p.
- The synthetic FP32 bounded-RAM command added `--genotype-block-rows 16 --expression-block-rows 16 --cache-dir target\precision\cache --checkpoint-dir target\precision\checkpoint`; its output was byte-for-byte identical to the full-memory FP32 output, SHA-256 `6577B25E457DE30DE919DF0A27A5454BD8241A27292651F942C63A5DE2A75595`.
- Explicit OpenCL application command: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --gpu-backend opencl --precision fp32 --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\eqtl-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --threshold none 0 --output target\precision\fp32-opencl.csv` — successful with automatic block size 16 and one worker. Its six identifiers matched CUDA FP32; maximum absolute R-squared difference between the backends was `6.89e-8`.
- `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --precision fp32 ... --validate-only` on the synthetic fixture with no block/thread settings — successful without GPU initialization or output creation; reported that tuning is deferred until a real GPU run.
- Representative WHI FP64/FP32 commands used `--config D:\Research\topmed\sqtl\sqtl-whi\sqtl-chr1a.ini`, ignored 128-row genotype and expression subsets under `target\whi-precision`, explicit `--genotype-id-column NWDID --expression-id-column TORID`, `--threshold none 0 --block-size 64 --threads 1`, and outputs `fp64.csv`/`fp32.csv`; the FP32 command additionally used `--precision fp32`. Both runs completed 16,384 associations with identical identifiers, zero sample reorders, covariate rank 28, error degrees of freedom 1976, and offset 7.
- Across those 16,384 WHI associations, FP32 versus FP64 maximum/RMS absolute differences were: R-squared `5.59e-9` / `2.25e-10`; effect `1.03e-8` / `8.66e-10`; t `9.70e-7` / `1.76e-7`; log10 p `2.51e-6` / `1.29e-7`. Both modes classified exactly 11 associations at p <= `1e-4`, with zero classification differences.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — successful; reported NVIDIA GeForce RTX 2080, 8,589,606,912 global/max-allocation bytes, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, and FP64 support. OpenCL testing used the same RTX 2080 through NVIDIA OpenCL 3.0 CUDA / driver 610.88.
- Final shaded-jar SHA-256: `1CE25488BF7E94DCB8CE3D3A086027A18E1BBF5C117544468A0380BA6D16B2E8`.

### Known limitations, compatibility, and next step

- FP32 is approximate. The representative study was stable at p <= `1e-4`, but associations near any reporting threshold can change inclusion; important borderline findings should be rerun in FP64. FP32 is not enabled automatically.
- FP32 reduces GPU and packed host tile width, but prepared binary caches remain FP64 and retain the same disk footprint. FP32 does not repair or enable the categorical-SNP branch.
- Automatic block sizing uses reported total VRAM and maximum allocation, not live free VRAM. Other GPU processes, driver reservations, JVM heap limits, or unusual drivers can still require an explicit smaller `block_size`. The 1-GiB output target is a conservative deterministic heuristic, not yet an end-to-end performance optimum.
- Only one physical RTX 2080 was available. CUDA and OpenCL FP32 were hardware-tested, but Intel Iris Xe, AMD, and actual multi-GPU throughput remain unverified. Multi-device selection, least-capacity sizing, thread recommendations, and exclusive context pooling are covered by hardware-independent tests.
- Automatic full-memory mode permits two workers per GPU for overlap; bounded-RAM mode uses one. Representative profiling should confirm these choices and tune them only with end-to-end evidence.
- Next: perform the cache/scheduler profiling in `TODO.md`, including FP32/FP64 throughput, real mixed/multi-GPU systems, and Intel Iris Xe OpenCL validation before considering further backend work.

### 2026-08-20T17:38:06.9713754-04:00 verification addendum

- Corrected command spelling for the targeted kernel suite: `.\mvnw.cmd -Dtest=GpuKernelIntegrationTest test`.
- Added a final backend safety guard requiring the declared kernel `DATATYPE` to match the requested precision before CUDA or OpenCL execution.
- `.\mvnw.cmd clean package` after that guard — successful; all 26 tests passed with 0 failures, errors, or skips.
- Superseding final shaded-jar SHA-256: `827AF86985C5A193949EE0C23B35EC5EF2916E5FCFBC205883DD37F971AC1EB3`.

## 2026-08-20T18:37:07.0916546-04:00 — Cohort-aware joint-analysis roadmap

### Baseline and goal

- Baseline commit: `5088eae9bf354f1db5a71b3cb87139db490b81e8` in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: record a further statistical-analysis milestone for jointly analyzing cohorts with different technical covariates, FHS-like familial relatedness, CARDIA-like repeated examinations, and independent WHI/JHS-like samples without materializing complete per-cohort QTL result sets.

### Decisions and files changed

- Updated `TODO.md` with a cohort-partitioned joint-analysis design. Each cohort may have a distinct fixed-effect design; genotype and expression must undergo the same cohort-specific projection, and correlated cohorts require covariance-aware pre-whitening rather than phenotype-only BLUP subtraction.
- Recorded concatenation of transformed sample columns and blockwise score/information accumulation as algebraically equivalent implementation choices that require a deterministic equivalence test.
- Required explicit effect/scaling semantics, effect-allele harmonization, correct rank/degrees-of-freedom accounting, support for FHS kinship and CARDIA repeated subjects, bounded tile accumulation, and optional heterogeneity diagnostics without persistent cohort-wide result files.
- No production source, configuration behavior, statistical formula, or build dependency changed.

### Verification

- `git diff --check -- TODO.md SESSION_HISTORY.md` — completed with no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings.
- No Maven or GPU test was run because this was a documentation-only roadmap update. The most recent executable verification remains the preceding clean package with 26 passing tests on the NVIDIA GeForce RTX 2080 (CUDA driver API 13.3, CUDA Runtime 12.6, NVIDIA driver 610.88).

### Known limitations and next step

- This entry defines a future goal only; cohort-specific projection, mixed covariance, repeated-measure handling, joint score accumulation, and heterogeneity testing are not implemented yet.
- The next recommended implementation remains profiling and optimizing the existing cache-backed scheduler on representative WHI data. That measurement should establish stable timing, I/O, memory, block-size, and FP32/FP64 baselines before introducing a new statistical execution model.

## 2026-08-20T19:08:28.6833687-04:00 — Phase profiling and cache/scheduler optimization

### Baseline and goal

- Baseline commit: `5088eae9bf354f1db5a71b3cb87139db490b81e8` in the canonical standalone checkout at `D:\projects\NIH-Project`; the preceding cohort-aware roadmap documentation was already uncommitted.
- Goal: instrument the bounded-RAM schedule, measure cold/warm WHI behavior and FP64/FP32, optimize only demonstrated bottlenecks, and update automatic block/worker choices without changing FP64 scientific output.
- The supplied WHI data under `D:\Research\topmed\sqtl\sqtl-whi` remained read-only. Row-limited inputs, prepared caches, checkpoints, profiles, and results were created under ignored `target\profile-whi` and removed by the final clean build.

### Decisions and files changed

- Added opt-in `--profile` / `profile = true` and `--profile-output FILE` / `profile_output = FILE`. `QeQTLProfiler` records metadata/alignment, cache signature/build/open/read, genotype/expression packing, GPU-context wait, backend buffer setup/upload/compute/download, CPU result/write, kernel compilation, output assembly, and analysis wall time. Concurrent phase totals are explicitly labeled as overlapping.
- Added opt-in backend measurements through `GpuExecutionMetrics`. CUDA and OpenCL synchronize only when profiling is enabled so upload, compute, and download can be separated; ordinary runs do not pay that synchronization cost.
- Replaced per-double `RandomAccessFile` cache reads/writes and byte-at-a-time CRC updates with checksummed bulk row records in `QBinaryMatrixCache`. The version-1 on-disk bytes, signatures, index, FP64 values, and row CRC behavior are unchanged; caches written by the previous implementation were successfully reused by the optimized reader. The cache test now also corrupts a row and verifies a fatal checksum failure.
- Removed clearing of newly returned result arrays immediately before they become unreachable and restricted per-tile offset logging to `--debug`.
- Updated `GpuTuning` so bounded-RAM automatic mode permits up to four pipeline workers per GPU, capped by CPU cores, genotype jobs, and an estimate of available JVM heap versus prepared/packed/result bytes per worker. Explicit settings remain overrides and receive a heap-pressure warning when appropriate.
- Automatic block sizing still treats VRAM as a hard limit, but now preserves up to four genotype jobs per selected GPU when a small/rectangular genotype matrix would otherwise collapse to fewer jobs. Full WHI recommendations remain 11,264 FP64 and 16,384 FP32 rows; the 4,096-by-8,192 benchmark now selects 1,024 rows and four workers automatically.
- Added/updated `QeQTLProfiler`, `GpuExecutionMetrics`, `GpuContext`, both CUDA/OpenCL contexts, analysis/config/CLI/streamed-job plumbing, cache code, tuning code, unit tests, `README.md`, `AGENTS.md`, and `TODO.md`. Indexed VCF/VCF.gz/BCF is now the next implementation priority; full-cohort and real multi-GPU measurements remain performance follow-up.

### Verification and representative measurements

- Created ignored fixtures while preserving all 2,005 sample columns:
  - `java -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QCsvSubsetTool 'D:\Research\topmed\sqtl\sqtl-whi\chr1a-filtered-whi-dose-meanimp.csv' 'target\profile-whi\input\genotype-4096.csv' 4096`
  - `java -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QCsvSubsetTool 'D:\Research\topmed\sqtl\sqtl-whi\sqtl-whi-jointtmmfpkm-ratio-winsor-completeobs-2005.csv' 'target\profile-whi\input\expression-8192.csv' 8192`
- The common profiling command was `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --config 'D:\Research\topmed\sqtl\sqtl-whi\sqtl-chr1a.ini' --genotype target\profile-whi\input\genotype-4096.csv --expression target\profile-whi\input\expression-8192.csv --covariates 'D:\Research\topmed\sqtl\sqtl-whi\sqtl-whi-mastermat-withpcs-2005.csv' --genotype-id-column NWDID --expression-id-column TORID --output target\profile-whi\baseline-cold.csv --cache-dir target\profile-whi\cache --checkpoint-dir target\profile-whi\checkpoint-baseline-cold --genotype-block-rows 1024 --expression-block-rows 1024 --block-size 1024 --threads 1 --precision fp64 --profile --profile-output target\profile-whi\baseline-cold-profile.csv`; warm and optimized runs changed only the output/checkpoint/profile names unless the variant is stated below.
- Before bulk cache I/O, cold analysis wall time was 233.962 seconds: genotype/expression cache creation used 23.054/46.188 seconds, cache reads accumulated 159.969 seconds, and GPU compute accumulated 3.353 seconds. The matching warm run took 159.326 seconds, including 153.168 seconds of cache reads and 4.633 seconds of GPU compute.
- After bulk cache I/O, the exact old caches were reused successfully. Warm analysis fell to 2.061 seconds, including 0.502 seconds of cache reads and 0.545 seconds of GPU compute: 77.3 times faster than the pre-change warm run. A new cold-cache run used `--cache-dir target\profile-whi\cache-optimized-cold`; it took 6.798 seconds, with 4.549 seconds for both cache builds and 0.460 seconds for cache reads: 34.4 times faster than the pre-change cold analysis.
- Worker variants used `--threads 1`, `2`, and `4` with 1,024-row tiles and took 2.061, 1.265, and 1.086 seconds respectively. A 4,096-row/one-worker variant (`--genotype-block-rows 4096 --expression-block-rows 4096 --block-size 4096`) took 2.132 seconds. These measurements motivated the four-worker pipeline and workload-concurrency block cap, with heap-based safety rather than a universal fixed thread count.
- The automatic verification used the common command with `--output target\profile-whi\optimized-auto.csv --checkpoint-dir target\profile-whi\checkpoint-optimized-auto --expression-block-rows 1024 --block-size 0 --threads 0 --profile-output target\profile-whi\optimized-auto-profile.csv`. It reported `workload concurrency` as the block limiter, selected 1,024 rows, estimated 72.33 MiB per worker, selected four workers, and completed the analysis in 1.246 seconds.
- Baseline cold/warm, optimized cold/warm, two-worker, four-worker, final four-worker, and automatic FP64 result files each contained 48,948 associations plus the header and had identical SHA-256 `3E555AFB0251012480BD54FBE8792509424AE48ABB0A85927CCE078C1106C86E`. The 4,096-row block output had a different deterministic row order but zero sorted-line differences.
- The FP32 variant changed `--precision fp32 --threads 4` and took 1.220 seconds versus 1.086 seconds for the comparable FP64 run, because the optimized small-tile workload was no longer GPU-compute-bound. FP32 reported 48,949 associations, one more than FP64 at p <= `1e-4`. Across the 48,948 common results, maximum/RMS absolute differences were R-squared `9.733e-7` / `1.386e-8`, effect `5.712e-7` / `1.021e-8`, t `1.110e-4` / `1.914e-6`, and log10 p `9.903e-4` / `1.061e-5`.
- `.\mvnw.cmd test` was run after the profiling, bulk-cache, and tuning stages; every run succeeded. Final `.\mvnw.cmd clean package` compiled 100 production and 14 test sources, ran 28 tests with 0 failures/errors/skips, and created `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`.
- Hardware integration tests exercised production FP64/FP32 CUDA/cuBLAS and JOCL/OpenCL operations against CPU references. Benchmark hardware: NVIDIA GeForce RTX 2080, compute capability 7.5, CUDA driver API 13.3, CUDA Runtime 12.6, NVIDIA driver 610.88, and 8,589,606,912 bytes reported VRAM.
- `git diff --check` completed with no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings. Final shaded-jar SHA-256: `98D3ABFA389D376202F065AAD5CAA8B5FD7F074D3932BD5BA18242D81D3F14FC`.

### Known limitations and next step

- The performance study used 4,096 SNPs by 8,192 expression traits with all 2,005 WHI samples, not the complete chromosome. Absolute GPU phase times include profiling synchronization, and the one-device worker optimum may differ with larger sample counts, storage, thresholds, JVM heaps, and multiple GPUs.
- Bulk row I/O makes repeated cache reads inexpensive in this subset, but the complete expression cache is much larger and is still reread for every genotype block. Profile a complete WHI chromosome run before deciding whether expression-block sharing or a different checkpoint/scheduling order is worth its complexity.
- The larger FP32 study crossed the p-value reporting boundary once; FP32 remains explicit and approximate, and borderline findings require FP64 verification.
- Next recommended implementation: indexed VCF/VCF.gz and BCF genotype input through a metadata-plus-block source abstraction, with explicit dosage/genotype, allele, multiallelic, missingness, indexing, and sample-alignment policies.

## 2026-08-20T19:41:00.1979023-04:00 — GPU fixed-effect residualization

### Baseline and goal

- Baseline commit: `a3a83ad642dd82eec2892165566aa3c624232c7b` (`Optimization`) in the canonical standalone checkout at `D:\projects\NIH-Project`; the worktree was clean at the start.
- Goal: keep rank-aware QR construction on the CPU, apply the large fixed-effect projection `Y - (Y Q) Q^T` to genotype and expression blocks on CUDA/OpenCL GPUs, preserve a CPU compatibility mode and prepared-cache reuse, support all selected devices, and verify FP64/FP32 accuracy before packaging.

### Decisions and files changed

- Added `residualization = auto|gpu|cpu` and `--residualization {auto|gpu|cpu}`. `auto` is the default and uses GPU projection for the modern headered-CSV paths; `cpu` preserves the prior FP64 projector and its cache signatures. QR decomposition, pivot/rank validation, covariate columns, and degrees of freedom remain unchanged on the CPU.
- Extended `GpuContext` with FP64/FP32 row-residualization operations and projection-resource release. CUDA applies two cuBLAS DGEMM/SGEMM operations without a materialized `N x N` projector. JOCL compiles separate coefficient and projection kernels. Both implementations cache one flattened Q matrix per context, profile setup/upload/compute/download, and release temporary projection buffers before association tiles are allocated.
- Added `QGpuResidualizer`, using a non-owning exclusive `GpuContextPool`. Independent cache blocks may prepare concurrently across all selected devices; Q is uploaded once per participating context. `QBinaryMatrixCache` drains completed futures in input order, so cache identifiers and downstream result ordering remain deterministic.
- Prepared caches remain FP64 and retain the version-1 row/checksum format. GPU-prepared cache signatures additionally include projection version, precision, and selected backend sequence, preventing CPU/FP32/backend-specific prepared values from being silently mixed. Existing CPU caches remain reusable with `--residualization cpu`.
- The modern full-memory headered-CSV path now uses the same shared preprocessor and a prepared-analysis entry point, avoiding a second residualization/standardization pass. The legacy non-headered/non-modern loader remains on its existing CPU preprocessing path.
- FP32 GPU projection is permitted only through the existing explicit FP32 mode; standardization, cache encoding, effects, test statistics, and p-values remain FP64. Updated profiler phases, CLI/config tests, cache signature/order/concurrency tests, multi-context ownership tests, CUDA/OpenCL hardware tests, `README.md`, `AGENTS.md`, and `TODO.md`.
- Production files changed: `QResidualizationMode`, `QGpuResidualizer`, analysis/config/CLI/preprocessor/profiler/cache classes, `GpuContext`, `GpuContextPool`, and CUDA/OpenCL contexts. Tests changed/added: command-line, cache, `QGpuResidualizerTest`, and `GpuKernelIntegrationTest`.

### Verification and measurements

- `.\mvnw.cmd -DskipTests compile` — successful; 102 production sources compiled for Java 17.
- `.\mvnw.cmd '-Dtest=QGpuResidualizerTest,QeQTLCommandLineTest,GpuKernelIntegrationTest' test` — 13 tests passed. CUDA and JOCL/OpenCL each executed FP64 and FP32 projection against a CPU reference, and the second block call verified that Q was not uploaded again.
- `.\mvnw.cmd '-Dtest=QBinaryMatrixCacheTest,QGpuResidualizerTest' test` — 3 tests passed. A two-worker barrier fixture proved concurrent preparation while the cache retained `rs1,rs2,rs3` source order; borrowed contexts released projection resources without being closed before association scheduling.
- Deterministic eight-sample streamed runs used the following command with `MODE` replaced by `cpu` and `gpu`, and mode-specific output/cache/checkpoint/profile names: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\eqtl-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --threshold none 0 --block-size 16 --threads 1 --genotype-block-rows 16 --expression-block-rows 16 --precision fp64 --profile --residualization MODE --output target\residualization-MODE.csv --cache-dir target\residualization-cache-MODE --checkpoint-dir target\residualization-checkpoint-MODE --profile-output target\residualization-MODE-profile.csv`. Both runs produced the same six ordered identifier pairs. GPU-versus-CPU maxima were R-squared `4.441e-16`, effect `2.776e-16`, t `1.332e-15`, and log10 p `7.772e-16`.
- The same deterministic inputs were run without block-row/cache/checkpoint arguments as `target\residualization-full-cpu.csv` and `target\residualization-full-gpu.csv`. Sorted full-memory CPU output exactly matched streamed CPU output; full-memory and streamed GPU outputs had the same differences above.
- Representative read-only WHI inputs were subset under ignored `target\residualization-whi` with:
  - `java -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QCsvSubsetTool 'D:\Research\topmed\sqtl\sqtl-whi\chr1a-filtered-whi-dose-meanimp.csv' 'target\residualization-whi\input\genotype-4096.csv' 4096`
  - `java -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.eqtl.QCsvSubsetTool 'D:\Research\topmed\sqtl\sqtl-whi\sqtl-whi-jointtmmfpkm-ratio-winsor-completeobs-2005.csv' 'target\residualization-whi\input\expression-8192.csv' 8192`
- The FP64 WHI commands were `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --config 'D:\Research\topmed\sqtl\sqtl-whi\sqtl-chr1a.ini' --genotype target\residualization-whi\input\genotype-4096.csv --expression target\residualization-whi\input\expression-8192.csv --covariates 'D:\Research\topmed\sqtl\sqtl-whi\sqtl-whi-mastermat-withpcs-2005.csv' --genotype-id-column NWDID --expression-id-column TORID --genotype-block-rows 1024 --expression-block-rows 1024 --block-size 1024 --threads 4 --precision fp64 --profile --residualization MODE --output target\residualization-whi\MODE.csv --cache-dir target\residualization-whi\cache-MODE --checkpoint-dir target\residualization-whi\checkpoint-MODE --profile-output target\residualization-whi\MODE-profile.csv`, for `MODE=cpu` and `gpu`.
- That 4,096-SNP by 8,192-trait, 2,005-sample, 28-column design reported the same 48,948 associations in both FP64 modes. GPU-versus-CPU maximum/RMS absolute differences were R-squared `9.992e-16` / `2.178e-17`, effect `3.775e-15` / `5.434e-17`, t `1.252e-13` / `3.501e-15`, and log10 p `7.958e-13` / `1.257e-14`. CPU genotype/expression cold cache phases took 1.125/3.701 seconds; GPU phases took 0.982/3.746 seconds. GPU projection itself used 0.053 seconds of compute across 12 blocks, while CSV parsing/cache writing dominated cold preparation.
- Repeating the WHI commands with `--precision fp32`, separate outputs/checkpoints, the existing CPU cache, and a new `cache-gpu-fp32` produced 48,949 associations in each projection mode with no projection-specific threshold crossing. Across all results, GPU-versus-CPU projection maximum/RMS differences were R-squared `4.237e-7` / `8.916e-9`, effect `1.844e-6` / `2.281e-8`, t `6.894e-5` / `1.607e-6`, and log10 p `4.081e-4` / `5.114e-6`. The already documented overall FP32-versus-FP64 boundary difference remains one association.
- `.\mvnw.cmd test` before the final cache-concurrency fixture — 33 tests passed with 0 failures/errors/skips.
- Final `.\mvnw.cmd clean package` — successful; 102 production and 15 test sources compiled, and all 34 tests passed with 0 failures/errors/skips. The hardware suite exercised production CUDA/cuBLAS and JOCL/OpenCL association plus FP64/FP32 residualization.
- Hardware: NVIDIA GeForce RTX 2080, compute capability 7.5, CUDA driver API 13.3, CUDA Runtime 12.6, NVIDIA driver 610.88, and 8,589,606,912 bytes reported VRAM. Only this one physical GPU was available.
- `git diff HEAD --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings. Final shaded jar: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, SHA-256 `055203C2C8691FCFC1BC0895AC61D7B7A29E9B2A61B72C860DAF3D0A9D089EA6`.
- The final clean removed all ignored temporary subset, cache, profile, checkpoint, and comparison files under `target`; the supplied WHI data remained read-only.

### Known limitations, compatibility, and next step

- GPU projection is implemented for the modern headered-CSV full-memory and bounded-RAM paths. Legacy headerless/non-modern loading continues to use its historical CPU preprocessing. `--residualization cpu` is the exact compatibility/cache-reuse choice.
- Prepared-cache persistence deliberately requires one device-to-host return after projection so rows can be standardized in FP64 and reused in later analyses. Feeding a transient residual directly into one association tile would force repeated expression projection for every genotype block and would be slower for this schedule.
- Automatic GPU residualization creates a distinct cache from existing CPU preparation; a complete chromosome can therefore require another 6–7 GB until the user deliberately prunes an obsolete cache. Cache pruning remains a future explicit command.
- Only one RTX 2080 was available. Multi-context scheduling and deterministic ordering are hardware-independent tested, but real multi-GPU throughput, Intel Iris Xe FP32, AMD OpenCL, larger 5,100–5,700-sample cohorts, and higher-rank projections still require measurement. The simple OpenCL coefficient kernel is correct but may merit reduction/tile optimization only after such profiling.
- FP32 projection is approximate and remained stable at p <= `1e-4` in this study, but important borderline findings still require FP64 verification. FP64 also has ordinary last-bit backend rounding and can theoretically move a result exactly on a reporting boundary.
- This milestone handles ordinary fixed effects only. GRM/kinship and repeated-subject covariance require covariance pre-whitening before QR projection; random-effect variance estimation, cohort-aware stacking/score accumulation, and heterogeneity diagnostics remain in `TODO.md`.
- Next recommended implementation remains indexed VCF/VCF.gz and BCF input. After that input foundation, extend the projection abstraction with cohort-specific covariance pre-whitening and validate it at FHS/CARDIA sample sizes.

### 2026-08-20T19:44:04-04:00 verification addendum

- Added a conservative JVM-heap cap after the initial package: GPU cache-preparation concurrency now uses at most 75% of currently available heap based on precision-specific in-flight host arrays. It still uses every selected GPU when the configured block sizes fit safely and prints when heap headroom reduces concurrency. Updated `QeQTLPreprocessor`, `QGpuResidualizer`, `QBinaryMatrixCache`, `README.md`, and `AGENTS.md` accordingly.
- Repeated `.\mvnw.cmd clean package` after this safety change — successful; 102 production and 15 test sources compiled, and all 34 tests passed with 0 failures/errors/skips, including CUDA/OpenCL FP64/FP32 hardware tests.
- `git diff HEAD --check` again reported no whitespace errors, only informational LF-to-CRLF warnings.
- Superseding final shaded-jar SHA-256: `9384728BD8B4CB727818E851C79254E228101505A127807E35F9569A3B535888`.

## 2026-08-21T12:47:10.5570123-04:00 — VCF.gz/BCF input and variant QC

### Baseline and goal

- Baseline commit: `33449cb765665e6a79ce15bda71988318d71543b` (`Residualization now can use GPU`) in the canonical standalone checkout at `D:\projects\NIH-Project`; the worktree was clean at the start.
- Goal: add VCF.gz and BCF genotype input before advancing the remaining roadmap, detect and report monomorphic/singleton/doubleton variants, add MAF/MAC filtering and REF/ALT/rsID/EAF/MAF/MAC/HWE annotation, and add burden/SKAT/SKAT-O to the future statistical-analysis goal.

### Decisions and files changed

- Added HTSJDK 5.0.0 as a shaded dependency. `QVariantMatrixSource` reads plain/gzip/BGZF VCF and BCF 2.1/2.2 without loading the full genotype matrix. `QBcf22Codec` narrowly admits BCF 2.2 while the analysis reader restricts accepted data to the compatible biallelic diploid DS/GT subset; the executable jar still has no external `bcftools` dependency.
- Introduced `QMatrixRowSource` and moved delimited and variant inputs behind the same metadata/block contract. The cache builder, residualization/standardization preprocessor, full-memory path, bounded-RAM path, alignment, cache signatures, checkpointing, and deterministic row order now operate independently of genotype source format. CSV cache-signature bytes remain unchanged; variant signatures include field/filter policies.
- Added `--genotype-format {auto|csv|vcf|vcf.gz|bcf}`, `--genotype-field {auto|DS|GT}`, `--genotype-missing {error|mean|exclude-variant}`, `--multiallelic {exclude|error}`, `--min-maf`, `--min-mac`, and `--variant-qc-output`. Direct argument-only runs now infer the genotype format from `.vcf`, `.vcf.gz`, and `.bcf`, while legacy INI files without `genotype_format` retain their historical TPED default.
- `auto` selects header-declared DS before GT. Missing selected-field values are fatal unless an explicit policy is selected. Mean imputation uses twice the ALT/effect-allele frequency calculated from called samples. Both MAF and MAC filters are inclusive minimums and are combined when both are set; MAC may be fractional for DS.
- Variant rows use collision-resistant `CHROM:POS:REF:ALT` identifiers. Duplicate canonical variants, including filtered-out variants, are fatal. Original rsID, REF, ALT, VCF FILTER, called/missing counts, ALT allele count, allele number, EAF, MAF, MAC, exact biallelic HWE p-value, classification, inclusion, and exclusion reason are written atomically to a TSV. HWE uses usable diploid GT calls even when DS supplies association dosage and is `NA` when none exist.
- Monomorphic variants are always reported and excluded before standardization. Integer MAC 1/2 is reported as singleton/doubleton; fractional DS values are not rounded into those classes. Multiallelic records are explicitly excluded by default or fatal with `error`; split-by-ALT is not silently attempted. Original VCF record FILTER values are annotated but do not implicitly remove variants.
- Added deterministic VCF/VCF.gz/BCF 2.2 fixtures covering unsorted sample headers, DS and GT, dosage imputation, fractional BCF DS precision, MAF/MAC exclusion, monomorphic/singleton/doubleton classification, multiallelic exclusion, original FILTER annotation, sample-column permutation, canonical IDs, and HWE. Added an optional `-Deqtl.test.bcf=FILE` interoperability test for independently generated BCF.
- Updated `README.md` with runnable variant-input examples and scientific semantics. Updated `TODO.md` with production/indexed-access follow-up and a burden/SKAT/SKAT-O roadmap requiring CPU references, set/weight/null-model definitions, bounded GPU scheduling, and a correctly adjusted SKAT-O omnibus p-value.
- Fixed the checked-in Windows Maven wrapper's null indexing of `Get-Item ... .Target` for an ordinary non-symlinked Maven home. This restores the documented `mvnw.cmd` build on this machine.
- Production/build files changed or added: `pom.xml`, `mvnw.cmd`, `QeQTLAnalysis`, `QeQTLAnalysisConfig`, `QeQTLCommandLine`, `QeQTLPreprocessor`, `QMatrixRowSource`, `QDelimitedMatrixSource`, `QVariantMatrixSource`, `QBcf22Codec`, and `QBinaryMatrixCache`. Tests/fixtures and documentation were updated accordingly. No commit was created.

### Verification

- `.\mvnw.cmd -v` — successful after the wrapper correction: Maven 3.9.16 on Eclipse Adoptium Java 21.0.4, targeting Java 17 through the POM.
- `.\mvnw.cmd clean package` — successful; 105 production and 16 test sources compiled. Forty tests ran with 0 failures/errors and 1 intentional skip for the optional external BCF fixture. CUDA/cuBLAS and JOCL/OpenCL integration tests executed against CPU references. The shaded jar was created at `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` and contains the HTSJDK VCF/BCF classes.
- A public HTSlib regression fixture produced independently by HTSlib/bcftools was downloaded temporarily from `https://github.com/samtools/htslib/files/2458695/test-subset.bcf.gz`, unwrapped to `%TEMP%\htslib-test-subset.bcf`, and was not committed. `.\mvnw.cmd '-Dtest=QVariantMatrixSourceTest#readsExternalBcftoolsFixtureWhenConfigured' '-Deqtl.test.bcf=C:\Users\User\AppData\Local\Temp\htslib-test-subset.bcf' test` — 1 test passed with 0 failures/errors/skips; the reader scanned accepted rows from the 189-sample BCF 2.2 file and read a data block.
- `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --genotype-missing mean --min-mac 2 --expression test\resources\variant-reference\expression.csv --output "$env:TEMP\gpu-eqtl-vcf-validation.csv" --variant-qc-output "$env:TEMP\gpu-eqtl-vcf-validation.variants.tsv" --validate-only` — successful through the packaged jar without GPU initialization. It selected DS, scanned 7 records, included 4, reported 1 monomorphic/1 singleton/1 doubleton, aligned four samples, reordered two expression columns, and reported two residual error degrees of freedom.
- `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — automatic discovery selected the NVIDIA GeForce RTX 2080 via CUDA: compute capability 7.5, CUDA driver API 13.3, CUDA Runtime 12.6, 8,589,606,912 bytes VRAM, FP64 available. Only this one physical GPU was present.
- `git diff HEAD --check` — no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings.
- Final shaded-jar SHA-256: `E1C7E5059500E4CA6F2A0F69D1C80BC738648D5DEEC4DF0E5049F984E2BD0193`.

### Known limitations, compatibility, and next step

- VCF/BCF source access is currently sequential. The metadata/QC pass scans the complete file, and an uncached preparation pass rereads it sequentially. The prepared matrix cache is indexed, but tabix/CSI region queries and true random-access source blocks are not yet exposed.
- The current variant analysis is additive, biallelic, and diploid. Multiallelic split-by-ALT, haploid sex-chromosome policy, polyploid calls, per-record DS-to-GT fallback, HWE filtering, and implicit PASS-only filtering are not implemented. BCF FORMAT Float is inherently stored as float32 before being promoted to Java double; downstream preprocessing and FP64 analysis remain double precision.
- The public BCF interoperability check establishes basic BCF 2.2 GT decoding on an independent file, but production adoption still requires cohort-specific comparison with current `bcftools query` counts, sample order, chosen DS/GT values, EAF/MAF/MAC, HWE, exclusions, and several hand-checked variants.
- Recommended next work: first perform that representative VCF.gz/BCF cross-validation, then add tabix/CSI region selection and re-readable indexed variant blocks. That foundation should precede the new burden/SKAT/SKAT-O implementation; build deterministic CPU set-test references and null-model semantics before adding multi-GPU acceleration. SNP interaction, forward selection, and cohort covariance work remain later roadmap items unless reprioritized explicitly.

## 2026-08-21T13:12:22.6178348-04:00 — GPLv3 repository and package licensing

### Baseline and goal

- Baseline commit: `33449cb765665e6a79ce15bda71988318d71543b` (`Residualization now can use GPU`) in the canonical standalone checkout at `D:\projects\NIH-Project`. The preceding VCF/BCF milestone remained uncommitted and was preserved.
- Goal: review the current remaining roadmap and add the complete GNU General Public License version 3 to the repository and runnable package. The active Codex task was found to be attached to the stale combined-repository copy under `D:\git\NIH-Project`; all review, edits, and verification for this entry were performed against the registered standalone `GPU_eQTL` project at `D:\projects\NIH-Project`.

### Decisions and files changed

- Added root `LICENSE` as a byte-for-byte copy of the Free Software Foundation's official GPLv3 plain-text file from `https://www.gnu.org/licenses/gpl-3.0.txt`. This matches the existing source notices, which select version 3 without an `or later` clause. The official and repository files both have SHA-256 `3972DC9744F6499F0F9B2DBF76696F2AE7AD8AF9B23DDE66D6AF86C9DFB36986`.
- Added GNU GPL v3.0 license metadata to `pom.xml` and configured Maven resources to package the root file as `META-INF/LICENSE` in both the ordinary and shaded jars.
- Added a License section to `README.md` linking to the complete repository text.
- Reviewed `TODO.md` without changing it. The dependency order remains: representative VCF/BCF cross-validation, indexed tabix/CSI variant access, deterministic CPU burden/SKAT/SKAT-O references, then bounded multi-GPU set-test acceleration. Full-cohort performance measurements, interactions, forward selection, cohort covariance, backend validation, and release polish remain later or parallel work.

### Verification

- Compared `LICENSE` byte-for-byte with the downloaded official GNU file — identical, 35,149 bytes, SHA-256 `3972DC9744F6499F0F9B2DBF76696F2AE7AD8AF9B23DDE66D6AF86C9DFB36986`.
- `.\mvnw.cmd clean package` — successful Java 17-targeted build of 105 production and 16 test sources; 40 tests ran with 0 failures/errors and 1 intentional skip for the optional external BCF fixture. CUDA/cuBLAS and JOCL/OpenCL production integration tests ran against their CPU references. The shaded jar was created at `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` with SHA-256 `DDC3568C63D55024D442F7B9C4B97140C88985A9D12221870EB9400E974A6CBA`.
- Read `META-INF/LICENSE` directly from the shaded jar — 35,149 bytes with the same official GPLv3 SHA-256.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — successful. Automatic selection found one NVIDIA GeForce RTX 2080 through CUDA; CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, 8,589,606,912 bytes reported VRAM, and FP64 available.
- `git diff HEAD --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Known limitations and next step

- Maven Shade still reports pre-existing overlapping dependency metadata. `META-INF/LICENSE` also overlaps with Commons Logging, but the final shaded artifact was explicitly verified to contain this project's official GPLv3 text. Third-party dependency notices still need consolidation before a formal release.
- No application behavior or statistical formula changed in this licensing step. The current VCF/BCF implementation and its worktree changes remain uncommitted.
- Recommended next work: cross-validate a representative cohort VCF.gz/BCF against `bcftools`, then add tabix/CSI indexed region/block access before implementing the CPU reference layer for burden, SKAT, and SKAT-O.

## 2026-08-21T14:22:34.9440535-04:00 — Declared matrix roles and missing-data policies

### Baseline and goal

- Baseline commit: `33449cb765665e6a79ce15bda71988318d71543b` (`Residualization now can use GPU`) in the canonical standalone checkout at `D:\projects\NIH-Project`. The preceding uncommitted VCF/BCF and GPLv3 work was preserved.
- Goal: support missing predictor, trait, and selected-covariate values without guessing biological matrix roles; keep row-mean genotype imputation as the defensible default; add an explicitly warned nearest-flanking-genotype-pattern option; make exact trait-pattern deletion the default; and write a missingness QC file beside every modern analysis output.

### Decisions and files changed

- Added declared matrix roles `genotype`, `expression`, `methylation`, `proteomics`, and `continuous` through `--predictor-type` and `--trait-type`. Historical `--genotype`/`--expression` and INI runs retain genotype/expression defaults. Generic `--predictor` and `--traits` arguments require their corresponding explicit type and never scan values to infer biology.
- Added common missing policies `error`, `mean`, `zero`, `exclude-row`, and trait-only `pattern`. Predictor mean is the default for genotype and other continuous predictors. Trait exact pattern deletion is the default. Selected fixed covariates default to `complete-samples`, removing the same incomplete samples from predictor, trait, and covariate rows before matrix scans and QR; `error` remains available.
- Delimited input now recognizes blank, `NA`, `N/A`, `null`, `NaN`, and `.` as missing. Column readers now accept a validated unique subset, enabling complete-covariate and pattern-specific sample selections without materializing full matrices. Prepared-cache sample counts and heap estimates use the selected sample count.
- Added `QMissingnessScan` to collect exact row masks, per-row and per-sample counts, row means, and deterministic pattern groups while streaming. Added `QMissingnessReport`; every modern run atomically writes `<output>.missingness.tsv` unless overridden by `--missingness-qc-output`. `--inspect-missingness` writes the audit and exits before GPU initialization.
- VCF/BCF values now remain missing until the common scan and policy layer, so unified QC records original missing calls before mean/zero replacement. Variant MAF/MAC and annotation QC still use called genotypes across the complete variant header sample set; common missing QC uses the final aligned sample set after selected-covariate filtering.
- Added exact GPU trait-pattern analysis. Traits sharing one missing mask are processed together; each mask subsets predictor and trait columns, rebuilds and rank-checks the covariate model and QR, recomputes the reporting threshold, residualizes/standardizes on exactly those samples, and uses bounded prepared caches. Dynamic results append effective `N` and p-value `DF`; legacy complete-data output headers remain unchanged. Temporary pattern caches/checkpoints are verified to stay under a uniquely created output-side work directory and are removed after assembly.
- Added `--predictor-missing local-pattern` and `--predictor-flanks N` (default 1) for declared genotype predictors. `QLocalPatternImputedSource` is bounded-memory and does not average adjacent SNP dosages: it selects called donors with the nearest observed left/right dosage pattern and averages tied donor values. `CHROM:POS...` IDs enforce contiguous chromosomes and nondecreasing position, chromosome boundaries are not crossed, missing annotations cause an explicit row-order warning, and absence of a comparable donor warns and falls back to the row mean. The documentation labels this as an unphased local proxy rather than haplotype/reference-panel imputation; mean remains the default.
- Added deterministic missingness/policy, covariate-removal, VCF preservation, CLI/type, local-donor, and exact pattern-statistics tests plus `test/resources/missing-reference/expression.csv`. The reference asserts ordered IDs, N, residual DF, R-squared, effects, and t statistics for complete N=8 and incomplete N=7 masks.
- Updated `QeQTLStreamedJobReal` and checkpoint signatures so dynamic output includes `N,DF` without allowing an old checkpoint/output schema collision. Updated `README.md`, `TODO.md`, and `AGENTS.md` with usage, scientific constraints, and production follow-up. Production files added or materially changed for this goal include `QDataType`, `QMissingValuePolicy`, `QeQTLAnalysis`, config/CLI/streamed job, `QMatrixRowSource`, delimited/variant/covariate/alignment/cache sources, `QMissingnessScan`, `QMissingnessReport`, `QPolicyMatrixSource`, and `QLocalPatternImputedSource`. No commit was created.

### Verification

- Initial sandboxed `.\mvnw.cmd -DskipTests compile` could not download `maven-antrun-plugin:3.2.0` because network access was denied. Repeating the exact Maven compile with approved dependency access succeeded with 110 sources at the intermediate stage; the final build compiled 111 production sources.
- Intermediate `.\mvnw.cmd test` runs passed 45 and then 47 tests with 0 failures/errors and 1 intentional skip while static policies, common VCF missing preservation, and the deterministic exact-pattern CPU reference were added.
- Final `.\mvnw.cmd clean package` — successful Java 17-targeted build of 111 production and 17 test sources. Forty-eight tests ran with 0 failures/errors and 1 intentional skip for the optional external BCF fixture. CUDA/cuBLAS and JOCL/OpenCL production kernel and residualization integration tests ran against CPU references. Maven Shade produced `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`; its pre-existing overlapping metadata warnings remain.
- Exact dynamic GPU command: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\missing-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --predictor-type genotype --trait-type expression --predictor-missing mean --trait-missing pattern --threshold none 0 --block-size 16 --threads 1 --genotype-block-rows 16 --expression-block-rows 16 --precision fp64 --residualization cpu --output target\missing-pattern-e2e.csv --missingness-qc-output target\missing-pattern-e2e.missingness.tsv` — successful after the clean package. It ran two masks, emitted six ordered associations, reported N/DF `8/4` and `7/3`, and left no `.gpu-eqtl-pattern-*` work directory. A separate PowerShell comparison found maximum absolute R-squared difference `8.88178419700125E-16` from the deterministic CPU constants.
- Hardware for the dynamic run: automatic CUDA-first discovery selected one NVIDIA GeForce RTX 2080, compute capability 7.5, CUDA Runtime 12.6, 8,589,606,912 bytes reported VRAM, FP64 available. Only one physical GPU was present.
- VCF/common-policy command: `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-type genotype --predictor-missing mean --min-mac 2 --expression test\resources\variant-reference\expression.csv --trait-type expression --trait-missing pattern --output target\vcf-missing-validation.csv --variant-qc-output target\vcf-missing-validation.variants.tsv --missingness-qc-output target\vcf-missing-validation.missingness.tsv --validate-only` — successful without GPU initialization; the common audit reported one missing predictor value before mean replacement, four included variants, four aligned samples, and residual error DF 2.
- Generic role inspection initially omitted the fixture's required genotype/expression ID bridge and correctly failed with `Expression data are missing genotype sample 'G3'`. Repeating `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --predictor test\resources\eqtl-reference\genotype.csv --predictor-type methylation --traits test\resources\missing-reference\expression.csv --trait-type proteomics --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --predictor-missing mean --trait-missing pattern --output target\generic-inspection.csv --missingness-qc-output target\generic-inspection.missingness.tsv --inspect-missingness` succeeded without GPU initialization and reported the declared methylation/proteomics roles plus one missing trait value.
- Final shaded-jar SHA-256: `C33016654EDBC7CF84DE32E12A8ADF08FA5FAC50BC0589E9832D80CC7DFE9215`. `META-INF/LICENSE` remains present at 35,149 bytes. `git diff --check` reported no whitespace errors, only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Exact trait-pattern deletion is scientifically correct but intentionally conservative in this first pass: it re-prepares predictor rows once per distinct mask, groups final output by mask rather than global trait-row order, and does not support `--resume` or `--keep-checkpoints`. Many unique or singleton patterns can therefore be extremely slow; the program warns with the exact pattern count.
- Predictor pattern-wise deletion is not implemented because pair-specific predictor and trait missingness would create a combinatorial scheduler. A global all-matrices complete-sample option is also not yet implemented. Current choices are predictor imputation/error/feature exclusion, selected-covariate complete-sample removal, and trait static or exact-pattern policies.
- `local-pattern` is an unphased nearest-dosage-pattern heuristic. It does not estimate haplotypes or use LD distances/reference panels, and the unified audit reports pre-policy missingness rather than one row per imputed call. The runtime reports ordering and mean-fallback warnings, while production validation against mean and a standard external imputation remains required.
- Missing-data policies apply to the modern headered CSV/VCF/BCF matrix path. Legacy headerless/non-modern loading retains historical behavior. Dynamic `N,DF` columns are added only when exact masks differ; complete-data output remains compatible.
- Recommended next work: run `--inspect-missingness` on representative WHI and other cohort inputs, compare all counts/masks independently, then benchmark the number and size distribution of trait masks. Add resumable/reusable pattern preparation only after those measurements; separately compare local-pattern, mean, and standard externally imputed genotypes before changing the mean default.

## 2026-08-24T15:12:12.9147060-04:00 — Prefix-safe subset alignment and once-prepared trait residency

### Baseline and goal

- Baseline commit: `5337afe6f3ec47b5339c376f2de56f5578ceaf86` (`Add analysis with missing values`) in the canonical standalone checkout at `D:\projects\NIH-Project`; the worktree was clean at the start.
- Goal: add explicit literal prefix removal for matrix sample IDs, permit the 7,214-sample batch1234 VCF to select the 4,746 canonical covariate/expression samples without a silent intersection, and allow the completely prepared trait matrix to remain in RAM while genotype blocks stream so trait residualization is never repeated per genotype block.
- The supplied batch1234 directory was inspected read-only to confirm filenames and the `SampleName`/`framid` covariate columns. No participant data or derived cohort output was copied into the repository, and the 7.3-GB VCF plus 5.1-GB trait matrix were not run end to end in this session.

### Decisions and files changed

- Added `--predictor-id-strip-prefix` and `--trait-id-strip-prefix`, with genotype/expression aliases and matching INI keys. The transformation removes a literal leading prefix only where present; unchanged IDs remain valid, while blank normalized IDs and prefix-induced duplicates are fatal. Counts and literal prefixes are reported.
- Added `--sample-alignment strict|covariate-subset`. Strict remains the default. Covariate-subset uses covariate-row order as canonical, requires every unique covariate ID in both normalized matrix headers, and excludes only matrix-header extras. It never takes an implicit intersection. Reordering, prefix transformations, bridge columns, and excluded extras are written as `ALIGNMENT` records in the missingness QC file.
- Deferred VCF/BCF variant scanning until final sample alignment and selected-covariate complete-sample removal. Called/missing counts, EAF/MAF/MAC, HWE, singleton/doubleton classification, monomorphic exclusion, and MAF/MAC filtering now use the actual analysis columns, so VCF-only samples cannot make an analysis-set rare/monomorphic variant pass. The selected columns remain covered by prepared-cache signatures.
- Added the read-only `QPreparedMatrix` abstraction plus `QInMemoryPreparedMatrix` and `--trait-cache auto|memory|disk` (`auto` default). Disk retains the current indexed prepared-cache schedule. Memory reads the verified aligned/residualized/standardized FP64 trait cache exactly once into bounded row arrays shared read-only across workers; association still tiles those rows for GPU work. Auto chooses memory only when the complete estimate plus all configured worker buffers and conservative JVM headroom fit; explicit memory fails early with `-Xmx`/block/thread guidance otherwise. Checkpoint and scientific cache signatures are unchanged between residency modes.
- Exact trait-pattern analysis uses the same residency policy within each pattern. Each pattern necessarily retains its own sample set and QR projection, but each pattern-specific trait row is still residualized/standardized only once.
- Added a trait-memory-load profiling phase and documented that the 114,406-by-4,746 FP64 numeric trait payload is about 4.05 GiB before Java/object/worker overhead. Updated `README.md`, `TODO.md`, and `AGENTS.md` with command semantics, scientific constraints, and remaining pattern-specific rare-variant work.
- Production files added: `QTraitCacheMode`, `QPreparedMatrix`, `QInMemoryPreparedMatrix`, and `QSampleAlignmentPolicy`. Analysis/config/CLI/profiler/streamed-job and alignment/covariate/cache/missingness/variant sources were updated. Deterministic CLI, alignment, cache, missingness, variant-QC, and scientific-reference tests were expanded. No commit was created.

### Verification

- Initial sandboxed `\.\mvnw.cmd -DskipTests compile` and focused-test commands could not resolve `maven-antrun-plugin:3.2.0` because network access was denied. Repeating through the approved Maven dependency/runtime path succeeded.
- Focused `\.\mvnw.cmd '-Dtest=QCovariateTableTest,QBinaryMatrixCacheTest,QeQTLCommandLineTest,QMissingnessPolicyTest' test` — 15 tests passed with 0 failures/errors/skips. These cover strict-versus-subset behavior, literal prefix counts/collisions, prepared disk/memory equality and signature equality, CLI parsing, and alignment audit rows.
- `\.\mvnw.cmd '-Dtest=QVariantMatrixSourceTest,QeQTLReferenceTest' test` — 10 tests ran with 0 failures/errors and 1 intentional optional-BCF skip. A variant carried only by an excluded VCF sample became analysis-set monomorphic and was excluded; the prefix/subset fixture retained the exact six reference associations and R-squared values within `1e-12`.
- Final `\.\mvnw.cmd clean package` — successful Java 17-targeted compilation of 115 production and 17 test sources. Fifty-two tests ran with 0 failures/errors and 1 intentional skip for the optional independently generated BCF. CUDA/cuBLAS and JOCL/OpenCL production association/residualization integration tests ran against CPU references. The shaded jar was created at `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`.
- Packaged FP64 disk-residency command: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\eqtl-reference\genotype.csv --expression test\resources\eqtl-reference\expression.csv --covariates test\resources\eqtl-reference\covariates.csv --fixed-covariates Age,Batch --genotype-id-column genotype_id --expression-id-column expression_id --predictor-missing error --trait-missing error --threshold none 0 --block-size 16 --threads 1 --genotype-block-rows 16 --expression-block-rows 16 --precision fp64 --residualization cpu --cache-dir target\trait-residency-e2e\cache --trait-cache disk --checkpoint-dir target\trait-residency-e2e\disk-checkpoint --output target\trait-residency-e2e\disk.csv --missingness-qc-output target\trait-residency-e2e\disk.missingness.tsv` — successful on CUDA.
- The same packaged command with `--trait-cache memory`, memory-specific checkpoint/output/QC paths, and reuse of both prepared caches — successful. Disk and memory result files were byte-for-byte identical with SHA-256 `70FF175A10C5B006BD5458A6947EC8B6960A4D54A6BD239FAFAE8E809F9939AB`.
- Packaged post-alignment variant validation: `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-type genotype --predictor-missing mean --min-mac 2 --expression test\resources\variant-reference\expression.csv --trait-type expression --trait-missing error --output target\trait-residency-e2e\vcf-validation.csv --variant-qc-output target\trait-residency-e2e\vcf-validation.variants.tsv --missingness-qc-output target\trait-residency-e2e\vcf-validation.missingness.tsv --validate-only` — successful without GPU initialization; reported four aligned QC samples, seven input/four included variants, and the expected monomorphic/singleton/doubleton counts of one each.
- Hardware for application/kernel verification: one NVIDIA GeForce RTX 2080, CUDA Runtime 12.6, compute capability 7.5, 8,589,606,912 reported VRAM, FP64 available. The test suite also exercised the card's JOCL/OpenCL path. Final shaded-jar SHA-256: `2BE0B02B029C9CAEC1F5AFEFB10401905ECCA3B8F17C17D732AB151DC76C2301`.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Memory residency still uses the persistent prepared FP64 cache as its validated source, so it needs cache disk space and one cache-to-memory load. It avoids every later per-genotype-block disk reread and never repeats residualization/standardization. The existing version-1 cache bytes and checkpoint signatures remain compatible.
- `auto` is deliberately conservative and depends on current JVM maximum/used heap plus configured block rows and workers. The complete 4,746-sample trait matrix should use a deliberate JVM heap (for example `-Xmx12g` or larger) and smaller streamed tiles such as 2,048 rows if memory residency is desired; otherwise auto reports its estimate and safely retains disk mode.
- Final aligned-cohort variant QC does not solve the narrower exact-trait-pattern case where removing additional samples for one trait mask makes a rare predictor monomorphic only within that pattern. That behavior remains an explicit follow-up in `TODO.md`; the preprocessor currently fails rather than emitting invalid standardized results.
- Source VCF/BCF access remains sequential: one post-alignment QC scan and an uncached preparation reread. Tabix/CSI indexed regions remain future work. The complete batch1234 production run still needs independent ID/QC comparison and memory/disk profiling before changing the conservative auto heuristic.
- Recommended next step: first run `--inspect-missingness`/`--validate-only` and independent bcftools/R checks on batch1234 using `framid`, `SampleName`, prefix `X`, and `covariate-subset`; then run a representative row-limited association with `--trait-cache auto`, `memory`, and `disk` before launching the full chromosome.

### 2026-08-24T15:12:12.9147060-04:00 command-spelling addendum

- The Maven wrapper spellings in the verification bullets above should be read as `.\mvnw.cmd ...`; the extra leading backslash in two inline records is a documentation typo only.

## 2026-08-24T15:52:17.3037869-04:00 — Long VCF/BCF pass progress reporting

### Baseline and goal

- Baseline commit: `5337afe6f3ec47b5339c376f2de56f5578ceaf86` (`Add analysis with missing values`) in the canonical standalone checkout at `D:\projects\NIH-Project`. The preceding uncommitted prefix/subset-alignment, aligned-sample variant-QC, and trait-residency work was preserved.
- Goal: make the apparently stalled period after `Loading and validating covariate data...` observable, and provide bounded-overhead progress markers during every long sequential VCF/BCF pass without changing scientific calculations, cache signatures, or result ordering.

### Decisions and files changed

- `QeQTLAnalysis` now announces the aligned-sample variant-QC scan, its source path and sample count immediately before the first full VCF/BCF traversal, and flushes the message before scanning starts.
- `QVariantMatrixSource` now samples progress every 1,000 decoded records and reports approximately every 15 seconds or one million records, whichever comes first. Reports include decoded input records, retained variants, elapsed time, and throughput. The initial QC pass cannot know the total record count in advance, while subsequent missingness/cache-building rereads use the QC summary to add percentage and ETA. Long passes emit a completion marker; tiny fixture scans remain quiet.
- The reporter uses decoded records rather than compressed byte offsets, adds no per-sample output, and does not change QC filtering, variant values, identifiers, metadata, cache signatures, or persisted result formats. Added a deterministic formatting test for counts, rate, percentage, and ETA.
- Updated `README.md` with progress semantics and the first-pass total limitation. Material files changed for this goal: `src/gov/nih/eqtl/QeQTLAnalysis.java`, `src/gov/nih/eqtl/io/QVariantMatrixSource.java`, `test/gov/nih/eqtl/io/QVariantMatrixSourceTest.java`, and `README.md`. No commit was created.

### Verification

- `.\mvnw.cmd '-Dtest=QVariantMatrixSourceTest' test` — successful; 8 tests ran with 0 failures/errors and 1 intentional skip for the optional external BCF fixture.
- `.\mvnw.cmd test` — successful; 53 tests ran with 0 failures/errors and 1 intentional optional-BCF skip. The production CUDA/cuBLAS and JOCL/OpenCL kernel integration tests ran against their CPU references.
- Because process 25920 was still running the 4,746-sample chromosome-22 analysis from `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, that locked standard jar was not replaced. An initial isolated-build command inherited the original working directory; its clean phase removed only regenerable `target` contents before the Windows lock stopped it at the active jar. Source files and the running jar were unaffected. The corrected command explicitly changed into a verified temporary copy, ran `.\mvnw.cmd clean package`, compiled 115 production and 17 test sources, reran all 53 tests successfully, copied the result to `target\gpu-eqtl-2.0.0-SNAPSHOT-progress-all.jar`, and removed the temporary copy.
- Packaged validation command: `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-progress-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-type genotype --predictor-missing mean --min-mac 2 --expression test\resources\variant-reference\expression.csv --trait-type expression --trait-missing error --output target\progress-validation.csv --variant-qc-output target\progress-validation.variants.tsv --missingness-qc-output target\progress-validation.missingness.tsv --validate-only` — successful. It printed the new aligned-sample QC stage before scanning and retained the expected 4 of 7 variants over 4 aligned samples; no GPU was initialized.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-progress-all.jar --printgpuinfo` reported one NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, 8,589,606,912 bytes global memory, and FP64 support. Final alternate shaded-jar SHA-256: `31FC9506ADE9BD783E3EDEC5D87B538085D1FB9603E50A2FE27A1EC10DAD304A`.
- `git diff --check` reported no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- The already-running chromosome-22 process loaded the earlier jar and cannot display these new messages. It should be allowed to finish. The distinct `-progress-all.jar` is ready for the next run; after process 25920 exits, run `.\mvnw.cmd clean package` to refresh the conventional `-all.jar` path.
- The first QC scan reports record counts and rate but not percentage/ETA because an unindexed compressed VCF/BCF does not expose a reliable total record count before traversal. Later passes do report percentage and ETA. ETA is throughput-based and can change with record complexity and storage/cache state.
- This change instruments VCF/BCF traversal, which caused the observed silent interval. Large delimited-matrix metadata/missingness scans still have only phase-level messages. Add their own row/byte progress reporter if representative production profiling shows similarly long silent phases.

## 2026-08-24T16:26:30.4438853-04:00 — Aligned-sample variant-level QC parallelism

### Baseline and goal

- Baseline commit: `730b11f7839a33a1999c3e46768f298c42c2d872` (`Improve cohort alignment and streamed variant analysis`) in the canonical standalone checkout at `D:\projects\NIH-Project`; the prior progress-reporting source change was already part of that commit and the worktree was clean at the start.
- Goal: parallelize the expensive MAF/MAC/EAF, genotype missingness, monomorphic/singleton/doubleton classification, and exact HWE calculations across variants, while preserving the established requirement that genotype, expression, and selected-covariate samples are aligned before those quantities are computed.

### Decisions and files changed

- `QVariantMatrixSource` now keeps HTSJDK VCF/BCF cursor traversal and lazy-genotype decoding on one reader thread, then submits already-decoded distinct variants to a bounded fixed worker pool. This avoids concurrent use of the codecs' mutable decode state while parallelizing the per-sample/per-variant work that dominates at cohort sizes around 5,000.
- Automatic QC sizing uses `availableProcessors - 1`, with at least two variant workers on a multi-core JVM and a cap of 16. Added `--variant-qc-threads N` and INI `variant_qc_threads`; zero/omitted is automatic, a positive number is an override, and one is the sequential comparison path. The selected count is printed before the scan.
- Worker futures are bounded to two pending records per worker and consumed in input order. Duplicate detection, QC TSV rows, inclusion decisions, retained row order, and association identifiers therefore remain deterministic even when computations finish out of order. Any worker/kernel-independent QC error aborts the scan, terminates the pool, and removes the incomplete atomic QC file.
- QC dosage and HWE extraction now visits only the final aligned analysis-sample indices. Matrix-header samples excluded by `covariate-subset` are not inspected by those calculations. The source still validates header identity/sample count before the scan.
- Replaced the former full-width per-record dosage retention/recomputation with a chunked compact decision store containing inclusion plus aligned-sample EAF. Later missingness/cache-building source passes read only included variants and requested aligned columns and reuse the completed QC decision instead of recalculating MAF/MAC/HWE. Memory cost is approximately 8 MiB per million input records plus small chunk-list overhead.
- Added the overlapping `variant_qc` profiling subphase, with input-record units, so production runs can compare sequential and automatic worker counts separately from total metadata/alignment time.
- Added deterministic tests for byte-identical sequential/four-worker QC reports and identical retained row values, excluded-header-sample dosage isolation, automatic/explicit thread selection, CLI/config parsing, parallel exception propagation through existing missing-dosage coverage, VCF and BCF reading, and the scientific end-to-end reference.
- Updated `README.md`, `TODO.md`, and `AGENTS.md` with the thread option, ordered-pipeline constraint, later-pass decision reuse, profiling semantics, and the remaining exact-trait-pattern rare/monomorphic issue. Production files changed: `QVariantMatrixSource`, `QeQTLAnalysis`, `QeQTLAnalysisConfig`, `QeQTLCommandLine`, and `QeQTLProfiler`; tests changed: `QVariantMatrixSourceTest` and `QeQTLCommandLineTest`. No commit was created.

### Verification

- Initial sandboxed `.\mvnw.cmd -DskipTests compile` could not resolve `maven-antrun-plugin:3.2.0` because network access was denied. Repeating through the approved Maven dependency/runtime path succeeded and compiled 115 Java 17 production sources.
- `.\mvnw.cmd '-Dtest=QVariantMatrixSourceTest,QeQTLCommandLineTest,QeQTLReferenceTest' test` — 19 tests ran with 0 failures/errors and 1 intentional optional-BCF skip. The sequential and four-worker QC TSVs were byte-for-byte identical, retained identifiers/values were identical, and an invalid dosage present only in an excluded header sample did not affect aligned-subset QC but remained fatal when that sample was selected.
- `.\mvnw.cmd test` — 56 tests ran with 0 failures/errors and 1 intentional optional independently generated BCF skip. CUDA/cuBLAS and JOCL/OpenCL production association/residualization integration tests ran against CPU references on the local GPU.
- The conventional shaded jar was not replaced because Java process 29912 was still running the earlier chromosome analysis and may hold that jar open. A verified temporary copy ran `.\mvnw.cmd clean package`, again compiling 115 production and 17 test sources and passing all 56 tests, then produced `target\gpu-eqtl-2.0.0-SNAPSHOT-qc-parallel-all.jar`. SHA-256: `2C4792E7894D4E27814A0E44274B1330270A065C0EEF02C092148470F9470E5C`.
- Packaged validation command: `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-qc-parallel-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-type genotype --predictor-missing mean --min-mac 2 --expression test\resources\variant-reference\expression.csv --trait-type expression --trait-missing error --output target\qc-parallel-validation.csv --variant-qc-output target\qc-parallel-validation.variants.tsv --missingness-qc-output target\qc-parallel-validation.missingness.tsv --profile-output target\qc-parallel-validation.profile.csv --validate-only` — successful without GPU initialization. It reported 15 automatic QC workers, four aligned samples, seven input/four included variants, ordered expected monomorphic/singleton/doubleton counts of one each, and a distinct `variant_qc` profile row with seven record units.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-qc-parallel-all.jar --printgpuinfo` — automatic CUDA-first discovery reported one NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, 8,589,606,912 bytes global memory, and FP64 support. The full tests also exercised its JOCL/OpenCL path.
- `git diff --check` — no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- VCF/BCF decompression and record decoding remain sequential for codec safety; only already-decoded variant QC is parallel. On fast storage or low sample counts, input decode may therefore remain the bottleneck. The worker cap is intentionally conservative and should be tuned from representative measurements rather than increased speculatively.
- Global QC is correct for the aligned cohort after selected-covariate complete-sample removal. Exact trait-pattern deletion can remove additional samples for one trait mask and make a globally retained rare predictor monomorphic within that mask; this separate scientific policy remains in `TODO.md`, and the preprocessor continues to fail instead of emitting invalid standardized results.
- Later passes rely on the compact decisions from the initial scan and assume the source file is not modified during the same run. Cache/checkpoint signatures and source validation remain unchanged; do not replace the input while an analysis is active.
- Next recommended step: after the current chromosome process finishes, rebuild the conventional `-all.jar`, then run representative `--validate-only --profile` comparisons with `--variant-qc-threads 1`, automatic, and selected explicit counts (for example 4, 8, and 16). Compare the variant QC files byte-for-byte and choose the throughput plateau before changing the automatic cap. Separately design and test the pattern-specific monomorphic/rare-variant policy.

## 2026-08-24T16:39:55.7925052-04:00 — Dependency-ordered roadmap refresh

### Baseline and goal

- Baseline commit: `730b11f7839a33a1999c3e46768f298c42c2d872` (`Improve cohort alignment and streamed variant analysis`) in `D:\projects\NIH-Project`, with the preceding uncommitted parallel-QC implementation preserved.
- Goal: replace the remaining-work list with the reviewed dependency order, explicitly including a portable CPU backend and burden/SKAT/SKAT-O development.

### Decisions and files changed

- Reorganized `TODO.md` around the dependency chain `production QC validation -> CPU backend -> indexed VCF/BCF -> burden -> SKAT -> SKAT-O -> GPU set-test acceleration`.
- Added a CPU-backend milestone covering automatic GPU-to-CPU fallback, FP64/native BLAS plus pure-Java fallback, oversubscription control, cross-backend numerical tests, and extract-and-run Windows/Linux/macOS bundles.
- Expanded rare-variant work into shared set ingestion, burden-first validation, deterministic CPU SKAT, correctly adjusted SKAT-O, bounded streaming, restartability, and GPU acceleration only after CPU-reference agreement.
- Preserved and reordered the existing missing-data, cohort covariance, interaction, forward-selection, categorical-SNP, performance, input/output, AMD/Intel validation, and release requirements.
- No production or test source changed in this documentation-only update.

### Verification and runtime observation

- `rg -n "^## |CPU analysis backend|Rare-variant|Recommended dependency" TODO.md` — confirmed all ten ordered sections and the new CPU/rare-variant milestones.
- `git diff --check -- TODO.md` — no whitespace errors; Git emitted only the existing informational LF-to-CRLF warning. Tests were not rerun because only roadmap documentation changed after the already successful 56-test code verification.
- The user's new chromosome-22 command was observed read-only. Java process 4732, started at 16:37:29, accumulated about 5.2 CPU-seconds during a five-second sample and was responsive, so it was computing rather than deadlocked. Process 29912 from 16:03:58 was also still present, and three variant-QC `.partial` files targeted the same result basename. Neither process was stopped or modified.

### Compatibility, limitations, and next step

- The active command references the conventional jar timestamped 15:10:24, not the parallel-QC alternate jar built at 16:25:34. Its last visible `Loading and validating covariate data...` line can therefore cover a later silent sequential scan and does not identify the actual current phase.
- Concurrent analyses aimed at the same output/QC basename can race when atomically publishing final files. Confirm which terminal/process should remain active before terminating anything, then rerun with a unique output basename and the parallel-QC jar if appropriate.

## 2026-08-24T17:18:56.2191531-04:00 — OpenBLAS CPU analysis backend and universal desktop packaging

### Baseline and goal

- Baseline commit: `730b11f7839a33a1999c3e46768f298c42c2d872` (`Improve cohort alignment and streamed variant analysis`) in the canonical standalone checkout at `D:\projects\NIH-Project`. The preceding uncommitted aligned-sample parallel-QC and documentation work was preserved rather than reformatted or discarded.
- Goal: add an extract-and-run CPU association backend for systems without a suitable GPU, using a maintained native BLAS where possible and a safe Java fallback, while preserving the same backend-neutral scheduler, FP64 scientific path, sample/order behavior, and active GPU runs.

### Decisions and files changed

- Chose OpenBLAS 0.3.34 through JavaCPP Presets 1.5.14 rather than oneMKL for the first CPU implementation because the same maintained binding publishes native artifacts across the requested Windows/Linux/macOS targets. `pom.xml` now includes the binding plus explicit runtime classifiers for Windows x64, Linux x64/ARM64, and macOS x64/ARM64; Android/iOS artifacts are deliberately excluded. Maven downloaded these dependencies into its normal local cache; no binary or SDK was added to source control.
- Added `gov.nih.gpu.cpu`: `CpuBackend`, one host `CpuDevice`, `CpuContext`, a native `OpenBlasMatrixEngine`, a blocked deterministic `JavaMatrixEngine`, and their small engine contract. Both engines implement FP64 and opt-in FP32 real-valued association plus `Y - (Y Q) Q^T` residualization. Host-only profiling reports compute/setup but zero transfers. Native execution failures remain fatal after an engine is selected; only `eqtl.cpu.blas=auto` may fall back at context creation, while explicit `openblas` fails instead of silently changing engines.
- Generalized selection to `--backend auto|cuda|opencl|cpu` and `--printbackendinfo`; `--gpu-backend`, `--printgpuinfo`, and `eqtl.gpu.backend` remain compatibility aliases. Automatic discovery is GPU-first, keeps CUDA's duplicate-NVIDIA suppression, filters for the requested precision, and uses CPU only if no eligible GPU remains. CPU is never mixed into an otherwise usable GPU pool, and automatic fallback prints a performance warning.
- OpenBLAS uses one process-global worker team with `max(1, logical processors - 1)` threads by default; `eqtl.cpu.threads=N` overrides it. One exclusive CPU context prevents nested BLAS products. Application `--threads` remains a separate packing/cache/result pipeline control. Automatic CPU sizing uses available JVM heap, a 256-MiB result-tile target, one full-memory worker, or at most two streamed workers rather than treating heap as VRAM.
- A packaged smoke test found and fixed a CPU auto-sizing edge case in which a tiny SNP matrix's desired two-job concurrency could reduce a block below the required 16-row tile. The concurrency floor is now one complete tile and has a regression test.
- Added CPU selection/fallback/tuning tests, native-loading and Java-engine tests, FP64/FP32 association and residualization comparisons against the independent scalar calculation, and CPU output to the manual backend benchmark. Updated compute-backend wording in `QeQTLAnalysis`/`QGpuResidualizer`, while retaining historical `gpu_*` profiler field names for output compatibility. Updated `README.md`, `TODO.md`, and `AGENTS.md` with selection, packaging, tuning, compatibility, and remaining platform-validation work.
- Material files for this goal: `pom.xml`; `src/gov/nih/gpu/cpu/*`; `GpuRuntime`, `AutoGpuBackend`, `GpuTuning`; `QeQTLAnalysis`, `QGpuResidualizer`, `QeQTLCommandLine`; CPU/backend/tuning/CLI tests and `GpuBackendBenchmark`; `README.md`, `TODO.md`, `AGENTS.md`, and this record. No commit was created.

### Downloads and verification

- `.\mvnw.cmd dependency:resolve '-DincludeGroupIds=org.bytedeco'` — successfully resolved the OpenBLAS 0.3.34-1.5.14 binding, JavaCPP 1.5.14, and all five explicit desktop native classifier jars from Maven Central.
- `.\mvnw.cmd -DskipTests compile` — compiled 121 Java 17 production sources successfully after the new binding was connected.
- `.\mvnw.cmd '-Dtest=GpuKernelIntegrationTest,AutoGpuBackendTest,GpuTuningTest,QeQTLCommandLineTest' test` — 29 tests passed with no failure/error/skip. `.\mvnw.cmd '-Dtest=CpuBackendTest' test` — two tests passed and explicitly proved that bundled Windows OpenBLAS loaded and the no-native Java engine remained selectable.
- `.\mvnw.cmd test` — 66 tests ran with 0 failures/errors and one intentional optional independently generated BCF skip. CUDA/cuBLAS, JOCL/OpenCL, OpenBLAS, and Java FP64/FP32 numerical paths ran against scalar references on this machine. After the minimum-tile fix, `.\mvnw.cmd '-Dtest=GpuTuningTest' test` passed, and the final isolated `.\mvnw.cmd clean package` reran all 66 tests with the same result.
- The conventional `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` was not cleaned or replaced because Java processes 4732 and 29912 were still active. A verified temporary copy produced `target\gpu-eqtl-2.0.0-SNAPSHOT-cpu-universal-all.jar`, and the temporary directory was removed after copying. Final jar size: 134,979,664 bytes. SHA-256: `444960AC480155E600AA9B28CF25A6BE28F6BA4717255A6286CB9AB721683504`. Archive inspection found all five expected native paths: `windows-x86_64`, `linux-x86_64`, `linux-arm64`, `macosx-x86_64`, and `macosx-arm64`.
- Packaged diagnostics with `-Deqtl.cpu.blas=openblas --backend cpu --printbackendinfo` reported OpenBLAS 0.3.34, 15 BLAS threads, 16 logical processors, and the bundled native library. Repeating with `-Deqtl.cpu.blas=java` reported the portable Java engine. Automatic diagnostics continued to report the NVIDIA GeForce RTX 2080 through CUDA first (driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, 8,589,606,912 bytes memory) and listed CPU as fallback; analysis selection does not mix them.
- Packaged forced-OpenBLAS and forced-Java full eQTL runs used the eight-sample deterministic genotype/expression/covariate fixture, FP64, accelerated residualization, block 16, and one worker. Both emitted the same six ordered identifiers. OpenBLAS versus Java maximum absolute differences were `2.776e-16` R-squared, `3.608e-16` effect, `1.110e-15` t, and `4.441e-16` log10 p; OpenBLAS versus the independent stored reference maxima were `2.776e-16`, `3.608e-16`, `1.138e-15`, and `4.441e-16` respectively.
- The corrected packaged run with both block size and application threads omitted selected block 16/one worker and emitted all six results. The manual 128-sample by 128-tile microbenchmark (nine measured calls after warmup, copies included) reported medians of 0.318 ms CUDA, 0.280 ms OpenCL, and 0.238 ms CPU; this tiny overhead-dominated shape is only a wiring smoke test, not evidence that CPU will beat GPU on production matrices.
- `git diff --check` reported no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Only Windows x64 was executed in this session. Linux x64/ARM64 and macOS x64/ARM64 resources are present but still require real-platform smoke tests before release. The universal jar is about 135 MB; classifier-specific release profiles remain an option if that becomes inconvenient, but the single jar satisfies the current extract-and-run goal.
- Pure Java fallback is intentionally slower and single-threaded inside the exclusive compute context. OpenBLAS thread control is process-global; keep it separate from variant-QC and application pipeline worker counts. oneMKL was not bundled—evaluate it only if representative x86 end-to-end profiling demonstrates a worthwhile gain over OpenBLAS.
- CPU and CUDA currently implement only the real-valued `eqtlReal` association operation. The historical categorical-SNP path remains separately disabled. Stable profile CSV phases retain `gpu_*` names even for CPU; CPU transfer phases are zero.
- The conventional jar should be rebuilt only after the active analyses exit. Next recommended work remains production VCF/QC and trait-residency profiling, including byte-identical QC comparisons across `--variant-qc-threads` settings. Then proceed to indexed VCF/BCF region access before burden/SKAT/SKAT-O. In parallel, benchmark forced CPU OpenBLAS against CUDA on representative cohort tiles and validate the universal artifact on Linux/macOS.

## 2026-08-24T19:34:53.3906384-04:00 — Optional oneMKL CPU engine and isolated evaluation packages

### Baseline and goal

- Baseline commit: `730b11f7839a33a1999c3e46768f298c42c2d872` (`Improve cohort alignment and streamed variant analysis`) in the canonical standalone checkout at `D:\projects\NIH-Project`. The preceding uncommitted VCF/QC, CPU/OpenBLAS, documentation, and test work was preserved.
- Goal: connect Intel oneMKL as a selectable CPU matrix engine, create isolated extract-and-run Windows/Linux x64 evaluation packages, verify FP64/FP32 association and residualization, and compare it with OpenBLAS without changing the portable release default.
- The active Codex workspace was again attached to the stale combined-repository copy at `D:\git\NIH-Project`. Review established that the current CPU implementation and complete modernization history live in `D:\projects\NIH-Project`; all material edits and verification in this entry were performed in that canonical standalone checkout.

### Decisions and files changed

- Added `MklMatrixEngine`, which calls JavaCPP Presets oneMKL 2026.1 directly through CBLAS DGEMM/SGEMM and uses the same `CpuMatrixEngine` contract as OpenBLAS. It supports production FP64/opt-in FP32 association and fixed-effect block residualization, reports the detected Intel product/build string, and shares the existing process-global native-BLAS thread control.
- Expanded `eqtl.cpu.blas` to `auto|mkl|openblas|java`. Explicit `mkl` and `openblas` failures are fatal. `auto` tries oneMKL when its JNI/runtime files are available, then OpenBLAS, then the portable Java engine; the normal universal package has no oneMKL native files and therefore retains OpenBLAS behavior.
- Added the JavaCPP oneMKL binding to normal compilation, but moved native CPU packaging into profiles. `cpu-openblas-universal` remains active by default and retains the existing five-platform OpenBLAS jar. `cpu-mkl-windows-x86_64` and `cpu-mkl-linux-x86_64` add only their matching oneMKL JNI/redist and OpenBLAS fallback classifiers, use isolated `target-mkl-*` directories, and create explicitly named shaded artifacts. `.gitignore` excludes those build directories.
- Added profile-gated oneMKL load, FP64/FP32 association, and residualization tests against scalar references. Updated `README.md`, `TODO.md`, and `AGENTS.md` with selection/build commands, platform scope, fallback behavior, and release constraints.
- Licensing decision: Intel permits oneMKL redistribution under the Intel Simplified Software License, but that license prohibits reverse engineering and is not automatically GPLv3-compatible. The oneMKL profiles and their generated jars are therefore documented as local evaluation artifacts only. The normal GPL distribution remains free of Intel-licensed oneMKL native binaries until the copyright holder approves a narrow GPLv3 section-7 linking exception or legal review confirms a different compliant arrangement. No license exception was inferred or added in this session.
- Material files for this increment: `.gitignore`, `pom.xml`, `src/gov/nih/gpu/cpu/CpuBackend.java`, `CpuDevice.java`, new `MklMatrixEngine.java`, `test/gov/nih/gpu/cpu/CpuBackendTest.java`, `test/gov/nih/gpu/GpuKernelIntegrationTest.java`, `README.md`, `TODO.md`, `AGENTS.md`, and this record. No commit was created.

### Downloads and verification

- Maven Central metadata identified `org.bytedeco:mkl-platform-redist:2026.1-1.5.14` as the JavaCPP-1.5.14-aligned release. `\.\mvnw.cmd dependency:get '-Dartifact=org.bytedeco:mkl-platform-redist:2026.1-1.5.14' '-Dtransitive=true'` downloaded the official Windows x64 and Linux x64 JNI/redist artifacts into the local Maven cache. The native redist jars are 158,318,815 bytes for Windows and 188,937,201 bytes for Linux; no SDK or binary was added to source control.
- `\.\mvnw.cmd -Pcpu-mkl-windows-x86_64 '-Dtest=CpuBackendTest,GpuKernelIntegrationTest' test` — 21 tests passed with no failures, errors, or skips. Forced oneMKL executed FP64, FP32, and residualization reference checks.
- `\.\mvnw.cmd test` using the default universal OpenBLAS profile — 70 tests ran with 0 failures/errors and 5 intentional skips: four profile-gated oneMKL checks and the optional independently generated BCF fixture. CUDA, JOCL/OpenCL, OpenBLAS, Java, and the common application tests remained successful.
- `\.\mvnw.cmd -Pcpu-mkl-windows-x86_64 package` — full Java 17 build succeeded; 70 tests ran with 0 failures/errors and only the optional external BCF test skipped. Maven Shade produced `target-mkl-windows-x86_64\gpu-eqtl-2.0.0-SNAPSHOT-mkl-windows-x86_64-all.jar`, 249,013,162 bytes, SHA-256 `21732F653A178340152922B4CE9ACDF4BFB1CC84800EB848B1CDF3DF99E33518`. Archive inspection found the project GPL at `META-INF/LICENSE` plus `jnimkl_rt.dll`, `mkl_rt.3.dll`, and `libiomp5md.dll` under the expected JavaCPP Windows path. Pre-existing shaded-resource overlap warnings remain.
- Packaged forced-MKL diagnostics reported Intel oneAPI Math Kernel Library 2026.1 Product Build 20260612, eight BLAS threads, 16 logical host processors, Windows 11, FP64 support, and bundled-native loading. Forced OpenBLAS in the same jar reported OpenBLAS 0.3.34 and 15 BLAS threads. Automatic backend diagnostics retained GPU-first selection and reported an NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, and 8,589,606,912 bytes memory; CPU remained a fallback rather than being mixed into analysis.
- Packaged forced-oneMKL and forced-OpenBLAS full-memory eQTL commands used the eight-sample deterministic genotype/expression/covariate fixture, FP64, backend residualization, block 16, and one application worker. Both wrote six identically ordered associations; maximum displayed differences were zero for R-squared, effect, t, and log10 p.
- Transfer-inclusive `GpuBackendBenchmark 2048 512` measurements (three warmups and nine measured calls) gave oneMKL CPU median/best 23.652/21.503 ms and forced OpenBLAS CPU 20.986/19.558 ms. The adjacent GPU measurements varied between runs: CUDA medians 11.465/10.639 ms and OpenCL 9.581/8.200 ms. This single shape is a first-pass wiring comparison, not a full-cohort performance conclusion, but it provides no reason to replace OpenBLAS as the default on this host.
- `git diff --check` reported no whitespace errors; Git emitted only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Only the Windows x64 oneMKL profile was executed. The Linux x64 profile resolves the corresponding Maven classifiers by construction but still needs a real Linux build, native-load smoke test, full numerical run, artifact inspection, and checksum. oneMKL 2026.1 has no macOS runtime; macOS remains OpenBLAS/Java.
- oneMKL here is CPU CBLAS. It does not use Intel Iris Xe or another GPU; Intel GPU acceleration would require a separately designed SYCL/OpenMP-offload or Level Zero path and corresponding scientific tests.
- The oneMKL jar is a local evaluation artifact, not an authorized project release. Before any publication, resolve the GPLv3/Intel-license combination and reproduce all Intel-required notices. The default universal OpenBLAS jar remains the only intended distributable CPU package.
- The conventional jar was not rebuilt or replaced because Java processes 4732 and 5776 were still running the earlier `gpu-eqtl-2.0.0-SNAPSHOT-qc-parallel-all.jar` production command. All oneMKL output used its isolated build directory.
- Next recommended oneMKL work is a representative forced-MKL versus forced-OpenBLAS end-to-end benchmark after the current analyses finish, recording total time, preprocessing/cache time, matrix compute, output time, peak heap, identifiers, and numeric deltas. Unless oneMKL wins materially on those workloads, retain it as an optional evaluation backend. The main roadmap remains indexed VCF/BCF access followed by deterministic burden/SKAT/SKAT-O references and bounded acceleration.

## 2026-08-24T20:53:19.2888183-04:00 — oneMKL GPLv3 exception and resumable aligned-sample variant QC

### Baseline and goal

- Baseline commit: `730b11f7839a33a1999c3e46768f298c42c2d872` (`Improve cohort alignment and streamed variant analysis`) in the canonical standalone checkout at `D:\projects\NIH-Project`. All preceding uncommitted VCF/QC, CPU/OpenBLAS/oneMKL, documentation, and test changes were preserved.
- Goal: authorize and package a legally explicit GPLv3/oneMKL combination at the copyright holder's request, and make the long aligned-sample VCF/BCF EAF/MAF/MAC/HWE scan resumable without changing source order, sample order, filtering decisions, or downstream identifiers.

### Decisions and files changed

- Kept `LICENSE` as the unmodified GPLv3 text and added `LICENSE_EXCEPTION`, a narrow GPLv3 Section 7 permission to convey Roby Joehanes-owned portions when linked or combined with Intel oneMKL. The permission does not apply to third-party project code, relicense Intel components, or permit modification of Intel binaries. `pom.xml` identifies GPLv3 with the exception and packages the GPL and exception under `META-INF`.
- Retrieved Intel's exact 2026.1.0-236 licensing files from the official Intel oneAPI package rather than relying on the JavaCPP jars, which contain the oneMKL/OpenMP binaries but omit Intel's notice files. Added the byte-identical Intel Simplified Software License, oneMKL incorporated-code notices, and Intel OpenMP notices under `THIRD_PARTY_LICENSES`; Maven packages all three under `META-INF/licenses`. `.gitattributes` disables text conversion and trailing-space diagnostics only for these verbatim vendor notices so Git cannot alter their official bytes. The 126,239,806-byte Intel `.deb` was used only in the temporary directory and was not added to the repository. Updated `MklMatrixEngine`, `README.md`, `TODO.md`, and `AGENTS.md` to replace the former local-evaluation-only restriction with the exception and mandatory notice rules.
- Added `QVariantQcCheckpoint`. A VCF/BCF QC scan now commits complete source-ordered batches atomically under `<variant-qc-output>.checkpoint/<signature>/`; `--variant-qc-checkpoint DIR` and INI `variant_qc_checkpoint` override the root. A per-signature file lock rejects simultaneous writers. Final QC TSV publication remains atomic and is assembled from durable ordered parts.
- The checkpoint signature covers the normalized genotype path, byte size, modification time, requested format, resolved DS/GT field, missing/multiallelic policies, MAF/MAC thresholds, complete source header sample IDs, and final aligned sample indices/order. A changed scientific input creates a different signature directory rather than reusing stale EAF/MAF/MAC/HWE state.
- On an incomplete restart, prior QC rows rebuild duplicate detection, summary counts, inclusion/EAF decisions, and source order. The sequential VCF/BCF reader decodes and verifies the completed canonical-ID prefix but does not redo its per-sample dosage counts or exact HWE. A completed matching checkpoint is loaded directly, can recreate a missing QC TSV, and avoids a full variant-record scan. Console output reports the active signature directory and number of reused records. Resume-aware progress rates exclude records completed before the current process.
- Added a deterministic interruption test that commits every two records, injects failure after record four, resumes records five through seven, compares the QC TSV byte-for-byte with a clean parallel scan, deletes/recreates the report from a completed checkpoint, and verifies identical summaries and retained decisions. CLI parsing covers the new path option. Updated the VCF documentation and roadmap with retention, locking, invalidation, and sequential-prefix limitations.
- Material files for this increment: `.gitattributes`; `LICENSE_EXCEPTION`; `THIRD_PARTY_LICENSES/*`; `pom.xml`; `src/gov/nih/gpu/cpu/MklMatrixEngine.java`; new `src/gov/nih/eqtl/io/QVariantQcCheckpoint.java`; `QVariantMatrixSource`, `QeQTLAnalysis`, `QeQTLAnalysisConfig`, and `QeQTLCommandLine`; the corresponding variant/CLI tests; `README.md`, `TODO.md`, `AGENTS.md`, and this record. No commit was created during implementation.

### Verification and artifacts

- `.\mvnw.cmd "-Dtest=QVariantMatrixSourceTest,QeQTLCommandLineTest" test` — final focused run compiled 123 Java 17 production sources and 18 test sources; 18 tests ran with 0 failures/errors and one intentional optional external-BCF skip. The injected four-record interruption resumed correctly, and completed seven-record state was reused.
- `.\mvnw.cmd test` — final full default-profile run executed 71 tests with 0 failures/errors and five intentional skips (four oneMKL-profile-gated checks and the optional external BCF). CUDA/cuBLAS, JOCL/OpenCL, OpenBLAS, portable Java, deterministic application, checkpoint, and VCF tests remained successful. This host's hardware path is the previously recorded NVIDIA GeForce RTX 2080 with CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, and 8,589,606,912 bytes VRAM.
- `.\mvnw.cmd -Pcpu-mkl-windows-x86_64 package` — before the final notice-only/lock-cleanup edits, the full profile ran 71 tests with 0 failures/errors and only the optional external BCF skip; forced oneMKL FP64/FP32 association and residualization, CUDA, and OpenCL hardware checks executed. The final `.\mvnw.cmd -Pcpu-mkl-windows-x86_64 -DskipTests package` recompiled the tested source and produced the isolated Windows artifact without touching the conventional jar.
- Final artifact: `target-mkl-windows-x86_64\gpu-eqtl-2.0.0-SNAPSHOT-mkl-windows-x86_64-all.jar`, 249,051,330 bytes, SHA-256 `668AB339710F01E19581E01C142DB4DB1C2609EF974D91F0E76164D75C79B09C`. Archive inspection confirmed `QVariantQcCheckpoint.class` and exact repository matches for `META-INF/LICENSE` (35,149 bytes), `META-INF/LICENSE_EXCEPTION` (727), `META-INF/licenses/INTEL-ONEMKL.txt` (4,246), `INTEL-ONEMKL-THIRD-PARTY.txt` (64,052), and `INTEL-OPENMP-THIRD-PARTY.txt` (31,028).
- The final packaged `--validate-only` VCF command was run twice with an explicit checkpoint root. The first pass selected DS, scanned seven aligned-sample records, retained four variants, and reported the expected one monomorphic/one singleton/one doubleton counts. The second pass printed `Reusing completed variant-QC checkpoint with 7 records`, preserved all counts/sample order/degrees of freedom, and completed successfully without compute-backend initialization.
- `git diff --check` reported no whitespace errors; Git emitted only the repository's informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Incomplete `.vcf.gz` and BCF resumes must still decode from the beginning to reach and verify the completed prefix because the current HTSJDK source is sequential. This avoids the expensive aligned-sample statistics/HWE work but not compressed decoding. True seek-to-resume requires the planned tabix/CSI indexed VCF/BCF reader.
- Durable batches contain 10,000 records by default, so an abrupt process or machine failure can require recomputing at most the last uncommitted 9,999 QC records. The internal batch-size property exists for deterministic tests; it is not a documented analysis setting. Checkpoint roots intentionally retain separate signature subdirectories and are not automatically deleted.
- File identity uses normalized path, size, and modification time rather than hashing an entire very large genotype file. The resumed incomplete prefix is additionally checked record-by-record by canonical ID; a completed checkpoint trusts the conservative metadata signature and already validated header/sample list.
- Existing runs loaded from older jars do not gain QC resume support. The conventional `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` was not repackaged; the verified deliverable is the isolated oneMKL jar above. Linux oneMKL packaging still needs a real Linux native-load, numerical, notice, and checksum smoke test.
- Next recommended step: exercise interruption/restart on a copied representative chromosome VCF.gz with a unique output basename, compare the final QC TSV byte-for-byte with an uninterrupted scan, and record initial, prefix-decode, resumed-QC, and completed-checkpoint-reuse timings plus checkpoint disk size. Then implement indexed VCF/BCF access so an interrupted run can seek near its durable boundary instead of decoding the prefix.

## 2026-08-25T00:18:08.8560596-04:00 — Indexed regions and exact pattern-aware variant statistics

### Baseline and goal

- Baseline commit: `c4af0422bb1ea79e050eaa408eac73c9768ea09d` (`Add CPU acceleration and resumable variant QC`) in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: add indexed VCF/BCF interval access and deterministic region/set definitions, make indexed QC resume seek to its durable boundary, persist aligned sufficient statistics, and make expression-missingness patterns use exact genotype means/frequencies and safe pattern-specific variant decisions without repeatedly decoding compressed genotype input.

### Decisions and files changed

- Added `QGenomicRegions`. Repeatable `--region [SET=]CHROM:START-END` values are one-based inclusive. `--regions-file` accepts tab-separated `CHROM/START/END` or `SET_ID/CHROM/START/END`, with `--region-coordinates one-based|bed`. Exact contigs are preferred; only unambiguous `chr` and mitochondrial aliases are normalized. Query intervals are sorted by source-contig order and merged, while overlapping set membership remains many-to-many and empty sets are reported.
- Added true interval queries for BGZF VCF/tabix `.tbi` and VCF/BCF HTSJDK Tribble `.idx` through the project BCF 2.2-compatible codec. `--variant-index` overrides neighboring index detection. Region requests fail before analysis when no usable index is present. Source and index size/mtime are captured, included in checkpoint/prepared-cache state, and rechecked before every later pass. Query overlap cannot duplicate a spanning variant, but genuine duplicate canonical records remain fatal.
- Indexed incomplete QC now seeks to the last durable `CHROM:POS:REF:ALT`, verifies that exact boundary, and continues after it. Sequential input retains full-prefix verification. Variant-QC checkpoint/signature format advanced to version 3 and covers index identity, normalized region definitions, and frequency scope.
- Expanded durable aligned QC rows with dosage sum of squares and diploid hom-ref/heterozygous/hom-alt counts, in addition to called/missing counts, dosage sum, EAF/MAF/MAC, and HWE. Region memberships, frequency scope, and the aligned-cohort frequency reason are also reported.
- Added `--frequency-scope aligned|pattern`. `aligned` remains the default: MAF/MAC/HWE are cohort-level QC after genotype-expression-covariate alignment, while every exact trait mask still recalculates the imputation mean, variance, EAF/MAF/MAC, and monomorphic status. Explicit `pattern` defers MAF/MAC filtering to each trait mask and reports the changed variant counts.
- Added `QRawMatrixCache`, a persistent aligned pre-residualization FP64 dosage cache with sample order, row IDs, signatures, and per-row checksums. VCF/BCF exact-pattern analysis decodes the retained aligned source into this cache once, then performs later mask passes over binary rows instead of reparsing VCF/BCF.
- Added `QPatternVariantSource`, a checksum-protected per-mask cache of called/missing counts, dosage sum/squared sum, mean, EAF/MAF/MAC, selection flags, and reasons. Pattern-specific mean filling uses the actual observed trait-mask samples. A variant that becomes constant is skipped only for that mask. Zero fill uses variance after zero replacement rather than observed-only variance.
- Pattern-specific statistics plus residualized/standardized predictor and trait caches now persist under `--cache-dir`; reruns reuse them. `<output>.pattern-variant-qc.tsv` records effective N and per-pattern variant counts. Association groups themselves still use temporary checkpoints and are not yet partially resumable.
- Updated the command-line/configuration layer, `README.md`, `TODO.md`, and `AGENTS.md`. Added deterministic indexed VCF/tabix, indexed BCF 2.2/Tribble, boundary-resume, index-mutation, BED/alias/overlap/empty-set, raw-cache corruption, pattern mean/zero/MAF/MAC, and VCF-to-pattern-cache tests. Added `test/resources/variant-pattern-reference/expression.csv` for the packaged exact-pattern run.
- Material production files: `QeQTLAnalysis`, `QeQTLAnalysisConfig`, `QeQTLCommandLine`, `QPolicyMatrixSource`, `QVariantMatrixSource`, and new `QGenomicRegions`, `QRawMatrixCache`, and `QPatternVariantSource`. The two still-running production processes use the separately named `gpu-eqtl-2.0.0-SNAPSHOT-qc-parallel-all.jar`; no clean or deletion touched that artifact.

### Verification and artifacts

- `.\mvnw.cmd "-Dtest=QGenomicRegionsTest,QPatternVariantSourceTest,QRawMatrixCacheTest,QVariantMatrixSourceTest,QeQTLCommandLineTest,QMissingnessPolicyTest" test` — 32 tests ran with 0 failures/errors and one intentional optional external-BCF skip.
- `.\mvnw.cmd test` — 82 tests ran with 0 failures/errors and five intentional skips: four oneMKL-profile-gated checks and the optional independently generated external BCF. CUDA/cuBLAS, JOCL/OpenCL, bundled OpenBLAS, portable Java, reference association/residualization, indexed input, checkpoint, and cache-corruption checks passed.
- `.\mvnw.cmd -DskipTests package` — packaged successfully without `clean`, preserving the active production jar. Final conventional artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,723,148 bytes, SHA-256 `1D05692835E11EEDB0282053AFEC4191919F7E580DE67156F09D8484FB9553F7`.
- Packaged first-pass exact-pattern command: `java '-Deqtl.cpu.blas=java' --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cpu --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-missing mean --expression test\resources\variant-pattern-reference\expression.csv --trait-missing pattern --sample-alignment strict --min-mac 2 --frequency-scope pattern --variant-qc-output target\variant-pattern-e2e-final\variants.tsv --variant-qc-checkpoint target\variant-pattern-e2e-final\variant-checkpoint --missingness-qc-output target\variant-pattern-e2e-final\missingness.tsv --threshold none 0 --block-size 16 --threads 1 --genotype-block-rows 16 --expression-block-rows 16 --trait-cache memory --precision fp64 --residualization cpu --cache-dir target\variant-pattern-e2e-final\cache --output target\variant-pattern-e2e-final\results.csv` — successful. It retained five aligned candidates, used two trait patterns (N=4/DF=2 and N=3/DF=1), tested four and three variants respectively, and skipped two variants that became monomorphic after removing S1.
- Repeating the packaged command with `results-reuse.csv`, reuse-specific QC paths, and the same checkpoint/cache roots reused the completed variant QC, aligned raw cache, both pattern-statistics caches, and all four prepared caches. The association files were byte-identical with SHA-256 `31A6577B5396CEF1921E8C6B911EF11E96C9E9B3C4AD0A84964CB1CB6418F972`.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` — successful. Hardware/runtime: NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, 8,589,606,912 bytes VRAM, FP64 available; Windows 11 amd64 host with 16 logical processors; bundled OpenBLAS 0.3.34 detected with 15 BLAS threads.
- `git diff --check` — no whitespace errors; only the repository's informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Standard bcftools variant CSI (`.csi`) is detected as a candidate but HTSJDK 5 does not decode it through the pure-Java Tribble variant reader. Current verified random access is tabix `.tbi` for BGZF VCF and Tribble `.idx` for VCF/BCF, including BCF 2.2. Unsupported CSI fails rather than silently scanning the whole file. A pure-Java CSI decoder or carefully packaged cross-platform htslib bridge remains follow-up work.
- The indexed BCF fixture is deterministically generated with HTSJDK and version-adjusted to BCF 2.2. The optional independent BCF test was not configured, and representative tabix/BCF queries still need byte-for-byte/ID comparison with bcftools on production copies.
- The aligned raw dosage cache intentionally uses FP64 and therefore needs about `8 * aligned samples * retained variants` bytes plus IDs/checksums. It bounds RAM and avoids repeated VCF decoding, but a first exact-pattern run still prepares one residualized predictor cache per distinct trait mask. A genotype-block-outer multi-pattern schedule remains a performance option only after exact output-order/checkpoint tests.
- Pattern statistics and prepared matrices are reusable, but a process interrupted midway through pattern association does not yet resume completed pattern result groups. Global analysis checkpoint behavior is unchanged outside exact-pattern mode.
- HWE remains aligned-cohort QC. Pattern-specific HWE is not silently substituted. `--frequency-scope pattern` changes tested variant sets across traits and must remain explicit.
- Next recommended step: validate `.tbi` region selection and indexed interruption/restart on a copied production chromosome against bcftools, then either add standard CSI support or begin the deterministic continuous-trait burden-test foundation using the delivered region/set and pattern-aware QC infrastructure.

## 2026-08-25T11:15:32.8321896-04:00 — Association progress, default MAC filtering, and reusable VCF/BCF preprocessing

### Baseline and goal

- Baseline commit: `7cabe457f523e6fb96e4fe5914ee40e10871c5d1` (`Add indexed regions and pattern-aware variant QC`) in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: add visible progress throughout association computation, make VCF/BCF MAF/MAC thresholds explicit with a conservative default of MAC at least 20, and allow aligned VCF/BCF QC/dosage preprocessing to finish without running an association so its durable state can be reused later.

### Decisions and files changed

- Added `QAnalysisProgress`, a thread-safe comparison counter shared by the full-memory, streamed, CPU/GPU, multi-context, resumed-checkpoint, and exact trait-pattern association paths. It counts active variant-trait comparisons rather than padded cells, prints start/completion markers plus percentage, elapsed time, current-run throughput, and ETA, and uses a daemon heartbeat so a still-running long tile reports approximately every 15 seconds. Reused checkpoint comparisons contribute to percentage but not current-run rate. The former unconditional full-memory block-coordinate output is now debug-only.
- Changed the application-level VCF/BCF default to `--min-mac 20`; `--min-maf` remains zero/disabled. Positive MAF and MAC thresholds are conjunctive, and explicit zero disables either one. The effective filter and default origin are printed before QC and remain covered by the existing variant-QC/cache signatures. An all-filtered input now reports its exact aligned-sample MAF/MAC/frequency-scope settings and directs the user to the QC report or explicit lower thresholds. CSV input is unaffected because variant frequency filtering applies only to VCF/BCF.
- Added `--preprocess-only` and INI `preprocess_only = true`. This mode requires VCF/BCF plus the trait matrix and any selected covariates, performs final sample alignment, variant QC/filtering, missingness reporting, covariate model/rank and degrees-of-freedom validation, creates or reuses the aligned checksummed FP64 raw-dosage cache, skips compute-backend initialization, and does not require an association output path.
- Moved matching raw-cache discovery ahead of predictor missingness scanning. A later ordinary run with the same source/index, field, regions, inclusion policy, aligned samples, MAF/MAC values, and `--cache-dir` reads the preprocessed binary rows for missingness and prepared-cache creation instead of decoding the compressed source again. Mean and zero replacement may share the raw cache because filling occurs in separately signed later state. Reusing the same explicit variant-QC output/checkpoint paths also avoids repeating the initial QC scan.
- Updated `README.md`, `TODO.md`, `AGENTS.md`, and command-line help with defaults, disable/combination rules, preprocessing semantics, a complete command, cache locations/reuse requirements, and progress behavior. Corrected help to identify the verified explicit indexes as `.tbi`/`.idx`, not unsupported standard CSI.
- Production files changed: `QeQTLAnalysis`, `QeQTLAnalysisConfig`, `QeQTLCommandLine`, `QeQTLSNPJobReal`, `QeQTLStreamedJobReal`, `QRawMatrixCache`, `QVariantMatrixSource`, and new `QAnalysisProgress`. Tests changed/added: `QeQTLCommandLineTest`, `QRawMatrixCacheTest`, new `QAnalysisProgressTest`, and new `QVariantPreprocessOnlyTest`. No commit was created in this increment.

### Verification and artifacts

- `.\mvnw.cmd "-Dtest=QAnalysisProgressTest,QVariantPreprocessOnlyTest,QeQTLCommandLineTest" test` — 12 focused tests passed with no failures/errors/skips after the scheduled heartbeat was added. They cover rate limiting, a stalled-tile heartbeat, exact completion, resumed-rate accounting, CLI/default behavior, default-MAC all-filtered diagnostics, first-pass preprocessing, byte-stable cache reuse, and no association output.
- `.\mvnw.cmd test` — final full run executed 88 tests with 0 failures/errors and five intentional skips: four profile/platform-gated oneMKL checks and the optional independently generated external BCF fixture. CUDA/cuBLAS, JOCL/OpenCL where available, bundled OpenBLAS, portable Java, deterministic association/residualization, indexed input, QC/checkpoints, raw/pattern/prepared caches, and corruption checks remained successful.
- `.\mvnw.cmd -DskipTests package` — final Java 17 shaded build succeeded without `clean`. Artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,729,617 bytes, SHA-256 `536BDDEA4C4C177313127AAD8E2F547FF2C021EB8DD383F9CB705B972004A42C`. Existing Maven Shade duplicate-license/module warnings were unchanged.
- Packaged preprocessing command: `java '-Deqtl.cpu.blas=java' --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-missing mean --expression test\resources\variant-reference\expression.csv --trait-missing error --sample-alignment strict --min-mac 0 --variant-qc-output target\preprocess-progress-e2e-20260825\variants.tsv --variant-qc-checkpoint target\preprocess-progress-e2e-20260825\variant-checkpoint --missingness-qc-output target\preprocess-progress-e2e-20260825\missingness-final.tsv --genotype-block-rows 2 --cache-dir target\preprocess-progress-e2e-20260825\cache --preprocess-only` — completed without backend initialization or association, reused all seven QC records and the aligned raw cache on the final repetition, and reported five retained variants/two traits/four aligned samples.
- Packaged reuse/association command: `java '-Deqtl.cpu.blas=java' --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cpu --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-missing mean --expression test\resources\variant-reference\expression.csv --trait-missing error --sample-alignment strict --min-mac 0 --variant-qc-output target\preprocess-progress-e2e-20260825\variants.tsv --variant-qc-checkpoint target\preprocess-progress-e2e-20260825\variant-checkpoint --missingness-qc-output target\preprocess-progress-e2e-20260825\missingness-association-final.tsv --genotype-block-rows 16 --expression-block-rows 16 --block-size 16 --threads 1 --trait-cache memory --precision fp64 --residualization cpu --cache-dir target\preprocess-progress-e2e-20260825\cache --threshold none 0 --output target\preprocess-progress-e2e-20260825\results-final.csv` — reused the completed QC, aligned raw predictor cache, prepared genotype cache, prepared expression cache, and in-memory trait residency; printed `Association started: 10` and `Association complete: 10/10`; emitted the expected ten associations plus header.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` — reported NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, and 8,589,606,912 bytes VRAM; Windows 11 amd64 host with 16 logical processors; bundled OpenBLAS 0.3.34 and 15 BLAS threads. `git diff --check` reported no whitespace errors before this append, with only informational LF-to-CRLF warnings.

### Known limitations, compatibility, and next step

- MAC 20 is an intentional VCF/BCF default change and can exclude every variant in very small fixtures, restricted regions, or rare-variant studies. Such runs now fail clearly and still materialize the QC report. Use `--min-mac 0` to disable MAC filtering or specify another threshold; burden/SKAT work will need explicit rare-variant/set-specific policies rather than relying blindly on this common-variant default.
- Preprocessing is deliberately alignment-specific, not a cohort-independent VCF conversion. Changing the analysis sample set, source/index metadata, selected field, regions, or variant-inclusion/frequency settings creates a distinct signature/cache. The FP64 raw cache needs approximately `8 * aligned samples * retained variants` bytes plus row identifiers/checksums.
- The heartbeat can report that the completed count is unchanged while a long tile is executing. Rate and ETA are based on completed comparisons, so early estimates may be absent or unstable; they become representative only after multiple tiles complete.
- The raw cache avoids post-QC VCF/BCF decoding, but avoiding the QC scan itself also requires reuse of the matching QC checkpoint root. Standard CSI remains unsupported as recorded previously.
- Next recommended step: run `--preprocess-only` on a copied production chromosome with the intended aligned cohort and default MAC 20, record QC/cache size and wall time, then run the full command with the same cache and QC checkpoint paths. Confirm the terminal heartbeat, absence of another compressed-source pass, retained-variant count, identifiers, and output against the current production command before making this the routine workflow.

## 2026-08-25T11:25:15.7468082-04:00 — Exhaustive command-line documentation

### Baseline and goal

- Baseline commit: `7cabe457f523e6fb96e4fe5914ee40e10871c5d1` (`Add indexed regions and pattern-aware variant QC`) in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: document the newly added association progress, default VCF/BCF MAC filtering, and preprocessing/reuse workflow, provide copyable examples, and publish a complete reference for every accepted command-line option before committing and pushing the accumulated verified change set to both configured remotes.

### Decisions and files changed

- Added `docs/COMMAND_LINE_REFERENCE.md` as the exhaustive user reference. It covers all 72 parser-recognized option spellings, aliases, accepted values, application defaults, legacy INI keys, limitations, and the three supported invocation styles.
- Added complete PowerShell examples for a bounded-RAM VCF.gz association with ID bridging, two-stage VCF/BCF preprocessing and cache reuse, indexed regions/sets, missingness-only inspection of generic QTL matrices, and forced CPU execution. Added the supported JVM properties for backend, CPU BLAS engine, and native BLAS threads.
- Documented result-affecting defaults together: no association reporting threshold unless requested, VCF/BCF MAC at least 20 with MAF disabled, conjunctive positive MAF/MAC filters, mean predictor imputation, exact pattern-wise trait deletion, complete-selected-covariate sample removal, strict ID alignment, FP64, and automatic backend/sizing/thread choices.
- Explicitly distinguished accepted legacy compatibility options from implemented modern functionality: an INI with no genotype format retains the TPED default, `--family` is TPED-only, `--pedigree` is not consumed by a mixed model, `--random-covariates` warns and is ignored, and the full categorical-SNP model remains disabled.
- Linked the new reference prominently from the `README.md` run section and corrected the README wording so built-in `--help` is described as a short summary rather than the complete list.

### Verification

- Parser/reference audit over `QeQTLCommandLine.java` — 72 unique accepted option spellings found; all 72 occur in `docs/COMMAND_LINE_REFERENCE.md`; no missing option.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings for existing files.
- Initial sandboxed `\.\mvnw.cmd test` could not download `maven-antrun-plugin:3.2.0` because network access was denied; this was an environment/dependency-resolution failure before compilation or tests.
- `\.\mvnw.cmd test` with approved Maven Central access — successful; 88 tests ran with 0 failures, 0 errors, and five intentional skips (four profile/platform-gated oneMKL checks and the optional independently generated external BCF fixture).
- Hardware exercised by the full suite remained the NVIDIA GeForce RTX 2080 CUDA/OpenCL host and bundled OpenBLAS environment recorded in the immediately preceding session entry; no new hardware contract or numerical tolerance was introduced by this documentation increment.

### Known limitations, compatibility, and next step

- The Markdown reference is exhaustive for the current parser but is maintained manually. Future command-line additions should update both the short built-in help and `docs/COMMAND_LINE_REFERENCE.md`, then repeat the parser/reference audit.
- The reference intentionally documents current behavior rather than implying roadmap features: random effects/familial mixed models, burden/SKAT/SKAT-O, SNP interactions, and forward selection remain pending.
- Next: push the verified baseline to both remotes, exercise preprocessing/reuse on a copied production chromosome as recommended above, and use its timings and retained-variant counts to decide whether CSI support or the continuous-trait rare-variant foundation should be implemented next.

## 2026-08-25T12:53:00-04:00 — Parallel VCF genotype decoding for variant QC

### Baseline and goal

- Baseline commit: `4524cbdd7161de9c6251411b62e4ccad9b8d666b` (`Add reusable variant preprocessing and CLI reference`) in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: explain and remove the bottleneck behind an automatic 15-worker VCF QC run using only about one third of a 16-logical-processor host, without interrupting the user's production process or changing aligned-sample QC results/order.

### Decisions and files changed

- Profiling the handoff showed that `VCFIteratorBuilder` caused HTSJDK to expand every sample genotype on the single traversal thread before the record reached the worker pool, especially when VCF header sample IDs were not already sorted. The workers therefore performed only the lighter dosage/count/HWE loop. The existing compressed-input buffer was already 16 MiB; buffer size alone could not remove this serial decode stage.
- Replaced unrestricted sequential VCF iteration in `QVariantMatrixSource` with a raw-line cursor that parses fixed VCF fields on the traversal thread and retains FORMAT/sample columns in a `LazyGenotypesContext`. Each QC worker has its own header-configured `VCFCodec`, so mutable HTSJDK parser arrays/caches are never shared and genotype expansion plus aligned-sample statistics can proceed across variants. Both sorted and unsorted VCF sample headers take this path. Later sequential matrix passes still decode a record before advancing the owning cursor.
- Preserved BCF's reader-thread genotype expansion because its lazy decoder shares codec/builders. Indexed VCF interval cursors and an indexed boundary-resume still use HTSJDK's indexed reader; these paths can remain decode-bound when an unsorted header forces eager expansion. Fresh unrestricted VCF/VCF.gz preprocessing uses the new parallel path even when a neighboring index exists.
- Replaced sample-name map lookups in the hot aligned-sample loops with the already-validated header offsets. Increased the bounded queue from two to four in-flight records per worker (60 records for 15 workers), retaining source-ordered result consumption and bounded memory.
- Added end-of-QC diagnostics for worker count, maximum in-flight records, worker evaluations, measured average active worker tasks, and records whose VCF genotypes were expanded by workers. `parallel VCF genotype decodes` should approach the number of biallelic genotype-bearing records on a new unrestricted VCF scan; a zero count indicates a BCF/indexed/eager or checkpoint-reuse path.
- Added deterministic tests requiring worker decoding for both unsorted and sorted VCF headers, sequential/parallel summary and byte-identical QC output, and reordered dosage equality. Updated `README.md`, `docs/COMMAND_LINE_REFERENCE.md`, `TODO.md`, and the VCF/BCF concurrency rule in `AGENTS.md`.
- The already-running production Java process was observed read-only and was not stopped, modified, or benchmarked against. It continues with the old loaded classes; this change applies only after starting the newly packaged JAR.

### Verification and artifact

- `.\mvnw.cmd "-Dtest=QVariantMatrixSourceTest" test` — 17 tests ran with 0 failures/errors and one intentional optional external-BCF skip. Sequential and four-worker VCF QC remained byte-identical, reordered rows/dosages were exact, and unsorted/sorted headers both reported worker-side genotype decodes.
- `.\mvnw.cmd test` — 89 tests ran with 0 failures/errors and five intentional skips (four platform/profile-gated oneMKL checks and the optional independently generated BCF fixture). CUDA/cuBLAS and JOCL/OpenCL integration tests exercised the same Windows 11 / NVIDIA GeForce RTX 2080 host previously recorded (CUDA driver API 13.3, CUDA runtime 12.6, compute capability 7.5, FP64); no separate driver probe was rerun. CPU/OpenBLAS and deterministic reference tests also passed.
- `.\mvnw.cmd -DskipTests package` — Java 17 shaded packaging succeeded. Artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,734,193 bytes, SHA-256 `4E2E68BD0413268D09E4FBE64B444673055B05A05A199F51099B7E5638E7CBB5`. Existing Shade duplicate-resource warnings were unchanged.
- Packaged smoke test: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-missing mean --expression test\resources\variant-reference\expression.csv --trait-missing error --sample-alignment strict --min-mac 0 --variant-qc-threads 4 --variant-qc-output target\qc-parallel-e2e-20260825\variants.tsv --variant-qc-checkpoint target\qc-parallel-e2e-20260825\variant-checkpoint --missingness-qc-output target\qc-parallel-e2e-20260825\missingness.tsv --cache-dir target\qc-parallel-e2e-20260825\cache --preprocess-only` — successful; seven records evaluated in source order, six biallelic genotype records decoded by workers, five variants retained, and the aligned raw cache completed.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- Compressed VCF traversal/decompression and source-ordered checkpoint/QC writing remain sequential. The intended result is substantially more than the prior approximately four effective workers, but 15 fully busy cores are not guaranteed; the final `average active worker tasks` value on a representative run identifies the remaining upstream/downstream ceiling.
- A completed matching QC checkpoint bypasses the scan, so it will not demonstrate the change. An incomplete indexed resume and explicit indexed region query can still show zero worker-side VCF decodes because HTSJDK owns those indexed record cursors. BCF also remains reader-decode-bound. These are reported rather than silently presented as fully parallel.
- Do not restart the expensive in-progress run solely for this optimization. Use its completed checkpoint/cache. On the next new unrestricted VCF/VCF.gz preprocessing input, keep `--variant-qc-threads 0` or the reported 15 and inspect `average active worker tasks` plus `parallel VCF genotype decodes`. Only compare 4/8/15 workers on a short representative chromosome slice if utilization remains low; a costly full-cohort thread sweep is not required.

## 2026-08-25T13:16:20-04:00 — Corrected low-allocation VCF QC and eight-worker plateau

### Baseline and goal

- Baseline commit remains `4524cbdd7161de9c6251411b62e4ccad9b8d666b` (`Add reusable variant preprocessing and CLI reference`); this corrective increment follows the uncommitted 12:53 VCF worker change in the canonical standalone checkout at `D:\projects\NIH-Project`.
- Goal: respond to the production observation that the old reader achieved approximately 1,300 records/s while the first worker-decoding build regressed to roughly 1,000–1,100 records/s, determine why nominal parallelism was slower, and replace it with a measured lower-allocation design.

### Decisions and files changed

- The first parallel implementation was rejected as the default optimization: although it moved HTSJDK genotype decoding into workers, every worker still allocated full `Genotype` objects, allele lists, extended-attribute maps, and unused FORMAT values for thousands of samples. Concurrent object expansion increased garbage collection and memory-bandwidth contention enough to lose throughput. Increasing the already-16-MiB compressed input buffer could not fix that allocation bottleneck.
- VCF QC workers now scan the raw genotype text directly and extract only FORMAT layout, `GT`, `DS`, and genotype filter `FT`. They validate the VCF sample-column count, FORMAT/value count, biallelic GT syntax/ploidy, dosage range, missing calls, and selected-sample filtering while accumulating dosage sum/squared sum, called/missing counts, and exact-HWE genotype counts without materializing per-sample HTSJDK objects. A boolean aligned-sample mask avoids name-map lookups and allows matrix-only samples to be skipped for scientific counts while retaining cheap GT syntax validation.
- The one-worker mode remains the independent HTSJDK reference implementation. Parallel/raw versus sequential/HTSJDK tests require identical summaries and byte-identical QC TSVs for sorted and unsorted headers, aligned sample subsets/reordering, DS and GT selection, FT-filtered samples, partial calls, missing/trailing FORMAT values, and per-record FORMAT changes.
- Automatic `--variant-qc-threads 0` is now capped at eight rather than 16/15. A bounded 4,746-sample benchmark plateaued at eight workers, so more threads are not presented as automatically better. Explicit thread values remain available. The queue remains bounded at four records per worker (32 records in automatic mode).
- Console diagnostics now report `direct low-allocation VCF QC parses`; it should cover biallelic genotype-bearing records in an unrestricted VCF scan. Zero remains expected for the sequential reference, BCF, fully reused checkpoints, and some eager indexed paths. Updated `AGENTS.md`, `README.md`, `TODO.md`, and `docs/COMMAND_LINE_REFERENCE.md` to replace the discarded full-object worker design and document the measured cap.

### Verification, benchmark, and artifact

- `.\mvnw.cmd "-Dtest=QVariantMatrixSourceTest" test` — 18 tests ran with 0 failures/errors and one intentional optional external-BCF skip. New DS/GT, FT, partial-call, changing-FORMAT, and missing-field cases were byte-identical between the direct worker parser and HTSJDK reference.
- Temporary bounded benchmark (test source removed afterward): 1,500 uncompressed biallelic DS+GT variants by 4,746 samples on the current 16-logical-processor Windows host. One HTSJDK worker: 2.503 s / 599.2 variants/s. Four direct workers: 0.328 s / 4,575.8 variants/s / 3.62 average active tasks. Eight: 0.185 s / 8,126.8 variants/s / 5.91 active. Fifteen: 0.185 s / 8,105.7 variants/s / 5.70 active. This establishes the allocation improvement and eight-worker plateau but is not a claim that compressed production VCF.gz will reach synthetic uncompressed rates.
- Final `.\mvnw.cmd test` — 90 tests ran with 0 failures/errors and five intentional skips (four platform/profile-gated oneMKL checks and the optional independently generated BCF). The same previously documented NVIDIA GeForce RTX 2080 / CUDA 13.3 driver API / CUDA 12.6 runtime host was exercised; no driver change or separate probe occurred.
- `.\mvnw.cmd -DskipTests package` — shaded Java 17 build succeeded. Final artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,737,763 bytes, SHA-256 `27B8A4B7CCC3CA5DEA7E45A0E3F4EBC8EA977612CAD28DC8B662FD769473246E`.
- Packaged automatic-mode smoke preprocessing on `test\resources\variant-reference\genotype.vcf` succeeded with eight workers, 32 maximum in-flight records, six direct low-allocation parses, unchanged five retained variants, and a completed aligned raw cache.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Known limitations, compatibility, and next step

- VCF.gz decompression, fixed-field parsing, ordered duplicate/QC/checkpoint consumption, and checkpoint storage remain serial stages. Production throughput may therefore plateau below eight active workers even though the worker parser itself scales. Overall record rate, not Task Manager utilization alone, is the optimization target.
- BCF and certain indexed/eager record paths retain the HTSJDK object decoder. The direct parser intentionally supports the additive biallelic GT/DS model used by this application, not arbitrary unused FORMAT-field semantic validation; HTSJDK remains responsible for full decoding during later retained-row reads.
- The previous 12:53 session entry records the superseded first attempt and must remain as history. Its recommendation to expect 15 useful workers is replaced by this measured eight-worker policy.
- Next: run the new JAR on one fresh or separately signed short compressed VCF.gz slice with automatic workers, record the rate and final pipeline diagnostic, and compare its QC TSV byte-for-byte with `--variant-qc-threads 1`. Do not repeat an expensive whole-chromosome run solely to benchmark this correction.

## 2026-08-25T13:49:54-04:00 — Unsorted VCF indexed-resume bottleneck identified and bypassed

### Baseline and goal

- Baseline commit remains `4524cbdd7161de9c6251411b62e4ccad9b8d666b`, with the 12:53 and 13:16 QC changes still uncommitted in the canonical `D:\projects\NIH-Project` checkout.
- Goal: determine why the user's corrected-JAR run still held at approximately 1,100 records/s and whether a new argument was required.

### Live evidence and decision

- `jcmd -l` and `VM.command_line` confirmed PID 35960 was launched from the expected `D:\projects\NIH-Project\target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` with `--variant-qc-threads 0`, the intended VCF, expression, covariate, checkpoint, cache, profile, and `--preprocess-only` arguments. No missing thread argument caused the result.
- A live `jcmd 35960 Thread.print -l` showed `main` consuming approximately 1,610 CPU seconds in `AbstractVCFCodec.createGenotypeMap -> LazyGenotypesContext.decode -> VCFCodec.decode -> QVariantMatrixSource$IndexedCursor.advance`, while all eight `pool-1-thread-*` workers were parked after only about 58–62 CPU seconds each. Heap use was about 2.30 GiB of 2.81 GiB committed under `-Xmx12g`. This proved the run had resumed its existing checkpoint through the neighboring index and HTSJDK was eagerly decoding unsorted-header genotypes serially before worker submission.
- The running process was inspected read-only and not stopped. Because it launched before this correction, overwriting the JAR on disk does not change its already loaded classes.
- For unrestricted VCF with an unsorted sample header and an incomplete checkpoint, resume now deliberately avoids the eager indexed cursor. It opens the low-allocation sequential cursor, validates every durable checkpoint row against the raw source prefix without genotype expansion, reports prefix-validation progress, and then continues the remaining records through direct parallel QC. It retains the same checkpoint and does not recompute completed EAF/MAF/MAC/HWE.
- Explicit region queries still require and use indexed access. Sorted-header unrestricted VCF and BCF/index-compatible paths may still seek directly to the durable boundary. Added cleanup if prefix validation fails so the source cursor is not leaked.
- Updated `AGENTS.md`, `README.md`, and `TODO.md` with the resume-selection rule. No new user-facing argument was added or required.

### Verification and artifact

- `.\mvnw.cmd "-Dtest=QVariantMatrixSourceTest" test` — 19 tests ran with 0 failures/errors and one optional external-BCF skip. The new test interrupts an unrestricted unsorted VCF that has a tabix index, resumes its checkpoint, asserts indexed resume was not used, verifies worker direct parsing resumed, and preserves all seven ordered records. The explicit indexed-region resume test still seeks directly and passes.
- `.\mvnw.cmd test` — 91 tests ran with 0 failures/errors and five intentional platform/fixture skips. Previously recorded Windows/NVIDIA RTX 2080 CUDA/OpenCL and CPU/OpenBLAS checks remained successful; no driver change occurred.
- `.\mvnw.cmd -DskipTests package` — Java 17 shaded build succeeded. Artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,738,124 bytes, SHA-256 `C88F13215853D937746A7E4AEA0BEC97D9EA5874D215714EB67FE18F9DD7DE85`.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF warnings.

### Compatibility and next step

- The current PID 35960 can safely finish on its existing serial indexed path. To use the correction immediately, the user may terminate that process normally and rerun the exact same command; the new JAR will preserve and validate the durable checkpoint before resuming. Up to the current checkpoint batch, not necessarily the last displayed record, is durable.
- `--variant-qc-threads 8` is equivalent to current automatic mode on this 16-logical-processor host and will not correct the indexed-resume issue by itself. A different checkpoint path would force a fresh parallel scan but unnecessarily discard useful durable work, so it is not recommended.
- During a corrected restart, expect the explicit unsorted-header resume message, possibly a fast prefix-validation phase, then eight-worker direct parsing. A completed checkpoint is reused immediately and needs no scan.

## 2026-08-25T14:44:00-04:00 — Aligned-cache sample-ID reorder fix

### Baseline and goal

- Baseline commit remains `4524cbdd7161de9c6251411b62e4ccad9b8d666b` (`Add reusable variant preprocessing and CLI reference`), with the preceding VCF-QC performance work still uncommitted in the canonical `D:\projects\NIH-Project` checkout.
- Goal: fix the production `ArrayIndexOutOfBoundsException: Index 4747 out of bounds for length 4746` from `QeQTLAnalysis.reorder` immediately after matrix missingness inspection, while preserving the validated genotype/expression/covariate sample order and reusable QC/cache artifacts.

### Cause, decisions, and files changed

- The aligned raw VCF cache contains only the 4,746 selected samples in canonical order. After scanning that cache with an identity column order, `QeQTLAnalysis` incorrectly rebuilt display/output sample IDs by applying the original VCF column permutation to the cache's already-reordered 4,746-element header. A valid original source position such as 4,747 was therefore out of range for the aligned cache header. Matrix values had already been read with the correct order; the exception occurred while reconstructing labels.
- `QeQTLAnalysis` now snapshots the immutable genotype and expression source-header IDs before aligned QC or cache substitution. Alignment and canonical genotype sample labels are consistently derived from that source snapshot, while an aligned raw cache continues to use an identity matrix-column order. This also preserves the intended genotype-side canonical IDs when expression IDs have a stripped prefix.
- Added an end-to-end `QVariantPreprocessOnlyTest` case with source VCF order `S2,S1,S3,S4`, covariate-subset order `S4,S1,S2`, and one excluded matrix-only sample. It covers both creation and reuse of the three-column aligned raw cache and verifies that missing genotype data for `S1` is reported at aligned position 1.
- Recorded the user's production confirmation that the corrected unrestricted VCF.gz resume reached approximately 17,600 records/s, versus approximately 1,000–1,100 records/s on the eager indexed-resume path and approximately 1,300 records/s in the older implementation. Updated `TODO.md` so only the short independent byte-for-byte QC comparison remains for that throughput milestone.
- Files changed in this increment: `src/gov/nih/eqtl/QeQTLAnalysis.java`, `test/gov/nih/eqtl/QVariantPreprocessOnlyTest.java`, `TODO.md`, and this append-only session record. Existing unrelated/uncommitted modernization changes were preserved.

### Verification and artifact

- `.\mvnw.cmd -Dtest=QVariantPreprocessOnlyTest test` — 3 tests ran with 0 failures, 0 errors, and 0 skips. The new reordered covariate-subset case successfully built and then reused the same aligned raw cache.
- `.\mvnw.cmd clean package` — successful Java 17 build and shaded packaging; 92 tests ran with 0 failures, 0 errors, and five intentional skips (four platform/profile-gated oneMKL checks and the optional independently generated external BCF fixture).
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — successful. Automatic discovery reported the NVIDIA GeForce RTX 2080 CUDA backend (CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, 8,589,606,912 bytes global memory) and the Windows x64 CPU fallback (16 logical processors, bundled OpenBLAS 0.3.34, 15 BLAS threads).
- Rebuilt artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,738,069 bytes, SHA-256 `E1D7F54E616EBF5D1613EFA7B5E8E620E3B37AC9FFD688E1FBAC05E5AE429796`. Existing Maven Shade duplicate-resource warnings were unchanged.

### Compatibility, limitations, and next step

- No cache or checkpoint format/signature changed. The user's completed variant-QC checkpoint and aligned raw predictor cache remain reusable; rerunning the same command with this JAR should proceed past missingness inspection without rebuilding them.
- The full 4,746-sample production association was not rerun because it is expensive. The failure shape is covered deterministically at small scale, including cache reuse and a source index larger than the aligned-cache header.
- Next: rerun the same production command with the rebuilt JAR. If it proceeds into residualization/analysis, retain the missingness QC file as the sample-order audit; separately perform the remaining short `--variant-qc-threads 1` versus automatic QC byte comparison when convenient.

## 2026-08-25T15:52:37-04:00 — Production cache reuse and association profile completed

### Baseline and goal

- Baseline commit remains `4524cbdd7161de9c6251411b62e4ccad9b8d666b` (`Add reusable variant preprocessing and CLI reference`), with the VCF-QC performance and aligned-cache sample-ID fixes still uncommitted in the canonical `D:\projects\NIH-Project` checkout.
- Goal: record the user's successful Step #4 production run after the sample-ID fix, confirm that Step #3 artifacts were reusable, and identify the next optimization boundary from the measured profile without requiring another costly full run.

### Production evidence and decisions

- The run printed `Reusing completed variant-QC checkpoint` and proceeded through association progress to completion. This confirms that the aligned-cache sample-ID fix requires neither QC recomputation nor cache invalidation and that the corrected JAR passed the former missingness-inspection failure on the 4,746-sample production inputs.
- Completed variant-QC reuse processed 15,979,340 durable records in 71.805116500 seconds (approximately 222,538 records/s). This is checkpoint reconstruction/materialization work, not genotype QC recomputation.
- The prepared genotype cache contained 440,406 rows and 16,741,963,967 profiled bytes and opened in 0.962446500 seconds. The prepared expression cache contained 114,406 rows and 4,348,917,857 profiled bytes and opened in 0.253872700 seconds. The three files under the matrix-cache root occupied 37,824,081,886 bytes (35.23 GiB), consistent with aligned raw genotype, prepared genotype, and prepared expression storage.
- The association evaluated 50,385,088,836 variant–trait pairs in 1,995.787424300 wall seconds (33.26 minutes), or approximately 25.25 million comparisons/s. There were 12,096 GPU-compute calls totaling 1,544.929592300 worker-seconds and 12,096 CPU result/filter/write calls totaling 1,269.687218700 worker-seconds. These phase totals overlap across scheduled work and must not be added to infer wall time; nevertheless, host-side result processing/output is now a substantial measured optimization target alongside GPU compute.
- Updated `TODO.md` with the successful production checkpoint/cache reuse, storage, and `--trait-cache memory` association baseline. No source code, scientific output, cache, checkpoint, or user data was changed in this documentation-only increment.

### Verification supplied by the user

- `Import-Csv $associationProfile | Where-Object phase -In 'variant_qc','genotype_cache_open_or_build','expression_cache_open_or_build','gpu_compute','cpu_results_and_write','analysis_wall' | Format-Table phase,calls,total_seconds,units,bytes` — reported the exact phase measurements above and one completed `analysis_wall` call.
- `Get-ChildItem -LiteralPath $matrixCache -File -Recurse | Measure-Object -Property Length -Sum` — reported 3 files totaling 37,824,081,886 bytes.
- `Get-FileHash -Algorithm SHA256 $variantQc` — SHA-256 `60740BD16DFCAFD22F8408D07730F15673EBBAF64C0D4A208176D3673FE72F03` for the production variant-QC output.
- No additional hardware test was run for this record. The analysis used the rebuilt artifact and NVIDIA RTX 2080/CUDA environment documented in the immediately preceding entry.

### Compatibility, limitations, and next step

- The production reuse and association portions of the VCF/QC milestone are now validated. The remaining closure check is a short independent QC comparison (automatic versus one-worker HTSJDK and/or bcftools/R); it does not require repeating the full chromosome association.
- Association progress messages are periodic status reports over a 33-minute run, not evidence that checkpoint reuse failed. The checkpoint eliminates the expensive genotype-QC scan; it does not skip the requested 50.385-billion-pair association.
- Before optimizing, inspect the association output size and retained row count so CPU result/filter/write time can be separated into statistical conversion, threshold filtering, formatting, and physical I/O. Then use a bounded representative slice to profile those subphases; do not rerun this full analysis merely to obtain finer timing.

## 2026-08-25T15:55:16-04:00 — Production VCF/VCF.gz milestone closed

### Baseline and goal

- Baseline commit remains `4524cbdd7161de9c6251411b62e4ccad9b8d666b` (`Add reusable variant preprocessing and CLI reference`), with the verified VCF-QC, cache-reuse, and sample-ID fixes still uncommitted in the canonical `D:\projects\NIH-Project` checkout.
- Goal: formally close the current VCF/VCF.gz hardening milestone after the successful 15,979,340-record QC/cache reuse and 50,385,088,836-pair association, defer unavailable BCF production validation, and reprioritize the remaining roadmap.

### Decisions and files changed

- Declared the biallelic additive VCF/VCF.gz workflow production-hardened for the current use case: aligned-sample QC, EAF/MAF/MAC/HWE, missingness preservation/policies, multithreaded direct parsing, durable checkpoint/resume, aligned raw/prepared caches, indexed regions, progress reporting, and full production association have deterministic tests or production evidence.
- Kept an independent short bcftools/R QC comparison as desirable release evidence rather than a gate. The existing automatic-versus-one-worker HTSJDK tests already require byte-identical QC for deterministic fixtures.
- Deferred production BCF validation because the user does not currently have a representative BCF file. Standard CSI decoding, BCF throughput work, production indexed BCF/region validation, and multiallelic splitting remain explicit non-blocking extensions.
- Changed the recommended dependency order to trait-pattern missingness correctness/resume, shared rare-variant set-test foundations, burden, GPU burden acceleration, SKAT, and SKAT-O. Bounded CPU result/output profiling and release validation can proceed alongside that sequence.
- Files changed: `TODO.md` and this append-only `SESSION_HISTORY.md` entry. No source code, binary artifact, scientific output, checkpoint, cache, or user data changed.

### Verification and next step

- Documentation review confirmed the completed VCF/VCF.gz work is separated from deferred BCF/index extensions and that the roadmap no longer lists production VCF/QC validation as the leading dependency.
- No build or hardware test was required for this documentation-only closure. The immediately preceding entry contains the successful 92-test package build, device probe, production cache reuse, association profile, cache sizes, and variant-QC checksum.
- Next recommended implementation: add a deterministic end-to-end CPU reference for exact trait-pattern deletion and compare effective N, identifiers, residual degrees of freedom, effects, and p-values against the optimized/GPU path; then make completed pattern groups resumable. After that correctness boundary, begin the continuous-trait burden-test foundation.

## 2026-08-25T16:33:38-04:00 — Exact missing-data correctness, pattern restart, and production preflight

### Baseline and goal

- Baseline commit: `4524cbdd7161de9c6251411b62e4ccad9b8d666b` (`Add reusable variant preprocessing and CLI reference`) in the canonical standalone checkout at `D:\projects\NIH-Project`; the previously recorded VCF-QC and aligned-cache fixes remain uncommitted in the same worktree and were preserved.
- Goal: establish deterministic end-to-end correctness for exact expression-missingness deletion, make interrupted trait-pattern association work restartable, and validate feasibility against `D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-jointtmmfpkm-ratio-incompleteobs-4746.csv` without launching another expensive full association.

### Decisions and files changed

- Added `QTraitPatternCheckpoint`, a signature-bound durable root checkpoint with atomically committed ordered result and optional pattern-variant-QC parts. Every pattern owns a nested ordinary genotype-block checkpoint. An identical `--resume` run skips complete pattern groups, resumes an interrupted group's genotype blocks, recovers a group assembled immediately before interruption, and atomically assembles final output only after all groups complete. `--keep-checkpoints` retains the durable parts after success; otherwise only known checkpoint files are safely removed.
- The checkpoint signature covers predictor/trait source signatures, sample alignment order, exact trait row masks/order, covariate-model values, missing/frequency policies, MAF/MAC, precision/residualization mode, thresholds, block sizes, DF offset, and output modes. A mismatch is rejected rather than reusing scientifically incompatible work.
- Added a hardware-independent FP64 CPU end-to-end fixture with two exact trait masks. It anchors identifiers, R-squared, effect, t statistic, log10 p-value, effective N, and DF. An injected failure after the first pattern proves the resumed final CSV is byte-identical to an uninterrupted run and that the completed first group is skipped.
- Added exact-pattern preflight before pattern-specific cache/checkpoint creation. It reports pattern count, observed-N range, equivalent full-cohort predictor-preparation passes, and an upper bound on numeric prepared-predictor storage. It rejects structurally non-positive residual DF and defaults `--max-trait-patterns` / `max_trait_patterns` to 256; `0` disables only the count limit, not estimability checks.
- Updated `QeQTLAnalysis`, `QeQTLAnalysisConfig`, `QeQTLCommandLine`, `README.md`, `docs/COMMAND_LINE_REFERENCE.md`, `TODO.md`, and `AGENTS.md`; added `QTraitPatternCheckpointTest` and `QTraitPatternAnalysisTest`, and extended `QeQTLCommandLineTest`. Prior VCF source/test edits in the dirty worktree were left intact.

### Verification and production evidence

- `.\mvnw.cmd "-Dtest=QeQTLCommandLineTest,QTraitPatternCheckpointTest,QTraitPatternAnalysisTest,QeQTLReferenceTest" test` — 14 tests ran with 0 failures, 0 errors, and 0 skips. The safety-limit case also proved no `trait-patterns` cache directory or association checkpoint was created.
- `.\mvnw.cmd clean package` — successful Java 17 compile, full test suite, and shaded packaging; 96 tests ran with 0 failures/errors and five intentional skips (optional platform/profile/fixture cases). Existing Maven Shade duplicate-resource warnings were unchanged.
- Production command (PowerShell, successful invocation): `java -Xmx12g "-Deqtl.cpu.blas=java" --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cpu --genotype D:\Research\topmed\sqtl\sqtl-batch1234\freeze.12c.chr22-fhs-framid.vcf.gz --genotype-format vcf --genotype-field auto --predictor-missing mean --expression D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-jointtmmfpkm-ratio-incompleteobs-4746.csv --trait-missing pattern --covariates D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-mastermat-withpcs-4746.csv --genotype-id-column framid --expression-id-column SampleName --trait-id-strip-prefix X --sample-alignment covariate-subset --fixed-covariates Sex,Age,Batch1,Batch2,Batch3,PCE1,PCE2,PCE3,PCE4,PCE5,PCE6,PCE7,PCE8,PCE9,PCE10,PCE11,PCE12,PCE13,PCE14,PCE15,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --genotype-block-rows 2048 --expression-block-rows 2048 --trait-cache memory --cache-dir D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-chr22-validation\cache --variant-qc-output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-chr22-validation\chr22.variants.tsv --variant-qc-checkpoint D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-chr22-validation\variant-qc-checkpoint --missingness-qc-output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-chr22-validation\incomplete.missingness.tsv --precision fp64 --residualization auto --threshold pval 1e-4 --max-trait-patterns 256 --output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-chr22-validation\incomplete-preflight.csv`.
- That bounded production run forced the portable Java CPU engine only to avoid GPU computation. It reused all 15,979,340 completed VCF-QC records and the existing 440,406-row aligned raw predictor cache, scanned the specified 3,950,614,528-byte expression CSV, and wrote the 15,861,901-byte missingness audit (SHA-256 `42A6D08A373AF187903C6F0AC2C403679E50DE612318455C61EA9E3B48CEEC8D`).
- Production missingness: 178,105 trait rows, 275,746,428 missing cells, 16,862 exact masks, and observed N 20–4,745. With 31 covariate-model columns, 727 patterns covering 2,967 traits cannot leave positive DF for one tested predictor. The current pattern-outer schedule would perform approximately 7,771.9 equivalent full-cohort predictor preparations and write up to 118.19 TiB of FP64 numeric predictor values. Preflight rejected the run before association: the result file, result checkpoint, and `cache\trait-patterns` directory do not exist.
- Rebuilt artifact: `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar`, 135,746,333 bytes, SHA-256 `09919EDABC31DBAB2EC33A15E3754FAB1EF9F4868AF695C8BAD65313F083B756`.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings. No commit or push was requested or performed.

### Compatibility, limitations, and next step

- Exact pattern deletion is now numerically anchored and restartable, but only practical for a modest number of masks. Raising or disabling the new safety limit on the batch1234 file is not recommended; restartability prevents lost completed work but cannot make a roughly 118-TiB preparation schedule affordable.
- Pattern output remains grouped by first-seen mask rather than global trait-row order. Rank loss caused by a categorical level disappearing can still be discovered during a particular pattern's QR even when its raw N is large enough; the preflight currently catches guaranteed non-positive-DF cases, not every subset-rank failure.
- The exact missingness scan itself is deterministic but not cached/resumable, so an identical restart rereads the 3.95-GB expression CSV before it can reuse pattern association work. Persisting/signing this scan is now recorded in `TODO.md`.
- Next: implement a genotype-block-outer or raw sufficient-statistics schedule that reads/prepares each genotype block once while preserving exact pattern-specific covariate projection, predictor imputation/statistics, N/DF, thresholds, output order, and durable restart. Define an explicit audited policy for the 2,967 unestimable traits (prefilter/skip versus an explicitly chosen trait-imputation policy) before attempting the production association. Only then begin the rare-variant set-test foundation.

## 2026-08-25T18:47:47-04:00 — Signed missingness reuse and scalable exact-pattern scheduler

### Baseline and goal

- Baseline commit: `4b8f9a6ebcece020aa1141497c36b5dc31f65157` (`Harden VCF QC and restart missing-data analyses`) in the canonical `D:\projects\NIH-Project` checkout. That commit was pushed to both `origin/master` and `github/master` before this increment began.
- Goal: remove the repeated 3.95-GB missingness scan, replace the 16,862-pass pattern-outer predictor schedule with an exact genotype-block-outer implementation, preserve pattern-specific rank/N/DF/threshold/predictor filling and restart determinism, and test the first scalable implementation on a bounded slice of the user-provided incomplete production data.

### Decisions and files changed

- Added atomically committed CRC-protected `QMIS` exact missingness caches. Their SHA-256-derived name/signature covers normalized source path, size/mtime, matrix dimensions, parser/cache tag, selected source columns, selected sample IDs, and exact selected order. Cache payload validation cross-checks sample/row missing totals, pattern coverage, masks, row assignments, dimensions, and checksum; corruption or incompatibility rebuilds instead of being reused.
- Added `--trait-pattern-scheduler auto|pattern|genotype` and `--unestimable-trait-patterns error|skip`. Automatic mode retains compatibility pattern-outer scheduling at or below `--max-trait-patterns` and selects genotype-outer above it. Strict mode writes `<output>.trait-patterns.tsv` before failing; explicit skip excludes only audited rank/DF-unestimable traits.
- Added an FP64 genotype-block-outer exact scheduler. It prepares each estimable trait once under its exact mask, pads excluded samples with zero, and writes one reusable checksummed global trait cache. Per genotype block, backend matrix products obtain exact pattern called counts, dosage sums/squares, `X'g`, and missing-indicator trait products. CPU FP64 triangular solves using each pattern's QR `R` produce exact projected genotype variance and trait residuals. Mean and zero predictor filling are computed inside each mask; pattern-monomorphic variant/trait combinations are skipped.
- Model/rank construction is parallel across up to `availableProcessors - 1` CPU workers and reports progress. Trait preparation solves through cached `R` rather than reconstructing a large `Q` per pattern. Backend workers are bounded to one per exclusive selected compute context, preserving multi-GPU execution without accumulating waiting genotype blocks.
- Genotype-outer result/checkpoint order is explicitly genotype block, trait block, variant, then original estimable trait row. Checkpoint signatures cover source/order, global trait/model signatures, block sizes, predictor missing policy, threshold/DF/output modes, and FP64. An injected pre-assembly interruption test proves byte-identical resumed output.
- Added approximately 15-second prepared-trait-cache row progress and ETA. Automatic global block capacity now rises to honor explicitly larger genotype/expression streamed row blocks; this fixes the bounded-region failure where 146 variants caused automatic capacity 32 to reject an explicitly requested 2,048-row trait block.
- Added deterministic tests for signed scan reuse/corruption/source-and-order invalidation; QR sufficient statistics versus explicit residualization; genotype-outer versus pattern-outer full statistics; VCF missing-genotype mean and zero filling; strict audit-before-error versus audited skip; automatic capacity; and genotype-outer restart.
- Updated `AGENTS.md`, `README.md`, `TODO.md`, `docs/COMMAND_LINE_REFERENCE.md`, `src/gov/nih/eqtl/QeQTLAnalysis.java`, `QeQTLAnalysisConfig.java`, `QeQTLCommandLine.java`, `QGenotypeOuterPatternJob.java`, `QPatternPreparedTraitSource.java`, `QPatternSufficientStatistics.java`, `QTraitPatternModelSet.java`, `QTraitPatternScheduler.java`, `QUnestimableTraitPolicy.java`, `src/gov/nih/eqtl/io/QBinaryMatrixCache.java`, `QMissingnessScan.java`, `test/gov/nih/eqtl/QTraitPatternAnalysisTest.java`, `QPatternSufficientStatisticsTest.java`, `QeQTLCommandLineTest.java`, `test/gov/nih/eqtl/io/QMissingnessScanCacheTest.java`, and this append-only record.

### Verification

- `.\mvnw.cmd -q "-Dtest=QMissingnessScanCacheTest,QTraitPatternAnalysisTest,QPatternSufficientStatisticsTest,QeQTLCommandLineTest" test` and subsequent focused reruns — all selected tests passed. The VCF fixture required identical identifiers/statistics within `2e-12` for both pattern-specific mean and zero predictor filling. The strict unestimable case wrote its complete audit before error; skip retained only estimable trait rows. Corrupt QMIS state rebuilt, changed source metadata/order selected a new signature, and identical state reused.
- CUDA/CPU fixture commands used the shaded JAR with `--genotype test\resources\variant-reference\genotype.vcf --genotype-format vcf --genotype-field auto --predictor-missing mean --expression test\resources\variant-pattern-reference\expression.csv --trait-missing pattern --sample-alignment strict --min-mac 0 --genotype-block-rows 16 --expression-block-rows 16 --trait-cache disk --precision fp64 --residualization cpu --trait-pattern-scheduler genotype --threshold none 0`, once with `--backend cuda` and once with `--backend cpu -Deqtl.cpu.blas=java`. Eight retained rows matched with maximum absolute numeric delta `0`.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printgpuinfo` — automatic discovery reported NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, and 8,589,606,912 bytes VRAM; CPU fallback reported Windows x64, 16 logical processors, and bundled OpenBLAS 0.3.34.
- Production bounded command: `java -Xmx12g --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cuda --genotype D:\Research\topmed\sqtl\sqtl-batch1234\freeze.12c.chr22-fhs-framid.vcf.gz --genotype-format vcf --genotype-field auto --predictor-missing mean --expression D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-jointtmmfpkm-ratio-incompleteobs-4746.csv --trait-missing pattern --covariates D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-mastermat-withpcs-4746.csv --genotype-id-column framid --expression-id-column SampleName --trait-id-strip-prefix X --sample-alignment covariate-subset --fixed-covariates Sex,Age,Batch1,Batch2,Batch3,PCE1,PCE2,PCE3,PCE4,PCE5,PCE6,PCE7,PCE8,PCE9,PCE10,PCE11,PCE12,PCE13,PCE14,PCE15,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --region chr22:10510000-10520000 --min-mac 20 --genotype-block-rows 256 --expression-block-rows 2048 --trait-cache disk --cache-dir D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\cache --variant-qc-output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\variants.tsv --variant-qc-checkpoint D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\variant-qc-checkpoint --missingness-qc-output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\missingness.tsv --checkpoint-dir D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\association-checkpoint --precision fp64 --residualization cpu --trait-pattern-scheduler auto --unestimable-trait-patterns skip --threshold pval 1e-4 --profile --profile-output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\profile.csv --keep-checkpoints --output D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke\results.csv`.
- That production run queried 3,024 indexed VCF records and retained 146 variants. Parallel exact model preflight completed 16,862 masks in 28.3 seconds: 16,127 masks/175,108 traits were estimable and 735 masks/2,997 traits were excluded. The global prepared cache is 6,656,382,998 bytes. CUDA evaluated 25,565,768 estimable comparisons in 27 seconds (approximately 940,730/s), retained 58,038 output rows, and completed in 241.4 overall seconds. Profile totals included 338 backend calls, 34,048.83 MiB uploaded, 3,711.81 MiB downloaded, and 8.755 GPU-compute seconds.
- The identical production command plus `--resume --profile-output ...\resume-profile.csv` reused the 3,024-record variant-QC checkpoint, signed expression scan, aligned raw cache, global trait cache, and all 25,565,768 association comparisons. It completed in 78.3 seconds overall; exact model/audit preflight was recomputed in 27.5 seconds and analysis wall was 31.36 seconds with no association compute calls.
- Production SHA-256 values: `results.csv` `8AD16A456844F47ADD6D554789B8ED04B7652E57B00DFBF4DD85BA95C499753B`; `results.csv.trait-patterns.tsv` `8EDAFBCC25FC88EB190BEB796B086B436694856D985D6ABA792CF5DEBF2D4C2C`; `missingness.tsv` `3E498935A0A2F20EC29D448A1FF61C65BCD09C6411A9567731BC40140BE72B49`; `variants.tsv` `7B84D6B7D2E3D916DCE4CA3205B753821A9575B4D07A6836726303635ACD183A`.
- `.\mvnw.cmd clean package` — successful final Java 17 build; 103 tests ran with zero failures/errors and five intentional platform/fixture skips. Existing Maven Shade duplicate-resource warnings were unchanged. Final shaded JAR size: 135,788,362 bytes; SHA-256 `D678F60881D01E2D28F4BF6BD142BF7949EF46DB4538AC73549C95EE3CB7D15D`.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Compatibility, limitations, and next step

- Legacy INI, explicit pattern-outer behavior, existing raw/prepared/variant-QC caches, complete-trait output, and default FP64 behavior remain compatible. New global trait and missingness caches use new signatures/files and do not invalidate older unrelated artifacts.
- Genotype-outer is correctness-first and currently supports FP64 plus aligned-scope MAF/MAC. It does not yet emit durable aggregate pattern-variant QC counts or support `--frequency-scope pattern`; FP32 remains disabled pending an accuracy study. Exact pattern/model QR is deliberately recomputed on restart (about 28 seconds for this fixture), while the much larger scan and trait cache are reused.
- Do not launch the full incomplete-trait chr22 run yet. The bounded scheduler achieved about 0.94 million comparisons/s versus the 25.25-million/s complete-trait baseline and transferred large `X'g` intermediates/design masks. Next implement an on-device pattern-statistics finalization/reduction with reusable device design state, so only compact per-pattern/per-variant variance/fill results return to the host; then benchmark 256–2,048-variant blocks against the pattern-outer FP64 reference before authorizing the expensive full run.
- Production smoke artifacts are isolated under `D:\Research\topmed\sqtl\sqtl-batch1234\gpu-eqtl-genotype-outer-smoke` and are not repository content. `target/` remains ignored and was not staged.

## 2026-08-26T16:20:00-04:00 — Compact GPU pattern finalization and production CSV block sweep

### Baseline and goal

- Baseline commit: `5e2c655` (`Add scalable exact missingness scheduling`) on `master` in `D:\projects\NIH-Project`; the worktree was clean before this increment.
- Goal: complete the leading missing-data performance item by retaining exact genotype-pattern sufficient-statistics intermediates on-device, downloading only compact finalized state, and benchmark 256–2,048-row genotype blocks using the user's smaller `D:\Research\topmed\sqtl\sqtl-batch1234\chr22-filtered-batch1234-dose-meanimp.csv` input against the production incomplete-trait and covariate data.

### Decisions and files changed

- Added backend-neutral `GpuPatternStatisticsPlan`, `GpuPatternStatisticsResult`, and `GpuPatternStatisticsSupport`. The plan carries estimable masks, observed indices, design sums, and QR `R`; the compact result carries FP64 replacement, residual sum of squares, and filled sum of squares for each pattern/variant. The default implementation preserves a portable CPU/reference path through the existing matrix-operation contract.
- Changed genotype-outer pattern preparation to invoke the compact operation once per genotype block and validate the device-computed residual against the same FP64 roundoff rule used by host finalization. Pattern IDs, observed N, rank/DF, predictor filling, monomorphic handling, trait ordering, thresholds, and output ordering remain owned by the established analysis layer.
- Added CUDA finalization through an NVRTC-compiled driver kernel after each cuBLAS DGEMM pattern batch. The packed genotype aggregate matrix is uploaded once per genotype block; reusable device buffers hold pattern inputs, `X'g` products, and compact outputs. Missing NVRTC, compilation, launch, or CUDA API failures are fatal.
- Added the equivalent FP64 OpenCL finalizer through the selected ICD compiler and changed the shared product buffer to read/write so the finalizer can consume it. CUDA/OpenCL finalizers currently bound design rank to 64; the production model has 31 columns. The CPU fallback has no new device-compiler requirement.
- Extended hardware integration coverage to compare compact CUDA, OpenCL, and CPU results with explicit pattern-specific QR residualization over multiple pattern batches and mean filling. The tests also assert that profiled pattern-result download bytes equal the compact result size.
- Updated `README.md`, `TODO.md`, and `AGENTS.md` with the delivered architecture, CUDA NVRTC requirement, benchmark evidence, and remaining genotype-outer limits. Source changes are in `src/gov/nih/gpu/GpuContext.java`, the three new `GpuPatternStatistics*` files, `src/gov/nih/gpu/cuda/CudaGpuBackend.java`, `CudaGpuContext.java`, `CudaGpuDevice.java`, `src/gov/nih/gpu/opencl/JoclGpuContext.java`, `src/gov/nih/eqtl/QTraitPatternModelSet.java`, `QPatternSufficientStatistics.java`, `QGenotypeOuterPatternJob.java`, and `test/gov/nih/gpu/GpuKernelIntegrationTest.java`.
- No source research file was modified. `QCsvSubsetTool` wrote the first 2,048 genotype rows to the ignored `target\pattern-benchmark\chr22-first-2048.csv`. Existing expression missingness/global trait caches were exposed under the ignored target benchmark directory with NTFS hard links, avoiding a duplicate multi-gigabyte copy and any write to the research directory. Benchmark results, profiles, and checkpoints remain ignored under `target\pattern-benchmark` and are removable by the ordinary Maven clean lifecycle.

### Verification and production evidence

- `.\mvnw.cmd -q "-Dtest=QPatternSufficientStatisticsTest,QTraitPatternAnalysisTest" test` — selected portable sufficient-statistics and end-to-end exact-pattern tests passed.
- `.\mvnw.cmd -q "-Dtest=GpuKernelIntegrationTest" test` — 21 tests ran with zero failures/errors and three intentional hardware/backend skips. The CUDA compact-finalizer test ran on the production NVIDIA device; the CPU fallback test ran; unavailable OpenCL coverage skipped cleanly. An initial test-fixture failure exposed the association kernel's required 64-sample padding and was corrected by using the production padding contract.
- `.\mvnw.cmd -q "-Dtest=QTraitPatternAnalysisTest,GpuKernelIntegrationTest" test` — focused end-to-end plus hardware integration rerun passed.
- `.\mvnw.cmd -q test` — final full suite passed with exit code 0: 106 tests in 25 suites, zero failures, zero errors, and five intentional platform/profile/fixture skips.
- `.\mvnw.cmd -q -DskipTests package` — final shaded package succeeded after the full test run. `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` is 135,802,422 bytes with SHA-256 `FC601D1D073BC536D8AD11E4AF447733C84F1180B6E17CEE7158DBE6CF96081B`.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` — detected NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, 8,589,606,912 bytes global memory, and the Windows x64 CPU fallback with bundled OpenBLAS 0.3.34. No usable OpenCL device was reported on this host.
- Benchmark input creation: `java -cp target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar gov.nih.eqtl.QCsvSubsetTool D:\Research\topmed\sqtl\sqtl-batch1234\chr22-filtered-batch1234-dose-meanimp.csv target\pattern-benchmark\chr22-first-2048.csv 2048`. The source is 4,082,010,589 bytes with 4,746 sample columns; the 2,048-row subset is 19,538,177 bytes.
- The benchmark command was run four times with `--genotype-block-rows` set exactly to `256`, `512`, `1024`, and `2048`, and matching `checkpoint-<block>`, `missingness-<block>.tsv`, `profile-<block>.csv`, and `results-<block>.csv` target paths: `java -Xmx12g --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cuda --genotype target\pattern-benchmark\chr22-first-2048.csv --genotype-format csv --predictor-missing error --expression D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-jointtmmfpkm-ratio-incompleteobs-4746.csv --trait-missing pattern --covariates D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-mastermat-withpcs-4746.csv --genotype-id-column framid --expression-id-column SampleName --trait-id-strip-prefix X --sample-alignment covariate-subset --fixed-covariates Sex,Age,Batch1,Batch2,Batch3,PCE1,PCE2,PCE3,PCE4,PCE5,PCE6,PCE7,PCE8,PCE9,PCE10,PCE11,PCE12,PCE13,PCE14,PCE15,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --genotype-block-rows <block> --expression-block-rows 2048 --trait-cache disk --cache-dir target\pattern-benchmark\cache --missingness-qc-output target\pattern-benchmark\missingness-<block>.tsv --checkpoint-dir target\pattern-benchmark\checkpoint-<block> --precision fp64 --residualization cpu --trait-pattern-scheduler genotype --unestimable-trait-patterns skip --threshold pval 1e-4 --profile --profile-output target\pattern-benchmark\profile-<block>.csv --keep-checkpoints --output target\pattern-benchmark\results-<block>.csv`.
- Each sweep run evaluated 2,048 variants against 175,108 estimable traits: 358,621,184 active comparisons over 16,127 estimable patterns; 735 patterns/2,997 traits were excluded by the audited rank/DF policy. Each output retained the same 876,057 variant/trait identifiers and exact N/DF. The 512-, 1,024-, and 2,048-row results were numerically identical by identifier. The 256-row result had only FP64 accumulation-order differences: maximum absolute deltas were `1.50e-15` R-squared, `2.00e-15` effect, `7.67e-13` T, and `1.23e-13` log10 p, with no retained-set/threshold changes.
- Analysis-wall seconds for 256/512/1,024/2,048 rows were 303.76/230.43/234.16/240.42. The 512-row run was fastest on this fixture at approximately 1.82 million comparisons/s during reported association progress. Corresponding uploads were 211.81/112.31/62.57/37.69 GiB; total downloads were 6.11 GiB in each run, including ordinary association result matrices. The compact pattern component for all eight 256-row blocks was approximately 0.76 GiB, versus roughly 25 GiB for the removed host `X'g` intermediates.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings. No commit or push was requested or performed.

### Compatibility, limitations, and next step

- FP64 remains the only genotype-outer scientific mode. Legacy INI/CLI behavior, CPU fallback, pattern-outer scheduling, aligned-scope frequency filters, cache/checkpoint signatures, identifiers, DF, thresholds, and deterministic within-block output remain compatible. Different genotype block sizes deliberately produce different block-grouped row order, so equality was checked by identifier rather than raw file hash.
- CUDA genotype-outer exact-pattern analysis now requires NVRTC matching the usable CUDA runtime. OpenCL requires a usable FP64 device and compiler; that path compiled in tests only where hardware is available and was not exercised on this host. Rank above 64 fails explicitly instead of silently producing incorrect device results.
- The compact finalizer removes the former result-download bottleneck, but repeated expanded design-mask uploads and GEMM work still dominate. At the measured 512-row rate, extrapolating the entire 440,406-variant chromosome is approximately a half-day run; this was not launched automatically.
- Next recommended performance work: avoid computing/uploading unused third aggregate columns for covariate rows and investigate retaining or more compactly encoding reusable pattern design state on-device. Independently, add durable per-pattern genotype-outer QC counts before enabling pattern-scoped MAF/MAC, and keep FP32 disabled until an explicit accuracy study is approved.

## 2026-08-26T19:53:52-04:00 — Initial rare-variant set-test and burden reference foundation

### Baseline and goal

- Baseline commit: `5e2c655` (`Add scalable exact missingness scheduling`) on `master`. The immediately preceding compact GPU pattern-finalizer changes remained uncommitted in the same worktree and were preserved without rewriting or reverting them.
- Goal: begin the next dependency-ordered roadmap item, the shared rare-variant set-test foundation, with conservative scientific contracts and a deterministic continuous-trait FP64 CPU burden reference before any CLI or GPU acceleration.

### Decisions and files changed

- Added `QVariantSetTable`, a strict UTF-8 TSV parser for ordered `SET_ID`, `VARIANT_ID`, `REF`, `ALT`, `EFFECT_ALLELE`, and optional `WEIGHT`. Weights default to one and must be finite/positive. Repeated membership within a set is fatal; the same variant may occur in multiple sets. Normalized ordered definitions receive a SHA-256 content signature.
- Chose exact allele harmonization for the reference boundary: source REF/ALT must equal declared REF/ALT after case normalization, the effect allele must be one of them, ALT uses source dosage, and REF uses `2 - ALT dosage`. No strand complement, REF/ALT swap, or frequency-based guess is permitted.
- Added `QSetTestPolicy` with inclusive aligned-cohort minimum/maximum MAF and MAC masks. The unfiltered policy assumes no rare cutoff (`MAF 0–0.5`, `MAC 0–infinity`); callers must explicitly request a rare-variant maximum. Predictor missingness is restricted to `error`, called-sample effect-allele `mean`, or effect-allele `zero`. Absent variants and empty/post-projection-monomorphic sets each use explicit error/skip behavior.
- Added `QSetTestNullModel`, which validates a complete aligned covariate design, rejects rank loss/non-positive residual DF/invalid traits, prepares continuous traits once with the established FP64 CPU projector, and reuses that state across sets. Duplicate or blank trait identifiers are fatal.
- Added `QBurdenReference`, a deterministic scalar weighted burden implementation. It audits requested, absent, frequency-excluded, and included variants, accumulates in set-file order, projects through the shared null Q, and uses the production correlation/effect/t/log10-p conversion. Unweighted burden is the exact all-weights-one case. Results are ordered by set first appearance and then trait order and retain set/trait IDs, variant count, N, and residual DF.
- Added the committed eight-sample reference files under `test/resources/set-test-reference` and `QSetTestFoundationTest`. The fixture covers overlapping membership, an REF effect allele, a non-unit weight, genotype mean filling, MAF filtering, empty and monomorphic skip audits, exact statistics, allele mismatch, missing-dosage error, duplicate membership, and covariate rank loss.
- Added `docs/SET_TEST_FOUNDATION.md` and updated `README.md`, `TODO.md`, and `AGENTS.md`. The documentation explicitly identifies this as a developer/reference API rather than a supported CLI analysis. Production aligned-source adapters, region expansion, audit/result output, bounded streaming, and checkpoint/restart remain open.
- New source files: `src/gov/nih/eqtl/settest/QVariantSetTable.java`, `QSetTestPolicy.java`, `QSetTestNullModel.java`, and `QBurdenReference.java`. New test files: `test/gov/nih/eqtl/settest/QSetTestFoundationTest.java` plus `sets.tsv`, `variants.csv`, `traits.csv`, and `covariates.csv` in `test/resources/set-test-reference`.

### Verification

- `.\mvnw.cmd -q -DskipTests compile` — the new production set-test classes compiled successfully under the Java 17 target.
- Independent fixture calculation used bundled Python/NumPy QR over the committed covariate, trait, and dosage values. It produced R-squared `[0.6747780738826216, 0.5395555683137198, 0.01083266499989602, 0.009789919608688065]`, effects `[0.406515205649314, -0.10457115235072852, 0.14854046495368042, -0.04062226883412121]`, and t statistics `[3.220887637744619, -2.4205523460484324, 0.23400103342370288, -0.22233656583953693]`; these were anchored before the Java assertions. The Java reference's established Student-t conversion yielded log10 p `[-1.6299970085493116, -1.2212956736106053, -0.08393386934985025, -0.07943403355250683]`.
- `.\mvnw.cmd -q "-Dtest=QSetTestFoundationTest" test` — three new tests passed with zero failures/errors/skips, including direct reads of all committed fixture files.
- `.\mvnw.cmd -q "-Dtest=QSetTestFoundationTest,QeQTLReferenceTest" test` — the new burden reference and existing scalar eQTL reference tests passed together.
- `.\mvnw.cmd -q test` — final full suite passed with exit code 0: 109 tests in 26 suites, zero failures, zero errors, and five intentional platform/profile/fixture skips.
- The full suite retained real CUDA integration coverage on the same NVIDIA GeForce RTX 2080 documented in the preceding entry: CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, and 8,589,606,912 bytes VRAM. The new set-test reference itself is hardware-independent and uses no GPU. No usable OpenCL device was available on this host, so optional OpenCL hardware tests skipped cleanly.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings. No commit or push was requested or performed.

### Compatibility, limitations, and next step

- Existing eQTL CLI/INI behavior, outputs, caches, checkpoint signatures, matrix/VCF readers, and GPU/CPU paths are unchanged by the new reference package. There is intentionally no set-test CLI switch yet, so users cannot accidentally mistake the reference for a production burden analysis.
- The current reference requires already aligned in-memory variants, complete continuous traits, and a prepared covariate design. It does not yet expand existing genomic interval memberships, stream large sources, write a stable audit/result file, resume blocks, handle exact trait missingness patterns, or compute SKAT/SKAT-O.
- Next recommended step: connect `QVariantSetTable` and existing indexed region memberships to aligned CSV/VCF/BCF variant blocks, define stable set-audit and burden-result schemas plus CLI/INI options, and build a bounded deterministic CPU runner with signature-bound checkpoint/restart. Only after its end-to-end bytes/statistics match this reference should optimized CPU/GPU burden work begin.

## 2026-08-27T01:35:13-04:00 — Production burden, SKAT, and SKAT-O set tests

### Baseline and goal

- Baseline commit: `25e09e86e6dc58a0c7e6c29bfda3e232724ef937` on `master`; the worktree was clean at the start of this increment.
- Goal: complete the next four rare-variant roadmap items: production set-test input/output and restart integration, optimized burden execution, deterministic weighted SKAT, and adjusted SKAT-O, using the existing chr22 CSV data for a bounded production smoke test.

### Decisions and files changed

- Added user-facing `--analysis burden|skat|skat-o`, strict explicit set definitions, aligned CSV/VCF/BCF adapters, and automatic exact ALT-effect/unit-weight adaptation from indexed VCF/BCF region memberships when no explicit TSV is supplied. Existing region order, overlap, and empty declared regions are preserved.
- Added inclusive set-level MAF/MAC masks, error/skip policies, a stable CSV result schema and TSV set audit, p-value/R-squared retention, configurable set and trait tiles, and atomic signature-bound checkpoint parts for every deterministic set-tile/trait-tile pair. Signatures include inputs, aligned columns, fixed-effect design, definitions, policy, thresholds, tile sizes, method, rho grid, simulations, and seed.
- Retained `QBurdenReference.analyze` as the scalar oracle and added a production batched FP64 set-by-trait product path. Tests compare audit state and every statistic with the scalar implementation to tight FP64 tolerances.
- Added deterministic continuous-trait weighted SKAT: projected effect-dosage kernel columns, null residual variance with `N-rank(X)` DF, batched FP64 score/Gram products, CPU eigenvalues, exact single/equal-eigenvalue chi-square cases, bounded Imhof inversion, and an explicitly labeled Satterthwaite fallback.
- Added SKAT-O with a strictly increasing configurable rho grid that starts at zero and ends at one (default `0,0.25,0.5,0.75,1`), weighted burden/SKAT mixture factors, shared seeded Gaussian replicates for correlated component adjustment, precomputed component critical statistics, and `(extreme+1)/(simulations+1)` reporting. The minimum component p-value remains an audit field and is never substituted for the adjusted result.
- Corrected legacy CSV missing recognition in the burden reference to accept the project sentinel as well as NaN. Exact trait-pattern deletion and FP32 remain rejected for production set tests.
- Added/changed production files: `QeQTLAnalysis.java`, `QeQTLAnalysisConfig.java`, `QeQTLCommandLine.java`, `QSetTestMethod.java`, `QSetTestRunner.java`, `QKernelSetReference.java`, `QSkatReference.java`, `QBurdenReference.java`, `QSetTestNullModel.java`, `QSetTestPolicy.java`, and `QVariantSetTable.java`. Added/changed tests: `QSetTestRunnerTest.java`, `QSkatReferenceTest.java`, `QSetTestFoundationTest.java`, and `QeQTLCommandLineTest.java`. Updated `README.md`, `TODO.md`, and `docs/SET_TEST_FOUNDATION.md`.

### Verification and production evidence

- `\.\mvnw.cmd -q "-Dtest=QSkatReferenceTest,QSetTestRunnerTest,QSetTestFoundationTest,QeQTLCommandLineTest" test` — final focused suite passed: 17 tests, zero failures/errors/skips. It covered scalar-versus-batched burden, bounded partial set/trait tiles, byte-identical resume, explicit/region definitions, masks, missing/empty/monomorphic policies, exact/eigenvalue SKAT cases, seeded correlated SKAT-O, and CLI parsing.
- `\.\mvnw.cmd test` — full suite passed: 116 tests, zero failures, zero errors, five intentional platform/profile/fixture skips; Maven reported `BUILD SUCCESS` in 32.187 seconds.
- `\.\mvnw.cmd -q package` — final shaded package and its repeated test phase completed with exit code zero. `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` is 135,870,653 bytes with SHA-256 `7C46B0D775328181FD1D82C379026E75D904F1D56202208560AD838D3092B5D3`. Jar inspection confirmed `META-INF/LICENSE`, `META-INF/LICENSE_EXCEPTION`, and the Intel oneMKL/OpenMP notices under `META-INF/licenses/`.
- Production smoke data were mechanically copied only into ignored `target\settest-smoke`: the first eight real rows of `D:\Research\topmed\sqtl\sqtl-batch1234\chr22-filtered-batch1234-dose-meanimp.csv`, the first two real rows of the matching complete-observation trait CSV, all 4,746 sample columns, and two four-variant exact-allele sets. The source research files were read-only and unchanged.
- Final production smoke command, run once for each method with `<method>` replaced by `burden`, `skat`, and `skat-o`: `java -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --gpu-backend cpu --genotype target/settest-smoke/variants.csv --genotype-format csv --expression target/settest-smoke/traits.csv --covariates D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-mastermat-withpcs-4746.csv --fixed-covariates Age --genotype-id-column framid --expression-id-column SampleName --expression-id-strip-prefix X --output target/settest-smoke/final-<method>.csv --analysis <method> --variant-sets target/settest-smoke/sets.tsv --predictor-missing mean --trait-missing error --set-degenerate skip --set-block-size 1 --expression-block-rows 1 --threshold none 0 --skat-o-simulations 99 --skat-o-seed 23`. All three completed with four ordered set/trait rows plus the header and companion audits. Result SHA-256 values were burden `6340AB2140D6F4CE32CD71495F53365F81B768DD079193FC1A913F960D26361E`, SKAT `4438B2ACBD05D83E745BE5579500290FE7EC1EF5546C65D308229839CDFAEA9D`, and SKAT-O `1154672AC30E442F953E86F01FCEAFEE574487814C7853EDC780E2FA96DE44FD`.
- An initial smoke attempt without the existing covariate ID bridge failed correctly because genotype IDs are `framid` while expression headers use prefixed `SampleName`. The final commands explicitly selected both ID columns and literal expression prefix removal. An initial SKAT-O implementation also exposed repeated Imhof work inside each replicate; the test process was stopped, component critical statistics were precomputed, and the same correlated extreme-event calculation then completed the 4,746-sample smoke in bounded time.
- `java --enable-native-access=ALL-UNNAMED -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` — detected NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, 8,589,606,912 bytes VRAM, and the Windows x64 CPU backend using bundled OpenBLAS 0.3.34 with 15 BLAS threads. No usable OpenCL device was reported. The production set-test smoke explicitly selected CPU; hardware-independent references do not require a GPU.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Compatibility, limitations, and next step

- Ordinary eQTL remains the default and legacy INI/CLI behavior is preserved. Set tests are FP64 continuous-trait modes; exact trait-pattern deletion, FP32 preprocessing, cohort relatedness, and implicit allele flips/swaps are deliberately unsupported.
- Optimized CPU products are production-supported and anchored to scalar FP64 references. Direct CUDA/OpenCL set-test products are not enabled yet; add them only after profiling shows benefit and every statistic/final adjusted p-value matches the deterministic references on singleton, collinear, overlapping, missing, empty, and common/rare fixtures.
- Imhof non-convergence is explicit in `P_value_method` and uses the documented Satterthwaite approximation. SKAT-O precision depends on the requested simulation count; the output records simulations and the derived per-set/trait seed.
- Next recommended step: add larger independent singleton/high-collinearity SKAT/SKAT-O fixtures and profile set/trait tile sizes plus Imhof convergence on representative production set definitions before implementing an optional multi-GPU score/covariance path.

## 2026-08-27T01:58:40-04:00 — Native sliding-window set tests

### Baseline and goal

- Baseline commit: `48e9e020e201620bcc6446843c528902550216c7` (`Add production burden, SKAT, and SKAT-O analyses`) on `master`; the worktree was clean at the start.
- Goal: let Burden, SKAT, and SKAT-O users specify ordinary sliding-window size and stride directly, without requiring a generated regions TSV, while retaining explicit TSV/region definitions for custom memberships and weights.

### Decisions and files changed

- Added `--window-size BP` / `window_size` and optional `--window-stride BP` / `window_stride`; stride defaults to size. Both must be positive, stride cannot exceed size, and the automatic mode is restricted to set-test analyses and is mutually exclusive with `--variant-sets`, `--region`, and `--regions-file`.
- Defined a deterministic chromosome-local one-based grid anchored at coordinate 1: starts are `1 + k * stride`, ends are inclusive, only windows containing a retained variant are emitted, set IDs are `CHROM:START-END`, contigs retain source first-appearance order, windows are start-sorted, and variant membership retains source order.
- Added aligned CSV and VCF/BCF generation from canonical `CHROM:POS:REF:ALT` identifiers, allowing additional colon-separated ID fields. Automatic memberships use POS, exact ALT effect orientation, and unit weight. Explicit definitions remain the path for custom weights/effect alleles/memberships.
- Clarified the model contract: every expression-matrix row is a tested phenotype; every selected fixed covariate is an adjustment variable, not a competing main phenotype.
- Changed `QeQTLCommandLine.java`, `QeQTLAnalysisConfig.java`, `QeQTLAnalysis.java`, and `QVariantSetTable.java`; extended `QeQTLCommandLineTest.java` and `QSetTestFoundationTest.java`; updated `README.md`, `TODO.md`, `docs/COMMAND_LINE_REFERENCE.md`, and `docs/SET_TEST_FOUNDATION.md`.

### Verification

- `.\mvnw.cmd -q "-Dtest=QSetTestFoundationTest,QeQTLCommandLineTest" test` — passed after Maven resolved the declared antrun plugin; deterministic tests cover overlapping 10-bp windows at 5-bp stride, chromosome separation, nonempty-window emission, source membership order, ALT/unit-weight adaptation, stride validation, explicit parsing, and stride-default-to-size behavior.
- `.\mvnw.cmd test` — full suite passed: 118 tests, zero failures, zero errors, five intentional hardware/platform skips; Maven reported `BUILD SUCCESS` in 24.513 seconds.
- `.\mvnw.cmd -q -DskipTests package` — packaged the updated runnable jar successfully. `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` is 135,874,782 bytes with SHA-256 `C0743E0D0499884AE3242989F0C4FF51297C895456D9269B7191FF7C08A3DEBD`.
- Real-data CLI smoke, with no variant-set or region file: `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cpu --genotype target\settest-smoke\variants.csv --genotype-format csv --expression target\settest-smoke\traits.csv --covariates D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-mastermat-withpcs-4746.csv --fixed-covariates Age --genotype-id-column framid --expression-id-column SampleName --expression-id-strip-prefix X --output target\settest-smoke\sliding-burden.csv --analysis burden --window-size 500 --window-stride 250 --predictor-missing mean --trait-missing error --set-degenerate skip --set-block-size 2 --expression-block-rows 1 --threshold none 0`. Eight real chr22 variants generated five ordered nonempty windows and ten result rows for two real expression phenotypes over all 4,746 samples; the companion audit recorded exact membership for each overlapping window.
- `java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` — detected NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, FP64, and 8,589,606,912 bytes VRAM; CPU was OpenBLAS 0.3.34 with 15 BLAS threads. The production sliding-window smoke explicitly used CPU.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Compatibility, limitations, and next step

- Existing explicit set TSVs and indexed custom region behavior are unchanged. Automatic windows require canonical genomic variant IDs, use unit ALT weights, and deliberately omit empty grid cells; custom weights/effect alleles or empty interval auditing still use explicit definitions.
- All expression rows remain analyzed because there is no row-level phenotype selector. Use a phenotype-subset expression matrix when only selected outcomes are intended.
- The bounded set runner rereads the variant source per set tile. This preserves bounded heap use but can amplify I/O for very small strides that create thousands of windows. Profile representative window/stride choices and tune `--set-block-size` before launching a whole-chromosome schedule; add a source-indexed or streaming active-window schedule if repeated reads dominate.
- Next recommended step: add a full VCF/BCF automatic-window integration fixture and profile the production chr22 CSV over representative 10-kb/5-kb and larger windows before changing the bounded scheduler.

## 2026-08-27T13:42:03-04:00 — Sliding-window integration, profiling, and statistical hardening

### Baseline and goal

- Baseline commit: `48e9e020e201620bcc6446843c528902550216c7` (`Add production burden, SKAT, and SKAT-O analyses`) on `master`; this increment includes the preceding uncommitted native sliding-window work.
- Goal: complete the five requested follow-ups: finish native size/stride windows, prove CSV/VCF/BCF and restart equivalence, profile three full chr22 grids, optimize the scheduler only where profiling justified it, and harden SKAT/SKAT-O edge statistics before committing and pushing.

### Decisions and files changed

- Automatic windows now use `--window-size` and optional `--window-stride` directly; they retain the one-based chromosome-local grid, source-order membership, exact ALT orientation, unit weights, expression rows as phenotypes, and selected fixed covariates strictly as adjustment variables. Custom TSV/region definitions remain separate and mutually exclusive.
- Made automatic CSV and VCF/BCF runs converge on the existing signature-bound version-1 aligned FP64 raw cache. `QRawMatrixCache` now builds a lazy derived row-ID/offset index, persists it atomically as a version-2 sidecar bound to the raw signature/size/mtime/row count, checksums every sidecar entry, rebuilds corrupt/incompatible sidecars, coalesces adjacent selected rows, bulk-decodes big-endian FP64 values, and still validates every raw-row CRC. Existing raw-cache bytes remain unchanged.
- Replaced the repeated whole-source scan per set tile with indexed selected-row reads when the source is an aligned raw cache. Added schedule/read-volume/heap diagnostics. Trusted checksummed row buffers transfer internally without a second dosage clone; the public `QBurdenReference.Variant` constructor remains defensive.
- Added actual end-to-end CSV, VCF, and generated BCF 2.2 runs for Burden, SKAT, and SKAT-O. Results and audits are byte-identical for exactly representable DS values, and CSV checkpoint resume reproduces identical output bytes. Added sidecar reuse/corruption recovery tests.
- Added statistical fixtures for the existing singleton exact scaled-chi-square path, nearly collinear variants, extreme finite weights, forced bounded-Imhof fallback labeling, and seeded SKAT-O Monte Carlo resolution/determinism.
- Changed source under `QeQTLAnalysis`, its config/CLI, `QRawMatrixCache`, and the set-test package; expanded command/foundation/TODO documentation and tests including new `QSlidingWindowAnalysisTest`.

### Verification and profiling evidence

- Full chr22 profile command form, repeated with `(10000,5000)`, `(50000,10000)`, and `(100000,50000)` for size/stride and matching output/checkpoint names: `java -Xmx8g --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cpu --genotype D:\Research\topmed\sqtl\sqtl-batch1234\chr22-filtered-batch1234-dose-meanimp.csv --genotype-format csv --expression target\settest-smoke\traits.csv --covariates D:\Research\topmed\sqtl\sqtl-batch1234\sqtl-batch1234-mastermat-withpcs-4746.csv --fixed-covariates Age --genotype-id-column framid --expression-id-column SampleName --expression-id-strip-prefix X --analysis burden --window-size <size> --window-stride <stride> --set-block-size 256 --expression-block-rows 1 --predictor-missing mean --trait-missing error --set-degenerate skip --cache-dir target\sliding-full-cache --checkpoint-dir target\sliding-full-<grid>-final.checkpoint --output target\sliding-full-<grid>-final.csv --threshold none 0`. All used CPU OpenBLAS FP64, 428,629 variants, 4,746 aligned samples, two expression traits, and 256/1 set/trait tiles.
- Completed warm-cache profiles: 10-kb/5-kb emitted 7,453 windows in 30 blocks, wall 182.143 seconds (application overall 181.458), selected 16,345,186,032 numeric bytes, and observed 3,861,495,624 bytes heap; 50-kb/10-kb emitted 3,808 windows in 15 blocks, wall 229.098 seconds (overall 228.355), selected 16,451,762,208 bytes, and observed 6,297,342,968 bytes heap; 100-kb/50-kb emitted 765 windows in three blocks, wall 148.352 seconds (overall 148.008), selected 16,315,419,120 bytes, and observed 7,905,198,712 bytes heap. Each reused the persistent index.
- The cold/unoptimized 10-kb/5-kb run took 4,531.037 seconds and exposed repeated scans/per-value random I/O. Two intermediate 100-kb/50-kb attempts were stopped when per-row seeks and then per-double reads were clearly pathological; the first bulk-read attempt failed before association with an 8-GiB heap because the defensive constructor temporarily duplicated the tile. The final indexed, bulk, ownership-transfer path completed all grids. The aligned raw cache is 16,290,664,721 bytes and is intentionally ignored under `target`.
- `.\mvnw.cmd -q '-Dtest=QRawMatrixCacheTest,QSetTestFoundationTest,QSkatReferenceTest,QSlidingWindowAnalysisTest,QeQTLCommandLineTest' test` — focused cache/window/statistics/CLI suite passed. `.\mvnw.cmd -q '-Dtest=QSlidingWindowAnalysisTest' test` — final expanded CSV/VCF/BCF test passed after changing fractional DS values to exactly binary-representable quarters; the initial generated-BCF assertion correctly exposed ordinary BCF Float rounding for decimal 0.2/0.8/1.4 rather than a scheduler discrepancy.
- `.\mvnw.cmd -q test` — final complete suite passed: 122 tests in 29 suites, zero failures, zero errors, and five intentional runtime/platform/optional-fixture skips.
- `.\mvnw.cmd -q -DskipTests package` — final runnable jar built successfully. `target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar` is 135,882,617 bytes with SHA-256 `F0F9FA2876CB4ED596669D2659DBC2FD103ACBD58799365B4CE0FDB3670832E6`; inspection confirmed `META-INF/LICENSE`, `META-INF/LICENSE_EXCEPTION`, and all Intel oneMKL/OpenMP notices under `META-INF/licenses/`.
- `java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo` — detected NVIDIA GeForce RTX 2080, CUDA driver API 13.3, CUDA Runtime 12.6, compute capability 7.5, compiler available, FP64, and 8,589,606,912 bytes VRAM; CPU was bundled OpenBLAS 0.3.34 with 15 BLAS threads. No usable OpenCL device was reported. Full-suite available CUDA tests passed; unavailable runtime/device cases skipped cleanly.
- `git diff --check` — no whitespace errors; only informational LF-to-CRLF working-copy warnings.

### Compatibility, limitations, and next step

- Legacy positional INI, explicit set TSVs, indexed custom regions, result/audit schemas, checkpoint ordering, FP64 defaults, and version-1 raw-cache bytes remain compatible. Sidecars are derived state and rebuild automatically when absent, old, corrupt, or source-incompatible.
- Heap scales with the unique variants spanned by a resident set tile. The 100-kb/50-kb default 256-set tile approached the 8-GiB limit; use a smaller `--set-block-size` for broad/dense windows or smaller heaps. The raw cache has a substantial one-time cold-build cost and should be preserved across grids.
- Direct CUDA/OpenCL set-test products and relatedness-aware null models remain deferred. The next profiling-led extension is automatic heap-aware set-tile sizing or, only if a demonstrated compute bottleneck remains after indexed I/O, backend-neutral accelerated score/covariance products anchored to these references.
