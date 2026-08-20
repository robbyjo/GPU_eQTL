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
