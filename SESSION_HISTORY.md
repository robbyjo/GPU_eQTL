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
