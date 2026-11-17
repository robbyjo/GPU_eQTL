# GPU eQTL 1.9.1.3

This patch fixes set-test checkpoint rejection when a job restarts on another cluster node that mounts the same shared files with a different filesystem device number. It includes [PR #7](https://github.com/robbyjo/GPU_eQTL/pull/7).

Covariate, pedigree and GRM files now contribute their normalized absolute path, size, modification time and SHA-256 content digest to the saved signature. Device/inode keys remain local cache-invalidation hints. Matching input bytes and metadata produce the same persistent identity across mounts; backend, device, runtime and analysis-setting checks remain in force.

Use a **fresh set-result checkpoint directory** when upgrading: the set signature is now version 6. Preserve previous results and checkpoints for provenance. Existing aligned raw dosage caches remain compatible. This release does not migrate old result checkpoints or change statistical formulas, FP64 defaults, dependencies or GPU operations.

Validation on the reviewed source passed 199 tests: 194 passed, five optional skips, zero failures or errors. This included Java/OpenBLAS references and CUDA/OpenCL on a Windows NVIDIA GeForce RTX 2080 (driver 621.14). The new regressions cover identical-byte replacement, changed content with preserved size/time, path/time binding and byte-identical related-SKAT result/audit resume. The required hosted Linux/Java 17 check also passed. Final artifact and publication evidence is recorded in `docs/release-1.9.1.3/verification.json` after publication.

Unchanged model inputs are hashed once per analysis when usable file keys are available; metadata changes trigger rehashing. Filesystems without keys are rehashed at every check. This metadata cache does not guarantee detection of deliberate in-place edits that restore all metadata. Shared-node portability still requires unchanged paths, input bytes, size/time, settings and compatible backend/device/runtime. This release verification did not rerun a production cohort or validate fresh multi-node, macOS, other-vendor or multi-GPU execution.

The universal JAR requires Java 17 or newer and includes OpenBLAS for Windows x64, Linux x64/ARM64 and macOS x64/ARM64. CUDA/OpenCL require compatible installed runtimes. Download the JAR with `SHA256SUMS` and `RELEASE-MANIFEST.txt` from this release, and verify the checksum before use.
