# Release verification checklist

GPU eQTL requires Java 17 or newer and retains FP64 as the scientific default. A release is not complete until the shaded jar has passed the hardware-independent Maven suite, its packaged notices have been checked, and the artifact has printed backend information on its target host.

## Artifacts

Build and verify the universal OpenBLAS jar on Windows PowerShell:

```powershell
.\scripts\verify-release.ps1 -Platform universal
```

Build a classifier-specific jar by selecting one of `windows-x86_64`, `linux-x86_64`, `linux-arm64`, `macosx-x86_64`, or `macosx-arm64`:

```powershell
.\scripts\verify-release.ps1 -Platform windows-x86_64
```

On Linux or macOS, run `./scripts/verify-release.sh linux-x86_64` (substitute the matching platform). Each invocation writes the jar, `SHA256SUMS`, and `RELEASE-MANIFEST.txt` below `target/release/PLATFORM`. The manifest records the source commit, Maven profile, Java runtime, build host, reproducible output timestamp, and verification gates. `target/` remains untracked.

Release verification requires a clean worktree so the manifest commit identifies the packaged source exactly. PowerShell `-AllowDirty` or POSIX `ALLOW_DIRTY=1` exists only for development validation and records that state in the manifest; do not publish such an artifact.

The universal artifact contains OpenBLAS 0.3.34 natives for Windows x64, Linux x64/ARM64, and macOS x64/ARM64. oneMKL evaluation artifacts remain separate under the existing `cpu-mkl-windows-x86_64` and `cpu-mkl-linux-x86_64` profiles because they carry Intel runtime licensing in addition to the project GPLv3 linking exception.

## Runtime matrix

| Host | Portable Java | OpenBLAS | oneMKL | CUDA/OpenCL |
|---|---:|---:|---:|---|
| Windows x64 | Supported | Bundled | Separate evaluation profile | Requires compatible installed driver/runtime |
| Linux x64 | Supported | Bundled | Separate evaluation profile | Requires compatible installed driver/runtime |
| Linux ARM64 | Supported | Bundled | Not packaged | OpenCL/CUDA require host validation |
| macOS x64/ARM64 | Supported | Bundled | Not packaged | OpenCL depends on the host ICD; CUDA is unavailable on current macOS |

Explicit native-engine selection is fatal if that engine cannot load. `--backend auto` remains GPU-first and falls back to the CPU backend only when no eligible GPU exists. `eqtl.cpu.blas=auto|mkl|openblas|java` selects the CPU engine; release smoke tests must record the selected engine printed by `--printbackendinfo`.

## Gates before publishing

1. Run the full Maven suite on Java 17 or newer.
2. Run the verification script on each published target host; do not infer Linux or macOS native loading from a Windows build.
3. Confirm `META-INF/LICENSE`, `META-INF/LICENSE_EXCEPTION`, and every notice under `META-INF/licenses/` are present.
4. Run deterministic FP64 resident, streamed, partial-tile, checkpoint/resume, exact-pattern, and set-test fixtures.
5. On hardware represented in the release notes, record GPU model, driver/runtime, precision, backend, and accelerated-versus-scalar tolerance results. Hardware not exercised is documented as unvalidated.
6. Publish the generated SHA-256 file beside the immutable artifact and retain the release manifest.
