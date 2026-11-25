# GPU eQTL

Java software for eQTL and related QTL analyses, including rare-variant burden, SKAT, and SKAT-O tests. It runs on CUDA, OpenCL, or CPU backends, with double precision (FP64) by default. **Java 17 or newer is required; a GPU is optional.**

## Start here

- **[Project website](https://robbyjo.github.io/GPU_eQTL/)** — an overview, workflow, and interactive command examples.
- **[Download the latest release](https://github.com/robbyjo/GPU_eQTL/releases/latest)** — the runnable universal JAR, checksum, and build manifest.
- **[Install and run](docs/USER_GUIDE.md#download-and-first-run)** — first steps, backend setup, and input validation.
- **[Updates and upgrade notes](https://robbyjo.github.io/GPU_eQTL/updates.html)** — what changed and what to check before upgrading.

## Find the right guide

| I want to… | Go to |
| --- | --- |
| Prepare inputs, handle missing values, use caches, or resume a run | [User guide](docs/USER_GUIDE.md) |
| Look up an option, INI setting, or complete command | [Command-line reference](docs/COMMAND_LINE_REFERENCE.md) |
| Run burden, SKAT, or SKAT-O | [Set-test guide](docs/SET_TEST_FOUNDATION.md) |
| Account for relatedness in set tests | [Related set tests](docs/RELATED_SET_TESTS.md) · [Spectral REML](docs/SPECTRAL_SET_NULL.md) |
| Tune set-test eigenvalues and CPU/GPU overlap | [Set-test performance](docs/SET_TEST_PERFORMANCE.md) |
| Fit single-variant mixed models | [Mixed-model guide](docs/MIXED_MODELS.md) |
| Build from source or contribute | [Build instructions](docs/USER_GUIDE.md#build-and-test) · [Contributor guide](AGENTS.md) |

For bugs and questions, [open a GitHub issue](https://github.com/robbyjo/GPU_eQTL/issues). Planned work is tracked in [TODO.md](TODO.md).

## License

GNU GPLv3, with a narrow [oneMKL linking permission](LICENSE_EXCEPTION). See [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY_LICENSES).
