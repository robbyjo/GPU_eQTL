# GPU eQTL 2.0.0-SNAPSHOT release notes

This development release keeps FP64 as the default scientific contract and adds GPU-first CUDA/OpenCL execution with a production CPU fallback (oneMKL evaluation builds, bundled OpenBLAS, or portable Java). The modern headered CSV and VCF/BCF paths enforce explicit sample alignment, covariate rank and degrees of freedom, missingness policies, signed preprocessing caches, and deterministic checkpoint assembly.

## Analysis additions

- Ordinary eQTL supports resident and bounded-RAM schedules, exact trait-pattern deletion, genotype-outer exact sufficient statistics, pattern-scoped MAF/MAC filtering, CPU/GPU fixed-effect residualization, and opt-in approximate FP32.
- Weighted burden, SKAT, and SKAT-O support explicit sets or native sliding windows with configurable window size and stride. Set tests remain FP64 and expression rows are the phenotypes.
- Cohort-aware ordinary eQTL supports cohort-specific block-diagonal fixed effects, independent cohorts, explicit aligned covariance pre-whitening, compound-symmetry repeated observations, and an audited one-row-per-subject contract.

## Operations and packaging

- Cache inspection reports types, sizes, timestamps, and SHA-256 readability. Stale pruning is dry-run by default, respects live file locks, and moves applied candidates to recoverable cache-local trash.
- Universal and single-platform OpenBLAS builds are reproducible from the source commit timestamp. Verification scripts run tests, validate packaged GPL/exception/third-party notices, smoke-test backend discovery, and emit SHA-256 and release manifests.
- Windows x64 verification on the development host detected an NVIDIA GeForce RTX 2080 through CUDA Runtime 12.6/driver API 13.3 and the bundled OpenBLAS 0.3.34 CPU engine. Linux and macOS artifacts still require smoke verification on their actual target hosts before publication.

## Known limitations

- Standard bcftools CSI is recognized but not decoded; indexed queries require tabix `.tbi` or HTSJDK Tribble `.idx` and never silently fall back to a sequential scan.
- Cohort covariance mode is not yet connected to burden/SKAT/SKAT-O, exact post-whitening trait-pattern deletion, heterogeneity statistics, or leave-one-cohort-out diagnostics.
- FP32 is opt-in and approximate. Genotype-outer exact trait-pattern analysis and set tests remain FP64.
- Direct GPU burden/SKAT products, relatedness-aware set tests, categorical-SNP repair, SNP interactions, and forward selection remain future work.

See [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for the host/runtime matrix and publication gates.
