# GPU eQTL

GPU-accelerated Java software for eQTL and related QTL analyses, with a portable CPU fallback. This codebase originated in 2011–2013 and is being modernized incrementally. FP64 remains the default scientific mode; FP32 is available as an explicit performance/compatibility choice.

## Current status

The application targets Java 17 and has selectable compute backends behind a vendor-neutral API:

| Setting | Implementation | Intended hardware |
| --- | --- | --- |
| `auto` (default) | CUDA-first discovery, OpenCL fallback, then CPU fallback | Mixed machines; avoids running the same NVIDIA card twice |
| `cuda` | JCuda 12.6.0 and cuBLAS DGEMM/SGEMM | NVIDIA GPUs |
| `opencl` | JOCL 2.0.6 and the production OpenCL C kernel | NVIDIA, AMD, or Intel GPUs exposed by an OpenCL driver |
| `cpu` | oneMKL 2026.1 or OpenBLAS 0.3.34 through JavaCPP 1.5.14, with a pure-Java fallback | Systems without a usable GPU, and reproducibility/diagnostic runs |

Automatic discovery uses CUDA for a usable NVIDIA device and also includes distinct OpenCL devices from other vendors. It filters those devices for the requested precision. If no eligible GPU remains, it selects one CPU context and prints a conspicuous performance warning. CPU is never mixed into an otherwise usable multi-GPU run. If CUDA cannot initialize, NVIDIA OpenCL remains available as a fallback. A native HIP/ROCm backend is not included yet; AMD cards currently use JOCL/OpenCL.

The real-valued eQTL calculation defaults to FP64 and is validated through CUDA, OpenCL, native oneMKL/OpenBLAS, and portable Java engines against the same scalar reference. Fixed-effect QR decomposition and rank validation remain on the CPU. By default, headered-CSV analyses apply the large block projection `Y - (Y Q) Q^T` on the selected compute backend; CUDA and the native CPU engines use two BLAS multiplications, while OpenCL uses two projection kernels. `--residualization cpu` retains the prior scalar FP64 projector. The legacy categorical-SNP path remains disabled; the CUDA and CPU backends currently implement `eqtlReal` only.

`precision = fp32` or `--precision fp32` uses SGEMM/a `float` OpenCL kernel and, when GPU residualization is selected, FP32 projection. Standardization, prepared-cache storage, effects, test statistics, and p-values remain FP64. Selecting CPU residualization keeps covariate adjustment FP64 even when the association product is FP32. This permits OpenCL execution on Intel Iris Xe devices that lack native FP64, provided their installed driver exposes a usable OpenCL GPU.

The legacy INI format remains supported, and every analysis setting now also has a command-line form. Headered CSV matrices are validated before computation: blank or duplicate row/sample IDs are fatal, genotype and expression samples are explicitly reordered, and a covariate table may bridge different ID namespaces. Text covariates are automatically one-hot encoded using the lexicographically first level as the reference; numeric-looking factors can be forced with `--factor-covariates`. Rank-deficient covariate models stop before analysis.

Headered CSV genotype and expression matrices can also be processed in bounded-RAM blocks. The first run creates indexed FP64 cache files after validation, reordering, covariate residualization, and standardization. GPU projection broadcasts the small Q matrix once to each context and prepares independent input blocks concurrently across selected devices while writing cache rows in source order. Later runs with the same source metadata, sample order, covariate projection, and preprocessing mode reuse those prepared rows instead of reparsing or reprocessing the CSV. Each cache row has an integrity checksum. SNP interactions, forward selection, and indexed VCF/BCF input remain on the roadmap.

## Build and test

Install a Java 17 or newer JDK. Maven itself does not have to be installed because the repository includes a checksum-pinned wrapper.

On Windows:

```powershell
.\mvnw.cmd clean package
```

On Linux or macOS:

```sh
./mvnw clean package
```

The runnable jar is `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar`. It includes the Java dependencies, JCuda/JOCL bindings, and OpenBLAS natives for Windows x64, Linux x64/ARM64, and macOS x64/ARM64. Users on those desktop platforms can copy the one jar and run it without installing BLAS. It does not bundle GPU drivers, CUDA, cuBLAS, or an OpenCL ICD. CUDA therefore still needs a compatible NVIDIA driver plus CUDA runtime/cuBLAS installation, and OpenCL needs the vendor's ICD. On an unbundled platform, CPU `auto` falls back to the slower pure-Java engine.

Hardware-independent CPU tests always run. CUDA and OpenCL FP64/FP32 numerical tests skip cleanly when their respective runtime/device is unavailable. A platform-gated test also verifies that the matching bundled OpenBLAS native loads on each supported desktop target.

### Optional oneMKL platform builds

The ordinary build remains the portable OpenBLAS distribution. Separate profiles create oneMKL-enabled jars for x64 Windows or Linux without modifying that default artifact:

```powershell
.\mvnw.cmd -Pcpu-mkl-windows-x86_64 package
java '-Deqtl.cpu.blas=mkl' --enable-native-access=ALL-UNNAMED `
  -jar target-mkl-windows-x86_64\gpu-eqtl-2.0.0-SNAPSHOT-mkl-windows-x86_64-all.jar `
  --backend cpu --printbackendinfo
```

```sh
./mvnw -Pcpu-mkl-linux-x86_64 package
java -Deqtl.cpu.blas=mkl --enable-native-access=ALL-UNNAMED \
  -jar target-mkl-linux-x86_64/gpu-eqtl-2.0.0-SNAPSHOT-mkl-linux-x86_64-all.jar \
  --backend cpu --printbackendinfo
```

Each profile includes only its matching oneMKL and OpenBLAS x64 natives. `auto` prefers oneMKL when its native runtime is present, then OpenBLAS, then portable Java. Intel permits binary oneMKL redistribution under the Intel Simplified Software License. GPU eQTL remains GPLv3 and adds a narrow GPLv3 Section 7 permission for combining Roby Joehanes-owned portions with oneMKL; see [`LICENSE_EXCEPTION`](LICENSE_EXCEPTION). Intel's terms remain in force for its binaries. The exact oneMKL 2026.1 license, oneMKL incorporated-code notices, and Intel OpenMP notices are reproduced in [`THIRD_PARTY_LICENSES`](THIRD_PARTY_LICENSES). All shaded jars package the GPL, exception, and vendor notices under `META-INF`. See Intel's [oneMKL license FAQ](https://www.intel.com/content/www/us/en/developer/articles/tool/onemkl-license-faq.html) and the GNU project's [GPL-incompatible library guidance](https://www.gnu.org/licenses/gpl-faq.html#GPLIncompatibleLibs).

## Eclipse

Use an Eclipse installation with Maven Integration for Eclipse (m2e) and a Java 17 or newer JDK. Import this directory with **File → Import → Maven → Existing Maven Projects**.

If the project was already open before the Maven metadata was added, select the project and run **Maven → Update Project…** (`Alt+F5`), enable **Force Update of Snapshots/Releases**, and then run **Project → Clean…**. Eclipse should show **Maven Dependencies** in the project tree and use the `JavaSE-17` execution environment.

The custom multi-delimiter JavaCSV jar remains attached directly with its source archive; JCuda, JOCL, OpenBLAS/oneMKL JavaCPP bindings, JDistlib, Commons Compress, and JUnit are supplied by Maven. The oneMKL native runtime is present only in an explicitly selected evaluation profile.

## Inspect and select compute backends

Automatic selection is the default:

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo
```

Force a backend with the application argument:

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cuda --printbackendinfo
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend opencl --printbackendinfo
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --backend cpu --printbackendinfo
```

The equivalent system property remains supported in PowerShell (quote the `-D` argument):

```powershell
java '-Deqtl.gpu.backend=cuda' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo
java '-Deqtl.gpu.backend=opencl' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo
java '-Deqtl.gpu.backend=cpu' -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo
```

In POSIX shells, the same property does not need PowerShell quoting:

```sh
java -Deqtl.gpu.backend=cpu -jar target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar --printbackendinfo
```

`--gpu-backend` and `--printgpuinfo` remain compatibility aliases. CPU `auto` tries oneMKL when its profile natives are present, then bundled OpenBLAS, then portable Java. Force its behavior for diagnostics with `-Deqtl.cpu.blas=mkl`, `openblas`, or `java`; explicitly selected native engines fail rather than silently changing implementations. Native BLAS uses `max(1, logical processors - 1)` threads by default; override that independent inner pool with `-Deqtl.cpu.threads=N`.

Recent Java runtimes may print a warning when a JNI binding loads. Add `--enable-native-access=ALL-UNNAMED` before `-jar` to suppress it. Java 17 remains the compilation target.

A manual, transfer-inclusive comparison is available after the build. Its results are indicative rather than a substitute for benchmarking a representative analysis:

```powershell
java --enable-native-access=ALL-UNNAMED -cp 'target\test-classes;target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar' gov.nih.gpu.GpuBackendBenchmark
```

## Precision, automatic sizing, and multiple GPUs

FP64 is the default. Select FP32 explicitly in an INI file or on the command line:

```ini
precision = fp32
```

```powershell
  --precision fp32
```

Fixed-effect residualization is selected independently:

```ini
residualization = auto
```

```powershell
  --residualization {auto|gpu|cpu}
```

`auto` and `gpu` use the selected compute contexts for headered-CSV preparation; this includes OpenBLAS when the analysis backend is CPU. Residualization mode `cpu` retains the previous scalar FP64 projector and its cache signatures. QR decomposition, pivot/rank decisions, and degrees-of-freedom accounting are identical in all three modes. The full `N x N` projection matrix is never formed.

FP32 changes results slightly and should not be substituted silently for an established FP64 pipeline. On the representative 128-SNP by 128-trait WHI subset (2,005 samples; 16,384 associations), CUDA FP32 versus FP64 had maximum absolute differences of `5.59e-9` in R-squared, `1.03e-8` in effect, `9.70e-7` in t, and `2.51e-6` in log10 p. Both modes selected the same 11 associations at p <= `1e-4`. A larger 4,096-SNP by 8,192-trait study (33,554,432 associations) found one FP32-only result at that reporting boundary; among the 48,948 common reported results, maximum absolute differences were `9.73e-7` in R-squared, `5.71e-7` in effect, `1.11e-4` in t, and `9.90e-4` in log10 p.

The same larger study directly compared the new projector. GPU FP64 and CPU FP64 residualization selected the same 48,948 associations; maximum absolute differences were `9.99e-16` in R-squared, `3.77e-15` in effect, `1.25e-13` in t, and `7.96e-13` in log10 p. Within FP32, GPU versus CPU projection selected the same 48,949 associations; maximum differences were `4.24e-7`, `1.84e-6`, `6.89e-5`, and `4.08e-4`, respectively. These are bounded accuracy studies, not guarantees for every matrix or threshold; revalidate important borderline associations in FP64.

When `block_size`/`--block-size` is omitted or zero, the application reads total VRAM and maximum-allocation limits from every selected GPU. It chooses one block size supported by the least-capable device, accounting for both input tiles, the quadratic result tile, FP32 versus FP64 element size, VRAM headroom, the per-allocation limit, and Java's array-size limit. It also caps the target result allocation at 1 GiB because consuming all available VRAM with one quadratic tile is usually slower and much harder on host memory. When the genotype matrix is smaller than the memory limit, the recommendation retains up to four genotype jobs per selected GPU so one oversized tile does not eliminate packing/post-processing overlap. The choice and its limiting reason are printed before computation; an explicit block size remains an override.

When `num_threads`/`--threads` is omitted or zero, bounded-RAM mode uses up to four workers per GPU to overlap cache reads, packing, GPU execution, and result processing; full-memory mode uses up to two. The bounded-RAM recommendation estimates prepared inputs, packed precision-specific inputs, and the result array per worker, retains JVM-heap headroom, is capped by required genotype jobs, and leaves one CPU core free. Explicit values remain supported, with warnings for idle GPUs, CPU oversubscription, or a streamed configuration estimated to consume more than 75% of the available JVM heap.

For CPU analysis, automatic block sizing uses available JVM heap and a 256-MiB result-tile target rather than pretending heap is VRAM. The application uses one exclusive CPU context: full-memory runs use one pipeline worker, and streamed runs may use two so cache preparation/result handling can overlap the native multiplication. oneMKL/OpenBLAS has its own inner thread pool as described above. Setting many application `--threads` does not create concurrent BLAS products and usually only adds memory pressure.

All distinct usable GPUs are opened, with one exclusive context per device. Association jobs and GPU cache preparation reserve contexts from shared exclusive pools, so a multi-GPU system uses multiple devices concurrently when JVM heap can safely hold their in-flight preparation blocks. Each device keeps one projection copy while preparation is active, and those temporary buffers are released before association tiles are allocated. Automatic discovery suppresses only the duplicate OpenCL representation of an NVIDIA device already selected through CUDA. On mixed-capacity systems, automatic block sizing follows the smallest selected GPU. Use `--backend cuda`, `opencl`, or `cpu` when you intentionally want only one backend family.

The old `lambda` INI setting no longer controls memory or result storage and is ignored with a compatibility notice. It should be removed from new configurations.

## Run an analysis

With automatic backend selection:

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar path\to\analysis.ini
```

The equivalent migration form is `--config`; command-line options override values read from the INI file:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --gpu-backend cuda `
  --threads 4
```

A run no longer needs an INI file:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype genotype.csv `
  --expression expression.csv `
  --covariates covariates.csv `
  --output results.csv `
  --fixed-covariates Age,Batch,PC1,PC2 `
  --threshold pval 1e-4 `
  --precision fp64 `
  --residualization auto
```

Omitting `--block-size` and `--threads` in this example enables the device-aware recommendations described above.

### Matrix roles and missing values

The same engine may analyze genotype, expression, methylation, proteomics, or another continuous matrix. Declare the roles instead of asking the program to infer them from the values:

```powershell
  --predictor-type genotype --trait-type expression
```

Accepted types are `genotype`, `expression`, `methylation`, `proteomics`, and `continuous`. The historical `--genotype`/`--expression` names and INI files retain genotype/expression defaults for compatibility. The generic `--predictor` and `--traits` spellings require `--predictor-type` and `--trait-type`, respectively. The declared type is reported in the missingness audit and controls whether genotype-only local-pattern imputation is allowed. It does not inspect a large matrix to guess its biological meaning.

Blank fields and the tokens `NA`, `N/A`, `null`, `NaN`, and `.` are missing in delimited matrices. Every modern headered-matrix run writes an atomic tab-separated missingness audit to `<output>.missingness.tsv`; change it with `--missingness-qc-output FILE`. The report includes matrix summaries, exact trait-row missingness patterns, affected row IDs, per-sample counts, and selected-covariate counts. Use `--inspect-missingness` to write this report and stop before GPU initialization.

The defaults reflect the common case of scant genotype missingness and expression-side missingness:

- Predictor missing values use row-mean imputation (`--predictor-missing mean`). The mean is computed over the samples actually used for the trait missingness pattern, rather than copied from the larger aligned cohort. For genotype dosage this is twice the pattern-specific EAF; it is also a defined, conservative policy for non-genotype continuous predictors.
- Trait missing values use exact pattern-wise deletion (`--trait-missing pattern`). Trait rows with the same complete-sample mask are residualized and analyzed together. This forces bounded-RAM processing, runs one compute pass per distinct pattern, reports the effective `N` and p-value `DF` on every result, and can still be slow when there are many patterns. For VCF/BCF, the aligned raw FP64 dosage rows are decoded once into a checksummed disk cache; exact pattern-specific `N`, dosage sum, squared-dosage sum, EAF/MAF/MAC, imputation mean, and variable-row decision are then cached by sample mask. Prepared predictor/trait caches are also persistent. Pattern mode does not yet resume partially completed association scheduling, and its output is grouped by missingness pattern rather than global trait-row order.
- Samples missing any selected fixed covariate are removed once from all three aligned inputs (`--covariate-missing complete-samples`). Use `error` to require selected covariates to be complete.

Predictor and trait alternatives are `error`, `mean`, `zero`, and `exclude-row`. `exclude-row` removes an entire predictor or trait row containing any missing value. Trait `pattern` performs dynamic complete-case selection; predictor `pattern` is intentionally rejected because it would create a separate sample set for every predictor-trait combination.

For sparse missing genotypes, `--predictor-missing local-pattern --predictor-flanks 1` enables a bounded-memory local proxy. It compares the sample's observed dosages at the requested number of variants on each side, finds called donors with the nearest flanking dosage pattern, and averages tied nearest donor dosages. It never averages the flanking SNP dosages themselves. `CHROM:POS...` row identifiers are validated for contiguous, nondecreasing genomic order and chromosome boundaries are respected. If location annotations are absent, the program explicitly warns that it is trusting input row order; if no comparable donor exists it warns and falls back to the row mean. This is inexpensive nearest-pattern imputation, not phasing, LD modeling, or reference-panel imputation, so important analyses should compare it with mean or a standard externally imputed genotype source.

### VCF.gz and BCF genotype input

Genotypes may be read from plain VCF, gzip/BGZF-compressed VCF, or BCF 2.1/2.2 while expression remains a headered delimited matrix. Variant input uses the same metadata/block interface, sample-ID alignment, prepared cache, checkpoint, full-memory path, and bounded-RAM path as CSV. For example:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype cohort.vcf.gz `
  --genotype-format vcf `
  --genotype-field auto `
  --genotype-missing mean `
  --min-maf 0.01 `
  --min-mac 5 `
  --variant-qc-output cohort.variant-qc.tsv `
  --variant-qc-checkpoint cohort.variant-qc.checkpoint `
  --variant-qc-threads 0 `
  --expression expression.csv `
  --covariates covariates.csv `
  --fixed-covariates Age,Batch,PC1,PC2 `
  --output results.csv `
  --genotype-block-rows 10240
```

`--genotype-format auto` infers VCF/BCF from `.vcf`, `.vcf.gz`, and `.bcf`; otherwise it retains CSV behavior. `--genotype-field auto` chooses header-declared `DS` first and falls back to `GT`. Force either field with `--genotype-field DS` or `GT`. VCF/BCF calls remain missing through the common QC scan and are then handled by `--predictor-missing`; the compatibility spelling `--genotype-missing` selects the same policy. Predictor mean imputation is the default. Use `error` to require complete calls, `exclude-row`/`exclude-variant` to remove the variant, `zero` for an explicit zero fill, or `local-pattern` for the documented flanking proxy.

The current variant model is biallelic, diploid, and additive. Multiallelic records are excluded by default and annotated as such; use `--multiallelic error` to make one fatal. Monomorphic variants are always detected, reported, and excluded because their standardized association row is undefined. `--min-maf` accepts values from 0 through 0.5, and `--min-mac` accepts a non-negative value; when both are present a variant must pass both. MAC may be fractional for imputed dosages. Singletons and doubletons are identified only when the computed MAC is within numerical tolerance of one or two.

By default, `--frequency-scope aligned` applies MAF/MAC once to the final genotype-expression-covariate aligned cohort; HWE is always cohort-level QC. Exact trait patterns still recompute their genotype mean, variance, EAF/MAF/MAC, and monomorphic status, and skip only variant-pattern combinations that become constant. Use `--frequency-scope pattern` with `--trait-missing pattern` to apply `--min-maf`/`--min-mac` independently after each trait sample mask. In that mode the cohort QC report retains frequency-filter candidates and records why they would have failed the aligned filter. The compact `<output>.pattern-variant-qc.tsv` summary records effective N and included/monomorphic/MAF/MAC counts per pattern. This explicit mode can test different variant sets for different traits, so use it only when that is the intended scientific contract.

Every VCF/BCF analysis performs an aligned-sample QC scan and writes a tab-separated variant report, defaulting to `<output>.variants.tsv`. `--variant-qc-output FILE` changes its location. It contains the canonical association identifier `CHROM:POS:REF:ALT`, rs ID, REF, ALT, original VCF FILTER, selected field, called/missing counts, effect-allele count and allele number, squared-dosage sum, EAF, MAF, MAC, exact biallelic HWE p-value, diploid hom-ref/heterozygous/hom-alt counts, classification, inclusion status, exclusion reason, region-set membership, frequency scope, and any aligned-cohort frequency reason. These sufficient statistics are persisted in the resumable checkpoint. HWE uses available diploid `GT` calls even when `DS` supplies the association dosage; it is `NA` when no usable GT calls exist. Original VCF FILTER values are annotated but are not an implicit analysis filter.

QC is parallelized across variants after sample alignment. `--variant-qc-threads 0` (the default) automatically selects a bounded worker count, normally leaving a CPU for input decoding and capping the pool at 16; multi-core machines use at least two variant workers. Set a positive value to override it, or `1` for a sequential comparison. The INI equivalent is `variant_qc_threads`. VCF/BCF decoding remains on the reader thread because the HTSJDK codecs are stateful, while dosage extraction, aligned-sample allele counts, MAF/MAC classification, and exact HWE run in the worker pool. Completed records are consumed in source order, so the QC report, duplicate detection, retained-variant order, and later association identifiers are deterministic. Matrix-header samples excluded by alignment are not inspected by dosage/QC calculations. The first scan retains only a compact inclusion/EAF decision per input record; later source passes reuse those decisions rather than recomputing frequency and HWE.

The aligned-sample QC scan is resumable. Ordered batches are committed atomically under `<variant-qc-output>.checkpoint/<signature>/` by default; `--variant-qc-checkpoint DIR` (INI: `variant_qc_checkpoint`) changes the root. The signature includes the normalized genotype path, source and index size/modification time, resolved VCF/BCF field and filtering policies, normalized regions, complete header sample list, and final aligned sample indices/order. A policy, file/index metadata, region, or alignment change therefore starts a separate state directory instead of silently reusing obsolete EAF/MAF/MAC/HWE decisions. Do not edit checkpoint parts. A process lock rejects simultaneous writers to the same signature. With an index, an interrupted scan seeks directly to the last durable canonical variant and verifies that boundary before continuing. Without an index it decodes and verifies the completed prefix while skipping its per-sample frequency/HWE work. A fully completed matching checkpoint is loaded directly and can recreate a missing QC TSV without opening a full variant-record pass. Checkpoints are retained intentionally for later reruns; remove a signature directory only when its analysis will not be resumed.

#### Indexed regions and sets

BGZF VCF plus tabix (`.tbi`) and BCF/VCF plus an HTSJDK Tribble index (`.idx`) support true interval queries. A neighboring index is detected automatically; use `--variant-index FILE` to specify it. Standard bcftools CSI (`.csi`) is recognized as an index candidate but is not yet decoded by the current pure-Java HTSJDK reader; region queries with an unsupported CSI fail before analysis rather than silently scanning the whole file. Plain gzip VCF remains valid for unrestricted sequential analysis but cannot service a region request.

Inline regions are one-based and inclusive, repeatable, and may have a set ID:

```powershell
  --region APOE=chr19:44900000-45000000 `
  --region nearby=chr19:44800000-45100000
```

`--regions-file FILE` accepts tab-separated `CHROM START END` or `SET_ID CHROM START END` rows, with an optional matching header. The default is one-based inclusive; `--region-coordinates bed` interprets file rows as zero-based half-open BED coordinates. Exact contig names are preferred; an unambiguous `chr` prefix or mitochondrial `M`/`MT` alias is normalized to the source. Query intervals are sorted in source-contig order and merged to prevent duplicate variant rows, while the QC report preserves all overlapping set IDs in definition order. Duplicate definitions are harmless, one variant may belong to several sets, and empty sets are warned and listed in the run summary. Canonical association IDs and within-source variant order do not change.

Large variant passes report progress approximately every 15 seconds (or every million input records, whichever comes first). The aligned-sample QC pass reports records scanned, variants retained, elapsed time, and throughput; its total is not known until that first pass ends. Later sequential missingness and cache-building rereads also report percentage and ETA using the input-record count learned during QC. These counts describe decoded VCF/BCF records rather than compressed bytes. The final QC report is assembled atomically from the durable checkpoint parts only after a complete scan.

Variant QC and aligned-scope frequency filters are computed over the final aligned analysis samples after validating unique header IDs and applying selected-covariate complete-sample removal. This is important under `covariate-subset`: VCF-only samples cannot make an analysis-set monomorphic or rare variant pass MAF/MAC filtering. The QC file's called/missing counts, EAF, MAF, MAC, HWE, singleton/doubleton classification, and monomorphic exclusion all describe that aligned set; the console reports its sample count. Strict alignment (the default) requires genotype, trait, and covariate sample sets to agree exactly. An explicitly requested covariate-subset alignment may instead select the canonical covariate rows from a larger genotype or trait header; all covariate samples must still be present and matrix-only extras are counted in the unified missingness/alignment audit. Association output uses the canonical variant identifier so missing or duplicated rs IDs cannot collide. Ordinary bounded-RAM analysis rereads source blocks without materializing the full genotype matrix. Exact trait-pattern VCF/BCF analysis instead creates a persistent aligned raw FP64 cache so the compressed source is decoded only once; budget approximately `8 * aligned_samples * retained_variants` bytes plus row IDs/checksums.

If genotype and expression headers use different identifiers, the program searches for covariate columns whose unique values exactly match each header. Ambiguous cases must be resolved explicitly:

```powershell
  --genotype-id-column NWDID --expression-id-column TORID
```

Use `--predictor-id-strip-prefix X` or `--trait-id-strip-prefix X` to remove a literal leading prefix before matching IDs. The prefix is removed only where it is present; IDs that do not start with it are unchanged. Blank IDs or collisions created by the transformation are fatal. The compatibility aliases are `--genotype-id-strip-prefix` and `--expression-id-strip-prefix`, and the INI keys are `predictor_id_strip_prefix` and `trait_id_strip_prefix`.

When a VCF contains a superset of the samples represented by the covariate rows, opt in explicitly:

```powershell
  --sample-alignment covariate-subset `
  --genotype-id-column framid `
  --expression-id-column SampleName `
  --trait-id-strip-prefix X
```

`strict` remains the default. `covariate-subset` never computes a silent intersection: every covariate ID must occur exactly once in each normalized matrix header. Matrix-only exclusions, transformed prefix counts, bridge columns, and reorder counts are written as `ALIGNMENT` records in the missingness QC file.

Use `--validate-only` to scan all row/sample IDs, resolve alignment, encode covariates, check model rank, and report degrees of freedom without running the association analysis.

### Bounded-RAM matrix mode

Set either block-row option to enable the streamed path; an omitted block-row value falls back to `--block-size`. Values are rounded up to the 16-row GPU tile boundary.

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --genotype-block-rows 10240 `
  --expression-block-rows 10240
```

The same settings can be added to an INI file as `genotype_block_rows` and `expression_block_rows`. ID bridge columns are `genotype_id_column` and `expression_id_column`. Plain CSV, `.csv.gz`, and `.csv.bz2` matrix sources are supported by the delimited metadata/block reader; VCF/VCF.gz/BCF uses the variant reader described above. CSV headers and a first-column row identifier are required for validation and streaming; legacy headerless genotype CSV files continue through the old full-memory loader.

By default, caches are placed in `.gpu-eqtl-cache` beside the output file. Choose another location, preferably fast local SSD storage, with:

```powershell
  --cache-dir D:\eqtl-cache
```

Use `--rebuild-cache` after an intentional forced rebuild. Cache signatures include the input path, size, modification time, row/sample metadata, column permutation, covariate projection, and GPU projection precision/backend family when applicable. A changed signature creates a different cache rather than silently reusing incompatible prepared values. Switching an existing run from CPU to automatic GPU residualization therefore creates a new cache; use `--residualization cpu` to reuse the prior CPU-prepared cache. Old cache versions are not automatically deleted. For the complete WHI chromosome example, allow roughly 6–7 GB of additional uncompressed cache storage.

Prepared rows are encoded and decoded in checksummed bulk records. This preserves the existing cache format and per-row CRC validation while avoiding a file operation and checksum update for every individual double. Caches created by the earlier version remain reusable.

Prepared traits can remain resident in heap memory while genotype blocks stream:

```powershell
  --genotype-block-rows 10240 `
  --expression-block-rows 10240 `
  --trait-cache memory
```

This does not residualize the trait matrix again. The source is aligned, residualized, standardized, and checksummed once while its reusable FP64 prepared cache is built; memory mode then reads that prepared cache once and retains its row blocks for every genotype block. GPU association still consumes bounded trait tiles. `--trait-cache auto` is the default and chooses memory only when the JVM heap can hold the complete prepared traits plus all configured worker buffers and conservative headroom; otherwise it reports the estimate and uses disk. `--trait-cache disk` always uses indexed disk reads. Explicit `memory` fails before association with `-Xmx`/block/thread guidance when the estimate is unsafe. INI files use `trait_cache = auto|memory|disk`.

For 114,406 traits by 4,746 samples, the numeric FP64 values alone require about 4.05 GiB, before Java row/identifier overhead and active worker buffers. Therefore “in memory” means once-prepared, chunk-built host residency—not one monolithic GPU projection or one giant GPU allocation. With exact trait-pattern deletion, each missingness pattern has a different sample set and covariate projection, so each pattern must still be prepared separately; within a pattern, its prepared trait rows are residualized only once and may use the same memory policy.

### Performance profiling

Add `--profile` to print phase totals for metadata/alignment, aligned-sample variant QC, cache creation and reads, accelerated residualization setup/upload/compute/download, host packing, association context wait, backend buffer setup/upload/compute/download, CPU statistics/result writing, and final assembly. The stable CSV phase names retain their historical `gpu_` prefix even in CPU mode; host-only CPU executions correctly report zero upload/download bytes and time. Variant QC is a subphase of metadata/alignment; worker-phase totals can overlap and therefore need not sum to wall time. Add `--profile-output FILE` to write the same measurements as CSV; specifying a profile output also enables profiling.

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --profile `
  --profile-output run-profile.csv
```

The INI equivalents are `profile = true` and `profile_output = run-profile.csv`. GPU profiling adds synchronization needed to separate upload, compute, and download time, so use the reported breakdown diagnostically and compare complete wall times for throughput decisions.

### Checkpoint and resume

Bounded-RAM analyses write one atomic result part per genotype block. Only a completely written part is considered finished, and final output is assembled in genotype-block order. On success the temporary checkpoint directory is removed. If execution stops, rerun the identical command with `--resume`:

```powershell
java --enable-native-access=ALL-UNNAMED -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --config path\to\analysis.ini `
  --genotype-block-rows 10240 `
  --expression-block-rows 10240 `
  --resume
```

The default directory is `<output-file>.checkpoint`. `--checkpoint-dir DIR` selects a different location, and `--keep-checkpoints` retains completed parts after successful assembly. Resume is rejected when input/cache signatures, precision, thresholds, degrees of freedom, block sizes, or output modes do not match. The corresponding INI keys are `cache_dir`, `rebuild_cache`, `checkpoint_dir`, `resume`, and `keep_checkpoints`.

Run `--help` for the complete argument list. See `TODO.md` for the remaining modernization work. The former `library_path` key is no longer used; runtime discovery is handled by the CUDA/OpenCL backends.

See `AGENTS.md` for contributor constraints and the ordered roadmap. See `SESSION_HISTORY.md` for timestamped changes, verification, and known limitations.

## License

GPU eQTL is distributed under the GNU General Public License version 3. See [`LICENSE`](LICENSE) for the complete license text. The copyright holder grants the narrow oneMKL linking permission in [`LICENSE_EXCEPTION`](LICENSE_EXCEPTION) for Roby Joehanes-owned portions of the program.

The default release does not contain Intel oneMKL native binaries. The platform-specific oneMKL artifacts also contain Intel-licensed components governed by the reproduced [`Intel and incorporated-code terms`](THIRD_PARTY_LICENSES). The exception permits conveying the combined work; it does not relicense or permit modification of Intel's binaries.
