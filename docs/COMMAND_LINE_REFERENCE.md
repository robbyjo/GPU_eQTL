# GPU eQTL command-line reference

This page is the reference for the command-line interface in GPU eQTL 2.0. It lists every option accepted by the current parser, its default, the corresponding legacy INI key where one exists, and important interactions. For background on input formats, missingness, caches, numerical precision, and backends, see the [main README](../README.md).

## Invocation

The runnable artifact is `target/gpu-eqtl-2.0.0-SNAPSHOT-all.jar`. On recent Java versions, put `--enable-native-access=ALL-UNNAMED` before `-jar` to suppress native-library access warnings.

```text
java [JVM options] -jar gpu-eqtl-2.0.0-SNAPSHOT-all.jar legacy.ini
java [JVM options] -jar gpu-eqtl-2.0.0-SNAPSHOT-all.jar --config legacy.ini [overrides]
java [JVM options] -jar gpu-eqtl-2.0.0-SNAPSHOT-all.jar --genotype FILE --expression FILE [options]
java [JVM options] -jar gpu-eqtl-2.0.0-SNAPSHOT-all.jar --predictor FILE --predictor-type TYPE --traits FILE --trait-type TYPE [options]
```

An association run requires predictor/genotype and trait/expression inputs plus `--output`. `--validate-only` and `--inspect-missingness` do not run association. `--preprocess-only` requires VCF/BCF plus the trait matrix and any selected covariates, but does not require `--output`.

Command-line paths are resolved against the current directory. Relative paths inside a legacy INI are resolved against the INI file's directory. Command-line values override the loaded INI. Do not combine the positional INI spelling with `--config`.

## Examples

### VCF.gz association with aligned IDs and bounded RAM

PowerShell:

```powershell
java -Xmx12g --enable-native-access=ALL-UNNAMED `
  -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype cohort.vcf.gz `
  --genotype-format vcf `
  --genotype-field auto `
  --predictor-missing mean `
  --expression traits.csv `
  --trait-missing pattern `
  --covariates covariates.csv `
  --genotype-id-column participant_id `
  --expression-id-column sample_id `
  --trait-id-strip-prefix X `
  --sample-alignment covariate-subset `
  --fixed-covariates Sex,Age,Batch,PC1,PC2 `
  --min-mac 20 `
  --genotype-block-rows 2048 `
  --expression-block-rows 2048 `
  --trait-cache memory `
  --precision fp64 `
  --residualization auto `
  --threshold pval 1e-4 `
  --output results.csv
```

The application aligns samples before computing VCF/BCF frequency statistics. Omitting `--block-size`, `--threads`, and `--variant-qc-threads` selects bounded automatic values. Trait pattern deletion is exact but can be slow when there are many distinct missingness patterns.

### Preprocess once and reuse VCF/BCF QC and dosage data

First create alignment-specific QC and the raw FP64 dosage cache without initializing a compute backend:

```powershell
java -Xmx12g --enable-native-access=ALL-UNNAMED `
  -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype cohort.vcf.gz `
  --genotype-format vcf `
  --genotype-field auto `
  --expression traits.csv `
  --covariates covariates.csv `
  --fixed-covariates Sex,Age,Batch,PC1,PC2 `
  --sample-alignment strict `
  --min-mac 20 `
  --variant-qc-output work\cohort.variants.tsv `
  --variant-qc-checkpoint work\variant-qc-checkpoint `
  --missingness-qc-output work\cohort.missingness.tsv `
  --cache-dir work\cache `
  --preprocess-only
```

Then run association with the same source, alignment, regions, field, frequency filters, QC paths, and cache directory:

```powershell
java -Xmx12g --enable-native-access=ALL-UNNAMED `
  -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype cohort.vcf.gz `
  --genotype-format vcf `
  --genotype-field auto `
  --expression traits.csv `
  --covariates covariates.csv `
  --fixed-covariates Sex,Age,Batch,PC1,PC2 `
  --sample-alignment strict `
  --min-mac 20 `
  --variant-qc-output work\cohort.variants.tsv `
  --variant-qc-checkpoint work\variant-qc-checkpoint `
  --cache-dir work\cache `
  --genotype-block-rows 2048 `
  --expression-block-rows 2048 `
  --threshold pval 1e-4 `
  --output results.csv
```

The matching QC checkpoint avoids repeating aligned MAF/MAC/HWE calculations. The matching raw cache avoids decoding the compressed genotype source during subsequent missingness and prepared-cache work. A changed sample set, source/index metadata, field, regions, or frequency/inclusion policy deliberately creates different state.

### Indexed regions or variant sets

```powershell
java --enable-native-access=ALL-UNNAMED `
  -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype cohort.vcf.gz `
  --genotype-format vcf `
  --variant-index cohort.vcf.gz.tbi `
  --region APOE=chr19:44900000-45000000 `
  --region nearby=chr19:44800000-45100000 `
  --expression traits.csv `
  --min-mac 0 `
  --threshold none 0 `
  --output region-results.csv
```

Inline coordinates are one-based inclusive. For a file of regions, use `--regions-file` and optionally `--region-coordinates bed`. Verified random-access indexes are tabix `.tbi` for BGZF VCF and HTSJDK Tribble `.idx` for VCF/BCF; standard bcftools CSI is not yet supported.

### Inspect missingness without association

```powershell
java -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --genotype predictors.csv `
  --predictor-type methylation `
  --expression proteins.csv `
  --trait-type proteomics `
  --predictor-missing mean `
  --trait-missing pattern `
  --missingness-qc-output missingness.tsv `
  --inspect-missingness
```

Generic `--predictor` and `--traits` spellings require their corresponding type arguments. The historical `--genotype` and `--expression` spellings retain genotype/expression role defaults.

### Force CPU execution

```powershell
java --enable-native-access=ALL-UNNAMED `
  -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --backend cpu `
  --genotype predictors.csv `
  --expression traits.csv `
  --threshold pval 1e-4 `
  --output results.csv
```

The standard jar tries bundled OpenBLAS and then portable Java. A platform-specific oneMKL build can be selected with the JVM property described under [JVM system properties](#jvm-system-properties).

## All accepted options

In the tables below, “none” means the option is absent unless supplied. Boolean INI keys use `true` or `false`. Application defaults are shown for the modern headered CSV/VCF/BCF path; the legacy TPED path is called out separately.

### Configuration and input

| Command-line option | Values and default | Legacy INI key | Meaning |
| --- | --- | --- | --- |
| positional `FILE` | One INI path; none | n/a | Loads a legacy INI. It cannot be combined with `--config`. |
| `--config FILE` | One INI path; none | n/a | Loads an INI, then applies command-line overrides. |
| `--genotype FILE` | Required predictor path | `genotype_file` | Compatibility spelling for a genotype predictor; defaults predictor type to `genotype`. |
| `--predictor FILE` | Required predictor path | `genotype_file` | Generic predictor spelling; requires `--predictor-type`. |
| `--expression FILE` | Required trait path | `expression_file` | Compatibility spelling for expression traits; defaults trait type to `expression`. |
| `--traits FILE` | Required trait path | `expression_file` | Generic trait spelling; requires `--trait-type`. |
| `--covariates FILE` | Optional path | `covariate_file` | Mixed numeric/categorical covariate table. At least one fixed covariate must be named when this file is supplied. |
| `--output FILE` | Required for association | `output_file` | Association CSV. Optional for validation, missingness inspection, and VCF/BCF preprocessing. |
| `--analysis METHOD` | `eqtl` (default), `burden`, `skat`, `skat-o` | `analysis` | Selects ordinary variant-by-trait eQTL or an FP64 set test. Expression rows are the tested phenotypes; fixed covariates are adjustment variables. |
| `--variant-sets FILE` | Optional TSV | `variant_sets` | Explicit set membership, exact REF/ALT/effect allele, and optional positive weight definitions. Required for custom CSV sets. |
| `--window-size BP` | Positive integer; none | `window_size` | Generates nonempty chromosome-local one-based sliding-window sets directly from canonical variant coordinates. Mutually exclusive with explicit/custom regions and sets. |
| `--window-stride BP` | Window size when omitted | `window_stride` | Distance between automatic window starts; positive and no larger than the window size. |
| `--set-audit-output FILE` | `<output>.sets.tsv` | `set_audit_output` | Set membership/filter/status audit. |
| `--set-min-maf V`, `--set-max-maf V` | `0`, `0.5` | `set_min_maf`, `set_max_maf` | Inclusive aligned-cohort set-member MAF mask. No rare cutoff is inferred. |
| `--set-min-mac V`, `--set-max-mac V` | `0`, unbounded | `set_min_mac`, `set_max_mac` | Inclusive aligned-cohort set-member MAC mask. |
| `--set-absent-variant POLICY` | `error` (default), `skip` | `set_absent_variant` | Handling for explicit definitions absent from the aligned source. |
| `--set-degenerate POLICY` | `error` (default), `skip` | `set_degenerate` | Handling for empty or post-projection monomorphic sets. |
| `--set-block-size N` | `0` automatic | `set_block_size` | Resident set tile size for bounded execution. Omitted/zero selects a heap-aware value from aligned sample count, method workspaces, and the actual set-membership density; a positive value remains an explicit override. |
| `--skat-o-rho-grid LIST` | `0,0.25,0.5,0.75,1` | `skat_o_rho_grid` | Strictly increasing rho grid beginning at zero and ending at one. |
| `--skat-o-simulations N` | `10000` | `skat_o_simulations` | Correlated parametric-null replicate count. |
| `--skat-o-seed N` | `20260827` | `skat_o_seed` | Deterministic SKAT-O simulation seed. |
| `--genotype-format FORMAT` | `auto` (CLI default), `csv`, `vcf`, `vcf.gz`, `bcf`; `tped` legacy | `genotype_format` | `auto` infers `.vcf`, `.vcf.gz`, and `.bcf`, otherwise headered CSV. An old INI that omits this key retains the historical `tped` default. |
| `--expression-format FORMAT` | `csv` | `expression_format` | Current modern trait input is a headered delimited matrix. |
| `--predictor-type TYPE` | `genotype` for compatibility spelling | `predictor_type` | `genotype`, `expression`, `methylation`, `proteomics`, or `continuous`. Aliases accepted for some types include `protein`, `proteomic`, `dna-methylation`, and `dnm`. |
| `--trait-type TYPE` | `expression` for compatibility spelling | `trait_type` | Same values as predictor type. Required with `--traits`. |
| `--family FILE` | Optional; none | `family_file` | PLINK family file used only by the legacy TPED loader. |
| `--pedigree FILE` | Optional; none | `pedigree_file` | Legacy placeholder: existence is checked and the path is reported, but no mixed/familial model currently consumes it. |
| `--no-genotype-header` | Flag; header is expected by default | `genotype_file_header = false` | Selects the legacy headerless delimited loader. It is not valid for the modern headered CSV alignment/streaming path. |

### Matrix model, missing values, and sample alignment

| Command-line option | Values and default | Legacy INI key | Meaning |
| --- | --- | --- | --- |
| `--genotype-model MODEL` | `additive` (default), `full` | `genotype_model` | VCF/BCF and current real-valued modern analysis support additive only. The old full categorical-SNP branch remains disabled/broken. |
| `--predictor-missing POLICY` | `mean` (default), `error`, `zero`, `local-pattern`, `exclude-row` | `predictor_missing` | Missing predictor handling. `local-pattern` is genotype-only. Predictor `pattern` is rejected. Aliases include `exclude`, `exclude-variant`, `flanking`, `neighbor-pattern`, and `local-haplotype`. |
| `--genotype-missing POLICY` | Same as predictor policy | `genotype_missing` | Compatibility alias. `predictor_missing` takes precedence if both keys exist. |
| `--predictor-flanks N` | Integer at least 1; default `1` | `predictor_flanks` | Variants on each side used by genotype-only local-pattern imputation. |
| `--trait-missing POLICY` | `pattern` (default), `mean`, `zero`, `error`, `exclude-row` | `trait_missing` | `pattern` groups exact complete-sample masks for dynamic deletion. Aliases include `dynamic`, `complete-case`, `exclude`, and `exclude-trait`. |
| `--max-trait-patterns N` | `256` (default); `0` disables | `max_trait_patterns` | Automatic scheduler switch threshold. Explicit pattern-outer scheduling treats it as a safety limit and reports repeated predictor-preparation work. It never bypasses rank/DF checks. |
| `--trait-pattern-scheduler MODE` | `auto` (default), `pattern`, `genotype` | `trait_pattern_scheduler` | `auto` uses pattern-outer at or below the threshold and scalable genotype-outer above it. Genotype-outer is FP64, aligned-frequency-scope only, and writes deterministic genotype-block/trait-block/variant/original-trait order. |
| `--unestimable-trait-patterns POLICY` | `error` (default), `skip` | `unestimable_trait_patterns` | Genotype-outer writes `<output>.trait-patterns.tsv`; `error` then stops if any mask is rank/DF-unestimable, while explicit `skip` excludes those audited traits. |
| `--covariate-missing POLICY` | `complete-samples` (default), `error` | `covariate_missing` | Removes samples missing selected fixed covariates across all aligned matrices, or fails. |
| `--inspect-missingness` | Flag; off | `inspect_missingness = true` | Writes the missingness/alignment QC report and stops before backend initialization. Cannot be combined with `--preprocess-only`. |
| `--missingness-qc-output FILE` | `<output>.missingness.tsv`; predictor-based name when no output | `missingness_qc_output` | Changes the missingness/alignment QC path. |
| `--fixed-covariates LIST` | Comma-separated; none | `covariate_fixed` | Covariates to include. Quote a whitespace-separated CLI list; INI lists are whitespace-separated. An intercept is added automatically. |
| `--factor-covariates LIST` | Comma-separated; none | `covariate_factor` | Forces named numeric-looking covariates to categorical factors. Text covariates are categorical automatically; the lexicographically first level is the reference. |
| `--random-covariates LIST` | Comma-separated; none | `covariate_random` | Accepted for compatibility but currently warns and is ignored; it does not fit random effects. |
| `--genotype-id-column NAME` | None; exact-set inference when possible | `genotype_id_column` | Covariate column containing genotype/predictor sample IDs. |
| `--expression-id-column NAME` | None; exact-set inference when possible | `expression_id_column` | Covariate column containing expression/trait sample IDs. |
| `--sample-alignment POLICY` | `strict` (default), `covariate-subset` | `sample_alignment` | `strict` requires equal normalized sample sets. `covariate-subset` uses covariate rows as the canonical subset; it never takes a silent intersection. |
| `--predictor-id-strip-prefix TEXT` | None | `predictor_id_strip_prefix` | Removes one literal leading prefix before ID validation/alignment. |
| `--genotype-id-strip-prefix TEXT` | None | `predictor_id_strip_prefix` | Compatibility alias for predictor prefix removal. The older INI key `genotype_id_strip_prefix` is also read. |
| `--trait-id-strip-prefix TEXT` | None | `trait_id_strip_prefix` | Removes one literal leading prefix before ID validation/alignment. |
| `--expression-id-strip-prefix TEXT` | None | `trait_id_strip_prefix` | Compatibility alias for trait prefix removal. The older INI key `expression_id_strip_prefix` is also read. |

### VCF/BCF fields, QC, filters, and regions

These options affect VCF/BCF genotype input. CSV predictors do not receive variant MAF/MAC/HWE filtering.

| Command-line option | Values and default | Legacy INI key | Meaning |
| --- | --- | --- | --- |
| `--genotype-field FIELD` | `auto` (default), `DS`, `GT` | `genotype_field` | `auto` prefers header-declared dosage and falls back to genotype calls. |
| `--multiallelic POLICY` | `exclude` (default), `error` | `multiallelic` | Current association model is biallelic. |
| `--min-maf VALUE` | `0` (disabled); range 0–0.5 | `min_maf` | Minimum aligned or pattern-specific minor-allele frequency. |
| `--min-mac VALUE` | `20` (default); `0` disables | `min_mac` | Minimum minor-allele count. When MAF and MAC are positive, both must pass. Use `--min-mac 0` for MAF-only or rare-variant preprocessing. |
| `--frequency-scope SCOPE` | `aligned` (default), `pattern` | `frequency_scope` | Applies MAF/MAC once after final sample alignment, or separately for each exact trait mask. HWE remains aligned-cohort QC. |
| `--variant-qc-output FILE` | `<output>.variants.tsv`; genotype-based name without output | `variant_qc_output` | Variant IDs/alleles, aligned counts, EAF/MAF/MAC, HWE, classifications, filters, and set membership. |
| `--variant-qc-threads N` | `0` automatic; `1` sequential | `variant_qc_threads` | Parallel low-allocation VCF FORMAT/GT/DS/FT QC workers. Automatic mode leaves one logical processor free and is capped at 8 based on the observed throughput plateau; `1` uses the HTSJDK reference path, while BCF genotype expansion remains reader-thread-bound. |
| `--variant-qc-checkpoint DIR` | `<variant-QC-output>.checkpoint` | `variant_qc_checkpoint` | Durable, signature-scoped QC state. Reuse resumes or avoids the aligned QC scan. |
| `--variant-index FILE` | Neighboring index auto-detected | `variant_index` | Explicit tabix `.tbi` or HTSJDK Tribble `.idx`. A standard `.csi` is recognized but not decoded; region requests fail clearly rather than falling back to a sequential scan. |
| `--region [SET=]CHROM:START-END` | Repeatable; none | `regions` | Indexed one-based inclusive interval, optionally assigned to a set. Repeated INI regions are semicolon-separated. |
| `--regions-file FILE` | Optional TSV | `regions_file` | `CHROM START END` or `SET_ID CHROM START END`, optionally with a header. |
| `--region-coordinates MODE` | `one-based` (default), `bed` | `region_coordinates` | Interpretation for region-file coordinates. Inline `--region` remains one-based inclusive. |
| `--preprocess-only` | Flag; off | `preprocess_only = true` | Aligns samples, validates covariates/DF, performs variant QC/filtering and missingness reporting, creates/reuses aligned raw FP64 dosage cache, then stops before association/backend initialization. VCF/BCF only. |

### Association, output, and compute selection

| Command-line option | Values and default | Legacy INI key | Meaning |
| --- | --- | --- | --- |
| `--threshold TYPE VALUE` | `none` when omitted | `threshold = TYPE VALUE` | `none`, `pval`, or `rsq`. Supply a numeric value even for `none`, for example `--threshold none 0`. No reporting cutoff is applied when omitted. |
| `--df-offset N` | Integer; default `0` | `df_offset` | Additional degrees of freedom subtracted from the p-value calculation. The final error DF must remain positive. |
| `--simplify-output` | Flag; off | `simplify_output = true` | Uses simplified identifier formatting retained from the legacy output path. |
| `--rsq-only` | Flag; off | `rsq_only = true` | Emits the reduced R-squared/direction output columns instead of full effect/t/log10-p columns. |
| `--precision MODE` | `fp64` (default), `fp32` | `precision` | Backend matrix-product precision. FP32 is explicit and changes results slightly; standardization and final statistics remain FP64. Aliases `double`/`64` and `float`/`single`/`32` are accepted. |
| `--residualization MODE` | `auto` (default), `gpu`, `cpu` | `residualization` | Where fixed-effect projection is applied during modern cache preparation. QR/rank decisions remain CPU FP64. |
| `--backend BACKEND` | `auto` (default), `cuda`, `opencl`, `cpu` | none | Selects the compute backend for this process. `auto` prefers usable CUDA, retains distinct OpenCL devices, then falls back to CPU. |
| `--gpu-backend BACKEND` | Same values | none | Compatibility alias for `--backend`. |
| `--block-size N` | `0`/omitted means automatic | `block_size` | Compute tile capacity. Automatic GPU sizing probes all selected devices; CPU sizing uses available JVM heap. |
| `--threads N` | `0`/omitted means automatic | `num_threads` | Application pipeline workers. Native CPU BLAS has a separate inner thread pool. |
| `--genotype-block-rows N` | `0`/omitted means full-memory unless required | `genotype_block_rows` | Predictor rows per prepared/cache block; a positive value enables bounded-RAM scheduling. |
| `--expression-block-rows N` | `0`/omitted means full-memory unless required | `expression_block_rows` | Trait rows per prepared/cache block; a positive value enables bounded-RAM scheduling. |
| `--trait-cache MODE` | `auto` (default), `memory`, `disk` | `trait_cache` | Keeps once-prepared traits in safe heap residency or reads indexed disk blocks. |
| `--expression-cache MODE` | Same values | `trait_cache` | Compatibility alias for `--trait-cache`. |

### Caches, checkpoints, diagnostics, and control flags

| Command-line option | Values and default | Legacy INI key | Meaning |
| --- | --- | --- | --- |
| `--cache-dir DIR` | `.gpu-eqtl-cache` beside output/input | `cache_dir` | Persistent signed missingness scans plus raw and prepared matrix caches, keyed by source/parser/sample-order and scientific preprocessing signatures. |
| `--rebuild-cache` | Flag; off | `rebuild_cache = true` | Replaces matching missingness/matrix caches instead of reusing them. It does not replace variant-QC checkpoint state. |
| `--checkpoint-dir DIR` | `<output>.checkpoint` | `checkpoint_dir` | Bounded-RAM association result parts. Pattern-outer nests blocks under ordered pattern groups; genotype-outer commits blocks spanning all estimable traits. |
| `--resume` | Flag; off | `resume = true` | Reuses completed signature-matching association blocks/groups and separately reuses matching signed missingness and prepared-trait caches. |
| `--keep-checkpoints` | Flag; off | `keep_checkpoints = true` | Retains completed block/pattern parts after final assembly. |
| `--profile` | Flag; off | `profile = true` | Prints phase timing and transfer/byte summaries. |
| `--profile-output FILE` | None | `profile_output` | Writes profiling CSV and implicitly enables profiling. |
| `--validate-only` | Flag; off | `validate_only = true` | Validates modern IDs/alignment, covariate model/rank, missingness policies, and degrees of freedom without association. Does not create the VCF/BCF raw cache. |
| `--debug` | Flag; off | none | Enables additional diagnostic output for this process. |
| `--printbackendinfo` | Flag | none | Prints detected compute and CPU engine information, then exits. |
| `--printgpuinfo` | Flag | none | Compatibility alias for `--printbackendinfo`; output also includes CPU fallback information. |
| `--help`, `-h` | Flag | none | Prints the built-in short help and exits. This page is the exhaustive reference. |

## JVM system properties

These are Java system properties, so they must appear before `-jar`. They are not ordinary application arguments or INI settings.

| Property | Values and default | Meaning |
| --- | --- | --- |
| `-Deqtl.gpu.backend=VALUE` | `auto` (default), `cuda`, `opencl`, `cpu` | Compatibility selection equivalent to `--backend`; the explicit application argument sets this property for the process. |
| `-Deqtl.cpu.blas=VALUE` | `auto` (default), `mkl`, `openblas`, `java` | Selects the CPU matrix engine. An explicitly forced unavailable native engine fails instead of silently falling back. |
| `-Deqtl.cpu.threads=N` | `max(1, logical processors - 1)` | Controls the native BLAS inner thread pool independently of application `--threads`. |

For example, in PowerShell the `-D` argument should be quoted:

```powershell
java '-Deqtl.cpu.blas=openblas' '-Deqtl.cpu.threads=8' `
  --enable-native-access=ALL-UNNAMED `
  -jar target\gpu-eqtl-2.0.0-SNAPSHOT-all.jar `
  --backend cpu --printbackendinfo
```

The internal `eqtl.variant.qc.checkpoint.records` property and `eqtl.test.*` properties exist only for deterministic development tests and are deliberately not supported user-facing analysis controls.

## Defaults that materially affect results

- Association reporting has no threshold by default. Use `--threshold pval 1e-4` when that is the intended cutoff.
- VCF/BCF filtering defaults to `MAC >= 20` and no MAF cutoff. Explicit `--min-mac 0` disables MAC filtering. Positive MAF and MAC cutoffs are conjunctive.
- Predictor missing values default to row-mean imputation; trait missing values default to exact pattern-wise deletion; selected-covariate missing values default to removing those samples from all aligned matrices.
- Sample alignment defaults to `strict`. Missing, duplicate, transformed-collision, or unmatched IDs are fatal unless the explicit `covariate-subset` policy applies.
- FP64 is the default. FP32 is never selected automatically merely because a device lacks FP64.
- Block size, application threads, variant-QC threads, trait residency, residualization location, and backend all use automatic selection when omitted.

## Legacy INI notes

The INI parser uses `key = value` lines. Command-line lists such as `Age,Batch,PC1` become whitespace-separated INI values such as `covariate_fixed = Age Batch PC1`. Boolean flags use `true`. A minimal modern INI can therefore look like:

```ini
genotype_file = cohort.vcf.gz
genotype_format = vcf
genotype_field = auto
expression_file = traits.csv
covariate_file = covariates.csv
covariate_fixed = Sex Age Batch PC1 PC2
sample_alignment = strict
min_mac = 20
threshold = pval 1e-4
output_file = results.csv
```

The obsolete `library_path` and `lambda` keys may still be read from old files but no longer control native-library discovery or memory saturation. Do not use them in new configurations.
