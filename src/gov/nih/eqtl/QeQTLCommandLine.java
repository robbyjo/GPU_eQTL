/*
 * Roby Joehanes
 *
 * Copyright 2011 Roby Joehanes
 * This file is distributed under the GNU General Public License version 3.0.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */
package gov.nih.eqtl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Command-line arguments with a compatibility path for the legacy INI file. */
public final class QeQTLCommandLine {
    public record Result(QeQTLAnalysisConfig config, boolean help, boolean printGpuInfo,
        boolean debug, String gpuBackend) { }

    private static final Map<String, String> VALUE_OPTIONS = new LinkedHashMap<>();
    static {
        VALUE_OPTIONS.put("--genotype", "genotype_file");
        VALUE_OPTIONS.put("--predictor", "genotype_file");
        VALUE_OPTIONS.put("--expression", "expression_file");
        VALUE_OPTIONS.put("--traits", "expression_file");
        VALUE_OPTIONS.put("--covariates", "covariate_file");
        VALUE_OPTIONS.put("--output", "output_file");
        VALUE_OPTIONS.put("--analysis", "analysis");
        VALUE_OPTIONS.put("--variant-sets", "variant_sets");
        VALUE_OPTIONS.put("--set-audit-output", "set_audit_output");
        VALUE_OPTIONS.put("--set-min-maf", "set_min_maf");
        VALUE_OPTIONS.put("--set-max-maf", "set_max_maf");
        VALUE_OPTIONS.put("--set-min-mac", "set_min_mac");
        VALUE_OPTIONS.put("--set-max-mac", "set_max_mac");
        VALUE_OPTIONS.put("--set-absent-variant", "set_absent_variant");
        VALUE_OPTIONS.put("--set-degenerate", "set_degenerate");
        VALUE_OPTIONS.put("--set-block-size", "set_block_size");
        VALUE_OPTIONS.put("--window-size", "window_size");
        VALUE_OPTIONS.put("--window-stride", "window_stride");
        VALUE_OPTIONS.put("--skat-o-rho-grid", "skat_o_rho_grid");
        VALUE_OPTIONS.put("--skat-o-simulations", "skat_o_simulations");
        VALUE_OPTIONS.put("--skat-o-seed", "skat_o_seed");
        VALUE_OPTIONS.put("--family", "family_file");
        VALUE_OPTIONS.put("--pedigree", "pedigree_file");
        VALUE_OPTIONS.put("--genotype-format", "genotype_format");
        VALUE_OPTIONS.put("--expression-format", "expression_format");
        VALUE_OPTIONS.put("--genotype-model", "genotype_model");
        VALUE_OPTIONS.put("--genotype-field", "genotype_field");
        VALUE_OPTIONS.put("--genotype-missing", "genotype_missing");
        VALUE_OPTIONS.put("--predictor-type", "predictor_type");
        VALUE_OPTIONS.put("--trait-type", "trait_type");
        VALUE_OPTIONS.put("--predictor-missing", "predictor_missing");
		VALUE_OPTIONS.put("--predictor-flanks", "predictor_flanks");
        VALUE_OPTIONS.put("--trait-missing", "trait_missing");
        VALUE_OPTIONS.put("--covariate-missing", "covariate_missing");
        VALUE_OPTIONS.put("--missingness-qc-output", "missingness_qc_output");
        VALUE_OPTIONS.put("--multiallelic", "multiallelic");
        VALUE_OPTIONS.put("--min-maf", "min_maf");
        VALUE_OPTIONS.put("--min-mac", "min_mac");
        VALUE_OPTIONS.put("--variant-qc-output", "variant_qc_output");
        VALUE_OPTIONS.put("--variant-qc-threads", "variant_qc_threads");
        VALUE_OPTIONS.put("--variant-qc-checkpoint", "variant_qc_checkpoint");
        VALUE_OPTIONS.put("--variant-index", "variant_index");
        VALUE_OPTIONS.put("--region", "regions");
        VALUE_OPTIONS.put("--regions-file", "regions_file");
        VALUE_OPTIONS.put("--region-coordinates", "region_coordinates");
        VALUE_OPTIONS.put("--frequency-scope", "frequency_scope");
        VALUE_OPTIONS.put("--fixed-covariates", "covariate_fixed");
        VALUE_OPTIONS.put("--random-covariates", "covariate_random");
        VALUE_OPTIONS.put("--factor-covariates", "covariate_factor");
        VALUE_OPTIONS.put("--genotype-id-column", "genotype_id_column");
        VALUE_OPTIONS.put("--expression-id-column", "expression_id_column");
        VALUE_OPTIONS.put("--sample-alignment", "sample_alignment");
        VALUE_OPTIONS.put("--predictor-id-strip-prefix", "predictor_id_strip_prefix");
        VALUE_OPTIONS.put("--genotype-id-strip-prefix", "predictor_id_strip_prefix");
        VALUE_OPTIONS.put("--trait-id-strip-prefix", "trait_id_strip_prefix");
        VALUE_OPTIONS.put("--expression-id-strip-prefix", "trait_id_strip_prefix");
        VALUE_OPTIONS.put("--df-offset", "df_offset");
        VALUE_OPTIONS.put("--block-size", "block_size");
        VALUE_OPTIONS.put("--threads", "num_threads");
        VALUE_OPTIONS.put("--precision", "precision");
		VALUE_OPTIONS.put("--residualization", "residualization");
        VALUE_OPTIONS.put("--genotype-block-rows", "genotype_block_rows");
        VALUE_OPTIONS.put("--expression-block-rows", "expression_block_rows");
        VALUE_OPTIONS.put("--trait-cache", "trait_cache");
        VALUE_OPTIONS.put("--expression-cache", "trait_cache");
        VALUE_OPTIONS.put("--cache-dir", "cache_dir");
        VALUE_OPTIONS.put("--checkpoint-dir", "checkpoint_dir");
		VALUE_OPTIONS.put("--max-trait-patterns", "max_trait_patterns");
		VALUE_OPTIONS.put("--trait-pattern-scheduler", "trait_pattern_scheduler");
		VALUE_OPTIONS.put("--unestimable-trait-patterns", "unestimable_trait_patterns");
        VALUE_OPTIONS.put("--profile-output", "profile_output");
    }

    private static final List<String> PATH_OPTIONS = List.of(
        "genotype_file", "expression_file", "covariate_file", "output_file", "family_file", "pedigree_file",
        "cache_dir", "checkpoint_dir", "profile_output", "variant_qc_output",
        "variant_qc_checkpoint", "variant_index", "regions_file", "missingness_qc_output",
        "variant_sets", "set_audit_output");

    private QeQTLCommandLine() { }

    public static Result parse(String[] args) throws IOException {
        if (args == null)
            args = new String[0];
        String configFile = null;
        List<String> positional = new ArrayList<>();
        for (int i = 0; i < args.length; i++) {
            if ("--config".equals(args[i])) {
                if (++i >= args.length)
                    throw new IllegalArgumentException("--config requires a file name");
                if (configFile != null)
                    throw new IllegalArgumentException("Only one --config file may be specified");
                configFile = args[i];
            } else if (!args[i].startsWith("-")) {
                positional.add(args[i]);
            } else if (VALUE_OPTIONS.containsKey(args[i]) || "--threshold".equals(args[i])
                || "--gpu-backend".equals(args[i]) || "--backend".equals(args[i])) {
                int values = "--threshold".equals(args[i]) ? 2 : 1;
                i += values;
                if (i >= args.length)
                    throw new IllegalArgumentException(args[i - values] + " requires " + values + " value(s)");
            }
        }
        if (configFile == null && !positional.isEmpty()) {
            if (positional.size() > 1)
                throw new IllegalArgumentException("Unexpected positional arguments: " + positional);
            configFile = positional.get(0);
        } else if (configFile != null && !positional.isEmpty()) {
            throw new IllegalArgumentException("Do not combine a positional INI file with --config");
        }

        QeQTLAnalysisConfig config;
        if (configFile != null) {
            config = QeQTLAnalysisConfig.loadConfig(configFile);
        } else {
            Map<String, String> values = new LinkedHashMap<>();
            String workingDirectory = new File(".").getAbsoluteFile().toPath().normalize().toString();
            if (!workingDirectory.endsWith(File.separator))
                workingDirectory += File.separator;
            values.put("ini.path", workingDirectory);
            values.put("genotype_format", "auto");
            values.put("expression_format", "csv");
            config = new QeQTLAnalysisConfig(values);
        }

        boolean help = false;
        boolean printGpuInfo = false;
        boolean debug = false;
        String gpuBackend = null;
		boolean genericPredictor = false;
		boolean genericTraits = false;
        for (int i = 0; i < args.length; i++) {
            String argument = args[i];
            if (!argument.startsWith("-"))
                continue;
            if ("--config".equals(argument)) {
                i++;
            } else if (VALUE_OPTIONS.containsKey(argument)) {
				if ("--predictor".equals(argument)) genericPredictor = true;
				if ("--traits".equals(argument)) genericTraits = true;
                String key = VALUE_OPTIONS.get(argument);
                String value = requireValue(args, ++i, argument);
                if (PATH_OPTIONS.contains(key))
                    value = new File(value).getAbsoluteFile().toPath().normalize().toString();
                if (key.startsWith("covariate_"))
                    value = value.replace(',', ' ').trim();
                if (key.equals("regions") && config.get(key) != null && !config.get(key).isBlank())
                    config.set(key, config.get(key) + ";" + value);
                else
                    config.set(key, value);
            } else if ("--threshold".equals(argument)) {
                String type = requireValue(args, ++i, argument);
                String value = requireValue(args, ++i, argument);
                config.set("threshold", type + " " + value);
            } else if ("--gpu-backend".equals(argument) || "--backend".equals(argument)) {
                gpuBackend = requireValue(args, ++i, argument).toLowerCase();
                if (!(gpuBackend.equals("auto") || gpuBackend.equals("cuda")
					|| gpuBackend.equals("opencl") || gpuBackend.equals("cpu")))
					throw new IllegalArgumentException(argument + " must be auto, cuda, opencl, or cpu");
            } else if ("--simplify-output".equals(argument)) {
                config.set("simplify_output", "true");
            } else if ("--rsq-only".equals(argument)) {
                config.set("rsq_only", "true");
            } else if ("--validate-only".equals(argument)) {
                config.set("validate_only", "true");
            } else if ("--preprocess-only".equals(argument)) {
                config.set("preprocess_only", "true");
            } else if ("--inspect-missingness".equals(argument)) {
                config.set("inspect_missingness", "true");
            } else if ("--rebuild-cache".equals(argument)) {
                config.set("rebuild_cache", "true");
            } else if ("--resume".equals(argument)) {
                config.set("resume", "true");
            } else if ("--keep-checkpoints".equals(argument)) {
                config.set("keep_checkpoints", "true");
            } else if ("--profile".equals(argument)) {
                config.set("profile", "true");
            } else if ("--no-genotype-header".equals(argument)) {
                config.set("genotype_file_header", "false");
            } else if ("--debug".equals(argument)) {
                debug = true;
            } else if ("--printgpuinfo".equals(argument) || "--printbackendinfo".equals(argument)) {
                printGpuInfo = true;
            } else if ("--help".equals(argument) || "-h".equals(argument)) {
                help = true;
            } else {
                throw new IllegalArgumentException("Unknown option: " + argument);
            }
        }
		if (genericPredictor && config.get("predictor_type") == null)
			throw new IllegalArgumentException("--predictor requires --predictor-type; use --genotype for the compatibility default");
		if (genericTraits && config.get("trait_type") == null)
			throw new IllegalArgumentException("--traits requires --trait-type; use --expression for the compatibility default");
        return new Result(config, help, printGpuInfo, debug, gpuBackend);
    }

    private static String requireValue(String[] args, int index, String option) {
        if (index >= args.length)
            throw new IllegalArgumentException(option + " requires a value");
        return args[index];
    }

    public static String usage() {
        return """
            Usage:
              java -jar gpu-eqtl-...-all.jar legacy.ini
              java -jar gpu-eqtl-...-all.jar --config legacy.ini [overrides]
              java -jar gpu-eqtl-...-all.jar --genotype FILE --expression FILE --output FILE [options]
			  java -jar gpu-eqtl-...-all.jar --predictor FILE --predictor-type TYPE --traits FILE --trait-type TYPE --output FILE [options]

            Main options:
              --genotype-format {auto|csv|vcf|vcf.gz|bcf}
              --analysis {eqtl|burden|skat|skat-o} (default: eqtl)
              --variant-sets FILE              Explicit set/variant/REF/ALT/effect-allele TSV
              --set-audit-output FILE          Default: OUTPUT.sets.tsv
              --set-min-maf V --set-max-maf V  Inclusive aligned-cohort set-test MAF mask
              --set-min-mac V --set-max-mac V  Inclusive aligned-cohort set-test MAC mask
              --set-absent-variant {error|skip} --set-degenerate {error|skip}
              --set-block-size N               Resident set tile size (default: 256)
              --window-size BP                 Automatic one-based sliding-window size
              --window-stride BP               Window-start stride (default: window size)
              --skat-o-rho-grid LIST           Default: 0,0.25,0.5,0.75,1
              --skat-o-simulations N --skat-o-seed N
              --predictor-type {genotype|expression|methylation|proteomics|continuous}
              --trait-type {genotype|expression|methylation|proteomics|continuous}
              --genotype-field {auto|DS|GT}   Prefer DS when auto and declared (default: auto)
              --predictor-missing {error|mean|zero|local-pattern|exclude-row}  (default: mean)
              --predictor-flanks N              Flanks per side for local-pattern (default: 1)
              --trait-missing {error|mean|zero|exclude-row|pattern} (default: pattern)
              --genotype-missing VALUE        Compatibility alias for --predictor-missing
              --covariate-missing {error|complete-samples} (default: complete-samples)
              --inspect-missingness            Write QC and stop before GPU initialization
              --missingness-qc-output FILE     Default: OUTPUT.missingness.tsv
              --multiallelic {exclude|error}  Current biallelic policy (default: exclude)
              --min-maf VALUE                    Minimum MAF (default: 0, disabled)
              --min-mac VALUE                    Minimum MAC (default: 20; use 0 to disable)
              --variant-qc-output FILE        Variant annotation/QC TSV (default: OUTPUT.variants.tsv)
              --variant-qc-threads N          Variant-level QC workers (default: auto; 1 is sequential)
              --variant-qc-checkpoint DIR     Resumable QC state root (default: QC_OUTPUT.checkpoint)
              --variant-index FILE            Explicit .tbi/.idx index (otherwise auto-detected)
              --region [SET=]CHROM:START-END  Indexed one-based inclusive region; repeatable
              --regions-file FILE             CHROM/START/END or SET_ID/CHROM/START/END TSV
              --region-coordinates {one-based|bed}  Region-file coordinates (default: one-based)
              --frequency-scope {aligned|pattern}   MAF/MAC filtering scope (default: aligned)
              --preprocess-only                 Align/QC/cache VCF or BCF, then stop before association
              --covariates FILE                 Mixed numeric/categorical covariate table
              --fixed-covariates LIST          Names separated by commas (or quote a space-separated list)
              --factor-covariates LIST         Force numeric-looking variables to be categorical
              --genotype-id-column NAME        Covariate column containing genotype sample IDs
              --expression-id-column NAME      Covariate column containing expression sample IDs
              --sample-alignment {strict|covariate-subset} (default: strict)
              --predictor-id-strip-prefix TEXT Remove a literal leading prefix before ID matching
              --trait-id-strip-prefix TEXT     Remove a literal leading prefix before ID matching
              --threshold {none|pval|rsq} VALUE
              --precision {fp64|fp32}         GPU matrix-product precision (default: fp64)
			  --residualization {auto|gpu|cpu} Fixed-effect projection location (default: auto)
              --df-offset N  --block-size N  --threads N
              --genotype-block-rows N          Enable bounded-RAM genotype analysis
              --expression-block-rows N        Enable bounded-RAM expression analysis
              --trait-cache {auto|memory|disk} Prepared trait residency (default: auto)
              --cache-dir DIR  --rebuild-cache
              --checkpoint-dir DIR  --resume  --keep-checkpoints
			  --max-trait-patterns N          Safety limit for exact deletion (default: 256; 0 disables)
			  --trait-pattern-scheduler {auto|pattern|genotype} (default: auto)
			  --unestimable-trait-patterns {error|skip} (default: error)
              --profile  --profile-output FILE  Phase timing summary and CSV
              --backend {auto|cuda|opencl|cpu}  Compute backend (default: auto)
			  --gpu-backend VALUE             Compatibility alias for --backend
              --simplify-output  --rsq-only  --validate-only  --debug
              --printbackendinfo  --printgpuinfo  --help
            """;
    }
}
