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
        VALUE_OPTIONS.put("--expression", "expression_file");
        VALUE_OPTIONS.put("--covariates", "covariate_file");
        VALUE_OPTIONS.put("--output", "output_file");
        VALUE_OPTIONS.put("--family", "family_file");
        VALUE_OPTIONS.put("--pedigree", "pedigree_file");
        VALUE_OPTIONS.put("--genotype-format", "genotype_format");
        VALUE_OPTIONS.put("--expression-format", "expression_format");
        VALUE_OPTIONS.put("--genotype-model", "genotype_model");
        VALUE_OPTIONS.put("--fixed-covariates", "covariate_fixed");
        VALUE_OPTIONS.put("--random-covariates", "covariate_random");
        VALUE_OPTIONS.put("--factor-covariates", "covariate_factor");
        VALUE_OPTIONS.put("--genotype-id-column", "genotype_id_column");
        VALUE_OPTIONS.put("--expression-id-column", "expression_id_column");
        VALUE_OPTIONS.put("--df-offset", "df_offset");
        VALUE_OPTIONS.put("--block-size", "block_size");
        VALUE_OPTIONS.put("--threads", "num_threads");
        VALUE_OPTIONS.put("--precision", "precision");
        VALUE_OPTIONS.put("--genotype-block-rows", "genotype_block_rows");
        VALUE_OPTIONS.put("--expression-block-rows", "expression_block_rows");
        VALUE_OPTIONS.put("--cache-dir", "cache_dir");
        VALUE_OPTIONS.put("--checkpoint-dir", "checkpoint_dir");
    }

    private static final List<String> PATH_OPTIONS = List.of(
        "genotype_file", "expression_file", "covariate_file", "output_file", "family_file", "pedigree_file",
        "cache_dir", "checkpoint_dir");

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
                || "--gpu-backend".equals(args[i])) {
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
            values.put("genotype_format", "csv");
            values.put("expression_format", "csv");
            config = new QeQTLAnalysisConfig(values);
        }

        boolean help = false;
        boolean printGpuInfo = false;
        boolean debug = false;
        String gpuBackend = null;
        for (int i = 0; i < args.length; i++) {
            String argument = args[i];
            if (!argument.startsWith("-"))
                continue;
            if ("--config".equals(argument)) {
                i++;
            } else if (VALUE_OPTIONS.containsKey(argument)) {
                String key = VALUE_OPTIONS.get(argument);
                String value = requireValue(args, ++i, argument);
                if (PATH_OPTIONS.contains(key))
                    value = new File(value).getAbsoluteFile().toPath().normalize().toString();
                if (key.startsWith("covariate_"))
                    value = value.replace(',', ' ').trim();
                config.set(key, value);
            } else if ("--threshold".equals(argument)) {
                String type = requireValue(args, ++i, argument);
                String value = requireValue(args, ++i, argument);
                config.set("threshold", type + " " + value);
            } else if ("--gpu-backend".equals(argument)) {
                gpuBackend = requireValue(args, ++i, argument).toLowerCase();
                if (!(gpuBackend.equals("auto") || gpuBackend.equals("cuda") || gpuBackend.equals("opencl")))
                    throw new IllegalArgumentException("--gpu-backend must be auto, cuda, or opencl");
            } else if ("--simplify-output".equals(argument)) {
                config.set("simplify_output", "true");
            } else if ("--rsq-only".equals(argument)) {
                config.set("rsq_only", "true");
            } else if ("--validate-only".equals(argument)) {
                config.set("validate_only", "true");
            } else if ("--rebuild-cache".equals(argument)) {
                config.set("rebuild_cache", "true");
            } else if ("--resume".equals(argument)) {
                config.set("resume", "true");
            } else if ("--keep-checkpoints".equals(argument)) {
                config.set("keep_checkpoints", "true");
            } else if ("--no-genotype-header".equals(argument)) {
                config.set("genotype_file_header", "false");
            } else if ("--debug".equals(argument)) {
                debug = true;
            } else if ("--printgpuinfo".equals(argument)) {
                printGpuInfo = true;
            } else if ("--help".equals(argument) || "-h".equals(argument)) {
                help = true;
            } else {
                throw new IllegalArgumentException("Unknown option: " + argument);
            }
        }
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

            Main options:
              --covariates FILE                 Mixed numeric/categorical covariate table
              --fixed-covariates LIST          Names separated by commas (or quote a space-separated list)
              --factor-covariates LIST         Force numeric-looking variables to be categorical
              --genotype-id-column NAME        Covariate column containing genotype sample IDs
              --expression-id-column NAME      Covariate column containing expression sample IDs
              --threshold {none|pval|rsq} VALUE
              --precision {fp64|fp32}         GPU matrix-product precision (default: fp64)
              --df-offset N  --block-size N  --threads N
              --genotype-block-rows N          Enable bounded-RAM CSV analysis
              --expression-block-rows N        Enable bounded-RAM CSV analysis
              --cache-dir DIR  --rebuild-cache
              --checkpoint-dir DIR  --resume  --keep-checkpoints
              --gpu-backend {auto|cuda|opencl}
              --simplify-output  --rsq-only  --validate-only  --debug
              --printgpuinfo  --help
            """;
    }
}
