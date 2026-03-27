#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 09A: Core ME-Trait Analysis
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load one or more module objects
#   - Load analysis-ready sample information
#   - Resolve sample IDs and align samples to module eigengenes (MEs)
#   - Compute ME-trait correlations (Pearson and/or bicor)
#   - Save full stats tables
#   - Save significant-only tables
#   - Optionally save outcome-only tables if an outcome trait file is provided
#   - Save Excel workbooks
#   - Save standard full heatmaps
#   - Save standard top-signals heatmaps
#   - Save run logs
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --sample_info  : path to sample information file (.xlsx, .csv, .tsv, .txt)
#   --modules_v1   : path to v1 Modules.rds
#
# OPTIONAL INPUTS
#   --sample_id_col         : column in sample_info containing sample IDs
#                             [default = use existing rownames]
#   --modules_v2            : optional path to v2 Modules.rds
#   --modules_v3            : optional path to v3 Modules.rds
#   --run_pearson           : TRUE/FALSE [default = FALSE]
#   --run_bicor             : TRUE/FALSE [default = TRUE]
#   --p_thresh              : significance threshold [default = 0.05]
#   --top_n                 : number of top associations for TOP heatmap [default = 250]
#   --module_dendro_distance: bicor or pearson [default = bicor]
#   --trait_exclude_file    : optional text file, one trait per line
#   --outcome_traits_file   : optional text file, one trait per line
#
# OUTPUTS
#   project_root/comethyl_output/09a_me_trait_analysis/<cpg_label>/<region_label>/<variant>/
#       ME_Trait_Correlation_Stats_Pearson.tsv
#       ME_Trait_Correlation_Stats_Pearson_significant.tsv
#       ME_Trait_Correlation_Stats_Pearson_outcome_only.tsv
#       ME_Trait_Correlation_Stats_Pearson_outcome_only_significant.tsv
#       ME_Trait_Correlations_Pearson.xlsx
#       ME_Trait_Correlation_Heatmap_Pearson_FULL.pdf
#       ME_Trait_Correlation_Heatmap_Pearson_TOP.pdf
#
#       ME_Trait_Correlation_Stats_Bicor.tsv
#       ME_Trait_Correlation_Stats_Bicor_significant.tsv
#       ME_Trait_Correlation_Stats_Bicor_outcome_only.tsv
#       ME_Trait_Correlation_Stats_Bicor_outcome_only_significant.tsv
#       ME_Trait_Correlations_Bicor.xlsx
#       ME_Trait_Correlation_Heatmap_Bicor_FULL.pdf
#       ME_Trait_Correlation_Heatmap_Bicor_TOP.pdf
#
#       sample_info_non_numeric_columns_removed.txt
#       sample_info_all_na_numeric_columns_removed.txt
#       sample_info_zero_variance_numeric_columns_removed.txt
#       traits_requested_exclude.txt
#       traits_found_exclude.txt
#       traits_missing_exclude.txt
#       traits_requested_outcome.txt
#       traits_found_outcome.txt
#       traits_missing_outcome.txt
#       run_parameters.txt

# ================================================================
# SCRIPT 09A: Core ME-Trait Analysis
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(openxlsx)
  library(AnnotationHub)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
})

# ------------------------------------------------------------
# Load helper.R
# ------------------------------------------------------------

# if 09a_me_trait_analysis.R is in scripts/

# then helper.R should also be in scripts/

# Example:

# project/
# └── scripts/
#     ├── helper.R
#     └── 09a_me_trait_analysis.R

# Then running:

# Rscript scripts/09a_me_trait_analysis.R ...

# will automatically load scripts/helper.R.

script_file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_file_arg) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) {
  stop("helper.R not found next to this script: ", helper_file)
}
source(helper_file)

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--sample_info", type = "character",
              help = "Path to sample information file (.xlsx, .csv, .tsv, .txt)"),

  make_option("--sample_id_col", type = "character", default = NULL,
              help = "Column in sample_info containing sample IDs [default = rownames]"),

  make_option("--modules_v1", type = "character",
              help = "Path to v1 Modules.rds"),

  make_option("--modules_v2", type = "character", default = NULL,
              help = "Optional path to v2 Modules.rds"),

  make_option("--modules_v3", type = "character", default = NULL,
              help = "Optional path to v3 Modules.rds"),

  make_option("--run_pearson", type = "character", default = "FALSE",
              help = "Run Pearson ME-trait correlations: TRUE or FALSE [default = FALSE]"),

  make_option("--run_bicor", type = "character", default = "TRUE",
              help = "Run bicor ME-trait correlations: TRUE or FALSE [default = TRUE]"),

  make_option("--p_thresh", type = "double", default = 0.05,
              help = "Significance threshold [default = 0.05]"),

  make_option("--top_n", type = "integer", default = 250,
              help = "Number of top associations for TOP heatmap [default = 250]"),

  make_option("--module_dendro_distance", type = "character", default = "bicor",
              help = "Distance for module dendrogram ordering: bicor or pearson [default = bicor]"),

  make_option("--trait_exclude_file", type = "character", default = NULL,
              help = "Optional text file with one trait to exclude per line"),

  make_option("--outcome_traits_file", type = "character", default = NULL,
              help = "Optional text file with one outcome trait per line")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$sample_info)) stop("--sample_info is required")
if (is.null(opt$modules_v1)) stop("--modules_v1 is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$sample_info)) stop("sample_info not found: ", opt$sample_info)
if (!file.exists(opt$modules_v1)) stop("modules_v1 not found: ", opt$modules_v1)
if (!is.null(opt$modules_v2) && !file.exists(opt$modules_v2)) stop("modules_v2 not found: ", opt$modules_v2)
if (!is.null(opt$modules_v3) && !file.exists(opt$modules_v3)) stop("modules_v3 not found: ", opt$modules_v3)
if (!is.null(opt$trait_exclude_file) && !file.exists(opt$trait_exclude_file)) {
  stop("trait_exclude_file not found: ", opt$trait_exclude_file)
}
if (!is.null(opt$outcome_traits_file) && !file.exists(opt$outcome_traits_file)) {
  stop("outcome_traits_file not found: ", opt$outcome_traits_file)
}

run_pearson <- parse_bool(opt$run_pearson, "--run_pearson")
run_bicor <- parse_bool(opt$run_bicor, "--run_bicor")

if (!run_pearson && !run_bicor) {
  stop("At least one of --run_pearson or --run_bicor must be TRUE.")
}

if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$top_n < 1) stop("--top_n must be >= 1")

module_dendro_distance <- tolower(opt$module_dendro_distance)
if (!module_dendro_distance %in% c("bicor", "pearson")) {
  stop("--module_dendro_distance must be 'bicor' or 'pearson'")
}

# ------------------------------------------------------------
# Configure cache and threads
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)
WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Derive lineage from modules_v1
# ------------------------------------------------------------
v1_variant_dir <- dirname(opt$modules_v1)
v1_region_dir <- dirname(v1_variant_dir)
region_label <- basename(v1_region_dir)
cpg_label <- basename(dirname(v1_region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "09a_me_trait_analysis")
out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Load sample info
# ------------------------------------------------------------
sample_info <- readSampleInfo(
  file = opt$sample_info,
  sample_id_col = opt$sample_id_col,
  verbose = TRUE
)

if (nrow(sample_info) < 2) stop("Sample info must contain at least 2 samples after loading.")
if (is.null(rownames(sample_info))) stop("Sample info must have rownames after ID resolution.")
if (anyDuplicated(rownames(sample_info))) {
  dup_ids <- unique(rownames(sample_info)[duplicated(rownames(sample_info))])
  stop("Sample info has duplicated sample IDs. Example: ",
       paste(head(dup_ids, 10), collapse = ", "))
}

# ------------------------------------------------------------
# Optional trait exclusion
# ------------------------------------------------------------
trait_exclude_requested <- readTraitFile(opt$trait_exclude_file, verbose = TRUE)
trait_exclude_resolved <- resolveTraits(
  requested_traits = trait_exclude_requested,
  available_traits = colnames(sample_info),
  label = "trait_exclude",
  verbose = TRUE
)

if (length(trait_exclude_resolved$found) > 0) {
  sample_info <- sample_info[, !colnames(sample_info) %in% trait_exclude_resolved$found, drop = FALSE]
}

# ------------------------------------------------------------
# Keep numeric usable traits only
# ------------------------------------------------------------
is_numeric_atomic <- vapply(sample_info, function(x) is.numeric(x) && is.atomic(x), logical(1))
removed_non_numeric <- names(sample_info)[!is_numeric_atomic]
sample_info_num <- sample_info[, is_numeric_atomic, drop = FALSE]

if (ncol(sample_info_num) == 0) {
  stop("No numeric trait columns available in sample_info after filtering.")
}

all_na <- vapply(sample_info_num, function(x) all(is.na(x)), logical(1))
zero_var <- vapply(sample_info_num, function(x) {
  vals <- x[!is.na(x)]
  if (length(vals) <= 1) return(TRUE)
  length(unique(vals)) == 1
}, logical(1))

removed_all_na <- names(sample_info_num)[all_na]
removed_zero_var <- names(sample_info_num)[zero_var]
sample_info_num <- sample_info_num[, !(all_na | zero_var), drop = FALSE]

if (ncol(sample_info_num) == 0) {
  stop("No usable numeric traits remain after removing all-NA and zero-variance columns.")
}

# ------------------------------------------------------------
# Optional outcome traits
# ------------------------------------------------------------
outcome_traits_requested <- readTraitFile(opt$outcome_traits_file, verbose = TRUE)
outcome_traits_resolved <- resolveTraits(
  requested_traits = outcome_traits_requested,
  available_traits = colnames(sample_info_num),
  label = "outcome_traits",
  verbose = TRUE
)

# ------------------------------------------------------------
# Build variant inputs
# ------------------------------------------------------------
variant_inputs <- list(
  v1_all_pcs = opt$modules_v1
)

if (!is.null(opt$modules_v2)) {
  variant_inputs[[basename(dirname(opt$modules_v2))]] <- opt$modules_v2
}

if (!is.null(opt$modules_v3)) {
  variant_inputs[[basename(dirname(opt$modules_v3))]] <- opt$modules_v3
}

# ------------------------------------------------------------
# Run per variant
# ------------------------------------------------------------
for (variant_name in names(variant_inputs)) {
  message("\n==============================")
  message("Running ME-trait analysis for variant: ", variant_name)
  message("==============================\n")

  modules <- validate_modules_object(
    readRDS(variant_inputs[[variant_name]]),
    paste0(variant_name, " modules object")
  )

  MEs <- modules$MEs
  message("[", variant_name, "] MEs dimensions: ", nrow(MEs), " samples x ", ncol(MEs), " modules")

  common_samples <- intersect(rownames(MEs), rownames(sample_info_num))
  if (length(common_samples) == 0) {
    stop("[", variant_name, "] No overlapping samples between MEs and sample_info.")
  }

  MEs_use <- MEs[common_samples, , drop = FALSE]
  traits_use <- sample_info_num[common_samples, , drop = FALSE]

  variant_out_dir <- file.path(out_dir, variant_name)
  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  write_vector_file(removed_non_numeric, file.path(variant_out_dir, "sample_info_non_numeric_columns_removed.txt"))
  write_vector_file(removed_all_na, file.path(variant_out_dir, "sample_info_all_na_numeric_columns_removed.txt"))
  write_vector_file(removed_zero_var, file.path(variant_out_dir, "sample_info_zero_variance_numeric_columns_removed.txt"))

  if (!is.null(opt$trait_exclude_file)) {
    write_vector_file(trait_exclude_resolved$requested, file.path(variant_out_dir, "traits_requested_exclude.txt"))
    write_vector_file(trait_exclude_resolved$found, file.path(variant_out_dir, "traits_found_exclude.txt"))
    write_vector_file(trait_exclude_resolved$missing, file.path(variant_out_dir, "traits_missing_exclude.txt"))
  }

  if (!is.null(opt$outcome_traits_file)) {
    write_vector_file(outcome_traits_resolved$requested, file.path(variant_out_dir, "traits_requested_outcome.txt"))
    write_vector_file(outcome_traits_resolved$found, file.path(variant_out_dir, "traits_found_outcome.txt"))
    write_vector_file(outcome_traits_resolved$missing, file.path(variant_out_dir, "traits_missing_outcome.txt"))
  }

  moduleDendro <- getDendro(MEs_use, distance = module_dendro_distance)

  if (run_pearson) {
    message("[", variant_name, "] Running Pearson ME-trait correlations")
    MEtraitCor_pearson <- getMEtraitCor(
      MEs_use,
      colData = traits_use,
      corType = "pearson",
      file = tempfile(fileext = ".tsv")
    )

    save_me_trait_method_outputs(
      MEtraitCor = MEtraitCor_pearson,
      method_name = "Pearson",
      out_dir = variant_out_dir,
      moduleDendro = moduleDendro,
      p_thresh = opt$p_thresh,
      top_n = opt$top_n,
      outcome_traits_found = outcome_traits_resolved$found
    )
  }

  if (run_bicor) {
    message("[", variant_name, "] Running bicor ME-trait correlations")
    MEtraitCor_bicor <- getMEtraitCor(
      MEs_use,
      colData = traits_use,
      corType = "bicor",
      file = tempfile(fileext = ".tsv")
    )

    save_me_trait_method_outputs(
      MEtraitCor = MEtraitCor_bicor,
      method_name = "Bicor",
      out_dir = variant_out_dir,
      moduleDendro = moduleDendro,
      p_thresh = opt$p_thresh,
      top_n = opt$top_n,
      outcome_traits_found = outcome_traits_resolved$found
    )
  }

  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("modules_file:", variant_inputs[[variant_name]]),
      paste("sample_info:", opt$sample_info),
      paste("sample_id_col:", ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
      paste("run_pearson:", run_pearson),
      paste("run_bicor:", run_bicor),
      paste("p_thresh:", opt$p_thresh),
      paste("top_n:", opt$top_n),
      paste("module_dendro_distance:", module_dendro_distance),
      paste("trait_exclude_file:", ifelse(is.null(opt$trait_exclude_file), "NULL", opt$trait_exclude_file)),
      paste("outcome_traits_file:", ifelse(is.null(opt$outcome_traits_file), "NULL", opt$outcome_traits_file)),
      paste("n_samples_overlap:", length(common_samples)),
      paste("n_modules:", ncol(MEs_use)),
      paste("n_traits_tested:", ncol(traits_use)),
      paste("n_outcome_traits_found:", length(outcome_traits_resolved$found))
    ),
    file.path(variant_out_dir, "run_parameters.txt")
  )

  message("✓ Finished variant: ", variant_name)
  message("  Outputs in: ", variant_out_dir)
}

# ------------------------------------------------------------
# Save top-level run parameters
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("sample_info:", opt$sample_info),
    paste("sample_id_col:", ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
    paste("modules_v1:", opt$modules_v1),
    paste("modules_v2:", ifelse(is.null(opt$modules_v2), "NULL", opt$modules_v2)),
    paste("modules_v3:", ifelse(is.null(opt$modules_v3), "NULL", opt$modules_v3)),
    paste("run_pearson:", run_pearson),
    paste("run_bicor:", run_bicor),
    paste("p_thresh:", opt$p_thresh),
    paste("top_n:", opt$top_n),
    paste("module_dendro_distance:", module_dendro_distance),
    paste("trait_exclude_file:", ifelse(is.null(opt$trait_exclude_file), "NULL", opt$trait_exclude_file)),
    paste("outcome_traits_file:", ifelse(is.null(opt$outcome_traits_file), "NULL", opt$outcome_traits_file)),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL CORE ME-TRAIT ANALYSES COMPLETE ✓")
message("Outputs saved under:\n  ", out_dir)