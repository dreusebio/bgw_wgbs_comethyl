#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 12A: Consensus ME-Trait Analysis
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Associates consensus module eigengenes with sample traits
#   per dataset using the consensus MEs from script 09:
#     1) Load Consensus_Modules.rds and extract per-dataset MEs
#     2) Load sample info per dataset and retain numeric traits
#     3) Compute ME-trait correlations (Pearson and/or bicor)
#     4) Save full and significant-only stats tables per dataset
#     5) Save full and top-N ME-trait heatmaps per dataset
#     6) Optionally filter to outcome traits of interest
#
# REQUIRED INPUTS
#   --project_root          : root directory of the analysis project
#   --consensus_modules_rds : path to Consensus_Modules.rds from script 09
#   --regions_file          : path to Filtered_Regions.txt from script 02
#   --dataset1_label        : label for dataset 1
#   --dataset1_sample_info  : path to sample info file for dataset 1
#
# OPTIONAL INPUTS
#   --dataset2_label         : label for dataset 2
#   --dataset2_sample_info   : sample info file for dataset 2
#   --dataset3_label         : label for dataset 3
#   --dataset3_sample_info   : sample info file for dataset 3
#   --sample_id_col          : sample ID column in sample info files
#   --adjustment_version     : label matching script 09 run [default = unadjusted]
#   --run_pearson            : run Pearson correlations TRUE/FALSE [default = FALSE]
#   --run_bicor              : run bicor correlations TRUE/FALSE [default = TRUE]
#   --p_thresh               : significance threshold [default = 0.05]
#   --top_n                  : top N associations for heatmap [default = 250]
#   --module_dendro_distance : bicor or pearson for module ordering [default = bicor]
#   --trait_exclude_file     : text file with traits to exclude
#   --outcome_traits_file    : text file with outcome traits of interest
#   --max_p_outliers         : bicor outlier fraction [default = 0.1]
#
# OUTPUTS
#   comethyl_output/consensus/12a_me_trait_analysis/
#       shared/<cpg_label>/<region_label>/<adjustment_version>/
#           run_parameters.txt
#           sessionInfo.txt
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           ME_Trait_Correlation_Stats_Bicor.tsv
#           ME_Trait_Correlation_Stats_Bicor_significant.tsv
#           ME_Trait_Correlation_Stats_Bicor_outcome_only.tsv       (if outcome traits provided)
#           ME_Trait_Correlations_Bicor.xlsx
#           ME_Trait_Correlation_Heatmap_Bicor_FULL.pdf
#           ME_Trait_Correlation_Heatmap_Bicor_TOP.pdf
#           (same for Pearson if --run_pearson TRUE)
#           sample_info_*_removed.txt
#           run_parameters.txt
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - The stats TSV outputs feed into scripts 12b, 14, and 16
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/12a_me_trait_analysis_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --dataset1_label Baseline \
#     --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/Baseline_sample_info.xlsx \
#     --dataset2_label 36_38wks \
#     --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/36_38wks_sample_info.xlsx \
#     --dataset3_label Postpartum \
#     --dataset3_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/Postpartum_sample_info.xlsx \
#     --adjustment_version v1_all_pcs \
#     --run_bicor TRUE \
#     --p_thresh 0.05 \
#     --outcome_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/outcome_traits.txt
# ================================================================
# SCRIPT 12A: Consensus ME-Trait Analysis
#
# PURPOSE
#   - Load Consensus_Modules.rds from script 09
#   - Load per-dataset sample info files
#   - Compute ME-trait correlations per dataset
#   - Save full stats tables, significant-only tables
#   - Save Excel workbooks
#   - Save full and top heatmaps per dataset
# ================================================================

message("Starting Script 12a ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(readr)
  library(dplyr)
})

# ------------------------------------------------------------
# Load helper.R
# ------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
]

if (length(script_file_arg) == 0) {
  stop("Could not determine script path.")
}

script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")

if (!file.exists(helper_file)) {
  stop("helper.R not found: ", helper_file)
}

source(helper_file)

required_helpers <- c(
  "readTraitFile",
  "resolveTraits",
  "readSampleInfo",
  "write_vector_file",
  "write_log_lines",
  "save_me_trait_method_outputs",
  "getDendro"
)

missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
if (length(missing_helpers) > 0) {
  stop("Missing helper functions in helper.R: ", paste(missing_helpers, collapse = ", "))
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE")
}

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",           type = "character"),
  make_option("--consensus_modules_rds",  type = "character"),
  make_option("--regions_file",           type = "character"),
  make_option("--adjustment_version",     type = "character", default = "unadjusted"),

  make_option("--dataset1_label",         type = "character"),
  make_option("--dataset1_sample_info",   type = "character"),
  make_option("--dataset2_label",         type = "character", default = NULL),
  make_option("--dataset2_sample_info",   type = "character", default = NULL),
  make_option("--dataset3_label",         type = "character", default = NULL),
  make_option("--dataset3_sample_info",   type = "character", default = NULL),

  make_option("--sample_id_col",          type = "character", default = NULL),
  make_option("--run_pearson",            type = "character", default = "FALSE"),
  make_option("--run_bicor",              type = "character", default = "TRUE"),
  make_option("--p_thresh",               type = "double",    default = 0.05),
  make_option("--top_n",                  type = "integer",   default = 250),
  make_option("--module_dendro_distance", type = "character", default = "bicor"),
  make_option("--trait_exclude_file",     type = "character", default = NULL),
  make_option("--outcome_traits_file",    type = "character", default = NULL),
  make_option("--max_p_outliers",         type = "double",    default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root))          stop("--project_root is required")
if (is.null(opt$consensus_modules_rds)) stop("--consensus_modules_rds is required")
if (is.null(opt$regions_file))          stop("--regions_file is required")
if (is.null(opt$dataset1_label))        stop("--dataset1_label is required")
if (!dir.exists(opt$project_root))      stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
if (!file.exists(opt$regions_file))     stop("regions_file not found: ", opt$regions_file)
if (is.null(opt$dataset1_sample_info))  stop("--dataset1_sample_info is required")

if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
if (!file.exists(opt$dataset1_sample_info)) stop("dataset1_sample_info not found: ", opt$dataset1_sample_info)

run_pearson <- parse_bool(opt$run_pearson, "--run_pearson")
run_bicor   <- parse_bool(opt$run_bicor, "--run_bicor")

if (!run_pearson && !run_bicor) {
  stop("At least one of --run_pearson or --run_bicor must be TRUE")
}

if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$top_n < 1) stop("--top_n must be >= 1")

module_dendro_distance <- tolower(opt$module_dendro_distance)
if (!module_dendro_distance %in% c("bicor", "pearson")) {
  stop("--module_dendro_distance must be bicor or pearson")
}

if (!is.null(opt$trait_exclude_file) && !file.exists(opt$trait_exclude_file)) {
  stop("trait_exclude_file not found: ", opt$trait_exclude_file)
}
if (!is.null(opt$outcome_traits_file) && !file.exists(opt$outcome_traits_file)) {
  stop("outcome_traits_file not found: ", opt$outcome_traits_file)
}

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_sample_info)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_sample_info)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_sample_info)) {
    stop("Provide both --dataset2_label and --dataset2_sample_info")
  }
  if (!file.exists(opt$dataset2_sample_info)) {
    stop("dataset2_sample_info not found: ", opt$dataset2_sample_info)
  }
}

if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_sample_info)) {
    stop("Provide both --dataset3_label and --dataset3_sample_info")
  }
  if (!file.exists(opt$dataset3_sample_info)) {
    stop("dataset3_sample_info not found: ", opt$dataset3_sample_info)
  }
}

dataset_inputs <- list(
  list(label = opt$dataset1_label, sample_info_file = opt$dataset1_sample_info)
)
if (dataset2_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(
    label = opt$dataset2_label,
    sample_info_file = opt$dataset2_sample_info
  )
}
if (dataset3_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(
    label = opt$dataset3_label,
    sample_info_file = opt$dataset3_sample_info
  )
}

WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Output dirs
# ------------------------------------------------------------

# ----------------------------------------------------------------
# Derive cpg_label / region_label from --regions_file path:
#   .../02_reference_region_filter/<ds>/<cpg_label>/<region_label>/Filtered_Regions.txt
# ----------------------------------------------------------------
{
  rfile_dir    <- dirname(opt$regions_file)
  region_label <- basename(rfile_dir)
  cpg_label    <- basename(dirname(rfile_dir))
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "12a_me_trait_analysis")
shared_dir    <- file.path(step_dir, "shared", cpg_label, region_label, opt$adjustment_version)
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
message("Output root: ", step_dir)

# ------------------------------------------------------------
# Load consensus modules
# ------------------------------------------------------------
consensusMods <- readRDS(opt$consensus_modules_rds)

if (is.null(consensusMods$colors)) {
  stop("Consensus object does not contain $colors")
}
if (is.null(consensusMods$multiMEs)) {
  stop("Consensus object does not contain $multiMEs")
}

all_colors    <- consensusMods$colors
module_colors <- unique(all_colors[all_colors != "grey"])
n_modules     <- length(module_colors)

message("Non-grey consensus modules: ", n_modules,
        " (", paste(module_colors, collapse = ", "), ")")

if (n_modules == 0) {
  stop("No non-grey modules found. Cannot run ME-trait analysis.")
}

# ------------------------------------------------------------
# Optional trait files
# ------------------------------------------------------------
trait_exclude_requested  <- readTraitFile(opt$trait_exclude_file, verbose = TRUE)
outcome_traits_requested <- readTraitFile(opt$outcome_traits_file, verbose = TRUE)

# ------------------------------------------------------------
# Per-dataset loop
# ------------------------------------------------------------
for (ds in dataset_inputs) {

  ds_label <- ds$label
  ds_dir   <- file.path(step_dir, ds_label, cpg_label, region_label, opt$adjustment_version)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n==============================")
  message("Dataset: ", ds_label)
  message("==============================")

  if (!ds_label %in% names(consensusMods$multiMEs)) {
    stop(
      ds_label, ": not found in consensusMods$multiMEs. Available datasets: ",
      paste(names(consensusMods$multiMEs), collapse = ", ")
    )
  }

  MEs_full <- consensusMods$multiMEs[[ds_label]]$data
  if (is.null(MEs_full)) {
    stop(ds_label, ": no eigengene data found in consensusMods$multiMEs[[ds_label]]$data")
  }

  MEs_full <- as.data.frame(MEs_full)

  if (is.null(rownames(MEs_full))) {
    stop(ds_label, ": MEs must have rownames as sample IDs")
  }

  grey_col <- grep("^MEgrey$", colnames(MEs_full), value = TRUE)
  MEs_full <- MEs_full[, !colnames(MEs_full) %in% grey_col, drop = FALSE]

  if (ncol(MEs_full) == 0) {
    message(ds_label, ": No real MEs after removing grey — skipping")
    next
  }

  message(ds_label, ": ", nrow(MEs_full), " samples x ", ncol(MEs_full), " modules")

  sample_info <- readSampleInfo(
    file          = ds$sample_info_file,
    sample_id_col = opt$sample_id_col,
    verbose       = TRUE
  )

  if (nrow(sample_info) < 2) {
    stop(ds_label, ": sample_info must have >= 2 samples")
  }
  if (is.null(rownames(sample_info))) {
    stop(ds_label, ": sample_info must have rownames as sample IDs")
  }

  trait_exclude_resolved <- resolveTraits(
    requested_traits = trait_exclude_requested,
    available_traits = colnames(sample_info),
    label            = paste0(ds_label, "_trait_exclude"),
    verbose          = TRUE
  )

  if (length(trait_exclude_resolved$found) > 0) {
    sample_info <- sample_info[
      , !colnames(sample_info) %in% trait_exclude_resolved$found, drop = FALSE
    ]
  }

  is_num <- vapply(sample_info, function(x) is.numeric(x) && is.atomic(x), logical(1))
  removed_non_numeric <- names(sample_info)[!is_num]
  sample_info_num <- sample_info[, is_num, drop = FALSE]

  if (ncol(sample_info_num) == 0) {
    stop(ds_label, ": No numeric trait columns available after filtering")
  }

  all_na <- vapply(sample_info_num, function(x) all(is.na(x)), logical(1))
  zero_var <- vapply(sample_info_num, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) <= 1) return(TRUE)
    length(unique(vals)) == 1
  }, logical(1))

  removed_all_na   <- names(sample_info_num)[all_na]
  removed_zero_var <- names(sample_info_num)[zero_var]

  sample_info_num <- sample_info_num[, !(all_na | zero_var), drop = FALSE]

  if (ncol(sample_info_num) == 0) {
    stop(ds_label, ": No usable numeric traits after removing all-NA and zero-variance columns")
  }

  outcome_traits_resolved <- resolveTraits(
    requested_traits = outcome_traits_requested,
    available_traits = colnames(sample_info_num),
    label            = paste0(ds_label, "_outcome_traits"),
    verbose          = TRUE
  )

  common_samples <- intersect(rownames(MEs_full), rownames(sample_info_num))
  if (length(common_samples) < 3) {
    stop(ds_label, ": fewer than 3 overlapping samples between MEs and sample_info")
  }

  MEs_use    <- MEs_full[common_samples, , drop = FALSE]
  traits_use <- sample_info_num[common_samples, , drop = FALSE]

  if (!identical(rownames(MEs_use), rownames(traits_use))) {
    stop(ds_label, ": sample alignment failed between MEs and traits")
  }

  message(ds_label, ": ", length(common_samples), " overlapping samples")
  message(ds_label, ": ", ncol(traits_use), " numeric traits retained")

  moduleDendro <- getDendro(MEs_use, distance = module_dendro_distance)

  write_vector_file(
    removed_non_numeric,
    file.path(ds_dir, paste0(ds_label, "_sample_info_non_numeric_columns_removed.txt"))
  )
  write_vector_file(
    removed_all_na,
    file.path(ds_dir, paste0(ds_label, "_sample_info_all_na_numeric_columns_removed.txt"))
  )
  write_vector_file(
    removed_zero_var,
    file.path(ds_dir, paste0(ds_label, "_sample_info_zero_variance_numeric_columns_removed.txt"))
  )

  if (!is.null(opt$trait_exclude_file)) {
    write_vector_file(
      trait_exclude_resolved$requested,
      file.path(ds_dir, paste0(ds_label, "_traits_requested_exclude.txt"))
    )
    write_vector_file(
      trait_exclude_resolved$found,
      file.path(ds_dir, paste0(ds_label, "_traits_found_exclude.txt"))
    )
    write_vector_file(
      trait_exclude_resolved$missing,
      file.path(ds_dir, paste0(ds_label, "_traits_missing_exclude.txt"))
    )
  }

  if (!is.null(opt$outcome_traits_file)) {
    write_vector_file(
      outcome_traits_resolved$requested,
      file.path(ds_dir, paste0(ds_label, "_traits_requested_outcome.txt"))
    )
    write_vector_file(
      outcome_traits_resolved$found,
      file.path(ds_dir, paste0(ds_label, "_traits_found_outcome.txt"))
    )
    write_vector_file(
      outcome_traits_resolved$missing,
      file.path(ds_dir, paste0(ds_label, "_traits_missing_outcome.txt"))
    )
  }

  for (method in c("pearson", "bicor")) {

    if (method == "pearson" && !run_pearson) next
    if (method == "bicor"   && !run_bicor)   next

    method_label <- if (method == "pearson") "Pearson" else "Bicor"
    message(ds_label, ": Running ", method_label, " ME-trait correlations")

    MEtraitCor <- getMEtraitCor(
      MEs_use,
      colData = traits_use,
      corType = method
    )

    save_me_trait_method_outputs(
      MEtraitCor           = MEtraitCor,
      method_name          = method_label,
      out_dir              = ds_dir,
      moduleDendro         = moduleDendro,
      p_thresh             = opt$p_thresh,
      top_n                = opt$top_n,
      outcome_traits_found = outcome_traits_resolved$found
    )

    message(ds_label, ": Saved ", method_label, " outputs")
  }

  write_log_lines(
    c(
      paste("dataset_label:", ds_label),
      paste("sample_info_file:", ds$sample_info_file),
      paste("consensus_modules_rds:", opt$consensus_modules_rds),
      paste("adjustment_version:", opt$adjustment_version),
      paste("run_pearson:", run_pearson),
      paste("run_bicor:", run_bicor),
      paste("p_thresh:", opt$p_thresh),
      paste("top_n:", opt$top_n),
      paste("module_dendro_distance:", module_dendro_distance),
      paste("n_samples_overlap:", length(common_samples)),
      paste("n_modules:", ncol(MEs_use)),
      paste("n_traits_tested:", ncol(traits_use)),
      paste("n_outcome_traits_found:", length(outcome_traits_resolved$found)),
      paste("date:", as.character(Sys.time()))
    ),
    file.path(ds_dir, "run_parameters.txt")
  )

  message("✓ Finished dataset: ", ds_label)
}

write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("consensus_modules_rds:", opt$consensus_modules_rds),
    paste("adjustment_version:", opt$adjustment_version),
    paste("datasets:", paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
    paste("run_pearson:", run_pearson),
    paste("run_bicor:", run_bicor),
    paste("p_thresh:", opt$p_thresh),
    paste("top_n:", opt$top_n),
    paste("module_dendro_distance:", module_dendro_distance),
    paste("n_real_modules:", n_modules),
    paste("module_colors:", paste(module_colors, collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(shared_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(shared_dir, "sessionInfo.txt")
)

message(" Script 12A complete: consensus ME-trait analysis finished")
message("Outputs saved under: ", step_dir)