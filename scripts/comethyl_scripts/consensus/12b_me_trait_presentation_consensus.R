#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 12B: Consensus ME-Trait Presentation Heatmaps
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Generates focused, publication-ready ME-trait heatmaps by
#   subsetting the full stats from script 12a to user-defined
#   trait sets:
#     1) Load ME-trait stats per dataset from script 12a
#     2) For each trait-set file, subset to traits of interest
#     3) Generate per-dataset full and top-N heatmaps
#     4) Generate cross-dataset comparison heatmaps (same trait
#        set, different datasets side-by-side in cross_dataset/)
#
# REQUIRED INPUTS
#   --project_root        : root directory of the analysis project
#   --regions_file        : path to Filtered_Regions.txt from script 02
#   --dataset1_label      : label for dataset 1
#   --dataset1_stats_file : path to ME-trait stats TSV from script 12a
#
# OPTIONAL INPUTS
#   --consensus_modules_rds  : path to Consensus_Modules.rds (for module ordering)
#   --dataset2_label         : label for dataset 2
#   --dataset2_stats_file    : stats TSV for dataset 2
#   --dataset3_label         : label for dataset 3
#   --dataset3_stats_file    : stats TSV for dataset 3
#   --adjustment_version     : label matching script 12a run [default = unadjusted]
#   --set_dir                : directory of .txt trait-set files (all processed)
#   --set_file               : path to single trait-set .txt file
#   --set_file2 to set_file5 : additional trait-set files
#   --module_dendro_distance : bicor or pearson [default = bicor]
#   --p_thresh               : significance threshold [default = 0.05]
#   --top_n                  : top N associations for top heatmap [default = 250]
#   --full_width/full_height : full heatmap dimensions in inches [default = 12 x 12]
#   --top_width/top_height   : top heatmap dimensions in inches [default = 9 x 6]
#
# OUTPUTS
#   comethyl_output/consensus/12b_me_trait_presentation/
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           <set_name>/
#               <dataset_label>_<set_name>_ME_Trait_Heatmap_FULL.pdf
#               <dataset_label>_<set_name>_ME_Trait_Heatmap_TOP.pdf
#               <dataset_label>_<set_name>_subset_stats.tsv
#               <dataset_label>_<set_name>_subset_stats_significant.tsv
#               <dataset_label>_<set_name>_traits_requested/found/missing.txt
#               run_parameters.txt
#       cross_dataset/<cpg_label>/<region_label>/<adjustment_version>/
#           <set_name>/
#               <set_name>_Comparison_<dataset_label>_FULL.pdf
#               <set_name>_Comparison_<dataset_label>_TOP.pdf
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - Trait sets are plain text files with one trait name per line
#   - If --set_dir is provided, all .txt files in that directory are processed
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/12b_me_trait_presentation_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --dataset1_label Baseline \
#     --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --dataset2_label 36_38wks \
#     --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/36_38wks/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --adjustment_version v1_all_pcs \
#     --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/trait_sets/
# ================================================================
# SCRIPT 12B: Consensus ME-Trait Presentation Heatmaps
#
# PURPOSE
#   - Load ME-trait stats from 12A for one or more datasets
#   - Generate focused heatmaps per dataset per trait set
#   - Generate comparison-ready heatmaps per trait set across datasets
# ================================================================

message("Starting Script 12b ✓")

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

required_helpers <- c("write_vector_file", "write_log_lines", "getDendro", "plotMEtraitCor")
missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]
if (length(missing_helpers) > 0) {
  stop("Missing helper functions in helper.R: ", paste(missing_helpers, collapse = ", "))
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
read_trait_set_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  unique(trimws(x[nzchar(trimws(x))]))
}

get_set_name <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

collect_set_files <- function(set_dir, set_file, set_file2, set_file3, set_file4, set_file5) {
  direct <- c(set_file, set_file2, set_file3, set_file4, set_file5)
  direct <- direct[!is.na(direct) & nzchar(direct)]

  if (length(direct) > 0) {
    missing <- direct[!file.exists(direct)]
    if (length(missing) > 0) {
      stop("Trait set files not found:\n  ", paste(missing, collapse = "\n  "))
    }
  }

  from_dir <- character(0)
  if (!is.null(set_dir) && nzchar(set_dir)) {
    if (!dir.exists(set_dir)) stop("set_dir not found: ", set_dir)
    from_dir <- list.files(set_dir, pattern = "\\.txt$", full.names = TRUE)
    if (length(from_dir) == 0) stop("No .txt files found in set_dir: ", set_dir)
  }

  all_files <- unique(normalizePath(c(direct, from_dir), mustWork = TRUE))

  if (length(all_files) == 0) {
    stop("No trait set files provided. Use --set_dir or --set_file.")
  }

  all_files
}

validate_stats_table <- function(df, label) {
  required <- c("module", "trait", "p")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "))
  }

  has_effect_col <- any(c("bicor", "cor") %in% colnames(df))
  if (!has_effect_col) {
    stop(label, " must contain either 'bicor' or 'cor' column")
  }

  df
}

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",           type = "character"),
  make_option("--consensus_modules_rds",  type = "character", default = NULL),
  make_option("--regions_file",           type = "character"),
  make_option("--adjustment_version",     type = "character", default = "unadjusted"),

  make_option("--dataset1_label",         type = "character"),
  make_option("--dataset1_stats_file",    type = "character"),
  make_option("--dataset2_label",         type = "character", default = NULL),
  make_option("--dataset2_stats_file",    type = "character", default = NULL),
  make_option("--dataset3_label",         type = "character", default = NULL),
  make_option("--dataset3_stats_file",    type = "character", default = NULL),

  make_option("--set_dir",                type = "character", default = NULL),
  make_option("--set_file",               type = "character", default = NULL),
  make_option("--set_file2",              type = "character", default = NULL),
  make_option("--set_file3",              type = "character", default = NULL),
  make_option("--set_file4",              type = "character", default = NULL),
  make_option("--set_file5",              type = "character", default = NULL),

  make_option("--module_dendro_distance", type = "character", default = "bicor"),
  make_option("--p_thresh",               type = "double",    default = 0.05),
  make_option("--top_n",                  type = "integer",   default = 250),
  make_option("--full_width",             type = "double",    default = 12),
  make_option("--full_height",            type = "double",    default = 12),
  make_option("--top_width",              type = "double",    default = 9),
  make_option("--top_height",             type = "double",    default = 6)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root))        stop("--project_root is required")
if (is.null(opt$regions_file))        stop("--regions_file is required")
if (is.null(opt$dataset1_label))      stop("--dataset1_label is required")
if (is.null(opt$dataset1_stats_file)) stop("--dataset1_stats_file is required")

if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$regions_file)) stop("regions_file not found: ", opt$regions_file)
if (!file.exists(opt$dataset1_stats_file)) stop("dataset1_stats_file not found: ", opt$dataset1_stats_file)

module_dendro_distance <- tolower(opt$module_dendro_distance)
if (!module_dendro_distance %in% c("bicor", "pearson")) {
  stop("--module_dendro_distance must be bicor or pearson")
}
if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$top_n < 1) stop("--top_n must be >= 1")

dataset_inputs <- list(
  list(label = opt$dataset1_label, stats_file = opt$dataset1_stats_file)
)

for (i in 2:3) {
  lbl <- opt[[paste0("dataset", i, "_label")]]
  sf  <- opt[[paste0("dataset", i, "_stats_file")]]

  if (!is.null(lbl) || !is.null(sf)) {
    if (is.null(lbl) || is.null(sf)) {
      stop("Provide both --dataset", i, "_label and --dataset", i, "_stats_file")
    }
    if (!file.exists(sf)) {
      stop("dataset", i, "_stats_file not found: ", sf)
    }
    dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = lbl, stats_file = sf)
  }
}

set_files <- collect_set_files(
  set_dir   = opt$set_dir,
  set_file  = opt$set_file,
  set_file2 = opt$set_file2,
  set_file3 = opt$set_file3,
  set_file4 = opt$set_file4,
  set_file5 = opt$set_file5
)

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
step_dir      <- file.path(pipeline_root, "12b_me_trait_presentation")
cross_dir     <- file.path(step_dir, "cross_dataset", cpg_label, region_label, opt$adjustment_version)

dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)
message("Output root: ", step_dir)

# ------------------------------------------------------------
# Optional module ordering from consensus RDS
# ------------------------------------------------------------
module_order <- NULL

if (!is.null(opt$consensus_modules_rds)) {
  if (!file.exists(opt$consensus_modules_rds)) {
    stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
  }

  consensusMods <- readRDS(opt$consensus_modules_rds)

  if (is.null(consensusMods$multiMEs)) {
    stop("Consensus object does not contain $multiMEs")
  }

  first_label <- dataset_inputs[[1]]$label
  if (!first_label %in% names(consensusMods$multiMEs)) {
    stop(
      "dataset1 label '", first_label,
      "' not found in consensusMods$multiMEs. Available: ",
      paste(names(consensusMods$multiMEs), collapse = ", ")
    )
  }

  MEs_ref <- consensusMods$multiMEs[[first_label]]$data
  MEs_ref <- as.data.frame(MEs_ref)

  grey_col <- grep("^MEgrey$", colnames(MEs_ref), value = TRUE)
  if (length(grey_col) > 0) {
    MEs_ref <- MEs_ref[, !colnames(MEs_ref) %in% grey_col, drop = FALSE]
  }

  if (ncol(MEs_ref) >= 2) {
    moduleDendro <- getDendro(MEs_ref, distance = module_dendro_distance)
    module_order <- moduleDendro$order
    message("Module ordering derived from consensus MEs using dataset: ", first_label)
  } else {
    message("Not enough non-grey modules to derive dendrogram-based module order")
  }
}

# ------------------------------------------------------------
# Load all stats tables
# ------------------------------------------------------------
stats_list <- list()

for (ds in dataset_inputs) {
  df <- readr::read_tsv(ds$stats_file, show_col_types = FALSE)
  df <- validate_stats_table(df, paste0(ds$label, " stats_file"))

  df$module <- as.character(df$module)
  df$trait  <- as.character(df$trait)

  stats_list[[ds$label]] <- df

  message(
    "Loaded stats for: ", ds$label,
    " (", nrow(df), " rows, ",
    length(unique(df$module)), " modules, ",
    length(unique(df$trait)), " traits)"
  )
}

# ------------------------------------------------------------
# Per-dataset presentation heatmaps + comparison-ready outputs
# ------------------------------------------------------------
for (set_file in set_files) {

  set_name <- get_set_name(set_file)
  requested_traits <- read_trait_set_file(set_file)

  message("\n--- Trait set: ", set_name, " ---")

  # ----------------------------------------------------------
  # Per-dataset outputs
  # ----------------------------------------------------------
  for (ds in dataset_inputs) {

    ds_label   <- ds$label
    ds_set_dir <- file.path(step_dir, ds_label, cpg_label, region_label, opt$adjustment_version, set_name)
    dir.create(ds_set_dir, recursive = TRUE, showWarnings = FALSE)

    stats_df <- stats_list[[ds_label]]
    available_traits <- unique(stats_df$trait)
    found_traits <- intersect(requested_traits, available_traits)
    missing_traits <- setdiff(requested_traits, available_traits)

    write_vector_file(
      requested_traits,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_traits_requested.txt"))
    )
    write_vector_file(
      found_traits,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_traits_found.txt"))
    )
    write_vector_file(
      missing_traits,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_traits_missing.txt"))
    )

    if (length(found_traits) == 0) {
      message(ds_label, "/", set_name, ": no traits found — skipping")
      next
    }

    subset_df <- dplyr::filter(stats_df, trait %in% found_traits)

    # Re-factor explicitly for plotMEtraitCor()
    subset_df$module <- factor(
      as.character(subset_df$module),
      levels = unique(as.character(stats_df$module))
    )

    subset_df$trait <- factor(
      as.character(subset_df$trait),
      levels = found_traits
    )

    subset_sig <- dplyr::filter(subset_df, !is.na(p), p < opt$p_thresh)

    readr::write_tsv(
      subset_df,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_subset_stats.tsv"))
    )
    readr::write_tsv(
      subset_sig,
      file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_subset_stats_significant.tsv"))
    )

    full_pdf <- file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_ME_Trait_Heatmap_FULL.pdf"))
    top_pdf  <- file.path(ds_set_dir, paste0(ds_label, "_", set_name, "_ME_Trait_Heatmap_TOP.pdf"))

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = module_order,
        p               = opt$p_thresh,
        topOnly         = FALSE,
        file            = full_pdf,
        width           = opt$full_width,
        height          = opt$full_height,
        colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
      )
      message(ds_label, "/", set_name, ": saved FULL heatmap")
    }, error = function(e) {
      message(ds_label, "/", set_name, ": FULL heatmap failed — ", conditionMessage(e))
    })

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = module_order,
        p               = opt$p_thresh,
        topOnly         = TRUE,
        nTop            = opt$top_n,
        label.type      = "p",
        label.size      = 4,
        label.nudge_y   = 0,
        legend.position = c(1.11, 0.795),
        colColorMargins = c(-1, 4.75, 0.5, 10.1),
        file            = top_pdf,
        width           = opt$top_width,
        height          = opt$top_height
      )
      message(ds_label, "/", set_name, ": saved TOP heatmap")
    }, error = function(e) {
      message(ds_label, "/", set_name, ": TOP heatmap failed — ", conditionMessage(e))
    })

    write_log_lines(
      c(
        paste("dataset_label:", ds_label),
        paste("set_name:", set_name),
        paste("set_file:", set_file),
        paste("stats_file:", ds$stats_file),
        paste("p_thresh:", opt$p_thresh),
        paste("top_n:", opt$top_n),
        paste("n_traits_requested:", length(requested_traits)),
        paste("n_traits_found:", length(found_traits)),
        paste("n_traits_missing:", length(missing_traits)),
        paste("n_rows_subset:", nrow(subset_df)),
        paste("n_rows_significant:", nrow(subset_sig)),
        paste("date:", as.character(Sys.time()))
      ),
      file.path(ds_set_dir, "run_parameters.txt")
    )
  }

  # ----------------------------------------------------------
  # Comparison-ready outputs in shared folder
  # ----------------------------------------------------------
  cross_set_dir <- file.path(cross_dir, set_name)
  dir.create(cross_set_dir, recursive = TRUE, showWarnings = FALSE)

  n_cross_saved <- 0

  for (ds in dataset_inputs) {
    ds_label <- ds$label
    stats_df <- stats_list[[ds_label]]

    found_traits <- intersect(requested_traits, unique(stats_df$trait))
    if (length(found_traits) == 0) {
      message("Cross/", set_name, "/", ds_label, ": no traits found — skipping")
      next
    }

    subset_df <- dplyr::filter(stats_df, trait %in% found_traits)

    # Re-factor explicitly for plotMEtraitCor()
    subset_df$module <- factor(
      as.character(subset_df$module),
      levels = unique(as.character(stats_df$module))
    )

    subset_df$trait <- factor(
      as.character(subset_df$trait),
      levels = found_traits
    )
    full_pdf <- file.path(cross_set_dir, paste0(set_name, "_Comparison_", ds_label, "_FULL.pdf"))
    top_pdf  <- file.path(cross_set_dir, paste0(set_name, "_Comparison_", ds_label, "_TOP.pdf"))

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = module_order,
        p               = opt$p_thresh,
        topOnly         = FALSE,
        file            = full_pdf,
        width           = opt$full_width,
        height          = opt$full_height,
        colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
      )
      message("Cross/", set_name, "/", ds_label, ": saved FULL")
    }, error = function(e) {
      message("Cross/", set_name, "/", ds_label, ": FULL failed — ", conditionMessage(e))
    })

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = module_order,
        p               = opt$p_thresh,
        topOnly         = TRUE,
        nTop            = opt$top_n,
        label.type      = "p",
        label.size      = 4,
        label.nudge_y   = 0,
        legend.position = c(1.11, 0.795),
        colColorMargins = c(-1, 4.75, 0.5, 10.1),
        file            = top_pdf,
        width           = opt$top_width,
        height          = opt$top_height
      )
      message("Cross/", set_name, "/", ds_label, ": saved TOP")
    }, error = function(e) {
      message("Cross/", set_name, "/", ds_label, ": TOP failed — ", conditionMessage(e))
    })

    n_cross_saved <- n_cross_saved + 1
  }

  write_log_lines(
    c(
      paste("set_name:", set_name),
      paste("set_file:", set_file),
      paste("datasets:", paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
      paste("n_traits_requested:", length(requested_traits)),
      paste("n_datasets_saved:", n_cross_saved),
      paste("p_thresh:", opt$p_thresh),
      paste("top_n:", opt$top_n),
      paste("date:", as.character(Sys.time()))
    ),
    file.path(cross_set_dir, "run_parameters.txt")
  )
}

# ------------------------------------------------------------
# Top-level log
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("consensus_modules_rds:", ifelse(is.null(opt$consensus_modules_rds), "NULL", opt$consensus_modules_rds)),
    paste("adjustment_version:", opt$adjustment_version),
    paste("datasets:", paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
    paste("set_files:", paste(set_files, collapse = ", ")),
    paste("module_dendro_distance:", module_dendro_distance),
    paste("p_thresh:", opt$p_thresh),
    paste("top_n:", opt$top_n),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(step_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(step_dir, "sessionInfo.txt")
)

message("\n✓ Script 12B complete: consensus ME-trait presentation heatmaps finished")
message("Outputs saved under: ", step_dir)