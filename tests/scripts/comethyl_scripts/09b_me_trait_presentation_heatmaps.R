#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 09B: ME-Trait Presentation Heatmaps
#
# PURPOSE
#   - Load ME-trait correlation stats from 09A
#   - Optionally load modules RDS for module ordering
#   - Load one or more trait-set text files, either directly or from a directory
#   - Generate focused presentation heatmaps for those trait sets
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --stats_file   : path to one ME-trait stats file from 09A
#
# TRAIT SET INPUTS
#   Provide either:
#     --set_dir  : directory containing one or more .txt trait set files
#   and/or:
#     --set_file, --set_file2, ..., --set_file5
#
# OPTIONAL INPUTS
#   --modules_rds             : optional Modules.rds file for module ordering
#   --module_dendro_distance  : bicor or pearson [default = bicor]
#   --p_thresh                : significance threshold [default = 0.05]
#   --top_n                   : number of top associations for TOP heatmap [default = 250]
#   --full_width              : full heatmap width [default = 12]
#   --full_height             : full heatmap height [default = 12]
#   --top_width               : top heatmap width [default = 9]
#   --top_height              : top heatmap height [default = 6]
#
# OUTPUTS
#   project_root/comethyl_output/09b_me_trait_presentation/<cpg_label>/<region_label>/<variant>/<set_name>/
#       traits_requested.txt
#       traits_found.txt
#       traits_missing.txt
#       subset_stats.tsv
#       subset_stats_significant.tsv
#       ME_Trait_Heatmap_FULL.pdf
#       ME_Trait_Heatmap_TOP.pdf
#       run_parameters.txt
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
})

# ------------------------------------------------------------
# Load helper.R
# ------------------------------------------------------------
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
# Extra helpers
# ------------------------------------------------------------
read_trait_set_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

get_set_name <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

validate_stats_table <- function(df, label = "stats_file") {
  required <- c("module", "trait", "p")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!("bicor" %in% colnames(df) || "cor" %in% colnames(df))) {
    stop(label, " must contain either 'bicor' or 'cor' column.")
  }
  df
}

collect_set_files <- function(set_dir = NULL,
                              set_file = NULL,
                              set_file2 = NULL,
                              set_file3 = NULL,
                              set_file4 = NULL,
                              set_file5 = NULL) {

  direct_files <- c(set_file, set_file2, set_file3, set_file4, set_file5)
  direct_files <- direct_files[!is.na(direct_files) & nzchar(direct_files)]

  if (length(direct_files) > 0) {
    missing_direct <- direct_files[!file.exists(direct_files)]
    if (length(missing_direct) > 0) {
      stop("The following trait set files do not exist:\n  ",
           paste(missing_direct, collapse = "\n  "))
    }
  }

  dir_files <- character(0)
  if (!is.null(set_dir) && nzchar(set_dir)) {
    if (!dir.exists(set_dir)) {
      stop("set_dir does not exist: ", set_dir)
    }

    dir_files <- list.files(
      set_dir,
      pattern = "\\.txt$",
      full.names = TRUE
    )

    if (length(dir_files) == 0) {
      stop("No .txt files found in set_dir: ", set_dir)
    }
  }

  all_files <- c(direct_files, dir_files)

  if (length(all_files) == 0) {
    stop("No trait set files provided. Use --set_dir or at least one --set_file.")
  }

  all_files <- unique(normalizePath(all_files, mustWork = TRUE))
  all_files
}

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--stats_file", type = "character",
              help = "Path to one ME-trait stats file from 09A"),

  make_option("--set_dir", type = "character", default = NULL,
              help = "Directory containing one or more .txt trait set files"),

  make_option("--set_file", type = "character", default = NULL,
              help = "Trait set text file, one trait per line"),

  make_option("--modules_rds", type = "character", default = NULL,
              help = "Optional Modules.rds file for module ordering"),

  make_option("--set_file2", type = "character", default = NULL),
  make_option("--set_file3", type = "character", default = NULL),
  make_option("--set_file4", type = "character", default = NULL),
  make_option("--set_file5", type = "character", default = NULL),

  make_option("--module_dendro_distance", type = "character", default = "bicor",
              help = "bicor or pearson [default = bicor]"),

  make_option("--p_thresh", type = "double", default = 0.05,
              help = "Significance threshold [default = 0.05]"),

  make_option("--top_n", type = "integer", default = 250,
              help = "Top N associations for TOP heatmap [default = 250]"),

  make_option("--full_width", type = "double", default = 12,
              help = "Full heatmap width [default = 12]"),

  make_option("--full_height", type = "double", default = 12,
              help = "Full heatmap height [default = 12]"),

  make_option("--top_width", type = "double", default = 9,
              help = "Top heatmap width [default = 9]"),

  make_option("--top_height", type = "double", default = 6,
              help = "Top heatmap height [default = 6]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$stats_file)) stop("--stats_file is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$stats_file)) stop("stats_file not found: ", opt$stats_file)

if (!is.null(opt$modules_rds) && !file.exists(opt$modules_rds)) {
  stop("modules_rds not found: ", opt$modules_rds)
}

module_dendro_distance <- tolower(opt$module_dendro_distance)
if (!module_dendro_distance %in% c("bicor", "pearson")) {
  stop("--module_dendro_distance must be 'bicor' or 'pearson'")
}
if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$top_n < 1) stop("--top_n must be >= 1")

set_files <- collect_set_files(
  set_dir   = opt$set_dir,
  set_file  = opt$set_file,
  set_file2 = opt$set_file2,
  set_file3 = opt$set_file3,
  set_file4 = opt$set_file4,
  set_file5 = opt$set_file5
)

# ------------------------------------------------------------
# Configure cache
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)
WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Derive lineage from stats_file
# Expected:
#   .../09a_me_trait_analysis/<cpg_label>/<region_label>/<variant>/ME_Trait_Correlation_Stats_*.tsv
# ------------------------------------------------------------
variant_dir <- dirname(opt$stats_file)
region_dir <- dirname(variant_dir)
variant_name <- basename(variant_dir)
region_label <- basename(region_dir)
cpg_label <- basename(dirname(region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "09b_me_trait_presentation")
out_dir <- file.path(step_dir, cpg_label, region_label, variant_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Load stats
# ------------------------------------------------------------
stats_df <- readr::read_tsv(opt$stats_file, show_col_types = FALSE)
stats_df <- validate_stats_table(stats_df, "stats_file")

stats_df$module <- factor(as.character(stats_df$module), levels = unique(as.character(stats_df$module)))
stats_df$trait  <- factor(as.character(stats_df$trait), levels = unique(as.character(stats_df$trait)))

# ------------------------------------------------------------
# Optional module ordering
# ------------------------------------------------------------
module_order <- seq_along(levels(stats_df$module))

if (!is.null(opt$modules_rds)) {
  modules <- validate_modules_object(readRDS(opt$modules_rds), "modules_rds")
  MEs <- modules$MEs
  moduleDendro <- getDendro(MEs, distance = module_dendro_distance)

  dendro_modules <- colnames(MEs)[moduleDendro$order]
  current_levels <- levels(stats_df$module)
  module_order <- match(intersect(dendro_modules, current_levels), current_levels)
  module_order <- module_order[!is.na(module_order)]

  if (length(module_order) == 0) {
    module_order <- seq_along(current_levels)
  }
}

# ------------------------------------------------------------
# Run each trait set
# ------------------------------------------------------------
for (set_file in set_files) {
  set_name <- get_set_name(set_file)
  set_out_dir <- file.path(out_dir, set_name)
  dir.create(set_out_dir, recursive = TRUE, showWarnings = FALSE)

  requested_traits <- read_trait_set_file(set_file)
  available_traits <- unique(as.character(stats_df$trait))

  found_traits <- intersect(requested_traits, available_traits)
  missing_traits <- setdiff(requested_traits, available_traits)

  write_vector_file(requested_traits, file.path(set_out_dir, "traits_requested.txt"))
  write_vector_file(found_traits, file.path(set_out_dir, "traits_found.txt"))
  write_vector_file(missing_traits, file.path(set_out_dir, "traits_missing.txt"))

  if (length(found_traits) == 0) {
    write_log_lines(
      c(
        paste("set_name:", set_name),
        paste("set_file:", set_file),
        paste("stats_file:", opt$stats_file),
        paste("modules_rds:", ifelse(is.null(opt$modules_rds), "NULL", opt$modules_rds)),
        paste("p_thresh:", opt$p_thresh),
        paste("top_n:", opt$top_n),
        "status: skipped",
        "reason: no requested traits found in stats table"
      ),
      file.path(set_out_dir, "run_parameters.txt")
    )
    message("Skipping set ", set_name, ": no requested traits found in stats table.")
    next
  }

  subset_df <- stats_df %>%
    dplyr::filter(as.character(trait) %in% found_traits)

  subset_sig_df <- subset_df %>%
    dplyr::filter(!is.na(p), p < opt$p_thresh)

  readr::write_tsv(subset_df, file.path(set_out_dir, "subset_stats.tsv"))
  readr::write_tsv(subset_sig_df, file.path(set_out_dir, "subset_stats_significant.tsv"))

  full_pdf <- file.path(set_out_dir, "ME_Trait_Heatmap_FULL.pdf")
  top_pdf <- file.path(set_out_dir, "ME_Trait_Heatmap_TOP.pdf")

  plotMEtraitCor(
    subset_df,
    moduleOrder = module_order,
    p = opt$p_thresh,
    topOnly = FALSE,
    file = full_pdf,
    width = opt$full_width,
    height = opt$full_height,
    colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
  )

  plotMEtraitCor(
    subset_df,
    moduleOrder = module_order,
    p = opt$p_thresh,
    topOnly = TRUE,
    nTop = opt$top_n,
    label.type = "p",
    label.size = 4,
    label.nudge_y = 0,
    legend.position = c(1.11, 0.795),
    colColorMargins = c(-1, 4.75, 0.5, 10.1),
    file = top_pdf,
    width = opt$top_width,
    height = opt$top_height
  )

  write_log_lines(
    c(
      paste("set_name:", set_name),
      paste("set_file:", set_file),
      paste("stats_file:", opt$stats_file),
      paste("modules_rds:", ifelse(is.null(opt$modules_rds), "NULL", opt$modules_rds)),
      paste("module_dendro_distance:", module_dendro_distance),
      paste("p_thresh:", opt$p_thresh),
      paste("top_n:", opt$top_n),
      paste("n_requested_traits:", length(requested_traits)),
      paste("n_found_traits:", length(found_traits)),
      paste("n_missing_traits:", length(missing_traits)),
      paste("n_rows_subset:", nrow(subset_df)),
      paste("n_rows_significant:", nrow(subset_sig_df)),
      paste("full_heatmap:", full_pdf),
      paste("top_heatmap:", top_pdf)
    ),
    file.path(set_out_dir, "run_parameters.txt")
  )

  message("✓ Done set: ", set_name)
  message("  Traits found: ", length(found_traits))
  message("  Output: ", set_out_dir)
}

# ------------------------------------------------------------
# Top-level run log
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("stats_file:", opt$stats_file),
    paste("modules_rds:", ifelse(is.null(opt$modules_rds), "NULL", opt$modules_rds)),
    paste("set_dir:", ifelse(is.null(opt$set_dir), "NULL", opt$set_dir)),
    paste("set_files:", paste(set_files, collapse = ", ")),
    paste("module_dendro_distance:", module_dendro_distance),
    paste("p_thresh:", opt$p_thresh),
    paste("top_n:", opt$top_n),
    paste("variant_name:", variant_name),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL 09B PRESENTATION HEATMAPS COMPLETE ✓")
message("Outputs saved under:\n  ", out_dir)