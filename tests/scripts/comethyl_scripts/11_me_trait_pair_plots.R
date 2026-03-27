#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 11: ME-Trait Pair Plots
#
# PURPOSE
#   - Load one modules object
#   - Load one ME-trait stats file from 09A
#   - Load one sample info file
#   - Load one or more trait-set text files
#   - Generate raw pair plots for module-trait pairs
#
# DESIGN
#   - stats_file decides which module-trait pairs exist and which are significant
#   - sample_info + modules_rds provide the raw values used for plotting
#   - trait sets are supplied externally via text files
#
# REQUIRED INPUTS
#   --project_root : root directory of the project
#   --modules_rds  : path to one Modules.rds file
#   --stats_file   : one stats file from 09A
#   --sample_info  : path to sample info file (.xlsx, .csv, .tsv, .txt)
#
# TRAIT SET INPUTS
#   Provide either:
#     --set_dir
#   and/or:
#     --set_file, --set_file2, ..., --set_file5
#
# OPTIONAL INPUTS
#   --sample_id_col     : column in sample_info containing sample IDs
#   --trait_exclude_file: optional text file of traits to remove before plotting
#   --p_thresh          : significance threshold [default = 0.05]
#   --only_significant  : TRUE/FALSE [default = TRUE]
#   --top_n_pairs       : optional integer; keep top N pairs per set by p-value
#   --cat_plot_style    : dot, violin, or both [default = violin]
#   --continuous_width  : width for continuous plots [default = 8]
#   --continuous_height : height for continuous plots [default = 8]
#   --categorical_width : width for categorical plots [default = 8]
#   --categorical_height: height for categorical plots [default = 8]
#
# OUTPUTS
#   project_root/comethyl_output/11_me_trait_pair_plots/<cpg_label>/<region_label>/<variant>/<set_name>/
#       traits_requested.txt
#       traits_found.txt
#       traits_missing.txt
#       pair_list.tsv
#       <module>_vs_<trait>.pdf
#       run_parameters.txt
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(openxlsx)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(ggplot2)
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
# Helpers
# ------------------------------------------------------------
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

sanitize_filename <- function(x) {
  x <- gsub("[/\\?<>\\:*|\"'` ]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

parse_optional_int <- function(x, arg_name) {
  if (is.null(x) || is.na(x) || !nzchar(as.character(x))) return(NULL)
  out <- suppressWarnings(as.integer(x))
  if (is.na(out) || out < 1) stop(arg_name, " must be a positive integer if provided.")
  out
}

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--modules_rds", type = "character",
              help = "Path to one Modules.rds file"),

  make_option("--stats_file", type = "character",
              help = "Path to one stats file from 09A"),

  make_option("--sample_info", type = "character",
              help = "Path to sample info file (.xlsx, .csv, .tsv, .txt)"),

  make_option("--sample_id_col", type = "character", default = NULL,
              help = "Column in sample_info containing sample IDs [default = rownames]"),

  make_option("--set_dir", type = "character", default = NULL,
              help = "Directory containing trait set .txt files"),

  make_option("--set_file", type = "character", default = NULL,
              help = "Trait set file, one trait per line"),

  make_option("--set_file2", type = "character", default = NULL),
  make_option("--set_file3", type = "character", default = NULL),
  make_option("--set_file4", type = "character", default = NULL),
  make_option("--set_file5", type = "character", default = NULL),

  make_option("--trait_exclude_file", type = "character", default = NULL,
              help = "Optional text file of traits to exclude before plotting"),

  make_option("--p_thresh", type = "double", default = 0.05,
              help = "Significance threshold [default = 0.05]"),

  make_option("--only_significant", type = "character", default = "TRUE",
              help = "Plot only significant pairs: TRUE/FALSE [default = TRUE]"),

  make_option("--top_n_pairs", type = "character", default = NULL,
              help = "Optional integer: keep top N pairs per set by p-value"),

  make_option("--cat_plot_style", type = "character", default = "violin",
              help = "dot, violin, or both [default = violin]"),

  make_option("--continuous_width", type = "double", default = 8,
              help = "Continuous plot width [default = 8]"),

  make_option("--continuous_height", type = "double", default = 8,
              help = "Continuous plot height [default = 8]"),

  make_option("--categorical_width", type = "double", default = 8,
              help = "Categorical plot width [default = 8]"),

  make_option("--categorical_height", type = "double", default = 8,
              help = "Categorical plot height [default = 8]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$modules_rds)) stop("--modules_rds is required")
if (is.null(opt$stats_file)) stop("--stats_file is required")
if (is.null(opt$sample_info)) stop("--sample_info is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$modules_rds)) stop("modules_rds not found: ", opt$modules_rds)
if (!file.exists(opt$stats_file)) stop("stats_file not found: ", opt$stats_file)
if (!file.exists(opt$sample_info)) stop("sample_info not found: ", opt$sample_info)

if (!is.null(opt$trait_exclude_file) && !file.exists(opt$trait_exclude_file)) {
  stop("trait_exclude_file not found: ", opt$trait_exclude_file)
}

only_significant <- parse_bool(opt$only_significant, "--only_significant")
top_n_pairs <- parse_optional_int(opt$top_n_pairs, "--top_n_pairs")

cat_plot_style <- tolower(opt$cat_plot_style)
if (!cat_plot_style %in% c("dot", "violin", "both")) {
  stop("--cat_plot_style must be one of: dot, violin, both")
}
if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")

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
# Derive lineage from modules_rds
# Expected:
#   .../07_module_detection/<cpg_label>/<region_label>/<variant>/Modules.rds
# ------------------------------------------------------------
variant_dir <- dirname(opt$modules_rds)
region_dir <- dirname(variant_dir)
variant_name <- basename(variant_dir)
region_label <- basename(region_dir)
cpg_label <- basename(dirname(region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "11_me_trait_pair_plots")
out_dir <- file.path(step_dir, cpg_label, region_label, variant_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Load modules
# ------------------------------------------------------------
modules <- validate_modules_object(readRDS(opt$modules_rds), "modules_rds")
MEs <- as.data.frame(modules$MEs)

# ------------------------------------------------------------
# Load stats
# ------------------------------------------------------------
stats_df <- readr::read_tsv(opt$stats_file, show_col_types = FALSE)
stats_df <- validate_stats_table(stats_df, "stats_file")

stats_df$module <- as.character(stats_df$module)
stats_df$trait  <- as.character(stats_df$trait)

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
# Align samples
# ------------------------------------------------------------
common_samples <- intersect(rownames(MEs), rownames(sample_info))
if (length(common_samples) == 0) {
  stop("No overlapping samples between module eigengenes and sample_info.")
}

MEs_use <- MEs[common_samples, , drop = FALSE]
sample_info_use <- sample_info[common_samples, , drop = FALSE]

message("Samples in common: ", length(common_samples))
message("Modules available: ", ncol(MEs_use))
message("Traits available: ", ncol(sample_info_use))

# ------------------------------------------------------------
# Run each trait set
# ------------------------------------------------------------
for (set_file in set_files) {
  set_name <- get_set_name(set_file)
  set_out_dir <- file.path(out_dir, set_name)
  dir.create(set_out_dir, recursive = TRUE, showWarnings = FALSE)

  requested_traits <- read_trait_set_file(set_file)
  available_traits <- colnames(sample_info_use)

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
        paste("modules_rds:", opt$modules_rds),
        paste("stats_file:", opt$stats_file),
        paste("sample_info:", opt$sample_info),
        "status: skipped",
        "reason: no requested traits found in sample_info"
      ),
      file.path(set_out_dir, "run_parameters.txt")
    )
    message("Skipping set ", set_name, ": no requested traits found in sample_info.")
    next
  }

  pair_df <- stats_df %>%
    dplyr::filter(.data$trait %in% found_traits) %>%
    dplyr::filter(.data$module %in% colnames(MEs_use))

  if (only_significant) {
    pair_df <- pair_df %>%
      dplyr::filter(!is.na(.data$p), .data$p < opt$p_thresh)
  }

  pair_df <- pair_df %>%
    dplyr::arrange(.data$p)

  if (!is.null(top_n_pairs) && nrow(pair_df) > top_n_pairs) {
    pair_df <- pair_df %>% dplyr::slice_head(n = top_n_pairs)
  }

  readr::write_tsv(pair_df, file.path(set_out_dir, "pair_list.tsv"))

  if (nrow(pair_df) == 0) {
    write_log_lines(
      c(
        paste("set_name:", set_name),
        paste("set_file:", set_file),
        paste("modules_rds:", opt$modules_rds),
        paste("stats_file:", opt$stats_file),
        paste("sample_info:", opt$sample_info),
        paste("only_significant:", only_significant),
        paste("p_thresh:", opt$p_thresh),
        paste("top_n_pairs:", ifelse(is.null(top_n_pairs), "NULL", top_n_pairs)),
        "status: skipped",
        "reason: no module-trait pairs left after filtering"
      ),
      file.path(set_out_dir, "run_parameters.txt")
    )
    message("Skipping set ", set_name, ": no module-trait pairs left after filtering.")
    next
  }

  message("[", set_name, "] Pairs to plot: ", nrow(pair_df))

  for (i in seq_len(nrow(pair_df))) {
    mod <- as.character(pair_df$module[i])
    tr  <- as.character(pair_df$trait[i])

    if (!(mod %in% colnames(MEs_use))) next
    if (!(tr  %in% colnames(sample_info_use))) next

    ME_vec <- MEs_use[[mod]]
    trait_vec <- sample_info_use[[tr]]
    yl <- auto_ylim(ME_vec)

    subtitle_text <- make_me_trait_subtitle(stats_df, module = mod, trait = tr)

    make_cat <- is.factor(trait_vec) || is.character(trait_vec) || is_binary_like(trait_vec)

    if (make_cat) {
      if (!is.factor(trait_vec)) trait_vec <- as.factor(trait_vec)
      if (nlevels(trait_vec) == 2) trait_vec <- relabel_binary_to_yesno(trait_vec)
    }

    file_prefix <- file.path(set_out_dir, paste0(sanitize_filename(mod), "_vs_", sanitize_filename(tr)))
    pdf_file <- paste0(file_prefix, ".pdf")

    if (!make_cat) {
      p <- plotMEtraitScatter(
        ME = ME_vec,
        trait = trait_vec,
        ylim = yl,
        xlab = tr,
        ylab = paste0(mod, " Module Eigengene"),
        save = FALSE,
        axis.text.size = 16,
        verbose = FALSE
      )
      p <- add_plot_subtitle(p, subtitle_text)

      ggplot2::ggsave(
        pdf_file,
        plot = p,
        dpi = 600,
        width = opt$continuous_width,
        height = opt$continuous_height,
        units = "in"
      )

    } else {
      if (cat_plot_style == "dot") {
        p_dot <- plotMEtraitDot(
          ME = ME_vec,
          trait = trait_vec,
          ylim = yl,
          xlab = tr,
          ylab = paste0(mod, " Module Eigengene"),
          save = FALSE,
          axis.text.size = 16,
          verbose = FALSE
        )
        p_dot <- add_plot_subtitle(p_dot, subtitle_text)

        ggplot2::ggsave(
          pdf_file,
          plot = p_dot,
          dpi = 600,
          width = opt$categorical_width,
          height = opt$categorical_height,
          units = "in"
        )

      } else if (cat_plot_style == "violin") {
        p_v <- plotMEtraitViolin(
          ME_vec = ME_vec,
          trait_fac = trait_vec,
          module_name = mod,
          trait_name = tr,
          ylim = yl
        )
        p_v <- add_plot_subtitle(p_v, subtitle_text)

        if (!is.null(p_v)) {
          ggplot2::ggsave(
            pdf_file,
            plot = p_v,
            dpi = 600,
            width = opt$categorical_width,
            height = opt$categorical_height,
            units = "in"
          )
        }

      } else if (cat_plot_style == "both") {
        grDevices::pdf(pdf_file, width = opt$categorical_width, height = opt$categorical_height)

        p_dot <- plotMEtraitDot(
          ME = ME_vec,
          trait = trait_vec,
          ylim = yl,
          xlab = tr,
          ylab = paste0(mod, " Module Eigengene"),
          save = FALSE,
          axis.text.size = 16,
          verbose = FALSE
        )
        p_dot <- add_plot_subtitle(p_dot, subtitle_text)
        print(p_dot)

        p_v <- plotMEtraitViolin(
          ME_vec = ME_vec,
          trait_fac = trait_vec,
          module_name = mod,
          trait_name = tr,
          ylim = yl
        )
        p_v <- add_plot_subtitle(p_v, subtitle_text)
        if (!is.null(p_v)) print(p_v)

        dev.off()
      }
    }
  }

  write_log_lines(
    c(
      paste("set_name:", set_name),
      paste("set_file:", set_file),
      paste("modules_rds:", opt$modules_rds),
      paste("stats_file:", opt$stats_file),
      paste("sample_info:", opt$sample_info),
      paste("sample_id_col:", ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
      paste("trait_exclude_file:", ifelse(is.null(opt$trait_exclude_file), "NULL", opt$trait_exclude_file)),
      paste("only_significant:", only_significant),
      paste("p_thresh:", opt$p_thresh),
      paste("top_n_pairs:", ifelse(is.null(top_n_pairs), "NULL", top_n_pairs)),
      paste("cat_plot_style:", cat_plot_style),
      paste("n_requested_traits:", length(requested_traits)),
      paste("n_found_traits:", length(found_traits)),
      paste("n_missing_traits:", length(missing_traits)),
      paste("n_pairs_plotted:", nrow(pair_df))
    ),
    file.path(set_out_dir, "run_parameters.txt")
  )

  message("✓ Done set: ", set_name)
  message("  Output: ", set_out_dir)
}

# ------------------------------------------------------------
# Top-level run log
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("modules_rds:", opt$modules_rds),
    paste("stats_file:", opt$stats_file),
    paste("sample_info:", opt$sample_info),
    paste("sample_id_col:", ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
    paste("set_dir:", ifelse(is.null(opt$set_dir), "NULL", opt$set_dir)),
    paste("set_files:", paste(set_files, collapse = ", ")),
    paste("trait_exclude_file:", ifelse(is.null(opt$trait_exclude_file), "NULL", opt$trait_exclude_file)),
    paste("only_significant:", only_significant),
    paste("p_thresh:", opt$p_thresh),
    paste("top_n_pairs:", ifelse(is.null(top_n_pairs), "NULL", top_n_pairs)),
    paste("cat_plot_style:", cat_plot_style),
    paste("variant_name:", variant_name),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL ME-TRAIT PAIR PLOTS COMPLETE ✓")
message("Outputs saved under:\n  ", out_dir)