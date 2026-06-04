#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 09B: ME-Trait Presentation Heatmaps
#
# PURPOSE
#   - Load ME-trait correlation stats from Script 09A (up to 3 variants)
#   - Load one or more trait-set text files, either directly or from a directory
#   - Generate focused presentation heatmaps for those trait sets
#   - Optionally order modules by biological association class:
#       1) exposure + outcome associated modules
#       2) exposure-only associated modules
#       3) outcome-only associated modules
#       4) other modules
#   - Optionally order traits/heatmap columns by the same biological classes:
#       1) traits linked to exposure + outcome modules
#       2) traits linked to exposure-only modules
#       3) traits linked to outcome-only modules
#       4) other traits
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --stats_file   : path to v1 ME-trait stats file from Script 09A
#
# TRAIT SET INPUTS
#   Provide either:
#     --set_dir  : directory containing one or more .txt trait-set files
#   and/or:
#     --set_file, --set_file2, ..., --set_file5
#
# OPTIONAL INPUTS
#   --stats_file_v2           : path to v2 ME-trait stats file from 09A
#   --stats_file_v3           : path to v3 ME-trait stats file from 09A
#   --modules_rds             : optional v1 Modules_for_downstream.rds for dendrogram ordering
#   --modules_rds_v2          : optional v2 Modules_for_downstream.rds for dendrogram ordering
#   --modules_rds_v3          : optional v3 Modules_for_downstream.rds for dendrogram ordering
#   --module_dendro_distance  : bicor or pearson [default = bicor]
#   --module_order_mode       : default, dendro, or exposure_outcome [default = dendro]
#   --module_order_trait_file : TSV/CSV file defining exposure/outcome traits
#   --module_order_p_thresh   : p-value cutoff for classifying modules [default = same as --p_thresh]
#   --module_order_score      : p or abs_cor, within-class ordering [default = p]
#   --trait_order_mode        : input or exposure_outcome [default = input]
#   --trait_order_p_thresh    : p-value cutoff for trait classification [default = module_order_p_thresh]
#   --trait_order_score       : p or abs_cor, within-class trait ordering [default = p]
#   --write_reordered_trait_files : TRUE/FALSE; save ordered trait text files [default = TRUE]
#   --drop_grey_modules       : TRUE/FALSE; remove grey from heatmaps [default = TRUE]
#   --p_thresh                : significance threshold for heatmap labels/top filtering [default = 0.05]
#   --top_n                   : number of top associations for TOP heatmap [default = 250]
#   --min_sig_traits          : minimum number of significant trait associations (p < p_thresh)
#                               a module must have to appear in the FOCUSED heatmap.
#                               Use this to reduce a large module set (e.g. 33 modules) down
#                               to only the most relevant modules for your trait set, making
#                               the heatmap easier to read in presentations.
#                               Set to 1 to include any module with at least one hit.
#                               Set to 2 or 3 to show only multi-trait associated modules.
#                               [default = 1]
#   --full_width              : full heatmap width [default = 12]
#   --full_height             : full heatmap height [default = 12]
#   --top_width               : top heatmap width [default = 9]
#   --top_height              : top heatmap height [default = 6]
#   --focused_width           : FOCUSED heatmap width in inches [default = same as --full_width]
#   --focused_height          : FOCUSED heatmap height in inches [default = 9]
#                               Set this independently from --full_height so that the filtered
#                               module subset has a fixed, presentation-ready size regardless
#                               of how many modules are in the full set.
#
# MODULE ORDER TRAIT FILE FORMAT
#   A tab-delimited file with two required columns:
#
#     category    trait
#     exposure    DDE
#     exposure    op_DDT
#     outcome     breast_cancer
#
#   category must be one of:
#     exposure
#     outcome
#
# OUTPUTS
#   project_root/comethyl_output/09b_me_trait_presentation/<cpg_label>/<region_label>/<variant>/<set_name>/
#       traits_requested.txt
#       traits_found.txt
#       traits_missing.txt
#       traits_ordered.txt
#       <set_name>_ORDERED.txt
#       trait_order_exposure_outcome_summary.tsv
#       subset_stats.tsv
#       subset_stats_significant.tsv
#       module_order_exposure_outcome_summary.tsv
#       ordered_modules.txt
#       ME_Trait_Heatmap_FULL.pdf
#       ME_Trait_Heatmap_TOP.pdf
#       ME_Trait_Heatmap_FOCUSED.pdf
#           A filtered version of the FULL heatmap showing only modules that meet
#           the --min_sig_traits threshold. Modules with fewer significant trait
#           associations than the threshold are dropped. Height scales automatically
#           with the number of retained modules. Useful for presentation slides
#           where text legibility matters. Replaced by _SKIPPED.txt if no modules
#           meet the threshold.
#       focused_modules.txt
#       run_parameters.txt
#
# NOTES
#   - This script changes module ordering only. It does not change the
#     underlying correlations or p-values.
#   - For downstream consistency, use Modules_for_downstream.rds from Script 08.
# ================================================================

message("Starting Script 09b")

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
# Load helper.R if present next to this script
# ------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_file_arg) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (file.exists(helper_file)) {
  source(helper_file)
} else {
  warning("helper.R not found next to this script: ", helper_file,
          ". Using internal fallback helper functions where possible.")
}

# ------------------------------------------------------------
# Fallback helpers, in case helper.R does not define them
# ------------------------------------------------------------
if (!exists("write_vector_file")) {
  write_vector_file <- function(x, file) {
    writeLines(as.character(x), con = file, useBytes = TRUE)
  }
}

if (!exists("write_log_lines")) {
  write_log_lines <- function(lines, file) {
    writeLines(as.character(lines), con = file, useBytes = TRUE)
  }
}

if (!exists("validate_modules_object")) {
  validate_modules_object <- function(x, label = "modules object") {
    if (!is.list(x)) stop(label, " must be a list-like object.")
    if (is.null(x$MEs)) stop(label, " is missing $MEs.")
    if (!(is.matrix(x$MEs) || is.data.frame(x$MEs))) {
      stop(label, "$MEs must be matrix-like.")
    }
    x$MEs <- as.matrix(x$MEs)
    if (is.null(colnames(x$MEs))) stop(label, "$MEs must have column names.")
    x
  }
}

# ------------------------------------------------------------
# General helpers
# ------------------------------------------------------------
str_to_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(as.character(x))) return(default)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

module_name_from_me_col <- function(x) {
  sub("^ME", "", as.character(x))
}

is_grey_module <- function(x) {
  tolower(module_name_from_me_col(x)) == "grey"
}

read_trait_set_file <- function(file) {
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]
  unique(x)
}

get_set_name <- function(file) {
  tools::file_path_sans_ext(basename(file))
}

get_cor_col <- function(stats_df) {
  if ("bicor" %in% colnames(stats_df)) return("bicor")
  if ("cor" %in% colnames(stats_df)) return("cor")
  if ("correlation" %in% colnames(stats_df)) return("correlation")
  stop("Stats table must contain one of: bicor, cor, correlation")
}

standardize_p_col <- function(stats_df) {
  if ("p" %in% colnames(stats_df)) return(stats_df)
  candidates <- c("pvalue", "p.value", "p_value", "P.Value", "P.value")
  hit <- candidates[candidates %in% colnames(stats_df)][1]
  if (is.na(hit)) return(stats_df)
  colnames(stats_df)[match(hit, colnames(stats_df))] <- "p"
  stats_df
}

validate_stats_table <- function(df, label = "stats_file") {
  df <- standardize_p_col(df)
  required <- c("module", "trait", "p")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
  get_cor_col(df)
  df
}

collect_set_files <- function(set_dir   = NULL,
                              set_file  = NULL,
                              set_file2 = NULL,
                              set_file3 = NULL,
                              set_file4 = NULL,
                              set_file5 = NULL,
                              recursive = FALSE) {

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
    if (!dir.exists(set_dir)) stop("set_dir does not exist: ", set_dir)
    dir_files <- list.files(
      set_dir,
      pattern = "\\.txt$",
      full.names = TRUE,
      recursive = recursive,
      ignore.case = TRUE
    )
    if (length(dir_files) == 0) stop("No .txt files found in set_dir: ", set_dir)
  }

  all_files <- c(direct_files, dir_files)
  if (length(all_files) == 0) {
    stop("No trait set files provided. Use --set_dir or at least one --set_file.")
  }

  unique(normalizePath(all_files, mustWork = TRUE))
}

read_module_order_trait_file <- function(file) {
  if (is.null(file) || !nzchar(file)) {
    stop("--module_order_trait_file is required when --module_order_mode exposure_outcome")
  }
  if (!file.exists(file)) stop("module_order_trait_file not found: ", file)

  ext <- tolower(tools::file_ext(file))
  x <- switch(ext,
    csv = read.csv(file, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
    txt = read.delim(file, stringsAsFactors = FALSE, check.names = FALSE),
    read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
  )

  required <- c("category", "trait")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop("module_order_trait_file missing required columns: ",
         paste(missing, collapse = ", "))
  }

  x$category <- tolower(trimws(as.character(x$category)))
  x$trait <- trimws(as.character(x$trait))

  x <- x[x$category %in% c("exposure", "outcome") &
         !is.na(x$trait) & x$trait != "", , drop = FALSE]

  if (nrow(x) == 0) {
    stop("module_order_trait_file has no valid exposure/outcome rows.")
  }

  unique(x[, c("category", "trait")])
}

# ------------------------------------------------------------
# Exposure/outcome module ordering
# ------------------------------------------------------------
build_exposure_outcome_module_order <- function(stats_df,
                                                trait_def_file,
                                                p_thresh = 0.05,
                                                score_by = "p",
                                                current_levels = NULL) {

  trait_def <- read_module_order_trait_file(trait_def_file)

  exposure_traits <- unique(trait_def$trait[trait_def$category == "exposure"])
  outcome_traits  <- unique(trait_def$trait[trait_def$category == "outcome"])

  if (length(exposure_traits) == 0) stop("No exposure traits defined in module_order_trait_file.")
  if (length(outcome_traits) == 0)  stop("No outcome traits defined in module_order_trait_file.")

  if (is.null(current_levels)) {
    current_levels <- unique(as.character(stats_df$module))
  }

  cor_col <- get_cor_col(stats_df)

  df <- stats_df %>%
    dplyr::mutate(
      module = as.character(module),
      trait = as.character(trait),
      p = suppressWarnings(as.numeric(p)),
      cor_value = suppressWarnings(as.numeric(.data[[cor_col]])),
      trait_class = dplyr::case_when(
        trait %in% exposure_traits ~ "exposure",
        trait %in% outcome_traits ~ "outcome",
        TRUE ~ "other"
      )
    ) %>%
    dplyr::filter(
      module %in% current_levels,
      trait_class %in% c("exposure", "outcome"),
      !is.na(p),
      p <= p_thresh
    )

  exposure_modules <- unique(df$module[df$trait_class == "exposure"])
  outcome_modules  <- unique(df$module[df$trait_class == "outcome"])

  module_summary <- data.frame(
    module = current_levels,
    module_class = dplyr::case_when(
      current_levels %in% exposure_modules & current_levels %in% outcome_modules ~ "exposure_and_outcome",
      current_levels %in% exposure_modules ~ "exposure_only",
      current_levels %in% outcome_modules ~ "outcome_only",
      TRUE ~ "other"
    ),
    stringsAsFactors = FALSE
  )

  score_df <- df %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      min_p = min(p, na.rm = TRUE),
      max_abs_cor = max(abs(cor_value), na.rm = TRUE),
      n_sig_traits = dplyr::n_distinct(trait),
      n_sig_exposure_traits = dplyr::n_distinct(trait[trait_class == "exposure"]),
      n_sig_outcome_traits = dplyr::n_distinct(trait[trait_class == "outcome"]),
      sig_exposure_traits = paste(sort(unique(trait[trait_class == "exposure"])), collapse = ";"),
      sig_outcome_traits = paste(sort(unique(trait[trait_class == "outcome"])), collapse = ";"),
      .groups = "drop"
    )

  module_summary <- module_summary %>%
    dplyr::left_join(score_df, by = "module") %>%
    dplyr::mutate(
      min_p = ifelse(is.na(min_p), Inf, min_p),
      max_abs_cor = ifelse(is.na(max_abs_cor), 0, max_abs_cor),
      n_sig_traits = ifelse(is.na(n_sig_traits), 0, n_sig_traits),
      n_sig_exposure_traits = ifelse(is.na(n_sig_exposure_traits), 0, n_sig_exposure_traits),
      n_sig_outcome_traits = ifelse(is.na(n_sig_outcome_traits), 0, n_sig_outcome_traits),
      sig_exposure_traits = ifelse(is.na(sig_exposure_traits), "", sig_exposure_traits),
      sig_outcome_traits = ifelse(is.na(sig_outcome_traits), "", sig_outcome_traits),
      class_rank = dplyr::case_when(
        module_class == "exposure_and_outcome" ~ 1L,
        module_class == "exposure_only" ~ 2L,
        module_class == "outcome_only" ~ 3L,
        TRUE ~ 4L
      )
    )

  if (score_by == "abs_cor") {
    module_summary <- module_summary %>%
      dplyr::arrange(class_rank, dplyr::desc(max_abs_cor), min_p, module)
  } else {
    module_summary <- module_summary %>%
      dplyr::arrange(class_rank, min_p, dplyr::desc(max_abs_cor), module)
  }

  ordered_modules <- module_summary$module
  module_order <- match(ordered_modules, current_levels)
  module_order <- module_order[!is.na(module_order)]

  list(
    module_order = module_order,
    ordered_modules = ordered_modules,
    module_summary = module_summary,
    exposure_traits = exposure_traits,
    outcome_traits = outcome_traits
  )
}


# ------------------------------------------------------------
# Exposure/outcome trait ordering
# ------------------------------------------------------------
build_exposure_outcome_trait_order <- function(stats_df,
                                               found_traits,
                                               module_order_summary,
                                               p_thresh = 0.05,
                                               score_by = "p") {

  found_traits <- unique(as.character(found_traits))

  if (length(found_traits) == 0) {
    return(list(
      trait_order = integer(0),
      ordered_traits = character(0),
      trait_summary = data.frame()
    ))
  }

  if (is.null(module_order_summary) || nrow(module_order_summary) == 0) {
    trait_summary <- data.frame(
      trait = found_traits,
      trait_class = "input_order",
      min_p = NA_real_,
      max_abs_cor = NA_real_,
      n_sig_modules = NA_integer_,
      sig_exposure_and_outcome_modules = "",
      sig_exposure_only_modules = "",
      sig_outcome_only_modules = "",
      sig_other_modules = "",
      original_rank = seq_along(found_traits),
      class_rank = seq_along(found_traits),
      stringsAsFactors = FALSE
    )
    return(list(
      trait_order = seq_along(found_traits),
      ordered_traits = found_traits,
      trait_summary = trait_summary
    ))
  }

  cor_col <- get_cor_col(stats_df)

  module_class_df <- module_order_summary %>%
    dplyr::select(module, module_class)

  sig_df <- stats_df %>%
    dplyr::mutate(
      module = as.character(module),
      trait = as.character(trait),
      p = suppressWarnings(as.numeric(p)),
      cor_value = suppressWarnings(as.numeric(.data[[cor_col]]))
    ) %>%
    dplyr::filter(
      trait %in% found_traits,
      !is.na(p),
      p <= p_thresh
    ) %>%
    dplyr::left_join(module_class_df, by = "module") %>%
    dplyr::mutate(
      module_class = ifelse(is.na(module_class), "other", module_class)
    )

  trait_summary <- data.frame(
    trait = found_traits,
    original_rank = seq_along(found_traits),
    stringsAsFactors = FALSE
  )

  score_df <- sig_df %>%
    dplyr::group_by(trait) %>%
    dplyr::summarise(
      has_exposure_and_outcome = any(module_class == "exposure_and_outcome"),
      has_exposure_only = any(module_class == "exposure_only"),
      has_outcome_only = any(module_class == "outcome_only"),
      has_other = any(module_class == "other"),
      min_p = min(p, na.rm = TRUE),
      max_abs_cor = max(abs(cor_value), na.rm = TRUE),
      n_sig_modules = dplyr::n_distinct(module),
      sig_exposure_and_outcome_modules = paste(sort(unique(module[module_class == "exposure_and_outcome"])), collapse = ";"),
      sig_exposure_only_modules = paste(sort(unique(module[module_class == "exposure_only"])), collapse = ";"),
      sig_outcome_only_modules = paste(sort(unique(module[module_class == "outcome_only"])), collapse = ";"),
      sig_other_modules = paste(sort(unique(module[module_class == "other"])), collapse = ";"),
      .groups = "drop"
    )

  trait_summary <- trait_summary %>%
    dplyr::left_join(score_df, by = "trait") %>%
    dplyr::mutate(
      has_exposure_and_outcome = ifelse(is.na(has_exposure_and_outcome), FALSE, has_exposure_and_outcome),
      has_exposure_only = ifelse(is.na(has_exposure_only), FALSE, has_exposure_only),
      has_outcome_only = ifelse(is.na(has_outcome_only), FALSE, has_outcome_only),
      has_other = ifelse(is.na(has_other), FALSE, has_other),
      trait_class = dplyr::case_when(
        has_exposure_and_outcome ~ "exposure_and_outcome",
        has_exposure_only ~ "exposure_only",
        has_outcome_only ~ "outcome_only",
        TRUE ~ "other"
      ),
      min_p = ifelse(is.na(min_p), Inf, min_p),
      max_abs_cor = ifelse(is.na(max_abs_cor), 0, max_abs_cor),
      n_sig_modules = ifelse(is.na(n_sig_modules), 0, n_sig_modules),
      sig_exposure_and_outcome_modules = ifelse(is.na(sig_exposure_and_outcome_modules), "", sig_exposure_and_outcome_modules),
      sig_exposure_only_modules = ifelse(is.na(sig_exposure_only_modules), "", sig_exposure_only_modules),
      sig_outcome_only_modules = ifelse(is.na(sig_outcome_only_modules), "", sig_outcome_only_modules),
      sig_other_modules = ifelse(is.na(sig_other_modules), "", sig_other_modules),
      class_rank = dplyr::case_when(
        trait_class == "exposure_and_outcome" ~ 1L,
        trait_class == "exposure_only" ~ 2L,
        trait_class == "outcome_only" ~ 3L,
        TRUE ~ 4L
      )
    )

  if (score_by == "abs_cor") {
    trait_summary <- trait_summary %>%
      dplyr::arrange(class_rank, dplyr::desc(max_abs_cor), min_p, original_rank)
  } else {
    trait_summary <- trait_summary %>%
      dplyr::arrange(class_rank, min_p, dplyr::desc(max_abs_cor), original_rank)
  }

  ordered_traits <- trait_summary$trait
  trait_order <- match(ordered_traits, found_traits)
  trait_order <- trait_order[!is.na(trait_order)]

  list(
    trait_order = trait_order,
    ordered_traits = ordered_traits,
    trait_summary = trait_summary
  )
}

build_dendro_module_order <- function(stats_df,
                                      modules_rds,
                                      module_dendro_distance = "bicor") {
  current_levels <- levels(stats_df$module)
  module_order <- seq_along(current_levels)

  if (is.null(modules_rds) || !nzchar(modules_rds)) return(module_order)
  if (!file.exists(modules_rds)) stop("modules_rds not found: ", modules_rds)

  modules <- validate_modules_object(readRDS(modules_rds), "modules_rds")
  MEs <- as.matrix(modules$MEs)

  # Remove grey and remove orphan/non-overlapping MEs for ordering only.
  me_modules <- module_name_from_me_col(colnames(MEs))
  keep <- !is_grey_module(colnames(MEs)) & me_modules %in% current_levels
  MEs <- MEs[, keep, drop = FALSE]

  if (ncol(MEs) < 2) return(module_order)

  moduleDendro <- getDendro(MEs, distance = module_dendro_distance)
  dendro_modules <- module_name_from_me_col(colnames(MEs)[moduleDendro$order])

  module_order <- match(intersect(dendro_modules, current_levels), current_levels)
  module_order <- module_order[!is.na(module_order)]

  remaining <- setdiff(seq_along(current_levels), module_order)
  module_order <- c(module_order, remaining)

  if (length(module_order) == 0) module_order <- seq_along(current_levels)
  module_order
}

# ------------------------------------------------------------
# Core per-variant runner
# ------------------------------------------------------------
run_presentation_for_variant <- function(stats_file,
                                          modules_rds,
                                          set_files,
                                          pipeline_root,
                                          module_dendro_distance,
                                          module_order_mode,
                                          module_order_trait_file,
                                          module_order_p_thresh,
                                          module_order_score,
                                          trait_order_mode,
                                          trait_order_p_thresh,
                                          trait_order_score,
                                          write_reordered_trait_files,
                                          drop_grey_modules,
                                          p_thresh,
                                          top_n,
                                          min_sig_traits,
                                          full_width,
                                          full_height,
                                          top_width,
                                          top_height,
                                          focused_width,
                                          focused_height) {

  # Derive output path from stats_file lineage
  # Expected: .../09a_me_trait_analysis/<cpg>/<region>/<variant>/ME_Trait_*.tsv
  variant_dir  <- dirname(stats_file)
  region_dir   <- dirname(variant_dir)
  variant_name <- basename(variant_dir)
  region_label <- basename(region_dir)
  cpg_label    <- basename(dirname(region_dir))

  step_dir <- file.path(pipeline_root, "09b_me_trait_presentation")
  out_dir  <- file.path(step_dir, cpg_label, region_label, variant_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("\n==============================")
  message("Running 09b for variant: ", variant_name)
  message("Output directory: ", out_dir)
  message("==============================\n")

  # Load and validate stats
  stats_df <- readr::read_tsv(stats_file, show_col_types = FALSE)
  stats_df <- validate_stats_table(stats_df, paste0("stats_file (", variant_name, ")"))

  stats_df <- stats_df %>%
    dplyr::mutate(
      module = as.character(module),
      trait = as.character(trait),
      p = suppressWarnings(as.numeric(p))
    )

  if (isTRUE(drop_grey_modules)) {
    n_before <- dplyr::n_distinct(stats_df$module)
    stats_df <- stats_df %>% dplyr::filter(!is_grey_module(module))
    n_after <- dplyr::n_distinct(stats_df$module)
    message("[", variant_name, "] Dropped grey module from stats if present: ",
            n_before, " -> ", n_after, " module(s)")
  }

  stats_df$module <- factor(stats_df$module, levels = unique(stats_df$module))
  stats_df$trait  <- factor(stats_df$trait,  levels = unique(stats_df$trait))

  current_levels <- levels(stats_df$module)

  # ----------------------------------------------------------
  # Module ordering
  # ----------------------------------------------------------
  module_order <- seq_along(current_levels)
  module_order_summary <- NULL
  ordered_modules <- current_levels

  if (module_order_mode == "exposure_outcome") {
    order_res <- build_exposure_outcome_module_order(
      stats_df = stats_df,
      trait_def_file = module_order_trait_file,
      p_thresh = module_order_p_thresh,
      score_by = module_order_score,
      current_levels = current_levels
    )

    module_order <- order_res$module_order
    module_order_summary <- order_res$module_summary
    ordered_modules <- order_res$ordered_modules

    if (length(module_order) == 0) {
      warning("[", variant_name, "] exposure_outcome ordering returned zero modules; using default order.")
      module_order <- seq_along(current_levels)
      ordered_modules <- current_levels
    }
  } else if (module_order_mode == "dendro") {
    module_order <- build_dendro_module_order(
      stats_df = stats_df,
      modules_rds = modules_rds,
      module_dendro_distance = module_dendro_distance
    )
    ordered_modules <- current_levels[module_order]
  } else {
    module_order <- seq_along(current_levels)
    ordered_modules <- current_levels
  }

  write_vector_file(ordered_modules, file.path(out_dir, "ordered_modules.txt"))
  if (!is.null(module_order_summary)) {
    readr::write_tsv(module_order_summary, file.path(out_dir, "module_order_exposure_outcome_summary.tsv"))
  }

  # ----------------------------------------------------------
  # Run each trait set
  # ----------------------------------------------------------
  for (set_file in set_files) {
    set_name    <- get_set_name(set_file)
    set_out_dir <- file.path(out_dir, set_name)
    dir.create(set_out_dir, recursive = TRUE, showWarnings = FALSE)

    requested_traits <- read_trait_set_file(set_file)
    available_traits <- unique(as.character(stats_df$trait))
    found_traits     <- intersect(requested_traits, available_traits)
    missing_traits   <- setdiff(requested_traits, available_traits)

    write_vector_file(requested_traits, file.path(set_out_dir, "traits_requested.txt"))
    write_vector_file(found_traits,     file.path(set_out_dir, "traits_found.txt"))
    write_vector_file(missing_traits,   file.path(set_out_dir, "traits_missing.txt"))
    write_vector_file(ordered_modules,  file.path(set_out_dir, "ordered_modules.txt"))

    if (!is.null(module_order_summary)) {
      readr::write_tsv(module_order_summary,
                       file.path(set_out_dir, "module_order_exposure_outcome_summary.tsv"))
    }

    if (length(found_traits) == 0) {
      write_log_lines(
        c(
          paste("set_name:", set_name),
          paste("variant_name:", variant_name),
          paste("stats_file:", stats_file),
          paste("modules_rds:", ifelse(is.null(modules_rds), "NULL", modules_rds)),
          paste("p_thresh:", p_thresh),
          paste("module_order_mode:", module_order_mode),
          "status: skipped",
          "reason: no requested traits found in stats table"
        ),
        file.path(set_out_dir, "run_parameters.txt")
      )
      message("Skipping set '", set_name, "' for variant '", variant_name,
              "': no requested traits found.")
      next
    }

    # ------------------------------------------------------
    # Trait ordering within this heatmap
    # ------------------------------------------------------
    ordered_traits <- found_traits
    trait_order_summary <- NULL

    if (trait_order_mode == "exposure_outcome") {
      if (is.null(module_order_summary)) {
        warning("[", variant_name, " / ", set_name,
                "] trait_order_mode=exposure_outcome requires module_order_mode=exposure_outcome. Using input trait order.")
      } else {
        # Order traits by the row order of the input TSV file, then append
        # any traits that are in the set but absent from the TSV at the end.
        if (!is.null(module_order_trait_file) && file.exists(module_order_trait_file)) {
          tsv_traits <- read_module_order_trait_file(module_order_trait_file)$trait
          tsv_traits <- unique(as.character(tsv_traits))
          # Keep only traits that are actually in this set, preserving TSV row order
          tsv_ordered <- tsv_traits[tsv_traits %in% found_traits]
          # Append any found traits not listed in the TSV (preserve original set order)
          remainder <- found_traits[!found_traits %in% tsv_ordered]
          ordered_traits <- c(tsv_ordered, remainder)
          message("[", variant_name, " / ", set_name, "] trait_order_mode=exposure_outcome: ",
                  "ordered ", length(tsv_ordered), " trait(s) by TSV row order, ",
                  length(remainder), " appended from set not in TSV.")
        } else {
          trait_order_res <- build_exposure_outcome_trait_order(
            stats_df = stats_df,
            found_traits = found_traits,
            module_order_summary = module_order_summary,
            p_thresh = trait_order_p_thresh,
            score_by = trait_order_score
          )
          ordered_traits <- trait_order_res$ordered_traits
          trait_order_summary <- trait_order_res$trait_summary
        }
      }
    }

    write_vector_file(ordered_traits, file.path(set_out_dir, "traits_ordered.txt"))
    if (!is.null(trait_order_summary)) {
      readr::write_tsv(trait_order_summary,
                       file.path(set_out_dir, "trait_order_exposure_outcome_summary.tsv"))
    }

    if (isTRUE(write_reordered_trait_files)) {
      write_vector_file(
        ordered_traits,
        file.path(set_out_dir, paste0(set_name, "_ORDERED.txt"))
      )
    }

    subset_df <- stats_df %>% dplyr::filter(as.character(trait) %in% found_traits)

    subset_df$module <- factor(
      as.character(subset_df$module),
      levels = current_levels
    )

    subset_df$trait <- factor(
      as.character(subset_df$trait),
      levels = ordered_traits
    )

    subset_sig_df <- subset_df %>% dplyr::filter(!is.na(p), p < p_thresh)

    # ----------------------------------------------------------
    # Focused module filtering (min_sig_traits)
    # Keep only modules with >= min_sig_traits significant trait
    # associations within this set. Used for the FOCUSED heatmap.
    # ----------------------------------------------------------
    focused_modules <- subset_sig_df %>%
      dplyr::group_by(module) %>%
      dplyr::summarise(n_sig = dplyr::n_distinct(trait), .groups = "drop") %>%
      dplyr::filter(n_sig >= min_sig_traits) %>%
      dplyr::pull(module) %>%
      as.character()

    focused_df <- subset_df %>%
      dplyr::filter(as.character(module) %in% focused_modules)

    focused_df$module <- factor(
      as.character(focused_df$module),
      levels = current_levels[current_levels %in% focused_modules]
    )
    focused_df$trait <- factor(
      as.character(focused_df$trait),
      levels = ordered_traits
    )
    focused_module_order <- seq_along(levels(focused_df$module))

    message("[", variant_name, " / ", set_name, "] Focused heatmap: ",
            length(focused_modules), " of ", length(current_levels),
            " modules have >= ", min_sig_traits, " significant trait association(s)")

    write_vector_file(focused_modules, file.path(set_out_dir, "focused_modules.txt"))

    readr::write_tsv(subset_df,     file.path(set_out_dir, "subset_stats.tsv"))
    readr::write_tsv(subset_sig_df, file.path(set_out_dir, "subset_stats_significant.tsv"))

    full_pdf    <- file.path(set_out_dir, "ME_Trait_Heatmap_FULL.pdf")
    top_pdf     <- file.path(set_out_dir, "ME_Trait_Heatmap_TOP.pdf")
    focused_pdf <- file.path(set_out_dir, "ME_Trait_Heatmap_FOCUSED.pdf")

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = module_order,
        p               = p_thresh,
        topOnly         = FALSE,
        file            = full_pdf,
        width           = full_width,
        height          = full_height,
        colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
      )
    }, error = function(e) {
      message("Failed FULL heatmap for set '", set_name, "' variant '", variant_name, "': ", conditionMessage(e))
      writeLines(paste("Failed FULL heatmap:", conditionMessage(e)),
                 con = file.path(set_out_dir, "ME_Trait_Heatmap_FULL_FAILED.txt"))
    })

    tryCatch({
      plotMEtraitCor(
        subset_df,
        moduleOrder     = module_order,
        p               = p_thresh,
        topOnly         = TRUE,
        nTop            = top_n,
        label.type      = "p",
        label.size      = 4,
        label.nudge_y   = 0,
        legend.position = c(1.11, 0.795),
        colColorMargins = c(-1, 4.75, 0.5, 10.1),
        file            = top_pdf,
        width           = top_width,
        height          = top_height
      )
    }, error = function(e) {
      message("Failed TOP heatmap for set '", set_name, "' variant '", variant_name, "': ", conditionMessage(e))
      writeLines(paste("Failed TOP heatmap:", conditionMessage(e)),
                 con = file.path(set_out_dir, "ME_Trait_Heatmap_TOP_FAILED.txt"))
    })

    # FOCUSED heatmap: only modules with >= min_sig_traits significant associations
    if (length(focused_modules) == 0) {
      message("[", variant_name, " / ", set_name, "] Skipping FOCUSED heatmap: no modules meet min_sig_traits >= ", min_sig_traits)
      writeLines(paste("Skipped: no modules with >= ", min_sig_traits, " significant trait associations"),
                 con = file.path(set_out_dir, "ME_Trait_Heatmap_FOCUSED_SKIPPED.txt"))
    } else {
      tryCatch({
        plotMEtraitCor(
          focused_df,
          moduleOrder     = focused_module_order,
          p               = p_thresh,
          topOnly         = FALSE,
          file            = focused_pdf,
          width           = focused_width,
          height          = focused_height,
          colColorMargins = c(-2.5, 4.21, 3.0, 12.07)
        )
      }, error = function(e) {
        message("Failed FOCUSED heatmap for set '", set_name, "' variant '", variant_name, "': ", conditionMessage(e))
        writeLines(paste("Failed FOCUSED heatmap:", conditionMessage(e)),
                   con = file.path(set_out_dir, "ME_Trait_Heatmap_FOCUSED_FAILED.txt"))
      })
    }

    write_log_lines(
      c(
        paste("set_name:", set_name),
        paste("variant_name:", variant_name),
        paste("set_file:", set_file),
        paste("stats_file:", stats_file),
        paste("modules_rds:", ifelse(is.null(modules_rds), "NULL", modules_rds)),
        paste("module_dendro_distance:", module_dendro_distance),
        paste("module_order_mode:", module_order_mode),
        paste("module_order_trait_file:", ifelse(is.null(module_order_trait_file), "NULL", module_order_trait_file)),
        paste("module_order_p_thresh:", module_order_p_thresh),
        paste("module_order_score:", module_order_score),
        paste("trait_order_mode:", trait_order_mode),
        paste("trait_order_p_thresh:", trait_order_p_thresh),
        paste("trait_order_score:", trait_order_score),
        paste("write_reordered_trait_files:", write_reordered_trait_files),
        paste("drop_grey_modules:", drop_grey_modules),
        paste("p_thresh:", p_thresh),
        paste("top_n:", top_n),
        paste("min_sig_traits:", min_sig_traits),
        paste("n_focused_modules:", length(focused_modules)),
        paste("n_requested_traits:", length(requested_traits)),
        paste("n_found_traits:", length(found_traits)),
        paste("n_missing_traits:", length(missing_traits)),
        paste("n_modules_in_heatmap:", length(current_levels)),
        paste("n_rows_subset:", nrow(subset_df)),
        paste("n_rows_significant:", nrow(subset_sig_df)),
        paste("full_heatmap:", full_pdf),
        paste("top_heatmap:", top_pdf),
        paste("date:", as.character(Sys.time()))
      ),
      file.path(set_out_dir, "run_parameters.txt")
    )

    message("Done set '", set_name, "' for variant '", variant_name, "'")
    message("  Traits found: ", length(found_traits), " / ", length(requested_traits))
    message("  Output: ", set_out_dir)
  }

  # Variant-level run log
  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("stats_file:", stats_file),
      paste("modules_rds:", ifelse(is.null(modules_rds), "NULL", modules_rds)),
      paste("module_dendro_distance:", module_dendro_distance),
      paste("module_order_mode:", module_order_mode),
      paste("module_order_trait_file:", ifelse(is.null(module_order_trait_file), "NULL", module_order_trait_file)),
      paste("module_order_p_thresh:", module_order_p_thresh),
      paste("module_order_score:", module_order_score),
      paste("trait_order_mode:", trait_order_mode),
      paste("trait_order_p_thresh:", trait_order_p_thresh),
      paste("trait_order_score:", trait_order_score),
      paste("write_reordered_trait_files:", write_reordered_trait_files),
      paste("drop_grey_modules:", drop_grey_modules),
      paste("p_thresh:", p_thresh),
      paste("top_n:", top_n),
      paste("set_files:", paste(set_files, collapse = ", ")),
      paste("cpg_label:", cpg_label),
      paste("region_label:", region_label),
      paste("ordered_modules_file:", file.path(out_dir, "ordered_modules.txt")),
      paste("date:", as.character(Sys.time()))
    ),
    file.path(out_dir, "run_parameters.txt")
  )

  message("Finished variant: ", variant_name)
}

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--stats_file", type = "character",
              help = "Path to v1 ME-trait stats file from 09A"),

  make_option("--stats_file_v2", type = "character", default = NULL,
              help = "Optional path to v2 ME-trait stats file from 09A"),

  make_option("--stats_file_v3", type = "character", default = NULL,
              help = "Optional path to v3 ME-trait stats file from 09A"),

  make_option("--set_dir", type = "character", default = NULL,
              help = "Directory containing one or more .txt trait set files"),

  make_option("--set_file", type = "character", default = NULL,
              help = "Trait set text file, one trait per line"),

  make_option("--set_file2", type = "character", default = NULL),
  make_option("--set_file3", type = "character", default = NULL),
  make_option("--set_file4", type = "character", default = NULL),
  make_option("--set_file5", type = "character", default = NULL),

  make_option("--set_dir_recursive", type = "logical", default = FALSE,
              help = "If TRUE, search --set_dir recursively for .txt files [default = FALSE]"),

  make_option("--modules_rds", type = "character", default = NULL,
              help = "Optional v1 Modules_for_downstream.rds for dendrogram ordering"),

  make_option("--modules_rds_v2", type = "character", default = NULL,
              help = "Optional v2 Modules_for_downstream.rds for dendrogram ordering"),

  make_option("--modules_rds_v3", type = "character", default = NULL,
              help = "Optional v3 Modules_for_downstream.rds for dendrogram ordering"),

  make_option("--module_dendro_distance", type = "character", default = "bicor",
              help = "bicor or pearson [default = bicor]"),

  make_option("--module_order_mode", type = "character", default = "dendro",
              help = "Module order mode: default, dendro, or exposure_outcome [default = dendro]"),

  make_option("--module_order_trait_file", type = "character", default = NULL,
              help = "TSV/CSV/TXT with columns category and trait. Required when --module_order_mode exposure_outcome."),

  make_option("--module_order_p_thresh", type = "double", default = NA,
              help = "P-value cutoff for exposure/outcome module classification. If not set, uses --p_thresh."),

  make_option("--module_order_score", type = "character", default = "p",
              help = "Within-class ordering score: p or abs_cor [default = p]"),

  make_option("--trait_order_mode", type = "character", default = "input",
              help = "Trait order mode inside each heatmap: input or exposure_outcome [default = input]"),

  make_option("--trait_order_p_thresh", type = "double", default = NA,
              help = "P-value cutoff for exposure/outcome trait ordering. If not set, uses --module_order_p_thresh."),

  make_option("--trait_order_score", type = "character", default = "p",
              help = "Within-class trait ordering score: p or abs_cor [default = p]"),

  make_option("--write_reordered_trait_files", type = "logical", default = TRUE,
              help = "Write a reordered trait text file inside each output set folder [default = TRUE]"),

  make_option("--drop_grey_modules", type = "logical", default = TRUE,
              help = "Remove grey module from presentation heatmaps [default = TRUE]"),

  make_option("--p_thresh", type = "double", default = 0.05,
              help = "Significance threshold [default = 0.05]"),

  make_option("--top_n", type = "integer", default = 250,
              help = "Top N associations for TOP heatmap [default = 250]"),

  make_option("--min_sig_traits", type = "integer", default = 1,
              help = paste("Minimum number of significant trait associations (p < p_thresh)",
                           "a module must have to appear in the FOCUSED heatmap.",
                           "Set to 1 to show any module with at least one hit.",
                           "Set higher (e.g. 2 or 3) to show only multi-trait modules.",
                           "[default = 1]")),

  make_option("--full_width", type = "double", default = 12,
              help = "Full heatmap width in inches [default = 12]"),

  make_option("--full_height", type = "double", default = 12,
              help = "Full heatmap height in inches [default = 12]"),

  make_option("--top_width", type = "double", default = 9,
              help = "Top heatmap width in inches [default = 9]"),

  make_option("--top_height", type = "double", default = 6,
              help = "Top heatmap height in inches [default = 6]"),

  # Dedicated dimensions for the FOCUSED heatmap so it can be sized
  # independently from the FULL heatmap (which may have many more rows).
  make_option("--focused_width", type = "double", default = NULL,
              help = paste("FOCUSED heatmap width in inches.",
                           "Defaults to --full_width if not set.")),

  make_option("--focused_height", type = "double", default = 9,
              help = paste("FOCUSED heatmap height in inches [default = 9].",
                           "Set independently from --full_height so the filtered",
                           "module subset has a fixed presentation-ready size."))
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$stats_file))   stop("--stats_file is required")

if (!dir.exists(opt$project_root))  stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$stats_file))   stop("stats_file not found: ", opt$stats_file)

if (!is.null(opt$stats_file_v2) && !file.exists(opt$stats_file_v2))
  stop("stats_file_v2 not found: ", opt$stats_file_v2)
if (!is.null(opt$stats_file_v3) && !file.exists(opt$stats_file_v3))
  stop("stats_file_v3 not found: ", opt$stats_file_v3)

if (!is.null(opt$modules_rds) && !file.exists(opt$modules_rds))
  stop("modules_rds not found: ", opt$modules_rds)
if (!is.null(opt$modules_rds_v2) && !file.exists(opt$modules_rds_v2))
  stop("modules_rds_v2 not found: ", opt$modules_rds_v2)
if (!is.null(opt$modules_rds_v3) && !file.exists(opt$modules_rds_v3))
  stop("modules_rds_v3 not found: ", opt$modules_rds_v3)

module_dendro_distance <- tolower(opt$module_dendro_distance)
if (!module_dendro_distance %in% c("bicor", "pearson"))
  stop("--module_dendro_distance must be 'bicor' or 'pearson'")

module_order_mode <- tolower(opt$module_order_mode)
if (!module_order_mode %in% c("default", "dendro", "exposure_outcome"))
  stop("--module_order_mode must be one of: default, dendro, exposure_outcome")

module_order_score <- tolower(opt$module_order_score)
if (!module_order_score %in% c("p", "abs_cor"))
  stop("--module_order_score must be one of: p, abs_cor")

trait_order_mode <- tolower(opt$trait_order_mode)
if (!trait_order_mode %in% c("input", "default", "exposure_outcome"))
  stop("--trait_order_mode must be one of: input, default, exposure_outcome")
if (trait_order_mode == "default") trait_order_mode <- "input"

trait_order_score <- tolower(opt$trait_order_score)
if (!trait_order_score %in% c("p", "abs_cor"))
  stop("--trait_order_score must be one of: p, abs_cor")

if (trait_order_mode == "exposure_outcome" && module_order_mode != "exposure_outcome") {
  stop("--trait_order_mode exposure_outcome requires --module_order_mode exposure_outcome")
}

if (module_order_mode == "exposure_outcome") {
  if (is.null(opt$module_order_trait_file) || !file.exists(opt$module_order_trait_file)) {
    stop("--module_order_trait_file is required and must exist when --module_order_mode exposure_outcome")
  }
}

if (opt$p_thresh <= 0 || opt$p_thresh >= 1) stop("--p_thresh must be > 0 and < 1")
if (opt$top_n < 1) stop("--top_n must be >= 1")
if (opt$min_sig_traits < 1) stop("--min_sig_traits must be >= 1")
if (!is.null(opt$focused_width) && opt$focused_width <= 0)  stop("--focused_width must be > 0")
if (opt$focused_height <= 0) stop("--focused_height must be > 0")

module_order_p_thresh <- opt$module_order_p_thresh
if (is.na(module_order_p_thresh)) module_order_p_thresh <- opt$p_thresh
if (module_order_p_thresh <= 0 || module_order_p_thresh >= 1) {
  stop("--module_order_p_thresh must be > 0 and < 1")
}

trait_order_p_thresh <- opt$trait_order_p_thresh
if (is.na(trait_order_p_thresh)) trait_order_p_thresh <- module_order_p_thresh
if (trait_order_p_thresh <= 0 || trait_order_p_thresh >= 1) {
  stop("--trait_order_p_thresh must be > 0 and < 1")
}

set_files <- collect_set_files(
  set_dir   = opt$set_dir,
  set_file  = opt$set_file,
  set_file2 = opt$set_file2,
  set_file3 = opt$set_file3,
  set_file4 = opt$set_file4,
  set_file5 = opt$set_file5,
  recursive = isTRUE(opt$set_dir_recursive)
)

# ------------------------------------------------------------
# Configure cache
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)
WGCNA::enableWGCNAThreads()

pipeline_root <- file.path(opt$project_root, "comethyl_output")

# ------------------------------------------------------------
# Build variant list and run each
# ------------------------------------------------------------
variant_list <- list(
  list(stats_file = opt$stats_file, modules_rds = opt$modules_rds)
)

if (!is.null(opt$stats_file_v2)) {
  variant_list[[length(variant_list) + 1]] <- list(
    stats_file  = opt$stats_file_v2,
    modules_rds = opt$modules_rds_v2
  )
}

if (!is.null(opt$stats_file_v3)) {
  variant_list[[length(variant_list) + 1]] <- list(
    stats_file  = opt$stats_file_v3,
    modules_rds = opt$modules_rds_v3
  )
}

message("Variants to process: ", length(variant_list))
message("Trait-set files to process: ", length(set_files))
message("Module order mode: ", module_order_mode)
message("Trait order mode: ", trait_order_mode)

for (v in variant_list) {
  run_presentation_for_variant(
    stats_file              = v$stats_file,
    modules_rds             = v$modules_rds,
    set_files               = set_files,
    pipeline_root           = pipeline_root,
    module_dendro_distance  = module_dendro_distance,
    module_order_mode       = module_order_mode,
    module_order_trait_file = opt$module_order_trait_file,
    module_order_p_thresh   = module_order_p_thresh,
    module_order_score      = module_order_score,
    trait_order_mode        = trait_order_mode,
    trait_order_p_thresh    = trait_order_p_thresh,
    trait_order_score       = trait_order_score,
    write_reordered_trait_files = isTRUE(opt$write_reordered_trait_files),
    drop_grey_modules       = isTRUE(opt$drop_grey_modules),
    p_thresh                = opt$p_thresh,
    top_n                   = opt$top_n,
    min_sig_traits          = opt$min_sig_traits,
    full_width              = opt$full_width,
    full_height             = opt$full_height,
    top_width               = opt$top_width,
    top_height              = opt$top_height,
    focused_width           = if (is.null(opt$focused_width)) opt$full_width else opt$focused_width,
    focused_height          = opt$focused_height
  )
}

message("Script 09b complete: Presentation ME-trait heatmaps finished")
message("Outputs saved under: ",
        file.path(pipeline_root, "09b_me_trait_presentation"))