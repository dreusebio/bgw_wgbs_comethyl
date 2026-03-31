#!/usr/bin/env Rscript
# ============================================================
# 14_module_region_methylation_heatmaps.R
#
# Purpose:
#   Generate module-region methylation heatmaps using comethyl::plotMethTrait()
#   from:
#     - region-level methylation matrix
#     - modules RDS
#     - sample info / trait table
#     - user-supplied plot specification file
#
# Inputs:
#   --project_root
#   --meth_rds
#   --modules_rds
#   --sample_info
#   --plot_spec
#   --sample_id_col           optional
#   --out_dir                 optional
#   --trait_exclude_file      optional
#
# Plot spec required columns:
#   module
#   trait
#
# Optional plot spec columns:
#   plot_type                 binary|continuous
#   trait_legend_title
#   trait_code                e.g. "No=0;Yes=1"
#   trait_colors              e.g. "No=#3366CC;Yes=#FF3366"
#   expandY                   numeric
#   trait_legend_position     e.g. "1.034,3.35"
#   file_stub
#
# Output:
#   one PDF per plot
#   run log
# ============================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(comethyl)
  library(Biobase)
})

# ============================================================
# 1) CLI helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

msg <- function(...) cat(sprintf(...), "\n")

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 2) parse helpers
# ============================================================
parse_named_mapping <- function(x) {
  # "No=0;Yes=1" -> named vector c(No="0", Yes="1")
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  parts <- trimws(strsplit(x, ";")[[1]])
  kv <- strsplit(parts, "=")
  keys <- vapply(kv, function(z) trimws(z[1]), character(1))
  vals <- vapply(kv, function(z) trimws(z[2]), character(1))
  stats::setNames(vals, keys)
}

parse_numeric_named_mapping <- function(x) {
  out <- parse_named_mapping(x)
  if (is.null(out)) return(NULL)
  out_num <- suppressWarnings(as.numeric(out))
  stats::setNames(out_num, names(out))
}

parse_position <- function(x) {
  # "1.034,3.35" -> c(1.034, 3.35)
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  vals <- trimws(strsplit(x, ",")[[1]])
  vals <- suppressWarnings(as.numeric(vals))
  if (length(vals) != 2 || any(is.na(vals))) return(NULL)
  vals
}

# ============================================================
# 3) trait loading / cleaning
# ============================================================
load_sample_info <- function(path, sample_id_col = NULL) {
  if (!file.exists(path)) stop("sample_info not found: ", path)

  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    df <- openxlsx::read.xlsx(path, rowNames = is.null(sample_id_col))
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported sample_info format: ", path)
  }

  if (!is.null(sample_id_col) && nzchar(sample_id_col)) {
    if (!(sample_id_col %in% colnames(df))) {
      stop("sample_id_col not found in sample_info: ", sample_id_col)
    }
    rownames(df) <- as.character(df[[sample_id_col]])
  }

  df
}

clean_sample_info <- function(df) {
  df[] <- lapply(df, function(x) {
    x_chr <- as.character(x)
    suppressWarnings({
      min_val <- min(as.numeric(x_chr[x_chr != "M"]), na.rm = TRUE)
      if (is.finite(min_val)) {
        x_chr[x_chr == "M"] <- as.character(min_val)
      }
      x_chr[x_chr %in% c("N", ".", "")] <- NA
      out_num <- suppressWarnings(as.numeric(x_chr))
      if (all(is.na(out_num) == is.na(x_chr))) {
        return(out_num)
      } else {
        return(x_chr)
      }
    })
  })

  df[, colSums(is.na(df)) < nrow(df), drop = FALSE]
}

apply_trait_exclusions <- function(df, trait_exclude_file = NULL) {
  if (is.null(trait_exclude_file) || is.na(trait_exclude_file) || !nzchar(trait_exclude_file)) {
    return(df)
  }
  if (!file.exists(trait_exclude_file)) {
    stop("trait_exclude_file not found: ", trait_exclude_file)
  }
  drops <- readLines(trait_exclude_file, warn = FALSE)
  drops <- trimws(drops)
  drops <- drops[drops != ""]
  df[, !colnames(df) %in% drops, drop = FALSE]
}

# ============================================================
# 4) load plot spec
# ============================================================
load_plot_spec <- function(path) {
  if (!file.exists(path)) stop("plot_spec not found: ", path)

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    spec <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    spec <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  required <- c("module", "trait")
  missing <- setdiff(required, colnames(spec))
  if (length(missing)) {
    stop("plot_spec missing required columns: ", paste(missing, collapse = ", "))
  }

  spec
}

# ============================================================
# 5) args
# ============================================================
project_root       <- get_arg("--project_root", getwd())
meth_rds           <- get_arg("--meth_rds", NULL)
modules_rds        <- get_arg("--modules_rds", NULL)
sample_info        <- get_arg("--sample_info", NULL)
plot_spec_file     <- get_arg("--plot_spec", NULL)
sample_id_col      <- get_arg("--sample_id_col", NULL)
trait_exclude_file <- get_arg("--trait_exclude_file", NULL)
out_dir            <- get_arg("--out_dir", file.path(project_root, "module_region_methylation_heatmaps"))

if (is.null(meth_rds) || is.null(modules_rds) || is.null(sample_info) || is.null(plot_spec_file)) {
  stop("Required: --meth_rds --modules_rds --sample_info --plot_spec")
}

setwd(project_root)
safe_dir_create(out_dir)

log_file <- file.path(out_dir, "run_log.txt")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

msg("project_root: %s", project_root)
msg("meth_rds: %s", meth_rds)
msg("modules_rds: %s", modules_rds)
msg("sample_info: %s", sample_info)
msg("plot_spec: %s", plot_spec_file)
msg("out_dir: %s", out_dir)

# ============================================================
# 6) load data
# ============================================================
meth <- readRDS(meth_rds)
modules <- readRDS(modules_rds)

if (!("regions" %in% names(modules))) {
  stop("modules_rds does not contain $regions")
}
regions <- modules$regions

colData_raw <- load_sample_info(sample_info, sample_id_col = sample_id_col)
colData <- clean_sample_info(colData_raw)
colData <- apply_trait_exclusions(colData, trait_exclude_file = trait_exclude_file)

spec <- load_plot_spec(plot_spec_file)

msg("meth dims: %d x %d", nrow(meth), ncol(meth))
msg("regions dims: %d x %d", nrow(regions), ncol(regions))
msg("sample info dims: %d x %d", nrow(colData), ncol(colData))
msg("plot spec rows: %d", nrow(spec))

# ============================================================
# 7) align samples
# ============================================================
common_samples <- intersect(colnames(meth), rownames(colData))
if (length(common_samples) == 0) {
  stop("No overlapping samples between meth colnames and sample_info rownames")
}

meth_use <- meth[, common_samples, drop = FALSE]
colData_use <- colData[common_samples, , drop = FALSE]

msg("aligned samples: %d", length(common_samples))

# ============================================================
# 8) run plots
# ============================================================
for (i in seq_len(nrow(spec))) {
  row <- spec[i, , drop = FALSE]

  module_name <- as.character(row$module[1])
  trait_name  <- as.character(row$trait[1])

  if (!(trait_name %in% colnames(colData_use))) {
    msg("[SKIP] trait not found: %s", trait_name)
    next
  }

  trait_vec <- colData_use[[trait_name]]

  plot_type <- if ("plot_type" %in% colnames(spec)) as.character(row$plot_type[1]) else NA_character_
  legend_title <- if ("trait_legend_title" %in% colnames(spec)) as.character(row$trait_legend_title[1]) else trait_name
  file_stub <- if ("file_stub" %in% colnames(spec) && nzchar(row$file_stub[1])) as.character(row$file_stub[1]) else paste0(module_name, "_", trait_name)

  out_file <- file.path(out_dir, paste0(file_stub, ".pdf"))

  trait_code <- if ("trait_code" %in% colnames(spec)) parse_numeric_named_mapping(row$trait_code[1]) else NULL
  trait_colors <- if ("trait_colors" %in% colnames(spec)) parse_named_mapping(row$trait_colors[1]) else NULL
  expandY <- if ("expandY" %in% colnames(spec)) suppressWarnings(as.numeric(row$expandY[1])) else NULL
  trait_legend_position <- if ("trait_legend_position" %in% colnames(spec)) parse_position(row$trait_legend_position[1]) else NULL

  msg("[RUN] module=%s | trait=%s | out=%s", module_name, trait_name, out_file)

  args <- list(
    module = module_name,
    regions = regions,
    meth = meth_use,
    trait = trait_vec,
    file = out_file
  )

  if (!is.na(plot_type) && plot_type == "binary") {
    if (!is.null(trait_code)) args$traitCode <- trait_code
    if (!is.null(trait_colors)) args$traitColors <- trait_colors
  }

  if (!is.null(legend_title) && nzchar(legend_title)) {
    args$trait.legend.title <- legend_title
  }

  if (!is.null(expandY) && is.finite(expandY)) {
    args$expandY <- expandY
  }

  if (!is.null(trait_legend_position)) {
    args$trait.legend.position <- trait_legend_position
  }

  tryCatch({
    do.call(plotMethTrait, args)
    msg("[OK] %s", out_file)
  }, error = function(e) {
    msg("[FAIL] module=%s trait=%s :: %s", module_name, trait_name, conditionMessage(e))
  })
}

msg("[✓] Done.")