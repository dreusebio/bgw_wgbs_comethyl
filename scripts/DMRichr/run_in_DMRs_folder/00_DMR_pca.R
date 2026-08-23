#!/usr/bin/env Rscript
# ==============================================================================
# smooth_meth_pca_plot.R
#
# Generates a PCA plot from a DMRichR/comethyl-style smoothed methylation
# matrix ("DMR_individual_smoothed_methylation_filtered.txt") joined against
# a sample_info.xlsx covariate sheet.
#
# Usage:
#   Rscript smooth_meth_pca_plot.R \
#somewhere     --file /path/to/DMR_individual_smoothed_methylation_filtered.txt \
#     --sample_info /path/to/sample_info.xlsx \
#     --output /path/to/output_PCA.pdf \
#     [--outcome Diagnosis] [--skip_cols 16] [--width 7] [--height 6]
#
# Note: --outcome is optional. If omitted, or if the named column isn't
# present in sample_info, the script falls back to an uncoloured PCA plot
# instead of failing.
#
# Example (batch, from bash):
#   while read -r dmr si out; do
#     Rscript smooth_meth_pca_plot.R --file "$dmr" --sample_info "$si" --output "$out"
#   done < runs.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ---- CLI arguments ----------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "Path to DMR_individual_smoothed_methylation_filtered.txt [required]"),
  make_option(c("-s", "--sample_info"), type = "character", default = NULL,
              help = "Path to sample_info.xlsx [required]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path for output PCA PDF [required]"),
  make_option(c("-c", "--outcome"), type = "character", default = NULL,
              help = "Column in sample_info holding the outcome/grouping variable (e.g. Diagnosis) to colour points by. Optional — if omitted or not found in sample_info, an uncoloured PCA plot is produced instead."),
  make_option(c("--skip_cols"), type = "integer", default = 16,
              help = "Number of leading metadata columns to drop from the methylation matrix [default: %default]"),
  make_option(c("--width"), type = "double", default = 7,
              help = "PDF width in inches [default: %default]"),
  make_option(c("--height"), type = "double", default = 6,
              help = "PDF height in inches [default: %default]"),
  make_option(c("--frame_type"), type = "character", default = "norm",
              help = "ggfortify frame.type for autoplot [default: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "Build a PCA plot from smoothed methylation data and a sample_info workbook."
)
opt <- parse_args(opt_parser)

required_args <- c("file", "sample_info", "output")
missing_args <- required_args[sapply(opt[required_args], is.null)]
if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop("Missing required argument(s): ", paste(missing_args, collapse = ", "), call. = FALSE)
}
if (!file.exists(opt$file)) stop("Input file not found: ", opt$file, call. = FALSE)
if (!file.exists(opt$sample_info)) stop("sample_info file not found: ", opt$sample_info, call. = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggfortify)
  library(openxlsx)
})

# ---- Core function ------------------------------------------------------------
SmoothMethPCAPlot <- function(file, sample_info, name, outcome = NULL,
                               skip_cols = 16, width = 7, height = 6,
                               frame_type = "norm") {

  message("Reading methylation matrix: ", file)
  x <- read.delim(file, header = TRUE, check.names = FALSE)

  if (ncol(x) <= skip_cols) {
    stop("skip_cols (", skip_cols, ") >= number of columns (", ncol(x),
         ") in input file; nothing left to transpose.", call. = FALSE)
  }
  x <- x[, -seq_len(skip_cols)]

  t_x <- data.table::transpose(x)
  colnames(t_x) <- rownames(x)
  rownames(t_x) <- colnames(x)

  message("Reading sample info: ", sample_info)
  sample_info_df <- read.xlsx(sample_info, rowNames = TRUE)

  has_outcome <- !is.null(outcome) && nzchar(outcome)
  if (has_outcome && !outcome %in% colnames(sample_info_df)) {
    warning("Outcome column '", outcome, "' not found in sample_info (available columns: ",
            paste(colnames(sample_info_df), collapse = ", "),
            "). Proceeding with an uncoloured PCA plot.", call. = FALSE, immediate. = TRUE)
    has_outcome <- FALSE
  }

  y <- merge(sample_info_df, t_x, by = "row.names", all = FALSE)
  if (nrow(y) == 0) {
    stop("No overlapping sample IDs between methylation matrix and sample_info after merge.",
         call. = FALSE)
  }

  n_meta_cols <- ncol(sample_info_df) + 1  # +1 for the Row.names column merge() adds
  PCA <- prcomp(y[, -seq_len(n_meta_cols)], center = TRUE, scale. = TRUE)

  if (has_outcome) {
    plot <- autoplot(PCA, data = y, colour = outcome, frame = TRUE, frame.type = frame_type)
  } else {
    message("No outcome column supplied/found — plotting PCA without colour grouping.")
    plot <- autoplot(PCA, data = y)
  }

  plot <- plot +
    theme_classic() +
    theme(
      plot.background = element_rect(fill = "white"),
      axis.text = element_text(color = "black", size = 12, family = "sans"),
      axis.title.x = element_text(color = "black", size = 12, family = "sans"),
      axis.title.y = element_text(color = "black", size = 12, family = "sans"),
      axis.line = element_line(color = "black", linewidth = 0.2),
      legend.position = "right"
    )

  dir.create(dirname(name), showWarnings = FALSE, recursive = TRUE)
  pdf(file = name, width = width, height = height)
  print(plot)
  dev.off()
  message("Wrote: ", name)

  return(plot)
}

# ---- Run ----------------------------------------------------------------------
result <- SmoothMethPCAPlot(
  file        = opt$file,
  sample_info = opt$sample_info,
  name        = opt$output,
  outcome     = opt$outcome,
  skip_cols   = opt$skip_cols,
  width       = opt$width,
  height      = opt$height,
  frame_type  = opt$frame_type
)