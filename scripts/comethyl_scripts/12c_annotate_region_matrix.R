#!/usr/bin/env Rscript

# ============================================================
# 12c_annotate_region_matrix.R
#
# Purpose:
#   Annotate ALL regions in Region_Methylation.rds with nearest
#   gene symbols -- independent of module assignment. Useful for
#   generating a full RegionID -> gene_symbol lookup table before
#   modules exist, or for QC/exploration alongside the matrix.
#
# Input:
#   --filtered_regions   Filtered_Regions.txt from step 02
#                         (columns: RegionID, chr, start, end, ...)
#   --methylation_matrix  (optional) Region_Methylation.rds, used
#                         only to sanity-check that RegionIDs match
#                         1:1 with filtered_regions before annotating
#
# Output (written to <project_root>/comethyl_output/03_region_methylation/
#         <cpg_label>/<region_label>/Region_Annotation/):
#   Region_Gene_Annotation.tsv   RegionID, chr, start, end, gene_symbol, gene_entrezID
#   Region_Gene_Annotation.xlsx
#   run_parameters.txt
#   run_log.txt
#
# Usage:
#   Rscript 12c_annotate_region_matrix.R \
#     --project_root /path/to/project \
#     --filtered_regions /path/to/Filtered_Regions.txt \
#     --methylation_matrix /path/to/Region_Methylation.rds \
#     --genome hg38 \
#     --annotation_mode auto \
#     [--helper_file /path/to/helper.R]   # optional; defaults to
#                                          # helper.R alongside this script
#
# Depends on: helper.R (sourced for annotate_regions_safe() and
# friends -- see the "ANNOTATION HELPERS" section at the bottom of
# helper.R)
# ============================================================

message("Starting Script 12c: Annotate Region Matrix")

suppressPackageStartupMessages({
  library(comethyl)
  library(dplyr)
  library(openxlsx)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(AnnotationDbi)
})

# ------------------------------------------------------------
# Shared helpers
# ------------------------------------------------------------
# get_arg()/trim_or_null() live inside helper.R itself, so we need a
# minimal standalone version here just to read --helper_file before
# helper.R has been sourced.
.get_arg_bootstrap <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) NULL
)
if (is.null(script_dir) || !nzchar(script_dir)) script_dir <- getwd()

# Default: helper.R lives alongside this script. Override with
# --helper_file /full/path/to/helper.R if it lives elsewhere.
helper_file <- .get_arg_bootstrap("--helper_file", file.path(script_dir, "helper.R"))
if (!file.exists(helper_file)) {
  stop("helper.R not found at: ", helper_file,
       "\nPass --helper_file /path/to/helper.R if it lives elsewhere.", call. = FALSE)
}
source(helper_file)

# ------------------------------------------------------------
# Read arguments
# ------------------------------------------------------------
project_root       <- trim_or_null(get_arg("--project_root"))
filtered_regions_f <- trim_or_null(get_arg("--filtered_regions"))
meth_matrix_f      <- trim_or_null(get_arg("--methylation_matrix"))
genome             <- trim_or_null(get_arg("--genome", "hg38"))
annotation_mode    <- trim_or_null(get_arg("--annotation_mode", "auto"))

stop_if_missing(project_root, "--project_root")
stop_if_missing(filtered_regions_f, "--filtered_regions")

if (!dir.exists(project_root)) stop("project_root not found: ", project_root, call. = FALSE)
validate_file_exists(filtered_regions_f, "filtered_regions")
if (!is.null(meth_matrix_f)) validate_file_exists(meth_matrix_f, "methylation_matrix")

annotation_mode <- match.arg(annotation_mode, choices = c("auto", "great", "offline"))

setwd(project_root)

AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(project_root, ".cache")
)

# ------------------------------------------------------------
# Derive output dir from filtered_regions path
#   .../03_region_methylation/<cpg_label>/<region_label>/Filtered... (or wherever it lives)
#   We mirror whatever <cpg_label>/<region_label> we can find; if the
#   path doesn't match that shape, fall back to a flat output dir.
# ------------------------------------------------------------
derive_out_dir <- function(filtered_regions_path, project_root) {
  region_dir <- dirname(filtered_regions_path)      # .../covMin4_methSD0p08 (or similar)
  cpg_dir    <- dirname(region_dir)                 # .../cov3_75pct
  region_label <- basename(region_dir)
  cpg_label    <- basename(cpg_dir)

  out_dir <- file.path(
    project_root, "comethyl_output", "03_region_methylation",
    cpg_label, region_label, "Region_Annotation"
  )
  safe_dir_create(out_dir)
  out_dir
}

out_dir <- derive_out_dir(filtered_regions_f, project_root)

log_file    <- file.path(out_dir, "run_log.txt")
params_file <- file.path(out_dir, "run_parameters.txt")

append_log(log_file, "project_root: ", project_root)
append_log(log_file, "filtered_regions: ", filtered_regions_f)
append_log(log_file, "methylation_matrix: ", ifelse(is.null(meth_matrix_f), "(not provided)", meth_matrix_f))
append_log(log_file, "genome: ", genome)
append_log(log_file, "annotation_mode: ", annotation_mode)
append_log(log_file, "out_dir: ", out_dir)

# ------------------------------------------------------------
# Load filtered regions (has RegionID, chr, start, end, ...)
# ------------------------------------------------------------
filtered_regions <- read.delim(filtered_regions_f, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("RegionID", "chr", "start", "end")
missing_cols <- setdiff(required_cols, colnames(filtered_regions))
if (length(missing_cols) > 0) {
  stop("filtered_regions is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

append_log(log_file, "Loaded filtered_regions rows: ", nrow(filtered_regions))

# ------------------------------------------------------------
# Optional sanity check: RegionIDs match the methylation matrix rows
# ------------------------------------------------------------
if (!is.null(meth_matrix_f)) {
  meth_mat <- readRDS(meth_matrix_f)
  meth_ids <- rownames(meth_mat)

  n_meth <- length(meth_ids)
  n_filt <- nrow(filtered_regions)

  if (n_meth != n_filt) {
    append_log(log_file, "WARNING: row count mismatch -- methylation_matrix has ", n_meth,
               " regions, filtered_regions has ", n_filt)
  }

  missing_from_filtered <- setdiff(meth_ids, filtered_regions$RegionID)
  if (length(missing_from_filtered) > 0) {
    append_log(log_file, "WARNING: ", length(missing_from_filtered),
               " RegionIDs in methylation_matrix are absent from filtered_regions (e.g. ",
               paste(head(missing_from_filtered, 5), collapse = ", "), ")")
  } else {
    append_log(log_file, "OK: all methylation_matrix RegionIDs found in filtered_regions")
  }
}

# ------------------------------------------------------------
# Build regions_df for annotate_regions_safe()
#   No "module" column exists yet at this stage, so we relax the
#   required-columns check via req_cols.
# ------------------------------------------------------------
regions_df <- filtered_regions %>%
  dplyr::select(RegionID, chr, start, end)

annotated_regions <- suppressWarnings(
  annotate_regions_safe(
    regions_df = regions_df,
    genome = genome,
    annotation_mode = annotation_mode,
    file_txt = NULL,
    verbose = TRUE,
    req_cols = c("RegionID", "chr", "start", "end")
  )
)

append_log(log_file, "Annotated regions rows: ", nrow(annotated_regions))
append_log(log_file, "Regions with a gene_symbol: ",
           sum(!is.na(annotated_regions$gene_symbol) & annotated_regions$gene_symbol != ""))

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------
out_tsv  <- file.path(out_dir, "Region_Gene_Annotation.tsv")
out_xlsx <- file.path(out_dir, "Region_Gene_Annotation.xlsx")

write.table(annotated_regions, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
openxlsx::write.xlsx(annotated_regions, out_xlsx, rowNames = FALSE)

params <- c(
  paste0("timestamp\t", timestamp_now()),
  paste0("project_root\t", project_root),
  paste0("filtered_regions\t", filtered_regions_f),
  paste0("methylation_matrix\t", ifelse(is.null(meth_matrix_f), "", meth_matrix_f)),
  paste0("genome\t", genome),
  paste0("annotation_mode\t", annotation_mode),
  paste0("out_dir\t", out_dir),
  paste0("n_regions\t", nrow(filtered_regions)),
  paste0("n_annotated\t", nrow(annotated_regions))
)
write_lines_safe(params, params_file)

append_log(log_file, "Wrote: ", out_tsv)
append_log(log_file, "Wrote: ", out_xlsx)
message("Script 12c complete: Region matrix annotation finished")