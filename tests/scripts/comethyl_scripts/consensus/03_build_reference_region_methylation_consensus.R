#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 03: Calculate Reference Region-Level Methylation
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Calculates region-level methylation values for the reference
#   dataset using:
#     1) the filtered BSseq object from consensus Script 01
#     2) the filtered region definitions from consensus Script 02
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --dataset_label  : reference dataset label (e.g. Baseline)
#   --input_bs       : path to Filtered_BSseq.rds from consensus Script 01
#   --input_regions  : path to Filtered_Regions.txt or Filtered_Regions.rds
#                      from consensus Script 02
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/03_region_methylation/<dataset_label>/<cpg_label>/<region_label>/
#       Region_Methylation.rds
#       run_parameters.txt
#       sessionInfo.txt
#
# EXAMPLE INPUT LINEAGE
#   input_bs:
#     .../consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds
#
#   input_regions:
#     .../consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt
#
#   output:
#     .../consensus/03_region_methylation/Baseline/cov3_75pct/covMin4_methSD0p08/Region_Methylation.rds
#
# NOTES
#   - This script should be run only on the reference dataset.
#   - The BSseq object and filtered regions should be derived from the
#     same upstream filtering lineage.
#   - This script assumes region IDs and coordinates in
#     Filtered_Regions.txt are compatible with the filtered BSseq object.

# Example
# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/03_build_reference_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
#   --dataset_label Baseline \
#   --input_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
#   --input_regions //quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p05/Filtered_Regions.rds

# ================================================================
message("Starting Script 03")
suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(data.table)
})

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--dataset_label", type = "character",
              help = "Reference dataset label (e.g. Baseline)"),
  make_option("--input_bs", type = "character",
              help = "Path to Filtered_BSseq.rds from consensus Script 01"),
  make_option("--input_regions", type = "character",
              help = "Path to Filtered_Regions.txt or Filtered_Regions.rds from consensus Script 02")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$dataset_label)) stop("--dataset_label is required")
if (is.null(opt$input_bs)) stop("--input_bs is required")
if (is.null(opt$input_regions)) stop("--input_regions is required")

opt$dataset_label <- trimws(opt$dataset_label)

if (!nzchar(opt$dataset_label)) {
  stop("--dataset_label cannot be empty")
}
if (!dir.exists(opt$project_root)) {
  stop("Project root does not exist: ", opt$project_root)
}
if (!file.exists(opt$input_bs)) {
  stop("Filtered BSseq file not found: ", opt$input_bs)
}
if (!file.exists(opt$input_regions)) {
  stop("Filtered regions file not found: ", opt$input_regions)
}

# ------------------------------------------------------------
# 3. Derive output directory automatically
# ------------------------------------------------------------
# input_bs:
#   .../consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds
#   -> cpg_label = cov3_75pct
#
# input_regions:
#   .../consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt
#   -> region_label = covMin4_methSD0p08

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "03_region_methylation")

cpg_label <- basename(dirname(opt$input_bs))
region_label <- basename(dirname(opt$input_regions))

out_dir <- file.path(step_dir, opt$dataset_label, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# 4. Define output file
# ------------------------------------------------------------

meth_out_path <- file.path(out_dir, "Region_Methylation.rds")

# ------------------------------------------------------------
# 5. Load inputs
# ------------------------------------------------------------

bs <- readRDS(opt$input_bs)

if (!inherits(bs, "BSseq")) {
  stop("Input object is not a BSseq object: ", opt$input_bs)
}

if (grepl("\\.rds$", opt$input_regions, ignore.case = TRUE)) {
  regions <- readRDS(opt$input_regions)
} else {
  regions <- data.table::fread(opt$input_regions, data.table = FALSE)
}

if (!is.data.frame(regions)) {
  stop("Filtered regions did not load as a data.frame.")
}

if (nrow(regions) == 0) {
  stop("Filtered regions file contains zero rows.")
}

required_cols <- c("RegionID")
missing_cols <- setdiff(required_cols, colnames(regions))
if (length(missing_cols) > 0) {
  stop("Filtered regions file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(regions$RegionID)) {
  dup_ids <- unique(regions$RegionID[duplicated(regions$RegionID)])
  stop("Duplicate RegionID values found. Example: ",
       paste(head(dup_ids, 10), collapse = ", "))
}

message("Loaded filtered regions: ", nrow(regions))
message("Loaded BSseq object successfully")

# ------------------------------------------------------------
# 6. Calculate region methylation
# ------------------------------------------------------------

meth <- getRegionMeth(
  regions,
  bs = bs,
  file = meth_out_path
)

if (!(is.matrix(meth) || is.data.frame(meth))) {
  stop("getRegionMeth() did not return a matrix-like object.")
}

meth <- as.matrix(meth)

if (!is.numeric(meth)) {
  stop("Region methylation matrix is not numeric.")
}
if (nrow(meth) == 0) {
  stop("Region methylation matrix has zero rows.")
}
if (ncol(meth) < 2) {
  stop("Region methylation matrix must have at least 2 samples.")
}
if (is.null(colnames(meth))) {
  stop("Region methylation matrix must have sample IDs in colnames.")
}
if (anyDuplicated(colnames(meth))) {
  dup_ids <- unique(colnames(meth)[duplicated(colnames(meth))])
  stop("Duplicate sample IDs found in region methylation matrix. Example: ",
       paste(head(dup_ids, 10), collapse = ", "))
}

message("Saved region methylation RDS: ", meth_out_path)
message("Matrix dimensions: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ------------------------------------------------------------
# 7. Save run parameters for reproducibility
# ------------------------------------------------------------

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("dataset_label:", opt$dataset_label),
    paste("input_bs:", opt$input_bs),
    paste("input_regions:", opt$input_regions),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("n_regions:", nrow(meth)),
    paste("n_samples:", ncol(meth)),
    paste("output_dir:", out_dir),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(out_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

message("Consensus Script 03 complete: reference region methylation calculated")