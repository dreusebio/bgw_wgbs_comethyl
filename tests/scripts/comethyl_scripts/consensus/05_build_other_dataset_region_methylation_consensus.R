#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 05: Build Other-Dataset Region Methylation for Consensus
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   For one or two non-reference datasets:
#     1) loads filtered BSseq objects aligned to reference CpGs
#     2) loads the canonical filtered reference regions
#     3) calculates region-level methylation using those same regions
#     4) saves matched region methylation matrices for downstream
#        consensus soft-power and consensus module analysis
#
# REQUIRED INPUTS
#   --project_root      : root directory of the analysis project
#   --reference_regions : path to reference Filtered_Regions.txt or .rds
#
# DATASET 1 INPUTS
#   --dataset1_label    : label for non-reference dataset 1
#   --dataset1_bs       : path to dataset 1 Filtered_BSseq.rds
#
# OPTIONAL DATASET 2 INPUTS
#   --dataset2_label
#   --dataset2_bs
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/03_region_methylation/<dataset_label>/<cpg_label>/<region_label>/
#       Region_Methylation.rds
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - This script is intended for non-reference datasets only.
#   - The same reference filtered regions are applied to all datasets.
#   - Slight CpG-level differences across datasets are acceptable as long
#     as region-level matrices are successfully generated over the same
#     reference region definitions.

# Example
# Rscript scripts/consensus/05_build_other_dataset_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
#   --reference_regions /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.rds \
#   --dataset1_label TP36_38weeks \
#   --dataset1_bs /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/04_align_to_reference_cpgs/TP36_38weeks/Filtered_BSseq.rds

# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/05_build_other_dataset_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --reference_regions /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p05/Filtered_Regions.rds \
#   --dataset1_label 36-38wk \
#   --dataset1_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/04_align_to_reference_cpgs/36-38wk/Filtered_BSseq.rds \
#   --dataset2_label Postpartum \
#   --dataset2_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/04_align_to_reference_cpgs/Postpartum/Filtered_BSseq.rds

# ================================================================
message("Starting Script 05")
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
  make_option("--reference_regions", type = "character",
              help = "Path to reference Filtered_Regions.txt or .rds"),

  make_option("--dataset1_label", type = "character",
              help = "Dataset 1 label"),
  make_option("--dataset1_bs", type = "character",
              help = "Path to dataset 1 Filtered_BSseq.rds"),

  make_option("--dataset2_label", type = "character", default = NULL,
              help = "Dataset 2 label"),
  make_option("--dataset2_bs", type = "character", default = NULL,
              help = "Path to dataset 2 Filtered_BSseq.rds")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$reference_regions)) stop("--reference_regions is required")

if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_bs)) stop("--dataset1_bs is required")

opt$dataset1_label <- trimws(opt$dataset1_label)
if (!nzchar(opt$dataset1_label)) stop("--dataset1_label cannot be empty")

if (!dir.exists(opt$project_root)) {
  stop("Project root does not exist: ", opt$project_root)
}
if (!file.exists(opt$reference_regions)) {
  stop("Reference filtered regions not found: ", opt$reference_regions)
}
if (!file.exists(opt$dataset1_bs)) {
  stop("dataset1_bs not found: ", opt$dataset1_bs)
}

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_bs)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_bs)) {
    stop("If dataset 2 is used, provide both --dataset2_label and --dataset2_bs.")
  }
  opt$dataset2_label <- trimws(opt$dataset2_label)
  if (!nzchar(opt$dataset2_label)) stop("--dataset2_label cannot be empty")
  if (!file.exists(opt$dataset2_bs)) {
    stop("dataset2_bs not found: ", opt$dataset2_bs)
  }
}

# ------------------------------------------------------------
# 3. Load reference regions
# ------------------------------------------------------------

if (grepl("\\.rds$", opt$reference_regions, ignore.case = TRUE)) {
  reference_regions <- readRDS(opt$reference_regions)
} else {
  reference_regions <- data.table::fread(opt$reference_regions, data.table = FALSE)
}

if (!is.data.frame(reference_regions)) {
  stop("Reference regions did not load as a data.frame.")
}
if (nrow(reference_regions) == 0) {
  stop("Reference regions contain zero rows.")
}

required_cols <- c("RegionID")
missing_cols <- setdiff(required_cols, colnames(reference_regions))
if (length(missing_cols) > 0) {
  stop("Reference regions are missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

if (anyDuplicated(reference_regions$RegionID)) {
  dup_ids <- unique(reference_regions$RegionID[duplicated(reference_regions$RegionID)])
  stop("Duplicate RegionID values in reference regions. Example: ",
       paste(head(dup_ids, 10), collapse = ", "))
}

message("Loaded reference filtered regions: ", nrow(reference_regions))

# Derive cpg_label and region_label from --reference_regions path:
#   .../02_reference_region_filter/<dataset_label>/<cpg_label>/<region_label>/Filtered_Regions.rds
# Walking up from the file:
#   dirname     = <region_label>/
#   dirname x2  = <cpg_label>/
{
  ref_dir      <- dirname(opt$reference_regions)
  region_label <- basename(ref_dir)
  cpg_label    <- basename(dirname(ref_dir))
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}

# ------------------------------------------------------------
# 4. Helper function
# ------------------------------------------------------------

process_dataset <- function(dataset_label, input_bs, reference_regions, cpg_label, region_label, project_root) {
  message("\n====================================================")
  message("Processing dataset: ", dataset_label)
  message("====================================================")

  bs <- readRDS(input_bs)

  if (!inherits(bs, "BSseq")) {
    stop("Input BS object is not a BSseq object for ", dataset_label, ": ", input_bs)
  }

  out_dir <- file.path(
    project_root,
    "comethyl_output",
    "consensus",
    "03_region_methylation",
    dataset_label,
    cpg_label,
    region_label
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  message("Output directory: ", out_dir)

  meth_out_path <- file.path(out_dir, "Region_Methylation.rds")

  meth <- getRegionMeth(
    reference_regions,
    bs = bs,
    file = meth_out_path
  )

  if (!(is.matrix(meth) || is.data.frame(meth))) {
    stop("getRegionMeth() did not return a matrix-like object for ", dataset_label)
  }

  meth <- as.matrix(meth)

  if (!is.numeric(meth)) {
    stop("Region methylation matrix is not numeric for ", dataset_label)
  }
  if (nrow(meth) == 0) {
    stop("Region methylation matrix has zero rows for ", dataset_label)
  }
  if (ncol(meth) < 2) {
    stop("Region methylation matrix must have at least 2 samples for ", dataset_label)
  }
  if (is.null(colnames(meth))) {
    stop("Region methylation matrix must have sample IDs in colnames for ", dataset_label)
  }
  if (anyDuplicated(colnames(meth))) {
    dup_ids <- unique(colnames(meth)[duplicated(colnames(meth))])
    stop("Duplicate sample IDs found in region methylation matrix for ", dataset_label,
         ". Example: ", paste(head(dup_ids, 10), collapse = ", "))
  }

  # Check whether region count matches reference regions
  region_count_matches_reference <- identical(nrow(meth), nrow(reference_regions))

  writeLines(
    c(
      paste("dataset_label:", dataset_label),
      paste("input_bs:", input_bs),
      paste("reference_regions:", opt$reference_regions),
      paste("cpg_label:", cpg_label),
      paste("region_label:", region_label),
      paste("n_reference_regions:", nrow(reference_regions)),
      paste("n_matrix_regions:", nrow(meth)),
      paste("n_samples:", ncol(meth)),
      paste("region_count_matches_reference:", region_count_matches_reference),
      paste("output_dir:", out_dir),
      paste("date:", as.character(Sys.time()))
    ),
    con = file.path(out_dir, "run_parameters.txt")
  )

  writeLines(
    capture.output(sessionInfo()),
    con = file.path(out_dir, "sessionInfo.txt")
  )

  message("✓ Completed dataset: ", dataset_label)
  message("  Matrix dimensions: ", nrow(meth), " regions x ", ncol(meth), " samples")

  invisible(meth)
}

# ------------------------------------------------------------
# 5. Run dataset 1
# ------------------------------------------------------------

process_dataset(
  dataset_label     = opt$dataset1_label,
  input_bs          = opt$dataset1_bs,
  reference_regions = reference_regions,
  cpg_label         = cpg_label,
  region_label      = region_label,
  project_root      = opt$project_root
)

# ------------------------------------------------------------
# 6. Run dataset 2 if provided
# ------------------------------------------------------------

if (dataset2_provided) {
  process_dataset(
    dataset_label = opt$dataset2_label,
    input_bs          = opt$dataset2_bs,
    reference_regions = reference_regions,
    cpg_label         = cpg_label,
    region_label      = region_label,
    project_root      = opt$project_root
  )
}

message("Consensus Script 05 complete")