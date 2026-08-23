#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 01: Filter CpGs and Define Reference Regions for Consensus
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Filters CpGs from the unfiltered BSseq object for the reference
#   dataset, then defines genomic regions that will be reused across
#   all datasets in the consensus workflow.
#
# REQUIRED INPUTS
#   --project_root  : root directory of the analysis project
#   --dataset_label : reference dataset label
#                     e.g. Baseline
#   --input_rds     : path to Unfiltered_BSseq.rds from consensus Script 00
#
# OPTIONAL INPUTS
#   --cov           : minimum read coverage per CpG [default = 3]
#   --per_sample    : fraction of samples meeting coverage threshold
#                     [default = 0.75]
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/01_reference_filter_regions/<dataset_label>/cov<cov>_<pct>pct>/
#       Filtered_BSseq.rds
#       Regions.txt
#       Regions.rds
#       Region_Plots.pdf
#       SD_Plots.pdf
#       Region_Totals.txt
#       Region_Totals.pdf
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - This script should be run only on the reference dataset.
#   - The resulting regions are the canonical region definitions for
#     downstream consensus analysis.
#   - All other datasets must reuse these same regions rather than
#     defining their own.

# Example

# Rscript /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/scripts/comethyl_scripts/consensus/01_filter_cpgs_call_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
#   --dataset_label Baseline \
#   --input_rds /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/00_cpg_extraction/Baseline/Unfiltered_BSseq.rds \
#   --cov 3 \
#   --per_sample 0.75

message("Starting Script 01 ")

# ================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
})

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--dataset_label", type = "character",
              help = "Reference dataset label (e.g. Baseline)"),
  make_option("--input_rds", type = "character",
              help = "Path to Unfiltered_BSseq.rds from consensus Script 00"),
  make_option("--cov", type = "integer", default = 3,
              help = "Minimum CpG coverage [default = 3]"),
  make_option("--per_sample", type = "double", default = 0.75,
              help = "Fraction of samples passing coverage [default = 0.75]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$dataset_label)) stop("--dataset_label is required")
if (is.null(opt$input_rds)) stop("--input_rds is required")

opt$dataset_label <- trimws(opt$dataset_label)

if (!nzchar(opt$dataset_label)) {
  stop("--dataset_label cannot be empty")
}
if (!dir.exists(opt$project_root)) {
  stop("Project root does not exist: ", opt$project_root)
}
if (!file.exists(opt$input_rds)) {
  stop("Input BSseq file not found: ", opt$input_rds)
}
if (opt$cov < 1) {
  stop("--cov must be >= 1")
}
if (opt$per_sample <= 0 || opt$per_sample > 1) {
  stop("--per_sample must be > 0 and <= 1")
}

# ------------------------------------------------------------
# 3. Define output directories
# ------------------------------------------------------------

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "01_reference_filter_regions")

percent_label <- round(opt$per_sample * 100)
version_label <- paste0("cov", opt$cov, "_", percent_label, "pct")

out_dir <- file.path(step_dir, opt$dataset_label, version_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# 4. Load unfiltered BSseq object
# ------------------------------------------------------------

unfilt_bs <- readRDS(opt$input_rds)

if (!inherits(unfilt_bs, "BSseq")) {
  stop("Input object is not a BSseq object: ", opt$input_rds)
}

# ------------------------------------------------------------
# 5. Filter CpGs
# ------------------------------------------------------------

filtered_rds <- file.path(out_dir, "Filtered_BSseq.rds")

bs <- filterCpGs(
  bs = unfilt_bs,
  cov = opt$cov,
  perSample = opt$per_sample,
  file = filtered_rds
)

if (!inherits(bs, "BSseq")) {
  stop("filterCpGs() did not return a BSseq object.")
}

if (nrow(bs) == 0) {
  stop("No CpGs retained after filtering.")
}

message("Filtered CpGs retained: ", nrow(bs))

# ------------------------------------------------------------
# 6. Call regions
# ------------------------------------------------------------

regions_txt <- file.path(out_dir, "Regions.txt")
regions_rds <- file.path(out_dir, "Regions.rds")

regions <- getRegions(
  bs = bs,
  file = regions_txt
)

if (!is.data.frame(regions)) {
  stop("getRegions() did not return a data.frame.")
}

if (nrow(regions) == 0) {
  stop("No regions were generated.")
}

if (!("RegionID" %in% colnames(regions))) {
  stop("Regions output does not contain RegionID column.")
}

if (anyDuplicated(regions$RegionID)) {
  dup_ids <- unique(regions$RegionID[duplicated(regions$RegionID)])
  stop("Duplicate RegionID values found. Example: ",
       paste(head(dup_ids, 10), collapse = ", "))
}

saveRDS(regions, regions_rds)

message("Regions generated: ", nrow(regions))

# ------------------------------------------------------------
# 7. Plot region and SD statistics
# ------------------------------------------------------------

plotRegionStats(
  regions,
  maxQuantile = 0.99,
  file = file.path(out_dir, "Region_Plots.pdf")
)

plotSDstats(
  regions,
  maxQuantile = 0.99,
  file = file.path(out_dir, "SD_Plots.pdf")
)

# ------------------------------------------------------------
# 8. Compute and plot region totals
# ------------------------------------------------------------

regionTotals <- getRegionTotals(
  regions,
  file = file.path(out_dir, "Region_Totals.txt")
)

plotRegionTotals(
  regionTotals,
  file = file.path(out_dir, "Region_Totals.pdf")
)

# ------------------------------------------------------------
# 9. Save reproducibility information
# ------------------------------------------------------------

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("dataset_label:", opt$dataset_label),
    paste("input_rds:", opt$input_rds),
    paste("cov:", opt$cov),
    paste("per_sample:", opt$per_sample),
    paste("n_filtered_cpgs:", nrow(bs)),
    paste("n_regions:", nrow(regions)),
    paste("output_dir:", out_dir),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(out_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

message("Consensus Script 01 complete for reference dataset: ", opt$dataset_label)