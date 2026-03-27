#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 03: Calculate Region-Level Methylation
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   Calculates region-level methylation values using:
#     1) the filtered BSseq object from Script 01
#     2) the filtered region definitions from Script 02
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --input_bs       : path to Filtered_BSseq.rds from Script 01
#   --input_regions  : path to Filtered_Regions.txt from Script 02
#
# OUTPUTS
#   project_root/comethyl_output/03_region_methylation/<cpg_label>/<region_label>/
#       Region_Methylation.rds
#       run_parameters.txt
#
# EXAMPLE INPUT LINEAGE
#   input_bs:
#     .../01_filter_regions/cov3_75pct/Filtered_BSseq.rds
#
#   input_regions:
#     .../02_filter_regions/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt
#
#   output:
#     .../03_region_methylation/cov3_75pct/covMin4_methSD0p08/Region_Methylation.rds
#
# NOTES
#   - The BSseq object and filtered regions should be derived from the
#     same upstream filtering lineage.
#   - This script assumes region IDs and coordinates in Filtered_Regions.txt
#     are compatible with the filtered BSseq object.
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(AnnotationHub)
})

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--input_bs", type = "character",
              help = "Path to Filtered_BSseq.rds from Script 01"),
  make_option("--input_regions", type = "character",
              help = "Path to Filtered_Regions.txt from Script 02")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$input_bs)) stop("--input_bs is required")
if (is.null(opt$input_regions)) stop("--input_regions is required")

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
# 3. Configure AnnotationHub cache
# ------------------------------------------------------------

AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

# ------------------------------------------------------------
# 4. Derive output directory automatically
# ------------------------------------------------------------
# We derive labels from the parent folders of the two inputs:
#
# input_bs:
#   .../01_filter_regions/cov3_75pct/Filtered_BSseq.rds
#   -> cpg_label = cov3_75pct
#
# input_regions:
#   .../02_filter_regions/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt
#   -> region_label = covMin4_methSD0p08

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "03_region_methylation")

cpg_label <- basename(dirname(opt$input_bs))
region_label <- basename(dirname(opt$input_regions))

out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# 5. Define output file
# ------------------------------------------------------------

meth_out_path <- file.path(out_dir, "Region_Methylation.rds")

# ------------------------------------------------------------
# 6. Load inputs
# ------------------------------------------------------------

bs <- readRDS(opt$input_bs)
regions <- read.delim(opt$input_regions, header = TRUE, sep = "\t")

# Basic sanity check
required_cols <- c("RegionID")
missing_cols <- setdiff(required_cols, colnames(regions))
if (length(missing_cols) > 0) {
  stop("Filtered regions file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

message("Loaded filtered regions: ", nrow(regions))
message("Loaded BSseq object successfully")

# ------------------------------------------------------------
# 7. Calculate region methylation
# ------------------------------------------------------------

meth <- getRegionMeth(
  regions,
  bs = bs,
  file = meth_out_path
)

message("Saved region methylation RDS: ", meth_out_path)

# ------------------------------------------------------------
# 8. Save run parameters for reproducibility
# ------------------------------------------------------------

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("input_bs:", opt$input_bs),
    paste("input_regions:", opt$input_regions),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("output_dir:", out_dir),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(out_dir, "run_parameters.txt")
)

message("✓ Script 03 complete: region methylation calculated")