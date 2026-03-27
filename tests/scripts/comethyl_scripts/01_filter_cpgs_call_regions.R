#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 01: Filter CpGs and Call Regions
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   Filters CpGs from the unfiltered BSseq object using user-defined
#   coverage and sample presence thresholds, then groups retained
#   CpGs into genomic regions for downstream network analysis.
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --input_rds    : path to Unfiltered_BSseq.rds from Script 00
#
# OPTIONAL INPUTS
#   --cov          : minimum read coverage per CpG [default = 3]
#   --per_sample   : fraction of samples meeting coverage threshold
#                    [default = 0.75]
#
# OUTPUTS
#   project_root/comethyl_output/01_filter_regions/cov<cov>_<pct>pct/
#       Filtered_BSseq.rds
#       Regions.txt
#       Region_Plots.pdf
#       SD_Plots.pdf
#       Region_Totals.txt
#       Region_Totals.pdf
#       run_parameters.txt
#
# EXAMPLES
#   cov = 3, per_sample = 0.75  -> cov3_75pct
#   cov = 5, per_sample = 0.85  -> cov5_85pct
#
# NOTES
#   - Parameter-based subfolders are created automatically so that
#     different filtering runs do not overwrite each other.
#   - This step is important for later consensus analysis because
#     region definitions must be handled consistently across datasets.
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
  make_option("--input_rds", type = "character",
              help = "Path to Unfiltered_BSseq.rds"),
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
if (is.null(opt$input_rds)) stop("--input_rds is required")

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
# 3. Configure AnnotationHub cache
# ------------------------------------------------------------

AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

# ------------------------------------------------------------
# 4. Define output directories automatically
# ------------------------------------------------------------
# All outputs for this step go under:
#   <project_root>/comethyl_output/01_filter_regions/<parameter_folder>/

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "01_filter_regions")

percent_label <- round(opt$per_sample * 100)
version_label <- paste0("cov", opt$cov, "_", percent_label, "pct")

variant_dir <- file.path(step_dir, version_label)
dir.create(variant_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", variant_dir)

# ------------------------------------------------------------
# 5. Load unfiltered BSseq object
# ------------------------------------------------------------

unfilt_bs <- readRDS(opt$input_rds)

# ------------------------------------------------------------
# 6. Filter CpGs
# ------------------------------------------------------------

filtered_rds <- file.path(variant_dir, "Filtered_BSseq.rds")

bs <- filterCpGs(
  bs = unfilt_bs,
  cov = opt$cov,
  perSample = opt$per_sample,
  file = filtered_rds
)

# ------------------------------------------------------------
# 7. Call regions
# ------------------------------------------------------------

regions <- getRegions(
  bs = bs,
  file = file.path(variant_dir, "Regions.txt")
)

# ------------------------------------------------------------
# 8. Plot region and SD statistics
# ------------------------------------------------------------

plotRegionStats(
  regions,
  maxQuantile = 0.99,
  file = file.path(variant_dir, "Region_Plots.pdf")
)

plotSDstats(
  regions,
  maxQuantile = 0.99,
  file = file.path(variant_dir, "SD_Plots.pdf")
)

# ------------------------------------------------------------
# 9. Compute and plot region totals
# ------------------------------------------------------------

regionTotals <- getRegionTotals(
  regions,
  file = file.path(variant_dir, "Region_Totals.txt")
)

plotRegionTotals(
  regionTotals,
  file = file.path(variant_dir, "Region_Totals.pdf")
)

# ------------------------------------------------------------
# 10. Save run parameters for reproducibility
# ------------------------------------------------------------

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("input_rds:", opt$input_rds),
    paste("cov:", opt$cov),
    paste("per_sample:", opt$per_sample),
    paste("output_dir:", variant_dir),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(variant_dir, "run_parameters.txt")
)

message("✓ Script 01 complete: CpG filtering and region calling complete")