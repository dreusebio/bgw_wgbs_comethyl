#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 00: Generate Unfiltered BSseq Object from Cytosine Reports
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   Reads Bismark CpG_report files and constructs an unfiltered
#   BSseq object containing all CpG sites across all samples.
#   This is the starting object for downstream filtering and
#   comethylation analysis.
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --meta_file    : Excel file containing sample metadata
#   --cyto_dir     : directory containing Bismark CpG_report files
#
# OPTIONAL INPUTS
#   --workers      : number of CPU cores for parallel processing
#                    [default = 4]
#
# OUTPUTS
#   project_root/comethyl_output/00_cpg_extraction/
#       Unfiltered_BSseq.rds
#       CpG_Totals.txt
#       CpG_Totals.pdf
#       run_parameters.txt
#
# NOTES
#   - Metadata row names should correspond to sample identifiers
#     expected by comethyl::getCpGs().
#   - getCpGs() reads files from the current working directory,
#     so this script temporarily switches into --cyto_dir only
#     for that function call.
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(comethyl)
  library(BiocParallel)
  library(AnnotationHub)
})

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--meta_file", type = "character",
              help = "Metadata Excel file"),
  make_option("--cyto_dir", type = "character",
              help = "Directory containing CpG_report files"),
  make_option("--workers", type = "integer", default = 4,
              help = "Number of parallel workers [default = 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$meta_file)) stop("--meta_file is required")
if (is.null(opt$cyto_dir)) stop("--cyto_dir is required")

if (!dir.exists(opt$project_root)) {
  stop("Project root does not exist: ", opt$project_root)
}
if (!file.exists(opt$meta_file)) {
  stop("Metadata file not found: ", opt$meta_file)
}
if (!dir.exists(opt$cyto_dir)) {
  stop("Cytosine report directory not found: ", opt$cyto_dir)
}
if (opt$workers < 1) {
  stop("--workers must be >= 1")
}

# ------------------------------------------------------------
# 3. Define output directories automatically
# ------------------------------------------------------------
# All outputs for this step go under:
#   <project_root>/comethyl_output/00_cpg_extraction/

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "00_cpg_extraction")

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", step_dir)

# ------------------------------------------------------------
# 4. Configure AnnotationHub cache
# ------------------------------------------------------------
# Cache is stored inside the project root so the workflow remains
# portable within the project.

AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

# ------------------------------------------------------------
# 5. Configure parallel processing
# ------------------------------------------------------------

param <- BiocParallel::MulticoreParam(workers = opt$workers)

# ------------------------------------------------------------
# 6. Load metadata
# ------------------------------------------------------------

data <- openxlsx::read.xlsx(
  xlsxFile = opt$meta_file,
  rowNames = TRUE
)

if (nrow(data) == 0) {
  stop("Metadata file appears empty or incorrectly formatted.")
}
# ------------------------------------------------------------
# 7. Construct unfiltered BSseq object
# ------------------------------------------------------------
# getCpGs() expects to read from the working directory, so we
# temporarily switch to the cytosine report folder.

out_unfilt <- file.path(step_dir, "Unfiltered_BSseq.rds")

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(opt$cyto_dir)

unfilt_bs <- getCpGs(
  colData = data,
  file = out_unfilt,
  BPPARAM = param
)

# ------------------------------------------------------------
# 8. Compute and plot CpG totals
# ------------------------------------------------------------

CpGtotals <- getCpGtotals(
  unfilt_bs,
  file = file.path(step_dir, "CpG_Totals.txt")
)

plotCpGtotals(
  CpGtotals,
  file = file.path(step_dir, "CpG_Totals.pdf")
)

# ------------------------------------------------------------
# 9. Save run parameters for reproducibility
# ------------------------------------------------------------

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("meta_file:", opt$meta_file),
    paste("cyto_dir:", opt$cyto_dir),
    paste("workers:", opt$workers),
    paste("output_dir:", step_dir),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(step_dir, "run_parameters.txt")
)

message("✓ Script 00 complete: unfiltered BSseq object created")