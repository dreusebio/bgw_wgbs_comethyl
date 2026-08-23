#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 00: Import CpG Reports for Consensus Workflow
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Reads Bismark CpG_report files for one dataset/timepoint and
#   constructs an unfiltered BSseq object containing all CpG sites
#   across all samples. This is the starting object for the
#   consensus workflow.
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --dataset_label: dataset/timepoint label
#                    e.g. Baseline, TP36_38weeks, Postpartum
#   --meta_file    : Excel file containing sample metadata
#   --cyto_dir     : directory containing Bismark CpG_report files
#
# OPTIONAL INPUTS
#   --workers      : number of CPU cores for parallel processing
#                    [default = 4]
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/00_cpg_extraction/<dataset_label>/
#       Unfiltered_BSseq.rds
#       CpG_Totals.txt
#       CpG_Totals.pdf
#       detected_cpg_report_files.txt
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - Metadata row names should correspond to sample identifiers
#     expected by comethyl::getCpGs().
#   - getCpGs() reads files from the current working directory,
#     so this script temporarily switches into --cyto_dir only
#     for that function call.
#   - This script is intended to be run separately for each dataset
#     used in the consensus workflow.
# Example

# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --dataset_label Baseline \
#   --meta_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/Baseline/cytosine_reports \
#   --workers 8

# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --dataset_label 36-38wk \
#   --meta_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/36-38wk/cytosine_reports \
#   --workers 8

# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --dataset_label Postpartum \
#   --meta_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/Postpartum/cytosine_reports \
#   --workers 8

message("Starting Script 00")
# ================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(comethyl)
  library(BiocParallel)
})

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--dataset_label", type = "character",
              help = "Dataset/timepoint label (e.g. Baseline, TP36_38weeks, Postpartum)"),
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
if (is.null(opt$dataset_label)) stop("--dataset_label is required")
if (is.null(opt$meta_file)) stop("--meta_file is required")
if (is.null(opt$cyto_dir)) stop("--cyto_dir is required")

opt$dataset_label <- trimws(opt$dataset_label)

if (!nzchar(opt$dataset_label)) {
  stop("--dataset_label cannot be empty")
}
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
# 3. Define output directories
# ------------------------------------------------------------

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "00_cpg_extraction")
dataset_dir <- file.path(step_dir, opt$dataset_label)

dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", dataset_dir)

# ------------------------------------------------------------
# 4. Configure parallel processing
# ------------------------------------------------------------

param <- BiocParallel::MulticoreParam(workers = opt$workers)

# ------------------------------------------------------------
# 5. Detect CpG report files
# ------------------------------------------------------------

cpg_files <- list.files(
  opt$cyto_dir,
  pattern = "CpG_report\\.txt\\.gz$",
  full.names = FALSE
)

if (length(cpg_files) == 0) {
  stop("No CpG_report.txt.gz files found in: ", opt$cyto_dir)
}

writeLines(
  sort(cpg_files),
  con = file.path(dataset_dir, "detected_cpg_report_files.txt")
)

message("Detected ", length(cpg_files), " CpG report files")

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

if (is.null(rownames(data)) || any(rownames(data) == "")) {
  stop("Metadata row names must contain sample IDs.")
}

if (anyDuplicated(rownames(data))) {
  dup_ids <- unique(rownames(data)[duplicated(rownames(data))])
  stop("Metadata row names contain duplicates. Example: ",
       paste(head(dup_ids, 10), collapse = ", "))
}

message("Loaded metadata for ", nrow(data), " samples")

# ------------------------------------------------------------
# 7. Construct unfiltered BSseq object
# ------------------------------------------------------------

out_unfilt <- file.path(dataset_dir, "Unfiltered_BSseq.rds")

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(opt$cyto_dir)

unfilt_bs <- getCpGs(
  colData = data,
  file = out_unfilt,
  BPPARAM = param
)

if (!inherits(unfilt_bs, "BSseq")) {
  stop("getCpGs() did not return a BSseq object.")
}

message("Unfiltered BSseq object created successfully")

# ------------------------------------------------------------
# 8. Compute and plot CpG totals
# ------------------------------------------------------------

CpGtotals <- getCpGtotals(
  unfilt_bs,
  file = file.path(dataset_dir, "CpG_Totals.txt")
)

plotCpGtotals(
  CpGtotals,
  file = file.path(dataset_dir, "CpG_Totals.pdf")
)

# ------------------------------------------------------------
# 9. Save reproducibility information
# ------------------------------------------------------------

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("dataset_label:", opt$dataset_label),
    paste("meta_file:", opt$meta_file),
    paste("cyto_dir:", opt$cyto_dir),
    paste("workers:", opt$workers),
    paste("n_metadata_samples:", nrow(data)),
    paste("n_detected_cpg_report_files:", length(cpg_files)),
    paste("output_dir:", dataset_dir),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(dataset_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(dataset_dir, "sessionInfo.txt")
)

message("Consensus Script 00 complete for dataset: ", opt$dataset_label)