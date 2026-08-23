#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 04: Align Other Datasets to Reference Filtered CpGs
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   For one or two non-reference datasets:
#     1) reads CpG_report files and constructs an unfiltered BSseq object
#     2) restricts CpGs to those overlapping the reference filtered BSseq
#     3) applies a permissive CpG filter
#     4) saves aligned/filtered BSseq objects for downstream region methylation
#
# REQUIRED INPUTS
#   --project_root      : root directory of the analysis project
#   --reference_bs      : path to reference Filtered_BSseq.rds
#
# DATASET 1 INPUTS
#   --dataset1_label    : label for non-reference dataset 1
#   --dataset1_meta     : metadata Excel file for dataset 1
#   --dataset1_cyto_dir : CpG_report directory for dataset 1
#
# OPTIONAL DATASET 2 INPUTS
#   --dataset2_label
#   --dataset2_meta
#   --dataset2_cyto_dir
#
# OPTIONAL INPUTS
#   --cov              : minimum read coverage per CpG after alignment
#                        [default = 1]
#   --per_sample       : fraction of samples meeting coverage threshold
#                        [default = 0.01]
#   --workers          : number of parallel workers [default = 4]
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/04_align_to_reference_cpgs/<dataset_label>/
#       Unfiltered_BSseq.rds
#       Aligned_To_Reference_BSseq.rds
#       Filtered_BSseq.rds
#       CpG_Counts_Summary.txt
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - This script is intended for non-reference datasets only.
#   - It uses the filtered reference BSseq as the canonical CpG universe.
#   - subsetByOverlaps() is used here, which is appropriate for CpG loci
#     when coordinates are consistent across datasets.
# Example 

# Rscript scripts/consensus/04_align_other_datasets_to_reference_cpgs_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
#   --reference_bs /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
#   --dataset1_label 36-38wk \
#   --dataset1_meta /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/data/demo/comethyl_minidata/36-38wk/sample_info.xlsx \
#   --dataset1_cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/data/demo/comethyl_minidata/36-38wk/cytosine_reports \
#   --cov 1 \
#   --per_sample 0.01 \
#   --workers 8

# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/04_align_other_datasets_to_reference_cpgs_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --reference_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
#   --dataset1_label 36-38wk \
#   --dataset1_meta /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
#   --dataset1_cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/36-38wk/cytosine_reports \
#   --dataset2_label Postpartum \
#   --dataset2_meta /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
#   --dataset2_cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/Postpartum/cytosine_reports \
#   --cov 1 \
#   --per_sample 0.01 \
#   --workers 8  

# ================================================================
message("Starting Script 04")
suppressPackageStartupMessages({
  library(optparse)
  library(openxlsx)
  library(comethyl)
  library(BiocParallel)
  library(IRanges)
  library(GenomicRanges)
})

# ------------------------------------------------------------
# 1. Parse command-line arguments
# ------------------------------------------------------------

option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),
  make_option("--reference_bs", type = "character",
              help = "Path to reference Filtered_BSseq.rds"),

  make_option("--dataset1_label", type = "character",
              help = "Dataset 1 label"),
  make_option("--dataset1_meta", type = "character",
              help = "Dataset 1 metadata Excel file"),
  make_option("--dataset1_cyto_dir", type = "character",
              help = "Dataset 1 CpG_report directory"),

  make_option("--dataset2_label", type = "character", default = NULL,
              help = "Dataset 2 label"),
  make_option("--dataset2_meta", type = "character", default = NULL,
              help = "Dataset 2 metadata Excel file"),
  make_option("--dataset2_cyto_dir", type = "character", default = NULL,
              help = "Dataset 2 CpG_report directory"),

  make_option("--cov", type = "integer", default = 1,
              help = "Minimum CpG coverage after alignment [default = 1]"),
  make_option("--per_sample", type = "double", default = 0.01,
              help = "Fraction of samples passing coverage [default = 0.01]"),
  make_option("--workers", type = "integer", default = 4,
              help = "Number of parallel workers [default = 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$reference_bs)) stop("--reference_bs is required")

if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meta)) stop("--dataset1_meta is required")
if (is.null(opt$dataset1_cyto_dir)) stop("--dataset1_cyto_dir is required")

if (!dir.exists(opt$project_root)) {
  stop("Project root does not exist: ", opt$project_root)
}
if (!file.exists(opt$reference_bs)) {
  stop("Reference filtered BSseq not found: ", opt$reference_bs)
}
if (!file.exists(opt$dataset1_meta)) {
  stop("dataset1_meta not found: ", opt$dataset1_meta)
}
if (!dir.exists(opt$dataset1_cyto_dir)) {
  stop("dataset1_cyto_dir not found: ", opt$dataset1_cyto_dir)
}

dataset2_provided <- !is.null(opt$dataset2_label) ||
  !is.null(opt$dataset2_meta) ||
  !is.null(opt$dataset2_cyto_dir)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meta) || is.null(opt$dataset2_cyto_dir)) {
    stop("If dataset 2 is used, provide --dataset2_label, --dataset2_meta, and --dataset2_cyto_dir together.")
  }
  if (!file.exists(opt$dataset2_meta)) {
    stop("dataset2_meta not found: ", opt$dataset2_meta)
  }
  if (!dir.exists(opt$dataset2_cyto_dir)) {
    stop("dataset2_cyto_dir not found: ", opt$dataset2_cyto_dir)
  }
}

if (opt$cov < 1) stop("--cov must be >= 1")
if (opt$per_sample <= 0 || opt$per_sample > 1) stop("--per_sample must be > 0 and <= 1")
if (opt$workers < 1) stop("--workers must be >= 1")

# ------------------------------------------------------------
# 3. Output root
# ------------------------------------------------------------

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "04_align_to_reference_cpgs")
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 4. Load reference filtered BSseq
# ------------------------------------------------------------

reference_bs <- readRDS(opt$reference_bs)

if (!inherits(reference_bs, "BSseq")) {
  stop("reference_bs is not a BSseq object: ", opt$reference_bs)
}

message("Reference filtered BSseq loaded")
message("Reference CpGs: ", nrow(reference_bs))

# ------------------------------------------------------------
# 5. Helper function
# ------------------------------------------------------------

process_dataset <- function(dataset_label, meta_file, cyto_dir, reference_bs, step_dir, workers, cov, per_sample) {
  message("\n====================================================")
  message("Processing dataset: ", dataset_label)
  message("====================================================")

  out_dir <- file.path(step_dir, dataset_label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_unfiltered <- file.path(out_dir, "Unfiltered_BSseq.rds")
  out_aligned    <- file.path(out_dir, "Aligned_To_Reference_BSseq.rds")
  out_filtered   <- file.path(out_dir, "Filtered_BSseq.rds")

  # Metadata
  colData_rep <- openxlsx::read.xlsx(meta_file, rowNames = TRUE)

  if (nrow(colData_rep) == 0) {
    stop("Metadata file appears empty or incorrectly formatted: ", meta_file)
  }
  if (is.null(rownames(colData_rep)) || any(rownames(colData_rep) == "")) {
    stop("Metadata row names must contain sample IDs: ", meta_file)
  }
  if (anyDuplicated(rownames(colData_rep))) {
    dup_ids <- unique(rownames(colData_rep)[duplicated(rownames(colData_rep))])
    stop("Duplicate metadata row names found for ", dataset_label,
         ". Example: ", paste(head(dup_ids, 10), collapse = ", "))
  }

  # Build unfiltered BSseq
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(cyto_dir)

  bs_rep <- getCpGs(
    colData = colData_rep,
    file = out_unfiltered,
    BPPARAM = BiocParallel::MulticoreParam(workers = workers)
  )

  if (!inherits(bs_rep, "BSseq")) {
    stop("getCpGs() did not return a BSseq object for ", dataset_label)
  }

  n_unfiltered <- nrow(bs_rep)

  # Align to reference CpGs using overlaps
  bs_rep_aligned <- IRanges::subsetByOverlaps(bs_rep, ranges = reference_bs)
  saveRDS(bs_rep_aligned, out_aligned)

  n_aligned <- nrow(bs_rep_aligned)

  if (n_aligned == 0) {
    stop("No CpGs remained after aligning ", dataset_label, " to the reference BSseq.")
  }

  # Apply permissive filtering
  bs_rep_filtered <- filterCpGs(
    bs = bs_rep_aligned,
    cov = cov,
    perSample = per_sample,
    file = out_filtered
  )

  if (!inherits(bs_rep_filtered, "BSseq")) {
    stop("filterCpGs() did not return a BSseq object for ", dataset_label)
  }

  n_filtered <- nrow(bs_rep_filtered)

  if (n_filtered == 0) {
    stop("No CpGs remained after filtering aligned BSseq for ", dataset_label)
  }

  # Save summary
  writeLines(
    c(
      paste("dataset_label:", dataset_label),
      paste("meta_file:", meta_file),
      paste("cyto_dir:", cyto_dir),
      paste("cov:", cov),
      paste("per_sample:", per_sample),
      paste("n_unfiltered_cpgs:", n_unfiltered),
      paste("n_aligned_to_reference_cpgs:", n_aligned),
      paste("n_filtered_cpgs:", n_filtered),
      paste("output_dir:", out_dir),
      paste("date:", as.character(Sys.time()))
    ),
    con = file.path(out_dir, "run_parameters.txt")
  )

  writeLines(
    c(
      paste("dataset_label", dataset_label, sep = "\t"),
      paste("n_unfiltered_cpgs", n_unfiltered, sep = "\t"),
      paste("n_aligned_to_reference_cpgs", n_aligned, sep = "\t"),
      paste("n_filtered_cpgs", n_filtered, sep = "\t")
    ),
    con = file.path(out_dir, "CpG_Counts_Summary.txt")
  )

  writeLines(
    capture.output(sessionInfo()),
    con = file.path(out_dir, "sessionInfo.txt")
  )

  message(" Completed dataset: ", dataset_label)
  message("  Unfiltered CpGs: ", n_unfiltered)
  message("  Aligned CpGs:    ", n_aligned)
  message("  Filtered CpGs:   ", n_filtered)
}

# ------------------------------------------------------------
# 6. Run dataset 1
# ------------------------------------------------------------

process_dataset(
  dataset_label = opt$dataset1_label,
  meta_file = opt$dataset1_meta,
  cyto_dir = opt$dataset1_cyto_dir,
  reference_bs = reference_bs,
  step_dir = step_dir,
  workers = opt$workers,
  cov = opt$cov,
  per_sample = opt$per_sample
)

# ------------------------------------------------------------
# 7. Run dataset 2 if provided
# ------------------------------------------------------------

if (dataset2_provided) {
  process_dataset(
    dataset_label = opt$dataset2_label,
    meta_file = opt$dataset2_meta,
    cyto_dir = opt$dataset2_cyto_dir,
    reference_bs = reference_bs,
    step_dir = step_dir,
    workers = opt$workers,
    cov = opt$cov,
    per_sample = opt$per_sample
  )
}

message("Consensus Script 04 complete")