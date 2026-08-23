#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 11: Consensus Module Membership
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Computes region-module membership (correlation of each region
#   with each module eigengene) per dataset using the consensus
#   module assignments from script 09:
#     1) Load Consensus_Modules.rds and extract shared module colors
#     2) For each dataset, load the methylation matrix and align
#        to the consensus region set
#     3) Compute bicor or Pearson correlation between each region
#        and each module eigengene
#     4) Save membership tables per dataset (.txt and .xlsx)
#
# REQUIRED INPUTS
#   --project_root          : root directory of the analysis project
#   --consensus_modules_rds : path to Consensus_Modules.rds from script 09
#   --regions_file          : path to Filtered_Regions.txt from script 02
#   --dataset1_label        : label for dataset 1
#   --dataset1_meth         : path to adjusted methylation RDS for dataset 1
#
# OPTIONAL INPUTS
#   --dataset2_label     : label for dataset 2
#   --dataset2_meth      : path to adjusted methylation RDS for dataset 2
#   --dataset3_label     : label for dataset 3
#   --dataset3_meth      : path to adjusted methylation RDS for dataset 3
#   --adjustment_version : label matching script 09 run [default = unadjusted]
#   --membership_cor     : bicor or pearson [default = bicor]
#   --max_p_outliers     : bicor outlier fraction [default = 0.1]
#
# OUTPUTS
#   comethyl_output/consensus/11_consensus_membership/
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           <dataset_label>_Consensus_Region_Module_Membership.txt
#           <dataset_label>_Consensus_Region_Module_Membership.xlsx
#           run_parameters.txt
#       shared/<cpg_label>/<region_label>/<adjustment_version>/
#           sessionInfo.txt
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - Membership scores are useful for identifying hub regions within
#     each module for downstream network analysis
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/11_membership_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --dataset1_label Baseline \
#     --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline_Adjusted_Region_Methylation_allPCs.rds \
#     --dataset2_label 36_38wks \
#     --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36_38wks/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/36_38wks_Adjusted_Region_Methylation_allPCs.rds \
#     --adjustment_version v1_all_pcs \
#     --membership_cor bicor
# ================================================================
# SCRIPT 11: Consensus Module Membership
#
# PURPOSE
#   - Load Consensus_Modules.rds
#   - Load dataset methylation matrices used in consensus analysis
#   - Compute region-module membership per dataset
#   - Save membership tables for each dataset
# ================================================================

message("Starting Script 11")
suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(openxlsx)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

option_list <- list(
  make_option("--project_root",          type = "character"),
  make_option("--consensus_modules_rds", type = "character"),
  make_option("--regions_file",          type = "character"),
  make_option("--dataset1_label",        type = "character"),
  make_option("--dataset1_meth",         type = "character"),
  make_option("--dataset2_label",        type = "character", default = NULL),
  make_option("--dataset2_meth",         type = "character", default = NULL),
  make_option("--dataset3_label",        type = "character", default = NULL),
  make_option("--dataset3_meth",         type = "character", default = NULL),
  make_option("--adjustment_version",    type = "character", default = "unadjusted"),
  make_option("--membership_cor",        type = "character", default = "bicor"),
  make_option("--max_p_outliers",        type = "double",    default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$project_root))
  stop("--project_root is required")
if (is.null(opt$consensus_modules_rds))
  stop("--consensus_modules_rds is required")
if (is.null(opt$regions_file))
  stop("--regions_file is required")
if (is.null(opt$dataset1_label))
  stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth))
  stop("--dataset1_meth is required")
if (!dir.exists(opt$project_root))
  stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$consensus_modules_rds))
  stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
if (!file.exists(opt$regions_file))
  stop("regions_file not found: ", opt$regions_file)
if (!file.exists(opt$dataset1_meth))
  stop("dataset1_meth not found: ", opt$dataset1_meth)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(label, " must be matrix-like.")
  }
  x <- as.matrix(x)

  if (!is.numeric(x)) {
    stop(label, " must be numeric.")
  }
  if (nrow(x) < 1) {
    stop(label, " has zero rows.")
  }
  if (ncol(x) < 2) {
    stop(label, " must have at least 2 columns.")
  }
  if (is.null(rownames(x))) {
    stop(label, " must have rownames (sample IDs).")
  }
  if (is.null(colnames(x))) {
    stop(label, " must have colnames (region IDs).")
  }
  if (anyDuplicated(rownames(x))) {
    stop(label, " has duplicated rownames.")
  }
  if (anyDuplicated(colnames(x))) {
    stop(label, " has duplicated colnames.")
  }

  x
}

validate_me_matrix <- function(x, label) {
  if (is.null(x)) {
    stop(label, ": eigengene matrix is NULL.")
  }
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(label, ": eigengene matrix must be matrix-like.")
  }

  x <- as.data.frame(x)

  if (nrow(x) < 2) {
    stop(label, ": eigengene matrix must have at least 2 samples.")
  }
  if (ncol(x) < 1) {
    stop(label, ": eigengene matrix has zero columns.")
  }
  if (is.null(rownames(x))) {
    stop(label, ": eigengene matrix must have rownames (sample IDs).")
  }
  if (is.null(colnames(x))) {
    stop(label, ": eigengene matrix must have colnames.")
  }
  if (anyDuplicated(rownames(x))) {
    stop(label, ": eigengene matrix has duplicated rownames.")
  }
  if (anyDuplicated(colnames(x))) {
    stop(label, ": eigengene matrix has duplicated colnames.")
  }

  x
}

write_tsv <- function(df, file) {
  write.table(
    df,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# ------------------------------------------------------------
# Validate args (continued)
# ------------------------------------------------------------
membership_cor <- tolower(opt$membership_cor)
if (!membership_cor %in% c("bicor", "pearson")) {
  stop("--membership_cor must be bicor or pearson")
}

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth)) {
    stop("Provide both dataset2_label and dataset2_meth")
  }
  if (!file.exists(opt$dataset2_meth)) {
    stop("dataset2_meth not found")
  }
}

if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth)) {
    stop("Provide both dataset3_label and dataset3_meth")
  }
  if (!file.exists(opt$dataset3_meth)) {
    stop("dataset3_meth not found")
  }
}

# ------------------------------------------------------------
# Load consensus modules
# ------------------------------------------------------------
consensusMods <- readRDS(opt$consensus_modules_rds)

if (is.null(consensusMods$colors)) {
  stop("Consensus object does not contain $colors")
}
if (is.null(consensusMods$multiMEs)) {
  stop("Consensus object does not contain $multiMEs")
}

all_colors <- consensusMods$colors
if (is.null(names(all_colors))) {
  stop("consensusMods$colors must be a named vector with region IDs as names.")
}

module_colors <- unique(all_colors[all_colors != "grey"])
n_modules <- length(module_colors)

message("Non-grey consensus modules detected: ", n_modules,
        " (", paste(module_colors, collapse = ", "), ")")

if (n_modules == 0) {
  stop("No non-grey modules found in consensus object. Cannot compute module membership.")
}

# ------------------------------------------------------------
# Build dataset list
# ------------------------------------------------------------
dataset_inputs <- list(
  list(label = opt$dataset1_label, meth_file = opt$dataset1_meth)
)

if (dataset2_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <-
    list(label = opt$dataset2_label, meth_file = opt$dataset2_meth)
}

if (dataset3_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <-
    list(label = opt$dataset3_label, meth_file = opt$dataset3_meth)
}

# ------------------------------------------------------------
# Output dirs
# ------------------------------------------------------------

# ----------------------------------------------------------------
# Derive cpg_label / region_label from --regions_file path:
#   .../02_reference_region_filter/<ds>/<cpg_label>/<region_label>/Filtered_Regions.txt
# ----------------------------------------------------------------
{
  rfile_dir    <- dirname(opt$regions_file)
  region_label <- basename(rfile_dir)
  cpg_label    <- basename(dirname(rfile_dir))
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "11_consensus_membership")
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", step_dir)

# ------------------------------------------------------------
# Per-dataset membership
# ------------------------------------------------------------
for (ds in dataset_inputs) {

  ds_label <- ds$label
  ds_dir   <- file.path(step_dir, ds_label, cpg_label, region_label, opt$adjustment_version)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  message("Processing dataset: ", ds_label)

  # ----------------------------------------------------------
  # Load methylation matrix
  # Expect: rows = samples, cols = regions
  # ----------------------------------------------------------
  meth <- validate_meth_matrix(readRDS(ds$meth_file), ds_label)
  expr_df <- as.data.frame(meth)

  message(ds_label, ": ", nrow(expr_df), " samples x ", ncol(expr_df), " regions")

  # ----------------------------------------------------------
  # Check region names against consensus colors
  # ----------------------------------------------------------
  missing_regions <- setdiff(names(all_colors), colnames(expr_df))
  extra_regions   <- setdiff(colnames(expr_df), names(all_colors))

  if (length(missing_regions) > 0) {
    stop(
      ds_label, ": methylation matrix is missing ",
      length(missing_regions), " consensus region(s). Example: ",
      paste(head(missing_regions, 10), collapse = ", ")
    )
  }

  if (!identical(colnames(expr_df), names(all_colors))) {
    expr_df <- expr_df[, names(all_colors), drop = FALSE]
    message(ds_label, ": regions reordered to match consensus module assignment")
  }

  if (length(extra_regions) > 0) {
    message(ds_label, ": methylation matrix contains ",
            length(extra_regions), " extra region(s) not in consensus; ignored")
  }

  # ----------------------------------------------------------
  # Get MEs for this dataset
  # ----------------------------------------------------------
  if (!ds_label %in% names(consensusMods$multiMEs)) {
    stop(
      ds_label, ": not found in consensusMods$multiMEs. Available datasets: ",
      paste(names(consensusMods$multiMEs), collapse = ", ")
    )
  }

  MEs <- validate_me_matrix(consensusMods$multiMEs[[ds_label]]$data, ds_label)

  grey_col <- grep("^MEgrey$", colnames(MEs), value = TRUE)
  if (length(grey_col) > 0) {
    MEs <- MEs[, !colnames(MEs) %in% grey_col, drop = FALSE]
    message(ds_label, ": MEgrey removed (", ncol(MEs), " real MEs used for membership)")
  }

  if (ncol(MEs) == 0) {
    message(ds_label, ": No real MEs after removing grey — skipping")
    writeLines("Skipped: no non-grey eigengenes available.",
               con = file.path(ds_dir, paste0(ds_label, "_membership_SKIPPED.txt")))
    next
  }

  # ----------------------------------------------------------
  # Align samples between methylation matrix and eigengenes
  # ----------------------------------------------------------
  common_samples <- intersect(rownames(expr_df), rownames(MEs))

  if (length(common_samples) < 3) {
    stop(ds_label, ": fewer than 3 overlapping samples between methylation matrix and MEs")
  }

  expr_df <- expr_df[common_samples, , drop = FALSE]
  MEs     <- MEs[common_samples, , drop = FALSE]

  if (!identical(rownames(expr_df), rownames(MEs))) {
    stop(ds_label, ": internal error: sample alignment failed")
  }

  message(ds_label, ": aligned to ", length(common_samples), " shared samples")

  # ----------------------------------------------------------
  # Compute membership correlations
  # expr_df: samples x regions
  # MEs:     samples x modules
  # Output:  regions x modules
  # ----------------------------------------------------------
  if (membership_cor == "bicor") {
    mm <- bicorAndPvalue(
      expr_df,
      MEs,
      use = "pairwise.complete.obs",
      maxPOutliers = opt$max_p_outliers
    )
    membership <- as.data.frame(mm$bicor)
    p_vals     <- as.data.frame(mm$p)
  } else {
    mm <- corAndPvalue(
      expr_df,
      MEs,
      use = "pairwise.complete.obs"
    )
    membership <- as.data.frame(mm$cor)
    p_vals     <- as.data.frame(mm$p)
  }

  # ----------------------------------------------------------
  # Clean column names
  # ----------------------------------------------------------
  colnames(membership) <- gsub("^ME", "", colnames(membership))
  colnames(p_vals)     <- paste0(gsub("^ME", "", colnames(p_vals)), "_p")

  # ----------------------------------------------------------
  # Build output table
  # Region rows should match the columns of expr_df
  # ----------------------------------------------------------
  out_df <- data.frame(
    RegionID = colnames(expr_df),
    Module   = unname(all_colors[colnames(expr_df)]),
    membership,
    p_vals,
    stringsAsFactors = FALSE
  )

  out_df <- out_df[order(out_df$Module == "grey", out_df$Module, out_df$RegionID), , drop = FALSE]

  txt_file  <- file.path(ds_dir, paste0(ds_label, "_Consensus_Region_Module_Membership.txt"))
  xlsx_file <- file.path(ds_dir, paste0(ds_label, "_Consensus_Region_Module_Membership.xlsx"))

  write_tsv(out_df, txt_file)

  openxlsx::write.xlsx(
    out_df,
    file = xlsx_file,
    rowNames = FALSE
  )

  message(ds_label, ": Saved membership table (",
          nrow(out_df), " regions, ", ncol(MEs), " modules)")

  writeLines(
    c(
      paste("project_root:", opt$project_root),
      paste("consensus_modules_rds:", opt$consensus_modules_rds),
      paste("dataset_label:", ds_label),
      paste("dataset_meth:", ds$meth_file),
      paste("adjustment_version:", opt$adjustment_version),
      paste("membership_cor:", membership_cor),
      paste("max_p_outliers:", opt$max_p_outliers),
      paste("n_real_modules:", n_modules),
      paste("module_colors:", paste(module_colors, collapse = ", ")),
      paste("n_regions:", nrow(out_df)),
      paste("n_samples_used:", length(common_samples)),
      paste("date:", as.character(Sys.time()))
    ),
    con = file.path(ds_dir, "run_parameters.txt")
  )
}

# ------------------------------------------------------------
# Session info
# ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), con = file.path(
  step_dir, "shared", cpg_label, region_label, opt$adjustment_version, "sessionInfo.txt"))

message("Script 11 complete: consensus membership finished")