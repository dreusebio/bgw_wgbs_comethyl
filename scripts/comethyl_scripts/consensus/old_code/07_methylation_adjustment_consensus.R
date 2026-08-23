#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 07: Adjust Region-Level Methylation Using Selected PC Sets
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   For a single dataset (run separately per dataset):
#     1) Load region-level methylation matrix from script 03/05/05b
#     2) Load PCs from script 06
#     3) Load PC-trait association statistics from script 06
#     4) Run methylation adjustment using one or more PC-selection modes:
#          v1) all PCs
#          v2) exclude PCs significantly associated with protected traits
#          v3) retain only PCs significantly associated with technical traits
#     5) Save adjusted region methylation matrices per mode
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --input_meth     : path to Region_Methylation.rds (script 03/05/05b)
#   --input_pcs      : path to PCs.rds from script 06
#   --pc_trait_stats : path to PC-trait stats TSV from script 06
#
# OPTIONAL INPUTS
#   --protected_traits_file : text file, one protected trait per line
#   --technical_traits_file : text file, one technical trait per line
#   --run_v1                : run all-PC adjustment TRUE/FALSE [default = TRUE]
#   --run_v2                : run protected-aware adjustment TRUE/FALSE [default = FALSE]
#   --run_v3                : run technical-PC-only adjustment TRUE/FALSE [default = FALSE]
#   --v2_cor_method         : bicor or pearson for v2 PC selection [default = bicor]
#   --v2_p_thresh           : p-value threshold for excluding protected PCs [default = 0.05]
#   --v3_p_thresh           : p-value threshold for selecting technical PCs [default = 0.05]
#
# OUTPUTS
#   comethyl_output/consensus/07_methylation_adjustment/<dataset_label>/<cpg_label>/<region_label>/
#       v1_all_pcs/
#           <dataset_label>_Adjusted_Region_Methylation_allPCs.rds
#           <dataset_label>_Used_PCs.txt
#           <dataset_label>_adjustment_log.txt
#       v2_exclude_protected_pcs/   (if --run_v2 TRUE)
#           <dataset_label>_Adjusted_Region_Methylation_excluding_protected_PCs_<method>.rds
#           <dataset_label>_Removed_Protected_PCs.txt
#           <dataset_label>_Used_PCs.txt
#           <dataset_label>_Removed_PCs_Associated_Traits.tsv
#           <dataset_label>_adjustment_log.txt
#       v3_technical_pcs_only/      (if --run_v3 TRUE)
#           <dataset_label>_Adjusted_Region_Methylation_technicalPCs_only.rds
#           <dataset_label>_Used_Technical_PCs.txt
#           <dataset_label>_Technical_PCs_Associated_Traits.tsv
#           <dataset_label>_adjustment_log.txt
#       traits_requested/found/missing_protected.txt (if provided)
#       traits_requested/found/missing_technical.txt (if provided)
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - Run this script once per dataset
#   - dataset_label, cpg_label, and region_label are derived automatically
#     from --input_meth; no manual labelling required
#   - The adjusted RDS outputs feed into script 08 (--dataset*_meth)
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/07_methylation_adjustment_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin4_methSD0p08/Baseline_Methylation_complete.rds \
#     --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/Baseline/cov3_75pct/covMin4_methSD0p08/PCs.rds \
#     --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/Baseline/cov3_75pct/covMin4_methSD0p08/PC_Trait_Correlation_Stats_Bicor.tsv \
#     --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/protected_traits.txt \
#     --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/technical_traits.txt \
#     --run_v1 TRUE \
#     --run_v2 TRUE \
#     --v2_cor_method bicor \
#     --v2_p_thresh 0.05
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/07_methylation_adjustment/
#     <dataset_label>/
#       <cpg_label>/
#         <region_label>/
#           v1_all_pcs/
#             <dataset_label>_Adjusted_Region_Methylation_allPCs.rds
#             <dataset_label>_Used_PCs.txt
#             <dataset_label>_adjustment_log.txt
#           v2_exclude_protected_pcs/
#             <dataset_label>_Adjusted_Region_Methylation_excluding_protected_PCs_<method>.rds
#             <dataset_label>_Removed_Protected_PCs.txt
#             <dataset_label>_Used_PCs.txt
#             <dataset_label>_Removed_PCs_Associated_Traits.tsv
#             <dataset_label>_adjustment_log.txt
#           v3_technical_pcs_only/
#             <dataset_label>_Adjusted_Region_Methylation_technicalPCs_only.rds
#             <dataset_label>_Used_Technical_PCs.txt
#             <dataset_label>_Technical_PCs_Associated_Traits.tsv
#             <dataset_label>_adjustment_log.txt
#           run_parameters.txt
#           sessionInfo.txt
#
# PATH DERIVATION
#   dataset_label, cpg_label, region_label derived from --input_meth path.
#   --input_pcs and --pc_trait_stats should point to the updated
#   06_pc_diagnostics/<dataset_label>/<cpg_label>/<region_label>/ outputs.
# ================================================================

message("Starting Script 07")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(readr)
  library(dplyr)
})

# ----------------------------------------------------------------
# Load helper.R
# ----------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
]
if (length(script_file_arg) == 0) stop("Could not determine script path.")
script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found: ", helper_file)
source(helper_file)

# ----------------------------------------------------------------
# Parse args
# ----------------------------------------------------------------
option_list <- list(
  make_option("--project_root",          type = "character"),
  make_option("--input_meth",            type = "character"),
  make_option("--input_pcs",             type = "character"),
  make_option("--pc_trait_stats",        type = "character"),
  make_option("--protected_traits_file", type = "character", default = NULL),
  make_option("--technical_traits_file", type = "character", default = NULL),
  make_option("--run_v1",                type = "character", default = "TRUE"),
  make_option("--run_v2",                type = "character", default = "FALSE"),
  make_option("--run_v3",                type = "character", default = "FALSE"),
  make_option("--v2_cor_method",         type = "character", default = "bicor"),
  make_option("--v2_p_thresh",           type = "double",    default = 0.05),
  make_option("--v3_p_thresh",           type = "double",    default = 0.05)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------------------------------------------
# Internal helpers (also in helper.R but kept local for safety)
# ----------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  if (is.logical(x)) return(x)
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

# ----------------------------------------------------------------
# Validate
# ----------------------------------------------------------------
if (is.null(opt$project_root))   stop("--project_root is required")
if (is.null(opt$input_meth))     stop("--input_meth is required")
if (is.null(opt$input_pcs))      stop("--input_pcs is required")
if (is.null(opt$pc_trait_stats)) stop("--pc_trait_stats is required")

if (!dir.exists(opt$project_root))   stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$input_meth))    stop("input_meth not found: ", opt$input_meth)
if (!file.exists(opt$input_pcs))     stop("input_pcs not found: ", opt$input_pcs)
if (!file.exists(opt$pc_trait_stats)) stop("pc_trait_stats not found: ", opt$pc_trait_stats)

if (!is.null(opt$protected_traits_file) && !file.exists(opt$protected_traits_file))
  stop("protected_traits_file not found: ", opt$protected_traits_file)
if (!is.null(opt$technical_traits_file) && !file.exists(opt$technical_traits_file))
  stop("technical_traits_file not found: ", opt$technical_traits_file)

run_v1 <- parse_bool(opt$run_v1, "--run_v1")
run_v2 <- parse_bool(opt$run_v2, "--run_v2")
run_v3 <- parse_bool(opt$run_v3, "--run_v3")

if (!run_v1 && !run_v2 && !run_v3)
  stop("At least one of --run_v1, --run_v2, --run_v3 must be TRUE.")

v2_cor_method <- tolower(opt$v2_cor_method)
if (!v2_cor_method %in% c("bicor", "pearson"))
  stop("--v2_cor_method must be bicor or pearson")

# ----------------------------------------------------------------
# Derive dataset_label / cpg_label / region_label from input_meth
# Same logic as script 06
# ----------------------------------------------------------------
input_filename <- basename(opt$input_meth)
input_dir      <- dirname(opt$input_meth)

variant_pattern <- "^(v[0-9]|unadjusted)"
if (grepl(variant_pattern, basename(input_dir))) {
  region_label  <- basename(dirname(input_dir))
  cpg_label     <- basename(dirname(dirname(input_dir)))
  dataset_label <- basename(dirname(dirname(dirname(input_dir))))
} else if (grepl("05b_shared_complete_regions", input_dir)) {
  dataset_label <- sub("_Methylation_complete\\.rds$", "", input_filename)
  dataset_label <- sub("\\.rds$", "", dataset_label)
  region_label  <- basename(input_dir)
  cpg_label     <- basename(dirname(input_dir))
} else {
  region_label  <- basename(input_dir)
  cpg_label     <- basename(dirname(input_dir))
  dataset_label <- basename(dirname(dirname(input_dir)))
}

message("Derived labels:")
message("  dataset_label : ", dataset_label)
message("  cpg_label     : ", cpg_label)
message("  region_label  : ", region_label)

# ----------------------------------------------------------------
# Output directory
# ----------------------------------------------------------------
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "07_methylation_adjustment")
out_dir       <- file.path(step_dir, dataset_label, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", out_dir)

# ----------------------------------------------------------------
# Load inputs
# ----------------------------------------------------------------
meth <- as.matrix(readRDS(opt$input_meth))
PCs  <- as.matrix(readRDS(opt$input_pcs))
pc_trait_stats_raw <- readr::read_tsv(opt$pc_trait_stats, show_col_types = FALSE)

# Validate methylation
if (!is.numeric(meth))            stop("Methylation matrix must be numeric.")
if (nrow(meth) < 1)               stop("Methylation matrix has zero rows.")
if (ncol(meth) < 2)               stop("Methylation matrix must have >= 2 samples.")
if (is.null(colnames(meth)))      stop("Methylation matrix must have sample IDs in colnames.")
if (anyDuplicated(colnames(meth))) stop("Duplicate methylation sample IDs.")

# Validate PCs
if (!is.numeric(PCs))             stop("PC matrix must be numeric.")
if (nrow(PCs) < 2)                stop("PC matrix must have >= 2 samples.")
if (ncol(PCs) < 1)                stop("PC matrix must have >= 1 PC.")
if (is.null(rownames(PCs)))       stop("PC matrix must have sample IDs in rownames.")
if (anyDuplicated(rownames(PCs))) stop("Duplicate PC sample IDs.")

# Align samples
common_samples <- intersect(colnames(meth), rownames(PCs))
if (length(common_samples) == 0)
  stop("No overlapping samples between methylation and PCs.")
meth <- meth[, common_samples, drop = FALSE]
PCs  <- PCs[common_samples,  , drop = FALSE]
message("After alignment: ", nrow(meth), " regions x ", ncol(meth), " samples; ",
        ncol(PCs), " PCs")

# ----------------------------------------------------------------
# Standardize PC-trait stats
# ----------------------------------------------------------------
has_bicor  <- "bicor" %in% colnames(pc_trait_stats_raw)
has_cor    <- "cor"   %in% colnames(pc_trait_stats_raw)

pc_trait_stats_bicor   <- if (has_bicor) .standardize_pc_trait_table(pc_trait_stats_raw, "bicor", "p") else NULL
pc_trait_stats_pearson <- if (has_cor)   .standardize_pc_trait_table(pc_trait_stats_raw, "cor",   "p") else NULL

if (run_v2 && v2_cor_method == "bicor"   && is.null(pc_trait_stats_bicor))
  stop("v2_cor_method = bicor but pc_trait_stats has no 'bicor' column.")
if (run_v2 && v2_cor_method == "pearson" && is.null(pc_trait_stats_pearson))
  stop("v2_cor_method = pearson but pc_trait_stats has no 'cor' column.")

pc_trait_stats_v2     <- if (v2_cor_method == "bicor") pc_trait_stats_bicor else pc_trait_stats_pearson
pc_trait_stats_for_v3 <- if (!is.null(pc_trait_stats_pearson)) pc_trait_stats_pearson else pc_trait_stats_bicor
if (run_v3 && is.null(pc_trait_stats_for_v3))
  stop("Could not find usable PC-trait stats for v3.")

# Restrict to PCs we actually have
pc_names <- colnames(PCs)
if (!is.null(pc_trait_stats_v2)) {
  pc_trait_stats_v2 <- dplyr::filter(pc_trait_stats_v2, PC %in% pc_names)
  if (nrow(pc_trait_stats_v2) == 0) stop("No overlapping PCs in pc_trait_stats for v2.")
}
if (!is.null(pc_trait_stats_for_v3)) {
  pc_trait_stats_for_v3 <- dplyr::filter(pc_trait_stats_for_v3, PC %in% pc_names)
  if (nrow(pc_trait_stats_for_v3) == 0) stop("No overlapping PCs in pc_trait_stats for v3.")
}

# ----------------------------------------------------------------
# Protected / technical traits
# ----------------------------------------------------------------
all_traits <- sort(unique(c(
  if (!is.null(pc_trait_stats_v2))    pc_trait_stats_v2$trait    else character(0),
  if (!is.null(pc_trait_stats_for_v3)) pc_trait_stats_for_v3$trait else character(0)
)))

protected_requested <- readTraitFile(opt$protected_traits_file, verbose = TRUE)
technical_requested <- readTraitFile(opt$technical_traits_file, verbose = TRUE)

protected_resolved <- resolveTraits(protected_requested, all_traits,
                                    label = "protected_traits", verbose = TRUE)
technical_resolved <- resolveTraits(technical_requested, all_traits,
                                    label = "technical_traits", verbose = TRUE)

if (!is.null(opt$protected_traits_file)) {
  write_vector_file(protected_resolved$requested, file.path(out_dir, "traits_requested_protected.txt"))
  write_vector_file(protected_resolved$found,     file.path(out_dir, "traits_found_protected.txt"))
  write_vector_file(protected_resolved$missing,   file.path(out_dir, "traits_missing_protected.txt"))
}
if (!is.null(opt$technical_traits_file)) {
  write_vector_file(technical_resolved$requested, file.path(out_dir, "traits_requested_technical.txt"))
  write_vector_file(technical_resolved$found,     file.path(out_dir, "traits_found_technical.txt"))
  write_vector_file(technical_resolved$missing,   file.path(out_dir, "traits_missing_technical.txt"))
}

# ----------------------------------------------------------------
# v1: all PCs
# ----------------------------------------------------------------
if (run_v1) {
  dir_v1 <- file.path(out_dir, "v1_all_pcs")
  dir.create(dir_v1, recursive = TRUE, showWarnings = FALSE)

  v1_out <- file.path(dir_v1, paste0(dataset_label, "_Adjusted_Region_Methylation_allPCs.rds"))
  adjustRegionMeth(meth, PCs = PCs, file = v1_out)
  write_vector_file(colnames(PCs), file.path(dir_v1, paste0(dataset_label, "_Used_PCs.txt")))
  write_log_lines(c(
    "mode: v1_all_pcs", "status: completed",
    paste("dataset_label:", dataset_label),
    paste("n_samples:", ncol(meth)), paste("n_regions:", nrow(meth)),
    paste("n_pcs_used:", ncol(PCs)),
    paste("pcs_used:", paste(colnames(PCs), collapse = ", ")),
    paste("output_file:", v1_out)
  ), file.path(dir_v1, paste0(dataset_label, "_adjustment_log.txt")))
  message("v1 complete [", dataset_label, "]: used all ", ncol(PCs), " PCs")
} else {
  message("Skipping v1 (run_v1 = FALSE)")
}

# ----------------------------------------------------------------
# v2: exclude PCs associated with protected traits
# ----------------------------------------------------------------
if (run_v2) {
  dir_v2 <- file.path(out_dir, "v2_exclude_protected_pcs")
  dir.create(dir_v2, recursive = TRUE, showWarnings = FALSE)
  log_v2 <- file.path(dir_v2, paste0(dataset_label, "_adjustment_log.txt"))

  if (length(protected_resolved$found) == 0) {
    warning("Skipping v2: no protected traits found.")
    write_log_lines(c("mode: v2_exclude_protected_pcs", "status: skipped",
                       paste("dataset_label:", dataset_label),
                       "reason: no protected traits found"), log_v2)
  } else {
    sig_protected <- dplyr::filter(pc_trait_stats_v2,
                                   trait %in% protected_resolved$found,
                                   !is.na(p), p <= opt$v2_p_thresh)
    pcs_to_remove <- unique(sig_protected$PC)
    PCs_v2        <- PCs[, !colnames(PCs) %in% pcs_to_remove, drop = FALSE]

    if (ncol(PCs_v2) == 0) {
      warning("Skipping v2: all PCs removed after protected-trait filtering.")
      write_log_lines(c("mode: v2_exclude_protected_pcs", "status: skipped",
                         paste("dataset_label:", dataset_label),
                         "reason: all PCs removed",
                         paste("n_pcs_removed:", length(pcs_to_remove))), log_v2)
    } else {
      v2_out <- file.path(dir_v2, paste0(dataset_label,
        "_Adjusted_Region_Methylation_excluding_protected_PCs_", v2_cor_method, ".rds"))
      adjustRegionMeth(meth, PCs = PCs_v2, file = v2_out)
      write_vector_file(pcs_to_remove, file.path(dir_v2, paste0(dataset_label, "_Removed_Protected_PCs.txt")))
      write_vector_file(colnames(PCs_v2), file.path(dir_v2, paste0(dataset_label, "_Used_PCs.txt")))
      readr::write_tsv(if (nrow(sig_protected) > 0) dplyr::arrange(sig_protected, PC, p)
                       else data.frame(PC = character(0), trait = character(0),
                                       effect = numeric(0), p = numeric(0)),
                       file.path(dir_v2, paste0(dataset_label, "_Removed_PCs_Associated_Traits.tsv")))
      write_log_lines(c(
        "mode: v2_exclude_protected_pcs", "status: completed",
        paste("dataset_label:", dataset_label),
        paste("v2_cor_method:", v2_cor_method),
        paste("v2_p_thresh:", opt$v2_p_thresh),
        paste("n_samples:", ncol(meth)), paste("n_regions:", nrow(meth)),
        paste("n_protected_traits_found:", length(protected_resolved$found)),
        paste("protected_traits_found:", paste(protected_resolved$found, collapse = ", ")),
        paste("n_pcs_total:", ncol(PCs)),
        paste("n_pcs_removed:", length(pcs_to_remove)),
        paste("pcs_removed:", paste(pcs_to_remove, collapse = ", ")),
        paste("n_pcs_used:", ncol(PCs_v2)),
        paste("pcs_used:", paste(colnames(PCs_v2), collapse = ", ")),
        paste("output_file:", v2_out)
      ), log_v2)
      message("v2 complete [", dataset_label, "]: removed ", length(pcs_to_remove),
              " PCs; used ", ncol(PCs_v2))
    }
  }
} else {
  message("Skipping v2 (run_v2 = FALSE)")
}

# ----------------------------------------------------------------
# v3: technical PCs only
# ----------------------------------------------------------------
if (run_v3) {
  dir_v3 <- file.path(out_dir, "v3_technical_pcs_only")
  dir.create(dir_v3, recursive = TRUE, showWarnings = FALSE)
  log_v3 <- file.path(dir_v3, paste0(dataset_label, "_adjustment_log.txt"))

  if (length(technical_resolved$found) == 0) {
    warning("Skipping v3: no technical traits found.")
    write_log_lines(c("mode: v3_technical_pcs_only", "status: skipped",
                       paste("dataset_label:", dataset_label),
                       "reason: no technical traits found"), log_v3)
  } else {
    sig_technical    <- dplyr::filter(pc_trait_stats_for_v3,
                                      trait %in% technical_resolved$found,
                                      !is.na(p), p <= opt$v3_p_thresh)
    tech_pcs_to_use  <- intersect(colnames(PCs), unique(sig_technical$PC))

    if (length(tech_pcs_to_use) == 0) {
      warning("Skipping v3: no PCs met the technical-trait threshold.")
      write_log_lines(c("mode: v3_technical_pcs_only", "status: skipped",
                         paste("dataset_label:", dataset_label),
                         paste("v3_p_thresh:", opt$v3_p_thresh),
                         "reason: no PCs selected"), log_v3)
    } else {
      PCs_v3  <- PCs[, tech_pcs_to_use, drop = FALSE]
      v3_out  <- file.path(dir_v3, paste0(dataset_label,
        "_Adjusted_Region_Methylation_technicalPCs_only.rds"))
      adjustRegionMeth(meth, PCs = PCs_v3, file = v3_out)
      write_vector_file(colnames(PCs_v3),
        file.path(dir_v3, paste0(dataset_label, "_Used_Technical_PCs.txt")))
      readr::write_tsv(if (nrow(sig_technical) > 0) dplyr::arrange(sig_technical, PC, p)
                       else data.frame(PC = character(0), trait = character(0),
                                       effect = numeric(0), p = numeric(0)),
                       file.path(dir_v3, paste0(dataset_label,
                                                  "_Technical_PCs_Associated_Traits.tsv")))
      write_log_lines(c(
        "mode: v3_technical_pcs_only", "status: completed",
        paste("dataset_label:", dataset_label),
        paste("v3_p_thresh:", opt$v3_p_thresh),
        paste("n_samples:", ncol(meth)), paste("n_regions:", nrow(meth)),
        paste("n_technical_traits_found:", length(technical_resolved$found)),
        paste("technical_traits_found:", paste(technical_resolved$found, collapse = ", ")),
        paste("n_pcs_total:", ncol(PCs)),
        paste("n_technical_pcs_used:", ncol(PCs_v3)),
        paste("technical_pcs_used:", paste(colnames(PCs_v3), collapse = ", ")),
        paste("output_file:", v3_out)
      ), log_v3)
      message("v3 complete [", dataset_label, "]: used ", ncol(PCs_v3), " technical PCs")
    }
  }
} else {
  message("Skipping v3 (run_v3 = FALSE)")
}

# ----------------------------------------------------------------
# Run parameters
# ----------------------------------------------------------------
writeLines(c(
  paste("project_root:",           opt$project_root),
  paste("input_meth:",             opt$input_meth),
  paste("input_pcs:",              opt$input_pcs),
  paste("pc_trait_stats:",         opt$pc_trait_stats),
  paste("protected_traits_file:",  ifelse(is.null(opt$protected_traits_file), "NULL", opt$protected_traits_file)),
  paste("technical_traits_file:",  ifelse(is.null(opt$technical_traits_file), "NULL", opt$technical_traits_file)),
  paste("dataset_label:",          dataset_label),
  paste("cpg_label:",              cpg_label),
  paste("region_label:",           region_label),
  paste("run_v1:", run_v1), paste("run_v2:", run_v2), paste("run_v3:", run_v3),
  paste("v2_cor_method:",          v2_cor_method),
  paste("v2_p_thresh:",            opt$v2_p_thresh),
  paste("v3_p_thresh:",            opt$v3_p_thresh),
  paste("n_samples:",              ncol(meth)),
  paste("n_regions:",              nrow(meth)),
  paste("n_total_pcs:",            ncol(PCs)),
  paste("n_protected_found:",      length(protected_resolved$found)),
  paste("n_technical_found:",      length(technical_resolved$found)),
  paste("output_dir:",             out_dir),
  paste("date:",                   as.character(Sys.time()))
), con = file.path(out_dir, "run_parameters.txt"))

writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

message("Script 07 complete: methylation adjustment finished")
message("Output: ", out_dir)