#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 05: Adjust Region-Level Methylation Using Selected PC Sets
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load region-level methylation from Script 03
#   - Load PCs from Script 04
#   - Load PC-trait association statistics from Script 04
#   - Run methylation adjustment using one or more PC-selection modes:
#       v1) all PCs
#       v2) exclude PCs associated with protected traits
#       v3) technical PCs only
#   - Save adjusted region methylation outputs and reproducibility logs
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --input_meth     : path to Region_Methylation.rds from Script 03
#   --input_pcs      : path to PCs.rds from Script 04
#   --pc_trait_stats : path to PC-trait association stats from Script 04
#
# OPTIONAL INPUTS
#   --protected_traits_file : text file with one protected trait per line
#   --technical_traits_file : text file with one technical trait per line
#   --run_v1                : TRUE/FALSE, run all-PC adjustment [default = TRUE]
#   --run_v2                : TRUE/FALSE, run protected-trait-aware adjustment [default = FALSE]
#   --run_v3                : TRUE/FALSE, run technical-PC-only adjustment [default = FALSE]
#   --v2_cor_method         : bicor or pearson [default = bicor]
#   --v2_p_thresh           : p-value threshold for excluding protected PCs [default = 0.05]
#   --v3_p_thresh           : p-value threshold for selecting technical PCs [default = 0.05]
#
# OUTPUTS
#   project_root/comethyl_output/05_methylation_adjustment/<cpg_label>/<region_label>/
#       v1_all_pcs/
#         Adjusted_Region_Methylation_allPCs.rds
#         adjustment_log.txt
#
#       v2_exclude_protected_pcs/
#         Adjusted_Region_Methylation_excluding_protected_PCs_<method>.rds
#         Removed_Protected_PCs.txt
#         Used_PCs.txt
#         Removed_PCs_Associated_Traits.tsv
#         adjustment_log.txt
#
#       v3_technical_pcs_only/
#         Adjusted_Region_Methylation_technicalPCs_only.rds
#         Used_Technical_PCs.txt
#         Technical_PCs_Associated_Traits.tsv
#         adjustment_log.txt
#
#       run_parameters.txt
#
# NOTES
#   - v1 runs by default.
#   - v2 runs only if requested and protected traits are available.
#   - v3 runs only if requested and technical traits are available.
#   - No project-specific trait cleaning is done here.
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(readr)
  library(dplyr)
  library(AnnotationHub)
})

# ------------------------------------------------------------
# Load helper functions
# ------------------------------------------------------------
script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) {
  stop("helper.R not found next to this script: ", helper_file)
}
source(helper_file)

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--input_meth", type = "character",
              help = "Path to Region_Methylation.rds from Script 03"),

  make_option("--input_pcs", type = "character",
              help = "Path to PCs.rds from Script 04"),

  make_option("--pc_trait_stats", type = "character",
              help = "Path to PC-trait stats TSV from Script 04"),

  make_option("--protected_traits_file", type = "character", default = NULL,
              help = "Optional text file with one protected trait per line"),

  make_option("--technical_traits_file", type = "character", default = NULL,
              help = "Optional text file with one technical trait per line"),

  make_option("--run_v1", type = "character", default = "TRUE",
              help = "Run v1 all-PC adjustment: TRUE or FALSE [default = TRUE]"),

  make_option("--run_v2", type = "character", default = "FALSE",
              help = "Run v2 protected-trait-aware adjustment: TRUE or FALSE [default = FALSE]"),

  make_option("--run_v3", type = "character", default = "FALSE",
              help = "Run v3 technical-PC-only adjustment: TRUE or FALSE [default = FALSE]"),

  make_option("--v2_cor_method", type = "character", default = "bicor",
              help = "Correlation method for v2 PC exclusion: bicor or pearson [default = bicor]"),

  make_option("--v2_p_thresh", type = "double", default = 0.05,
              help = "P-value threshold for excluding protected PCs in v2 [default = 0.05]"),

  make_option("--v3_p_thresh", type = "double", default = 0.05,
              help = "P-value threshold for selecting technical PCs in v3 [default = 0.05]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  if (is.logical(x)) return(x)
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

write_vector_file <- function(x, file) {
  x <- unique(as.character(x))
  if (length(x) == 0) x <- character(0)
  writeLines(x, con = file)
}

write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$input_meth)) stop("--input_meth is required")
if (is.null(opt$input_pcs)) stop("--input_pcs is required")
if (is.null(opt$pc_trait_stats)) stop("--pc_trait_stats is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$input_meth)) stop("input_meth not found: ", opt$input_meth)
if (!file.exists(opt$input_pcs)) stop("input_pcs not found: ", opt$input_pcs)
if (!file.exists(opt$pc_trait_stats)) stop("pc_trait_stats not found: ", opt$pc_trait_stats)

if (!is.null(opt$protected_traits_file) && !file.exists(opt$protected_traits_file)) {
  stop("protected_traits_file not found: ", opt$protected_traits_file)
}
if (!is.null(opt$technical_traits_file) && !file.exists(opt$technical_traits_file)) {
  stop("technical_traits_file not found: ", opt$technical_traits_file)
}

run_v1 <- parse_bool(opt$run_v1, "--run_v1")
run_v2 <- parse_bool(opt$run_v2, "--run_v2")
run_v3 <- parse_bool(opt$run_v3, "--run_v3")

if (!run_v1 && !run_v2 && !run_v3) {
  stop("At least one of --run_v1, --run_v2, or --run_v3 must be TRUE.")
}

v2_cor_method <- tolower(opt$v2_cor_method)
if (!v2_cor_method %in% c("bicor", "pearson")) {
  stop("--v2_cor_method must be 'bicor' or 'pearson'")
}
if (opt$v2_p_thresh < 0 || opt$v2_p_thresh > 1) {
  stop("--v2_p_thresh must be between 0 and 1")
}
if (opt$v3_p_thresh < 0 || opt$v3_p_thresh > 1) {
  stop("--v3_p_thresh must be between 0 and 1")
}

# ------------------------------------------------------------
# Configure cache
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

# ------------------------------------------------------------
# Derive output directory from methylation lineage
# Expected input:
#   .../03_region_methylation/<cpg_label>/<region_label>/Region_Methylation.rds
# ------------------------------------------------------------
input_dir <- dirname(opt$input_meth)
region_label <- basename(input_dir)
cpg_label <- basename(dirname(input_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "05_methylation_adjustment")
out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------
meth <- readRDS(opt$input_meth)
PCs <- readRDS(opt$input_pcs)
pc_trait_stats_raw <- readr::read_tsv(opt$pc_trait_stats, show_col_types = FALSE)

# Validate methylation
if (!(is.matrix(meth) || is.data.frame(meth))) {
  stop("input_meth must be a matrix-like object.")
}
meth <- as.matrix(meth)
if (!is.numeric(meth)) stop("Region methylation matrix must be numeric.")
if (nrow(meth) < 1) stop("Region methylation matrix has zero rows.")
if (ncol(meth) < 2) stop("Region methylation matrix must have at least 2 samples.")
if (is.null(colnames(meth))) stop("Region methylation matrix must have sample IDs in colnames.")
if (anyDuplicated(colnames(meth))) {
  dup_ids <- unique(colnames(meth)[duplicated(colnames(meth))])
  stop("Duplicate methylation sample IDs found. Example: ", paste(head(dup_ids, 10), collapse = ", "))
}

# Validate PCs
if (!(is.matrix(PCs) || is.data.frame(PCs))) {
  stop("input_pcs must be a matrix-like object.")
}
PCs <- as.matrix(PCs)
if (!is.numeric(PCs)) stop("PC matrix must be numeric.")
if (nrow(PCs) < 2) stop("PC matrix must have at least 2 samples.")
if (ncol(PCs) < 1) stop("PC matrix must have at least 1 PC.")
if (is.null(rownames(PCs))) stop("PC matrix must have sample IDs in rownames.")
if (is.null(colnames(PCs))) stop("PC matrix must have PC names in colnames.")
if (anyDuplicated(rownames(PCs))) {
  dup_ids <- unique(rownames(PCs)[duplicated(rownames(PCs))])
  stop("Duplicate PC sample IDs found. Example: ", paste(head(dup_ids, 10), collapse = ", "))
}
if (anyDuplicated(colnames(PCs))) {
  dup_pcs <- unique(colnames(PCs)[duplicated(colnames(PCs))])
  stop("Duplicate PC names found. Example: ", paste(head(dup_pcs, 10), collapse = ", "))
}

# Align methylation and PCs
common_samples <- intersect(colnames(meth), rownames(PCs))
if (length(common_samples) == 0) {
  stop("No overlapping samples between region methylation and PCs.")
}
meth <- meth[, common_samples, drop = FALSE]
PCs <- PCs[common_samples, , drop = FALSE]

message("After alignment: meth = ", nrow(meth), " regions x ", ncol(meth),
        " samples; PCs = ", nrow(PCs), " samples x ", ncol(PCs), " PCs")

# ------------------------------------------------------------
# Standardize PC-trait table
# ------------------------------------------------------------
pc_trait_stats_bicor <- NULL
pc_trait_stats_pearson <- NULL

available_cols <- colnames(pc_trait_stats_raw)
has_bicor <- "bicor" %in% available_cols
has_cor <- "cor" %in% available_cols

if (has_bicor) {
  pc_trait_stats_bicor <- .standardize_pc_trait_table(
    pc_trait_stats = pc_trait_stats_raw,
    cor_column = "bicor",
    p_column = "p"
  )
}
if (has_cor) {
  pc_trait_stats_pearson <- .standardize_pc_trait_table(
    pc_trait_stats = pc_trait_stats_raw,
    cor_column = "cor",
    p_column = "p"
  )
}

if (v2_cor_method == "bicor" && is.null(pc_trait_stats_bicor)) {
  stop("Requested v2_cor_method = bicor, but pc_trait_stats file does not contain a 'bicor' column.")
}
if (v2_cor_method == "pearson" && is.null(pc_trait_stats_pearson)) {
  stop("Requested v2_cor_method = pearson, but pc_trait_stats file does not contain a 'cor' column.")
}

pc_trait_stats_v2 <- if (v2_cor_method == "bicor") pc_trait_stats_bicor else pc_trait_stats_pearson
pc_trait_stats_for_v3 <- if (!is.null(pc_trait_stats_pearson)) pc_trait_stats_pearson else pc_trait_stats_bicor

if (is.null(pc_trait_stats_for_v3)) {
  stop("Could not identify usable PC-trait stats for v3.")
}

# Restrict stats to PCs we actually have
pc_names <- colnames(PCs)
pc_trait_stats_v2 <- pc_trait_stats_v2 %>% dplyr::filter(PC %in% pc_names)

if (nrow(pc_trait_stats_v2) == 0) {
  stop("No overlapping PCs between pc_trait_stats and input_pcs for v2.")
}

pc_trait_stats_for_v3 <- pc_trait_stats_for_v3 %>% dplyr::filter(PC %in% pc_names)

if (nrow(pc_trait_stats_for_v3) == 0) {
  stop("No overlapping PCs between pc_trait_stats and input_pcs for v3.")
}

# ------------------------------------------------------------
# Load and resolve protected / technical traits
# Trait files are authoritative inputs for v2 and v3
# ------------------------------------------------------------
all_traits_in_stats <- sort(unique(c(pc_trait_stats_v2$trait, pc_trait_stats_for_v3$trait)))

protected_traits_requested <- readTraitFile(opt$protected_traits_file, verbose = TRUE)
technical_traits_requested <- readTraitFile(opt$technical_traits_file, verbose = TRUE)

protected_resolved <- resolveTraits(
  requested_traits = protected_traits_requested,
  available_traits = all_traits_in_stats,
  label = "protected_traits",
  verbose = TRUE
)

technical_resolved <- resolveTraits(
  requested_traits = technical_traits_requested,
  available_traits = all_traits_in_stats,
  label = "technical_traits",
  verbose = TRUE
)

if (!is.null(opt$protected_traits_file)) {
  write_vector_file(protected_resolved$requested, file.path(out_dir, "traits_requested_protected.txt"))
  write_vector_file(protected_resolved$found, file.path(out_dir, "traits_found_protected.txt"))
  write_vector_file(protected_resolved$missing, file.path(out_dir, "traits_missing_protected.txt"))
}

if (!is.null(opt$technical_traits_file)) {
  write_vector_file(technical_resolved$requested, file.path(out_dir, "traits_requested_technical.txt"))
  write_vector_file(technical_resolved$found, file.path(out_dir, "traits_found_technical.txt"))
  write_vector_file(technical_resolved$missing, file.path(out_dir, "traits_missing_technical.txt"))
}

# ------------------------------------------------------------
# v1: all PCs
# ------------------------------------------------------------
if (run_v1) {
  dir_v1 <- file.path(out_dir, "v1_all_pcs")
  dir.create(dir_v1, recursive = TRUE, showWarnings = FALSE)

  v1_out <- file.path(dir_v1, "Adjusted_Region_Methylation_allPCs.rds")

  methAdj_all <- adjustRegionMeth(
    meth,
    PCs = PCs,
    file = v1_out
  )

  write_log_lines(
    c(
      "mode: v1_all_pcs",
      paste("n_samples:", ncol(meth)),
      paste("n_regions:", nrow(meth)),
      paste("n_pcs_used:", ncol(PCs)),
      paste("pcs_used:", paste(colnames(PCs), collapse = ", ")),
      paste("output_file:", v1_out)
    ),
    file.path(dir_v1, "adjustment_log.txt")
  )

  message("✓ v1 complete: adjusted using all PCs (n = ", ncol(PCs), ")")
} else {
  message("Skipping v1 (run_v1 = FALSE)")
}

# ------------------------------------------------------------
# v2: exclude PCs associated with protected traits
# ------------------------------------------------------------
if (run_v2) {
  dir_v2 <- file.path(out_dir, "v2_exclude_protected_pcs")
  dir.create(dir_v2, recursive = TRUE, showWarnings = FALSE)

  if (length(protected_resolved$found) == 0) {
    warning("Skipping v2: no protected traits were found.")
    write_log_lines(
      c(
        "mode: v2_exclude_protected_pcs",
        "status: skipped",
        "reason: no protected traits found"
      ),
      file.path(dir_v2, "adjustment_log.txt")
    )
  } else {
    sig_protected <- pc_trait_stats_v2 %>%
      dplyr::filter(trait %in% protected_resolved$found, !is.na(p), p <= opt$v2_p_thresh) %>%
      dplyr::arrange(PC, p)

    protected_pcs_to_remove <- unique(sig_protected$PC)
    PCs_v2 <- PCs[, !colnames(PCs) %in% protected_pcs_to_remove, drop = FALSE]

    if (ncol(PCs_v2) == 0) {
      warning("Skipping v2: all PCs were flagged as protected-associated and none remain for adjustment.")
      write_log_lines(
        c(
          "mode: v2_exclude_protected_pcs",
          "status: skipped",
          "reason: all PCs removed after protected-trait filtering",
          paste("n_protected_pcs_removed:", length(protected_pcs_to_remove))
        ),
        file.path(dir_v2, "adjustment_log.txt")
      )
    } else {
      v2_out <- file.path(
        dir_v2,
        paste0("Adjusted_Region_Methylation_excluding_protected_PCs_", v2_cor_method, ".rds")
      )

      methAdj_v2 <- adjustRegionMeth(
        meth,
        PCs = PCs_v2,
        file = v2_out
      )

      write_vector_file(
        protected_pcs_to_remove,
        file.path(dir_v2, "Removed_Protected_PCs.txt")
      )

      write_vector_file(
        colnames(PCs_v2),
        file.path(dir_v2, "Used_PCs.txt")
      )

      if (nrow(sig_protected) > 0) {
        readr::write_tsv(
          sig_protected %>% dplyr::arrange(PC, p),
          file.path(dir_v2, "Removed_PCs_Associated_Traits.tsv")
        )
      } else {
        readr::write_tsv(
          data.frame(PC = character(0), trait = character(0), effect = numeric(0), p = numeric(0)),
          file.path(dir_v2, "Removed_PCs_Associated_Traits.tsv")
        )
      }

      write_log_lines(
        c(
          "mode: v2_exclude_protected_pcs",
          "status: completed",
          paste("v2_cor_method:", v2_cor_method),
          paste("v2_p_thresh:", opt$v2_p_thresh),
          paste("n_samples:", ncol(meth)),
          paste("n_regions:", nrow(meth)),
          paste("n_protected_traits_found:", length(protected_resolved$found)),
          paste("protected_traits_found:", paste(protected_resolved$found, collapse = ", ")),
          paste("n_pcs_total:", ncol(PCs)),
          paste("n_protected_pcs_removed:", length(protected_pcs_to_remove)),
          paste("protected_pcs_removed:", paste(protected_pcs_to_remove, collapse = ", ")),
          paste("n_pcs_used:", ncol(PCs_v2)),
          paste("pcs_used:", paste(colnames(PCs_v2), collapse = ", ")),
          paste("output_file:", v2_out)
        ),
        file.path(dir_v2, "adjustment_log.txt")
      )

      message("✓ v2 complete: removed ", length(protected_pcs_to_remove),
              " protected-associated PCs; used ", ncol(PCs_v2), " PCs")
    }
  }
} else {
  message("Skipping v2 (run_v2 = FALSE)")
}

# ------------------------------------------------------------
# v3: technical PCs only
# ------------------------------------------------------------
if (run_v3) {
  dir_v3 <- file.path(out_dir, "v3_technical_pcs_only")
  dir.create(dir_v3, recursive = TRUE, showWarnings = FALSE)

  if (length(technical_resolved$found) == 0) {
    warning("Skipping v3: no technical traits were found.")
    write_log_lines(
      c(
        "mode: v3_technical_pcs_only",
        "status: skipped",
        "reason: no technical traits found"
      ),
      file.path(dir_v3, "adjustment_log.txt")
    )
  } else {
    sig_technical <- pc_trait_stats_for_v3 %>%
      dplyr::filter(trait %in% technical_resolved$found, !is.na(p), p <= opt$v3_p_thresh) %>%
      dplyr::arrange(PC, p)

    technical_pcs_to_use <- unique(sig_technical$PC)
    technical_pcs_to_use <- intersect(colnames(PCs), technical_pcs_to_use)

    if (length(technical_pcs_to_use) == 0) {
      warning("Skipping v3: no PCs met the technical-trait association threshold.")
      write_log_lines(
        c(
          "mode: v3_technical_pcs_only",
          "status: skipped",
          paste("v3_p_thresh:", opt$v3_p_thresh),
          "reason: no PCs selected by technical-trait filtering"
        ),
        file.path(dir_v3, "adjustment_log.txt")
      )
    } else {
      PCs_v3 <- PCs[, technical_pcs_to_use, drop = FALSE]
      v3_out <- file.path(dir_v3, "Adjusted_Region_Methylation_technicalPCs_only.rds")

      methAdj_v3 <- adjustRegionMeth(
        meth,
        PCs = PCs_v3,
        file = v3_out
      )

      write_vector_file(
        colnames(PCs_v3),
        file.path(dir_v3, "Used_Technical_PCs.txt")
      )

      if (nrow(sig_technical) > 0) {
        readr::write_tsv(
          sig_technical %>% dplyr::arrange(PC, p),
          file.path(dir_v3, "Technical_PCs_Associated_Traits.tsv")
        )
      } else {
        readr::write_tsv(
          data.frame(PC = character(0), trait = character(0), effect = numeric(0), p = numeric(0)),
          file.path(dir_v3, "Technical_PCs_Associated_Traits.tsv")
        )
      }

      write_log_lines(
        c(
          "mode: v3_technical_pcs_only",
          "status: completed",
          paste("v3_p_thresh:", opt$v3_p_thresh),
          paste("n_samples:", ncol(meth)),
          paste("n_regions:", nrow(meth)),
          paste("n_technical_traits_found:", length(technical_resolved$found)),
          paste("technical_traits_found:", paste(technical_resolved$found, collapse = ", ")),
          paste("n_pcs_total:", ncol(PCs)),
          paste("n_technical_pcs_used:", ncol(PCs_v3)),
          paste("technical_pcs_used:", paste(colnames(PCs_v3), collapse = ", ")),
          paste("output_file:", v3_out)
        ),
        file.path(dir_v3, "adjustment_log.txt")
      )

      message("✓ v3 complete: used ", ncol(PCs_v3), " technical-associated PCs")
    }
  }
} else {
  message("Skipping v3 (run_v3 = FALSE)")
}

# ------------------------------------------------------------
# Save run parameters
# ------------------------------------------------------------
run_params <- c(
  paste("project_root:", opt$project_root),
  paste("input_meth:", opt$input_meth),
  paste("input_pcs:", opt$input_pcs),
  paste("pc_trait_stats:", opt$pc_trait_stats),
  paste("protected_traits_file:", ifelse(is.null(opt$protected_traits_file), "NULL", opt$protected_traits_file)),
  paste("technical_traits_file:", ifelse(is.null(opt$technical_traits_file), "NULL", opt$technical_traits_file)),
  paste("run_v1:", run_v1),
  paste("run_v2:", run_v2),
  paste("run_v3:", run_v3),
  paste("v2_cor_method:", v2_cor_method),
  paste("v2_p_thresh:", opt$v2_p_thresh),
  paste("v3_p_thresh:", opt$v3_p_thresh),
  paste("cpg_label:", cpg_label),
  paste("region_label:", region_label),
  paste("n_samples_after_alignment:", ncol(meth)),
  paste("n_regions:", nrow(meth)),
  paste("n_total_pcs:", ncol(PCs)),
  paste("n_protected_traits_found:", length(protected_resolved$found)),
  paste("n_technical_traits_found:", length(technical_resolved$found)),
  paste("date:", as.character(Sys.time()))
)

writeLines(run_params, con = file.path(out_dir, "run_parameters.txt"))

message("✓ Script 05 complete: methylation adjustment finished")