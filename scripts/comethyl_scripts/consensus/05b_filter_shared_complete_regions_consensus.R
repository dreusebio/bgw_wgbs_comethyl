#!/usr/bin/env Rscript

# ============================================================================
# SCRIPT 05b: Build a Shared, Complete and Optionally Joint-Dynamic Region Set
# ============================================================================
#
# WHY THIS SCRIPT EXISTS
#   Consensus WGCNA requires every dataset to use the same genomic regions in
#   the same order. PCA/SVA and WGCNA also require a complete numeric matrix.
#
#   The historical pipeline selected regions using coverage and methylation SD
#   in one reference dataset, then projected those regions into the other
#   datasets. That mode remains available as `complete_only` for backward
#   compatibility.
#
#   For comparisons such as females versus males, regions dynamic in the
#   reference group may be nearly invariant in the other group. This script
#   therefore adds an OPTIONAL joint-variability filter, calculated only after
#   region methylation has been constructed for every dataset.
#
# RECOMMENDED SEX-CONSENSUS DESIGN
#   1. Script 02: use a defensible region coverage threshold, but set
#      --meth_sd 0 (or a deliberately permissive value). This prevents the
#      reference sex from removing regions that may be dynamic in the other sex.
#   2. Scripts 03/05: calculate methylation for the same region definitions.
#   3. This script: use --selection_mode joint_sd_all and --joint_meth_sd 0.07.
#   4. Use the matrices and region table written here in Scripts 06 onward.
#
# IMPORTANT INTERPRETATION
#   `joint_sd_all` does NOT force male and female correlation structures to be
#   identical. It only ensures every retained region has enough variability to
#   contribute meaningful correlations in every dataset.
#
# REQUIRED ARGUMENTS
#   --project_root
#   --dataset1_label --dataset1_meth
#   --dataset2_label --dataset2_meth
#   --regions_file        Region table used to build the input matrices
#
# OPTIONAL ARGUMENTS
#   --dataset3_label --dataset3_meth
#   --selection_mode      complete_only | joint_sd_all | joint_sd_any |
#                         joint_sd_min_n [default: complete_only]
#   --joint_meth_sd       within-dataset SD threshold [default: 0]
#   --min_datasets_sd     used only by joint_sd_min_n [default: 2]
#   --min_regions         stop if fewer regions survive [default: 1000]
#   --write_sd_table      TRUE/FALSE; write per-region SD diagnostics
#                         [default: TRUE]
#
# SELECTION MODES
#   complete_only : historical behavior; present and NA-free in every dataset.
#   joint_sd_all  : complete and SD >= threshold in every dataset. Recommended
#                   when the aim is a shared male/female consensus network.
#   joint_sd_any  : complete and SD >= threshold in at least one dataset.
#                   Useful for discovery, but may retain invariant regions in
#                   another group and is less suitable for strict consensus.
#   joint_sd_min_n: complete and passes SD in at least --min_datasets_sd sets.
#
# OUTPUT
#   .../05b_shared_complete_regions/<cpg_label>/<effective_region_label>/
#
#   The effective label includes the joint-selection rule, preventing an old
#   complete-only run from being silently overwritten by a joint-SD run.
# ============================================================================

message("Starting Script 05b")

suppressPackageStartupMessages({
  library(optparse)
})

# ----------------------------------------------------------------------------
# 1. Command-line arguments
# ----------------------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_meth", type = "character"),
  make_option("--dataset2_label", type = "character"),
  make_option("--dataset2_meth", type = "character"),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_meth", type = "character", default = NULL),
  make_option("--regions_file", type = "character"),
  make_option("--selection_mode", type = "character", default = "complete_only"),
  make_option("--joint_meth_sd", type = "double", default = 0),
  make_option("--min_datasets_sd", type = "integer", default = 2),
  make_option("--min_regions", type = "integer", default = 1000),
  make_option("--write_sd_table", type = "character", default = "TRUE")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------------------------------------------------------
# 2. Small reusable validation helpers
# ----------------------------------------------------------------------------
parse_bool <- function(x, name) {
  value <- tolower(trimws(as.character(x)))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(name, " must be TRUE or FALSE. Received: ", x)
}

safe_label <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Use stable labels such as 0p07 instead of machine-dependent scientific text.
number_label <- function(x) {
  if (identical(as.numeric(x), 0)) return("0")
  txt <- format(x, scientific = FALSE, trim = TRUE, digits = 12)
  txt <- sub("0+$", "", txt)
  txt <- sub("\\.$", "", txt)
  gsub("\\.", "p", txt)
}

required <- c("project_root", "dataset1_label", "dataset1_meth",
              "dataset2_label", "dataset2_meth", "regions_file")
for (arg in required) {
  if (is.null(opt[[arg]]) || !nzchar(opt[[arg]])) stop("--", arg, " is required")
}

if (!dir.exists(opt$project_root)) stop("Project root not found: ", opt$project_root)
for (path in c(opt$dataset1_meth, opt$dataset2_meth, opt$regions_file)) {
  if (!file.exists(path)) stop("Required input not found: ", path)
}

dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)
if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth)) {
    stop("Provide both --dataset3_label and --dataset3_meth.")
  }
  if (!file.exists(opt$dataset3_meth)) stop("Dataset 3 file not found: ", opt$dataset3_meth)
}

selection_mode <- tolower(opt$selection_mode)
valid_modes <- c("complete_only", "joint_sd_all", "joint_sd_any", "joint_sd_min_n")
if (!selection_mode %in% valid_modes) {
  stop("--selection_mode must be one of: ", paste(valid_modes, collapse = ", "))
}
if (!is.finite(opt$joint_meth_sd) || opt$joint_meth_sd < 0) {
  stop("--joint_meth_sd must be a finite number >= 0")
}
if (opt$min_regions < 1) stop("--min_regions must be >= 1")
write_sd_table <- parse_bool(opt$write_sd_table, "--write_sd_table")

# ----------------------------------------------------------------------------
# 3. Derive lineage labels and a collision-safe effective region label
# ----------------------------------------------------------------------------
# Expected regions_file lineage:
#   .../<cpg_label>/<reference_region_label>/Filtered_Regions.txt
rfile_dir <- dirname(opt$regions_file)
reference_region_label <- basename(rfile_dir)
cpg_label <- basename(dirname(rfile_dir))

if (selection_mode == "complete_only") {
  effective_region_label <- paste0(reference_region_label, "_completeOnly")
} else {
  mode_tag <- switch(
    selection_mode,
    joint_sd_all = "jointSDall",
    joint_sd_any = "jointSDany",
    joint_sd_min_n = paste0("jointSDmin", opt$min_datasets_sd)
  )
  effective_region_label <- paste0(
    reference_region_label, "_", mode_tag, number_label(opt$joint_meth_sd)
  )
}

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
out_dir <- file.path(pipeline_root, "05b_shared_complete_regions",
                     cpg_label, effective_region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Derived cpg_label             : ", cpg_label)
message("Reference region label        : ", reference_region_label)
message("Effective downstream label    : ", effective_region_label)
message("Selection mode                : ", selection_mode)
message("Output directory              : ", out_dir)

# ----------------------------------------------------------------------------
# 4. Load and validate matrices
# ----------------------------------------------------------------------------
load_meth <- function(path, label) {
  message("Loading ", label, ": ", path)
  x <- as.matrix(readRDS(path))
  if (!is.numeric(x)) stop(label, " matrix is not numeric")
  if (nrow(x) < 1 || ncol(x) < 2) stop(label, " matrix has invalid dimensions")
  if (is.null(rownames(x))) stop(label, " matrix needs RegionIDs as rownames")
  if (is.null(colnames(x))) stop(label, " matrix needs sample IDs as colnames")
  if (anyDuplicated(rownames(x))) stop(label, " has duplicate RegionIDs")
  if (anyDuplicated(colnames(x))) stop(label, " has duplicate sample IDs")
  message("  Loaded: ", nrow(x), " regions x ", ncol(x), " samples")
  x
}

input_specs <- list(
  list(label = opt$dataset1_label, file = opt$dataset1_meth),
  list(label = opt$dataset2_label, file = opt$dataset2_meth)
)
if (dataset3_provided) {
  input_specs[[3]] <- list(label = opt$dataset3_label, file = opt$dataset3_meth)
}

dataset_labels <- vapply(input_specs, `[[`, character(1), "label")
if (anyDuplicated(dataset_labels)) stop("Dataset labels must be unique")

meth_list <- lapply(input_specs, function(x) load_meth(x$file, x$label))
names(meth_list) <- dataset_labels
n_datasets <- length(meth_list)

if (selection_mode == "joint_sd_min_n" &&
    (opt$min_datasets_sd < 1 || opt$min_datasets_sd > n_datasets)) {
  stop("--min_datasets_sd must be between 1 and ", n_datasets)
}

# ----------------------------------------------------------------------------
# 5. Restrict to region IDs present in every matrix
# ----------------------------------------------------------------------------
# Preserve the row order of dataset 1. Reduce(intersect, ...) can otherwise
# make ordering less obvious and ordering must be identical downstream.
shared_ids <- rownames(meth_list[[1]])
for (i in seq_along(meth_list)[-1]) {
  shared_ids <- shared_ids[shared_ids %in% rownames(meth_list[[i]])]
}

message("Regions present in every dataset: ", length(shared_ids))
if (length(shared_ids) == 0) stop("No shared RegionIDs across datasets")

shared_meth <- lapply(meth_list, function(x) x[shared_ids, , drop = FALSE])

# ----------------------------------------------------------------------------
# 6. Require finite, complete methylation values in every dataset
# ----------------------------------------------------------------------------
# complete.cases removes NA/NaN. is.finite additionally rejects Inf/-Inf.
complete_matrix <- vapply(shared_meth, function(x) {
  complete.cases(x) & rowSums(!is.finite(x)) == 0
}, logical(length(shared_ids)))

colnames(complete_matrix) <- dataset_labels
complete_all <- rowSums(complete_matrix) == n_datasets
complete_ids <- shared_ids[complete_all]

message("Regions complete in every dataset: ", length(complete_ids))
if (length(complete_ids) == 0) stop("No complete shared regions remain")

complete_meth <- lapply(shared_meth, function(x) x[complete_ids, , drop = FALSE])

# ----------------------------------------------------------------------------
# 7. Calculate within-dataset methylation SD and apply joint rule
# ----------------------------------------------------------------------------
# SD is intentionally calculated BEFORE PC adjustment. It represents observed
# biological/technical variability used for feature eligibility. Adjustment is
# performed later and independently within each dataset.
sd_matrix <- vapply(complete_meth, function(x) {
  apply(x, 1, stats::sd)
}, numeric(length(complete_ids)))

if (is.null(dim(sd_matrix))) {
  sd_matrix <- matrix(sd_matrix, ncol = 1,
                      dimnames = list(complete_ids, dataset_labels))
} else {
  rownames(sd_matrix) <- complete_ids
  colnames(sd_matrix) <- dataset_labels
}

sd_pass <- sd_matrix >= opt$joint_meth_sd
n_datasets_passing <- rowSums(sd_pass)

keep_dynamic <- switch(
  selection_mode,
  complete_only = rep(TRUE, length(complete_ids)),
  joint_sd_all = n_datasets_passing == n_datasets,
  joint_sd_any = n_datasets_passing >= 1,
  joint_sd_min_n = n_datasets_passing >= opt$min_datasets_sd
)

final_ids <- complete_ids[keep_dynamic]
message("Final regions retained: ", length(final_ids))

if (length(final_ids) < opt$min_regions) {
  stop("Only ", length(final_ids), " regions passed. This is below --min_regions=",
       opt$min_regions, ". Review the SD threshold and QC table.")
}

# ----------------------------------------------------------------------------
# 8. Create diagnostic table BEFORE discarding rejected regions
# ----------------------------------------------------------------------------
diagnostic <- data.frame(
  RegionID = complete_ids,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (label in dataset_labels) {
  safe <- safe_label(label)
  diagnostic[[paste0("SD_", safe)]] <- sd_matrix[, label]
  diagnostic[[paste0("PassSD_", safe)]] <- sd_pass[, label]
}
diagnostic$n_datasets_passing_sd <- n_datasets_passing
diagnostic$retained <- keep_dynamic

if (write_sd_table) {
  write.table(diagnostic,
              file = file.path(out_dir, "joint_variability_diagnostics.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Compact dataset-level summary for rapid review.
sd_summary <- do.call(rbind, lapply(dataset_labels, function(label) {
  values <- sd_matrix[, label]
  data.frame(
    Dataset = label,
    N_samples = ncol(complete_meth[[label]]),
    N_complete_shared = length(complete_ids),
    N_passing_SD = sum(values >= opt$joint_meth_sd),
    Proportion_passing_SD = mean(values >= opt$joint_meth_sd),
    Median_SD = stats::median(values),
    Mean_SD = mean(values),
    stringsAsFactors = FALSE
  )
}))

write.table(sd_summary, file.path(out_dir, "joint_variability_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------------
# 9. Save final, identically ordered methylation matrices
# ----------------------------------------------------------------------------
final_meth <- lapply(complete_meth, function(x) x[final_ids, , drop = FALSE])

for (label in dataset_labels) {
  out_file <- file.path(out_dir,
                        paste0(safe_label(label), "_Methylation_jointEligible.rds"))
  saveRDS(final_meth[[label]], out_file)
  message("Saved matrix for ", label, ": ", out_file)
}

saveRDS(final_ids, file.path(out_dir, "joint_eligible_regionIDs.rds"))
writeLines(final_ids, file.path(out_dir, "joint_eligible_regionIDs.txt"))

# ----------------------------------------------------------------------------
# 10. Save matching region annotation tables
# ----------------------------------------------------------------------------
regions_df <- read.delim(opt$regions_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("RegionID", "chr", "start", "end")
missing_cols <- setdiff(needed, colnames(regions_df))
if (length(missing_cols) > 0) {
  stop("regions_file lacks: ", paste(missing_cols, collapse = ", "))
}
if (anyDuplicated(regions_df$RegionID)) stop("regions_file has duplicate RegionIDs")

idx <- match(final_ids, regions_df$RegionID)
if (anyNA(idx)) stop("Some retained RegionIDs were absent from regions_file")
regions_final <- regions_df[idx, , drop = FALSE]

write.table(regions_final[, needed, drop = FALSE],
            file.path(out_dir, "joint_eligible_regions_4col.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(regions_final,
            file.path(out_dir, "joint_eligible_regions_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ----------------------------------------------------------------------------
# 11. Reproducibility record and explicit handoff
# ----------------------------------------------------------------------------
run_params <- c(
  paste("project_root:", opt$project_root),
  paste("regions_file:", opt$regions_file),
  paste("cpg_label:", cpg_label),
  paste("reference_region_label:", reference_region_label),
  paste("effective_region_label:", effective_region_label),
  paste("selection_mode:", selection_mode),
  paste("joint_meth_sd:", opt$joint_meth_sd),
  paste("min_datasets_sd:", opt$min_datasets_sd),
  paste("n_datasets:", n_datasets),
  paste("datasets:", paste(dataset_labels, collapse = ", ")),
  paste("n_shared_ids:", length(shared_ids)),
  paste("n_complete_ids:", length(complete_ids)),
  paste("n_final_ids:", length(final_ids)),
  paste("output_dir:", out_dir),
  paste("date:", as.character(Sys.time()))
)
writeLines(run_params, file.path(out_dir, "run_parameters.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("\nScript 05b complete")
message("Use these files in Script 06:")
for (label in dataset_labels) {
  message("  ", label, ": ", file.path(out_dir,
          paste0(safe_label(label), "_Methylation_jointEligible.rds")))
}
message("Use this region table in Scripts 09-11:")
message("  ", file.path(out_dir, "joint_eligible_regions_4col.tsv"))
