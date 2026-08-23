#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08b: Soft-Threshold Power Benchmark (Subsampled)
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Fast, low-memory companion to Script 08. Instead of running
#   pickSoftThreshold() on the full region x sample matrix (which
#   OOMs on large region sets), this repeatedly draws a random
#   subsample of regions (same region indices across all datasets,
#   so datasets stay comparable) and runs comethyl::getSoftPower()
#   on each subsample. Results across seeds are averaged so you get
#   a stable fit estimate in minutes instead of an hour+ full run.
#
#   Use this to compare candidate covMin/methSD region filters
#   BEFORE committing to a full Script 08 run on the complete
#   region set. It is a screening tool, not a replacement for
#   Script 08 -- once you've picked a region filter, confirm with
#   the full run.
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --dataset1_label : label for dataset 1 (e.g. females)
#   --dataset1_meth  : path to adjusted Region_Methylation.rds for dataset 1
#
# OPTIONAL INPUTS
#   --dataset2_label     : label for dataset 2
#   --dataset2_meth      : path to adjusted methylation RDS for dataset 2
#   --dataset3_label     : label for dataset 3
#   --dataset3_meth      : path to adjusted methylation RDS for dataset 3
#   --adjustment_version : label for this run (e.g. v1_all_pcs) [default = unadjusted]
#   --softpower_cor      : bicor or pearson [default = pearson]
#   --power_min           : minimum soft power to test [default = 1]
#   --power_max           : maximum soft power to test [default = 20]
#   --subsample_size      : number of regions to sample per seed [default = 20000]
#   --n_seeds             : number of random subsamples to run [default = 5]
#   --seed_start           : first random seed used [default = 1001]
#                             (seeds used are seed_start .. seed_start + n_seeds - 1)
#   --block_size          : WGCNA block size [default = 5000]
#   --gc_interval         : garbage collection interval [default = 4999]
#   --threads             : WGCNA threads [default = 4]
#   --scale_free_cutoff   : R^2 cutoff for scale-free topology [default = 0.8]
#
# OUTPUTS
#   comethyl_output/consensus/08b_softpower_benchmark/
#       shared/<cpg_label>/<region_label>/<adjustment_version>/subsample<n>/
#           combined_softpower_benchmark.tsv       (every seed x power x dataset)
#           softpower_benchmark_summary.tsv        (mean/sd across seeds, per power x dataset)
#           SoftPower_Benchmark_Combined.pdf        (R^2 / connectivity vs power, ribbon = seed spread)
#           chosen_power_estimate.txt              (screening estimate only, see NOTES)
#           run_parameters.txt
#           sessionInfo.txt
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/subsample<n>/
#           SoftPower_seed<seed>.tsv
#           Selected_regions_seed<seed>.txt
#
# NOTES
#   - cpg_label and region_label are derived from --dataset1_meth path,
#     exactly as in Script 08.
#   - IMPORTANT ON ORIENTATION: comethyl::adjustRegionMeth() transposes
#     the matrix during adjustment. Its input (from getRegionMeth(),
#     e.g. Script 05b's output) is regions-in-rows / samples-in-columns,
#     but its OUTPUT (the files this script reads) is
#     samples-in-rows / regions-in-columns -- the orientation
#     getSoftPower()/getModules() expect directly. This script
#     subsamples along COLUMNS (regions), not rows (samples).
#   - Regions are subsampled from the INTERSECTION of region IDs
#     (colnames) present in all provided datasets, so the same
#     regions are scored for every dataset within a seed. Sample IDs
#     (rownames) are expected to differ across datasets and are not
#     compared.
#   - chosen_power_estimate.txt uses the same "smallest power where
#     all datasets meet cutoff" logic as Script 08, but computed on
#     the seed-averaged fit. Treat it as a fast screening estimate --
#     confirm with a full Script 08 run before finalizing a power for
#     Script 09 module detection.
#   - Lower --block_size/--gc_interval than Script 08's defaults are
#     used because getSoftPower() is being run repeatedly per call;
#     this keeps peak memory low, which is the point of this script.
#
# EXAMPLE
#   Rscript 08b_softpower_benchmark_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#     --dataset1_label females \
#     --dataset1_meth  .../07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
#     --dataset2_label males \
#     --dataset2_meth  .../07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
#     --adjustment_version v1_all_pcs \
#     --subsample_size 20000 \
#     --n_seeds 5 \
#     --threads 8
# ================================================================

message("Starting Script 08b")

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(comethyl)
})

# ----------------------------------------------------------------
# 1. Parse command-line arguments
# ----------------------------------------------------------------
option_list <- list(
  make_option("--project_root",       type = "character"),
  make_option("--dataset1_label",     type = "character"),
  make_option("--dataset1_meth",      type = "character"),
  make_option("--dataset2_label",     type = "character", default = NULL),
  make_option("--dataset2_meth",      type = "character", default = NULL),
  make_option("--dataset3_label",     type = "character", default = NULL),
  make_option("--dataset3_meth",      type = "character", default = NULL),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--softpower_cor",      type = "character", default = "pearson"),
  make_option("--power_min",          type = "integer",   default = 1),
  make_option("--power_max",          type = "integer",   default = 20),
  make_option("--subsample_size",     type = "integer",   default = 20000),
  make_option("--n_seeds",            type = "integer",   default = 5),
  make_option("--seed_start",         type = "integer",   default = 1001),
  make_option("--block_size",         type = "integer",   default = 5000),
  make_option("--gc_interval",        type = "integer",   default = 4999),
  make_option("--threads",            type = "integer",   default = 4),
  make_option("--scale_free_cutoff",  type = "double",    default = 0.8)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------------------------------------------
# 2. Helpers (mirrors Script 08's helpers so outputs stay comparable)
# ----------------------------------------------------------------
# Adjusted methylation matrices (output of comethyl::adjustRegionMeth(), the
# --dataset*_meth inputs to this script) are samples-in-rows, regions-in-columns.
validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) stop(label, " must be matrix-like.")
  x <- as.matrix(x)
  if (!is.numeric(x))          stop(label, " must be numeric.")
  if (nrow(x) < 2)             stop(label, " must have >= 2 samples.")
  if (ncol(x) < 1)             stop(label, " has zero region columns.")
  if (is.null(rownames(x)))    stop(label, " must have sample IDs in rownames.")
  if (is.null(colnames(x)))    stop(label, " must have region IDs in colnames.")
  if (anyDuplicated(rownames(x))) stop(label, " has duplicated sample IDs.")
  if (anyDuplicated(colnames(x))) stop(label, " has duplicated region IDs.")
  x
}

extract_sft_table <- function(sft) {
  if (is.data.frame(sft))      return(sft)
  if (is.list(sft) && !is.null(sft$fitIndices)) return(as.data.frame(sft$fitIndices))
  if (is.list(sft) && length(sft) >= 1)         return(as.data.frame(sft[[1]]))
  stop("Cannot extract SFT table from getSoftPower() output.")
}

get_sft_columns <- function(sft_table) {
  cn_lower <- tolower(colnames(sft_table))
  get_col <- function(candidates) {
    for (c in candidates) {
      hit <- match(c, cn_lower)
      if (!is.na(hit)) return(colnames(sft_table)[hit])
    }
    NA_character_
  }
  list(
    power  = get_col(c("power")),
    fit    = get_col(c("sft.r.sq", "sft.rsq")),
    mean   = get_col(c("mean.k.", "mean.k")),
    median = get_col(c("median.k.", "median.k")),
    max    = get_col(c("max.k.", "max.k"))
  )
}

choose_consensus_power <- function(summary_df, cutoff = 0.8) {
  power_summary <- summary_df %>%
    dplyr::group_by(Power) %>%
    dplyr::summarise(
      n_datasets      = dplyr::n_distinct(Dataset),
      n_meet_cutoff   = dplyr::n_distinct(Dataset[Mean_SFT_R2 >= cutoff]),
      min_SFT_R2      = min(Mean_SFT_R2, na.rm = TRUE),
      median_SFT_R2   = median(Mean_SFT_R2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Power)

  all_ok <- dplyr::filter(power_summary, n_meet_cutoff == n_datasets)
  if (nrow(all_ok) > 0) {
    return(list(power = all_ok$Power[1],
                reason = "smallest power where all datasets meet scale-free cutoff (seed-averaged)"))
  }
  fallback <- dplyr::arrange(power_summary, dplyr::desc(n_meet_cutoff), Power)
  list(power  = fallback$Power[1],
       reason = "fallback: smallest power maximizing datasets meeting cutoff (seed-averaged)")
}

# ----------------------------------------------------------------
# 3. Validate inputs
# ----------------------------------------------------------------
if (is.null(opt$project_root))   stop("--project_root is required")
if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth))  stop("--dataset1_meth is required")
if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$dataset1_meth)) stop("dataset1_meth not found: ", opt$dataset1_meth)

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth))
    stop("Provide both --dataset2_label and --dataset2_meth.")
  if (!file.exists(opt$dataset2_meth)) stop("dataset2_meth not found: ", opt$dataset2_meth)
}
if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth))
    stop("Provide both --dataset3_label and --dataset3_meth.")
  if (!file.exists(opt$dataset3_meth)) stop("dataset3_meth not found: ", opt$dataset3_meth)
}

softpower_cor <- tolower(opt$softpower_cor)
if (!softpower_cor %in% c("bicor", "pearson")) stop("--softpower_cor must be bicor or pearson")
if (opt$subsample_size < 500) stop("--subsample_size should be at least 500 for a meaningful fit estimate")
if (opt$n_seeds < 1) stop("--n_seeds must be >= 1")

powerVector <- seq(opt$power_min, opt$power_max)
seeds <- seq(opt$seed_start, opt$seed_start + opt$n_seeds - 1)

WGCNA::disableWGCNAThreads()
WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ----------------------------------------------------------------
# 4. Derive cpg_label / region_label (same logic as Script 08)
# ----------------------------------------------------------------
d1_dir <- dirname(opt$dataset1_meth)
variant_pattern <- "^(v[0-9]|unadjusted)"
if (grepl(variant_pattern, basename(d1_dir))) {
  region_label <- basename(dirname(d1_dir))
  cpg_label    <- basename(dirname(dirname(d1_dir)))
} else {
  region_label <- basename(d1_dir)
  cpg_label    <- basename(dirname(d1_dir))
}

message("Derived labels:")
message("  cpg_label    : ", cpg_label)
message("  region_label : ", region_label)

subsample_tag <- paste0("subsample", opt$subsample_size)

# ----------------------------------------------------------------
# 5. Output directories
# ----------------------------------------------------------------
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "08b_softpower_benchmark")
shared_dir    <- file.path(step_dir, "shared", cpg_label, region_label,
                            opt$adjustment_version, subsample_tag)
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
message("Shared output directory: ", shared_dir)

# ----------------------------------------------------------------
# 6. Load and align datasets
# ----------------------------------------------------------------
dataset_inputs <- list(list(label = opt$dataset1_label, file = opt$dataset1_meth))
if (dataset2_provided)
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset2_label, file = opt$dataset2_meth)
if (dataset3_provided)
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset3_label, file = opt$dataset3_meth)

meth_list <- lapply(dataset_inputs, function(ds) {
  validate_meth_matrix(readRDS(ds$file), ds$label)
})
names(meth_list) <- sapply(dataset_inputs, `[[`, "label")

common_regions <- Reduce(intersect, lapply(meth_list, colnames))
message("Regions common to all provided datasets: ", length(common_regions))
message("(Sample counts per dataset: ",
        paste(sprintf("%s=%d", names(meth_list), sapply(meth_list, nrow)), collapse = ", "), ")")

if (length(common_regions) == 0) {
  stop("No regions are shared across the provided datasets -- check that inputs come from the same region_label.")
}
if (length(common_regions) < opt$subsample_size) {
  warning("--subsample_size (", opt$subsample_size,
          ") exceeds the number of common regions (", length(common_regions),
          "). Using all common regions instead.")
  opt$subsample_size <- length(common_regions)
}

# ----------------------------------------------------------------
# 7. Run subsampled soft power across seeds
# ----------------------------------------------------------------
all_seed_results <- list()

for (current_seed in seeds) {
  message("\n============================== Seed: ", current_seed)
  set.seed(current_seed)
  selected_region_ids <- sample(common_regions, size = opt$subsample_size, replace = FALSE)

  for (ds in dataset_inputs) {
    ds_label <- ds$label
    meth_sub <- meth_list[[ds_label]][, selected_region_ids, drop = FALSE]
    message("  ", ds_label, " subsample: ", nrow(meth_sub), " samples x ", ncol(meth_sub), " regions")

    ds_out_dir <- file.path(step_dir, ds_label, cpg_label, region_label,
                             opt$adjustment_version, subsample_tag)
    dir.create(ds_out_dir, recursive = TRUE, showWarnings = FALSE)

    sft <- getSoftPower(
      meth_sub,
      powerVector = powerVector,
      corType     = softpower_cor,
      file        = file.path(ds_out_dir, paste0("SoftPower_seed", current_seed, ".rds")),
      blockSize   = opt$block_size,
      gcInterval  = opt$gc_interval
    )

    sft_table <- extract_sft_table(sft)
    cols <- get_sft_columns(sft_table)
    if (any(is.na(unlist(cols))))
      stop("Could not identify expected columns in SFT table for dataset: ", ds_label)

    result_df <- tibble::tibble(
      Dataset           = ds_label,
      Seed              = current_seed,
      Power             = sft_table[[cols$power]],
      SFT_R2            = abs(sft_table[[cols$fit]]),
      Mean_Connectivity = sft_table[[cols$mean]]
    )

    readr::write_tsv(result_df, file.path(ds_out_dir, paste0("SoftPower_seed", current_seed, ".tsv")))
    writeLines(selected_region_ids, file.path(ds_out_dir, paste0("Selected_regions_seed", current_seed, ".txt")))

    all_seed_results[[paste(ds_label, current_seed, sep = "_")]] <- result_df
    message("  Done: ", ds_label, " (seed ", current_seed, ")")

    rm(meth_sub, sft, sft_table)
    gc()
  }
}

# ----------------------------------------------------------------
# 8. Combine seeds and summarize
# ----------------------------------------------------------------
combined <- dplyr::bind_rows(all_seed_results)
readr::write_tsv(combined, file.path(shared_dir, "combined_softpower_benchmark.tsv"))

summary_df <- combined %>%
  dplyr::group_by(Dataset, Power) %>%
  dplyr::summarise(
    n_seeds              = dplyr::n(),
    Mean_SFT_R2          = mean(SFT_R2, na.rm = TRUE),
    SD_SFT_R2            = stats::sd(SFT_R2, na.rm = TRUE),
    Mean_Connectivity    = mean(Mean_Connectivity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Dataset, Power)

readr::write_tsv(summary_df, file.path(shared_dir, "softpower_benchmark_summary.tsv"))

power_result <- choose_consensus_power(summary_df, cutoff = opt$scale_free_cutoff)

writeLines(
  c(
    "# SCREENING ESTIMATE ONLY -- subsampled data, not the full region set.",
    "# Confirm with a full Script 08 run before using for Script 09.",
    paste("chosen_power:", power_result$power),
    paste("reason:", power_result$reason),
    paste("scale_free_cutoff:", opt$scale_free_cutoff),
    paste("subsample_size:", opt$subsample_size),
    paste("n_seeds:", opt$n_seeds),
    paste("adjustment_version:", opt$adjustment_version),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("datasets:", paste(names(meth_list), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(shared_dir, "chosen_power_estimate.txt")
)
message("Estimated consensus power (subsampled): ", power_result$power, " (", power_result$reason, ")")

# ----------------------------------------------------------------
# 9. Plot: fit and connectivity vs power, ribbon = spread across seeds
# ----------------------------------------------------------------
pdf(file.path(shared_dir, "SoftPower_Benchmark_Combined.pdf"), width = 12, height = 6)

p1 <- ggplot(summary_df, aes(x = Power, y = Mean_SFT_R2, color = Dataset, fill = Dataset)) +
  geom_ribbon(aes(ymin = Mean_SFT_R2 - SD_SFT_R2, ymax = Mean_SFT_R2 + SD_SFT_R2),
              alpha = 0.15, color = NA) +
  geom_line() + geom_point() +
  geom_hline(yintercept = opt$scale_free_cutoff, linetype = "dashed") +
  geom_vline(xintercept = power_result$power, linetype = "dotted") +
  scale_x_continuous(breaks = powerVector) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank()) +
  labs(x = "Soft Threshold (Power)", y = "Scale-Free R^2 (mean +/- SD across seeds)",
       title = paste0(region_label, " | subsample = ", opt$subsample_size,
                       " | seeds = ", opt$n_seeds))

p2 <- ggplot(summary_df, aes(x = Power, y = Mean_Connectivity, color = Dataset)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = powerVector) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank()) +
  labs(x = "Soft Threshold (Power)", y = "Mean Connectivity (mean across seeds)")

print(p1)
print(p2)
dev.off()

message("Saved: SoftPower_Benchmark_Combined.pdf")

# ----------------------------------------------------------------
# 10. Run parameters and session info
# ----------------------------------------------------------------
writeLines(c(
  paste("project_root:",        opt$project_root),
  paste("adjustment_version:",  opt$adjustment_version),
  paste("cpg_label:",           cpg_label),
  paste("region_label:",        region_label),
  paste("datasets:",            paste(names(meth_list), collapse = ", ")),
  paste("softpower_cor:",       softpower_cor),
  paste("power_min:",           opt$power_min),
  paste("power_max:",           opt$power_max),
  paste("subsample_size:",      opt$subsample_size),
  paste("n_seeds:",             opt$n_seeds),
  paste("seeds:",               paste(seeds, collapse = ", ")),
  paste("scale_free_cutoff:",   opt$scale_free_cutoff),
  paste("block_size:",          opt$block_size),
  paste("gc_interval:",         opt$gc_interval),
  paste("threads:",             opt$threads),
  paste("estimated_power:",     power_result$power),
  paste("estimate_reason:",     power_result$reason),
  paste("shared_dir:",          shared_dir),
  paste("date:",                as.character(Sys.time()))
), con = file.path(shared_dir, "run_parameters.txt"))

writeLines(capture.output(sessionInfo()), con = file.path(shared_dir, "sessionInfo.txt"))

message("\nScript 08b complete: subsampled soft power benchmark finished")
message("Shared output: ", shared_dir)
message("This is a screening estimate. Confirm with a full Script 08 run before finalizing power for Script 09.")