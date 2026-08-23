#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08c: Compare Soft-Power Benchmarks Across Region Filters
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Aggregates multiple 08b_softpower_benchmark_consensus.R runs
#   (each representing a different covMin/methSD region_label) into
#   one ranked comparison, so you can answer "which region filter
#   gives the best soft-power fit" without eyeballing separate PDFs
#   and logs for each combo.
#
#   Run 08b once per candidate covMin/methSD combo first (each is a
#   fast subsampled run), then point this script at the shared
#   cpg_label folder that contains all of them.
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --cpg_label    : the CpG-calling label shared by all combos being
#                     compared (e.g. cov3_75pct) -- this is the folder
#                     directly under 08b_softpower_benchmark/shared/
#
# OPTIONAL INPUTS
#   --adjustment_version : restrict to one adjustment version
#                           [default = all versions found]
#   --scale_free_cutoff  : R^2 cutoff used to rank combos [default = 0.8]
#   --region_totals_glob : optional glob pattern (relative to
#                           02_reference_region_filter/) used to look
#                           up the actual filtered region count for
#                           each region_label, e.g.
#                           "*/*/Region_Totals.txt" is not needed --
#                           by default this script reads
#                           "Filtered_Regions.txt" line counts from
#                           02_reference_region_filter/<dataset>/<cpg_label>/<region_label>/
#                           if present, and leaves n_regions_full = NA
#                           otherwise.
#   --reference_dataset_label : dataset subfolder to use when looking
#                           up n_regions_full under
#                           02_reference_region_filter/ [default =
#                           first dataset label found in the benchmark
#                           data]
#
# OUTPUTS
#   comethyl_output/consensus/08c_softpower_benchmark_compare/<cpg_label>/
#       region_filter_comparison.tsv     (one row per region_label x adjustment_version,
#                                          ranked by fit)
#       Region_Filter_Comparison.pdf     (fit curves for every combo, faceted by dataset)
#       run_parameters.txt
#
# NOTES
#   - This only compares combos that have already been run through
#     08b. It does not run any new soft-power calculations itself.
#   - Ranking logic: for each region_label x adjustment_version,
#     find the smallest power where ALL datasets' seed-averaged R^2
#     meets --scale_free_cutoff. If none reach cutoff within the
#     tested power range, fall back to the best R^2 achieved at the
#     max tested power. Combos are ranked by (1) whether cutoff was
#     reached, (2) power needed to reach it (lower is better -- it
#     leaves more headroom), (3) region count (fewer is better, all
#     else equal, since smaller matrices are cheaper downstream).
#
# EXAMPLE
#   Rscript 08c_softpower_benchmark_compare_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#     --cpg_label cov3_75pct \
#     --adjustment_version v1_all_pcs
# ================================================================

message("Starting Script 08c")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(stringr)
})

# ----------------------------------------------------------------
# 1. Parse command-line arguments
# ----------------------------------------------------------------
option_list <- list(
  make_option("--project_root",             type = "character"),
  make_option("--cpg_label",                type = "character"),
  make_option("--adjustment_version",       type = "character", default = NULL),
  make_option("--scale_free_cutoff",        type = "double",    default = 0.8),
  make_option("--reference_dataset_label",  type = "character", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$cpg_label))    stop("--cpg_label is required")
if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
benchmark_root <- file.path(pipeline_root, "08b_softpower_benchmark", "shared", opt$cpg_label)

if (!dir.exists(benchmark_root)) {
  stop("No 08b benchmark output found for cpg_label '", opt$cpg_label,
       "' at: ", benchmark_root,
       "\nRun 08b_softpower_benchmark_consensus.R for at least one region_label first.")
}

# ----------------------------------------------------------------
# 2. Find all benchmark summary files under this cpg_label
#
# Directory shape:
#   .../08b_softpower_benchmark/shared/<cpg_label>/<region_label>/
#       <adjustment_version>/<subsample_tag>/softpower_benchmark_summary.tsv
# ----------------------------------------------------------------
summary_files <- list.files(
  benchmark_root,
  pattern = "^softpower_benchmark_summary\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No softpower_benchmark_summary.tsv files found under: ", benchmark_root)
}

message("Found ", length(summary_files), " benchmark summary file(s)")

# ----------------------------------------------------------------
# 3. Parse region_label / adjustment_version / subsample_tag from each path
#    and optionally filter to one adjustment_version
# ----------------------------------------------------------------
parse_combo <- function(path) {
  # path relative to benchmark_root:
  #   <region_label>/<adjustment_version>/<subsample_tag>/softpower_benchmark_summary.tsv
  rel <- sub(paste0("^", benchmark_root, "/?"), "", path)
  parts <- strsplit(rel, "/")[[1]]
  if (length(parts) < 4) {
    warning("Unexpected path shape, skipping: ", path)
    return(NULL)
  }
  list(
    region_label       = parts[1],
    adjustment_version = parts[2],
    subsample_tag      = parts[3],
    path                = path
  )
}

combos <- purrr::map(summary_files, parse_combo)
combos <- purrr::compact(combos)

if (!is.null(opt$adjustment_version)) {
  keep <- purrr::map_lgl(combos, ~ .x$adjustment_version == opt$adjustment_version)
  combos <- combos[keep]
  if (length(combos) == 0)
    stop("No benchmark combos found for adjustment_version = ", opt$adjustment_version)
}

message("Comparing ", length(combos), " region_label x adjustment_version combo(s)")

# ----------------------------------------------------------------
# 4. Load and tag each summary table
# ----------------------------------------------------------------
combo_id <- function(c) paste(c$region_label, c$adjustment_version, c$subsample_tag, sep = " | ")

all_summaries <- purrr::map_dfr(combos, function(c) {
  df <- readr::read_tsv(c$path, show_col_types = FALSE)
  df$region_label       <- c$region_label
  df$adjustment_version <- c$adjustment_version
  df$subsample_tag       <- c$subsample_tag
  df$combo_id            <- combo_id(c)
  df
})

if (!"Mean_Slope" %in% colnames(all_summaries)) {
  warning("These benchmark files predate slope-aware screening. Mean_Slope is ",
          "unavailable, so no combination will be marked as reaching the ",
          "revised topology criterion. Rerun Script 08b for final ranking.")
  all_summaries$Mean_Slope <- NA_real_
}

# ----------------------------------------------------------------
# 5. Optional: look up actual filtered region count per region_label
#    from Script 02's Filtered_Regions.txt (fast line count, if present)
# ----------------------------------------------------------------
reference_dataset_label <- opt$reference_dataset_label
if (is.null(reference_dataset_label)) {
  reference_dataset_label <- unique(all_summaries$Dataset)[1]
  message("No --reference_dataset_label given, using: ", reference_dataset_label)
}

get_n_regions_full <- function(region_label) {
  # Revised joint-dynamic candidates live under Script 05b. Prefer the exact
  # region IDs used by the benchmark, then fall back to legacy Script 02.
  joint_ids <- file.path(pipeline_root, "05b_shared_complete_regions",
                         opt$cpg_label, region_label,
                         "joint_eligible_regionIDs.txt")
  if (file.exists(joint_ids)) {
    return(tryCatch(length(readLines(joint_ids)),
                    error = function(e) NA_integer_))
  }
  filt_path <- file.path(pipeline_root, "02_reference_region_filter",
                          reference_dataset_label, opt$cpg_label, region_label,
                          "Filtered_Regions.txt")
  if (!file.exists(filt_path)) return(NA_integer_)
  # subtract 1 for header
  n <- tryCatch(length(readLines(filt_path)) - 1L, error = function(e) NA_integer_)
  n
}

region_counts <- tibble::tibble(
  region_label = unique(all_summaries$region_label)
) %>%
  dplyr::mutate(n_regions_full = purrr::map_int(region_label, get_n_regions_full))

# ----------------------------------------------------------------
# 6. Rank combos
# ----------------------------------------------------------------
rank_table <- all_summaries %>%
  dplyr::group_by(combo_id, region_label, adjustment_version, Power) %>%
  dplyr::summarise(
    n_datasets    = dplyr::n_distinct(Dataset),
    n_meet_cutoff = dplyr::n_distinct(Dataset[
      Mean_SFT_R2 >= opt$scale_free_cutoff & !is.na(Mean_Slope) & Mean_Slope < 0
    ]),
    min_R2        = min(Mean_SFT_R2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(combo_id, Power)

best_per_combo <- rank_table %>%
  dplyr::group_by(combo_id, region_label, adjustment_version) %>%
  dplyr::group_modify(function(df, key) {
    all_ok <- dplyr::filter(df, n_meet_cutoff == n_datasets)
    if (nrow(all_ok) > 0) {
      row <- all_ok[1, ]
      tibble::tibble(
        cutoff_reached    = TRUE,
        power_at_cutoff   = row$Power,
        min_R2_at_cutoff  = row$min_R2
      )
    } else {
      row <- df[which.max(df$min_R2), ]
      tibble::tibble(
        cutoff_reached    = FALSE,
        power_at_cutoff   = row$Power,
        min_R2_at_cutoff  = row$min_R2
      )
    }
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(region_counts, by = "region_label") %>%
  dplyr::arrange(dplyr::desc(cutoff_reached), power_at_cutoff, n_regions_full)

# ----------------------------------------------------------------
# 7. Output directory and files
# ----------------------------------------------------------------
out_dir <- file.path(pipeline_root, "08c_softpower_benchmark_compare", opt$cpg_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_tsv(best_per_combo, file.path(out_dir, "region_filter_comparison.tsv"))
message("Saved: ", file.path(out_dir, "region_filter_comparison.tsv"))
message("\nRanked combos (best first):")
print(as.data.frame(best_per_combo))

# ----------------------------------------------------------------
# 8. Comparison plot: fit curve per combo, faceted by dataset
# ----------------------------------------------------------------
pdf(file.path(out_dir, "Region_Filter_Comparison.pdf"), width = 12, height = 7)

p <- ggplot(all_summaries, aes(x = Power, y = Mean_SFT_R2, color = combo_id)) +
  geom_line() + geom_point(size = 1) +
  geom_hline(yintercept = opt$scale_free_cutoff, linetype = "dashed") +
  facet_wrap(~ Dataset) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(x = "Soft Threshold (Power)", y = "Scale-Free R^2 (seed-averaged)",
       color = "region_label | adjustment_version | subsample",
       title = paste0("Region filter comparison: ", opt$cpg_label))

print(p)
dev.off()

message("Saved: ", file.path(out_dir, "Region_Filter_Comparison.pdf"))

# ----------------------------------------------------------------
# 9. Run parameters
# ----------------------------------------------------------------
writeLines(c(
  paste("project_root:",            opt$project_root),
  paste("cpg_label:",                opt$cpg_label),
  paste("adjustment_version_filter:", ifelse(is.null(opt$adjustment_version), "ALL", opt$adjustment_version)),
  paste("scale_free_cutoff:",        opt$scale_free_cutoff),
  paste("reference_dataset_label:",  reference_dataset_label),
  paste("n_combos_compared:",        length(combos)),
  paste("out_dir:",                  out_dir),
  paste("date:",                     as.character(Sys.time()))
), con = file.path(out_dir, "run_parameters.txt"))

message("\nScript 08c complete: region filter comparison finished")
message("Output: ", out_dir)
