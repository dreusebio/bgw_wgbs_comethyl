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
#   --subsample_tag      : restrict to one benchmark size, e.g. subsample30000
#                           [default = all sizes found]
#   --region_label_regex : optional regular expression used to retain selected
#                           joint/legacy region labels [default = all]
#   --scale_free_cutoff  : R^2 cutoff used to rank combos [default = 0.8]
#   --reference_dataset_label : dataset subfolder to use when looking
#                           up n_regions_full under
#                           02_reference_region_filter/. Set this
#                           explicitly to your reference dataset (the
#                           one Script 02 was actually run on) --
#                           otherwise this script guesses by taking
#                           the first dataset label it encounters,
#                           which may not have a 02_reference_region_filter
#                           folder and will leave n_regions_full as NA.
#
# OUTPUTS
#   comethyl_output/consensus/08c_softpower_benchmark_compare/<cpg_label>/
#       region_filter_comparison.tsv     (one row per region_label x
#                                          adjustment_version x subsample_tag,
#                                          ranked by fit)
#       region_filter_plot_key.tsv       (short plot ID -> complete settings)
#       Region_Filter_Comparison.pdf     (page 1: fit curves with short IDs;
#                                          page 2: complete ID key)
#       run_parameters.txt
#
# NOTES
#   - This only compares combos that have already been run through
#     08b. It does not run any new soft-power calculations itself.
#   - Ranking logic: for each region_label x adjustment_version x subsample,
#     find the smallest power where ALL datasets have seed-averaged R2
#     meeting --scale_free_cutoff AND negative mean slope. If none pass,
#     report the power with the strongest worst-dataset R2 but mark the
#     combination cutoff_reached=FALSE. Combos are ranked by (1) whether cutoff was
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
  make_option("--subsample_tag",            type = "character", default = NULL),
  make_option("--region_label_regex",       type = "character", default = NULL),
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
  # Both paths are produced by list.files(full.names=TRUE), so a literal prefix
  # check is safer than building a regex from a filesystem path.
  if (!startsWith(path, benchmark_root)) {
    warning("Benchmark file is outside expected root, skipping: ", path)
    return(NULL)
  }
  rel <- substring(path, nchar(benchmark_root) + 1L)
  rel <- sub("^/+", "", rel)
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

if (!is.null(opt$subsample_tag)) {
  keep <- purrr::map_lgl(combos, ~ .x$subsample_tag == opt$subsample_tag)
  combos <- combos[keep]
  if (length(combos) == 0)
    stop("No benchmark combinations found for subsample_tag = ", opt$subsample_tag)
}

if (!is.null(opt$region_label_regex)) {
  keep <- purrr::map_lgl(
    combos,
    ~ grepl(opt$region_label_regex, .x$region_label, perl = TRUE)
  )
  combos <- combos[keep]
  if (length(combos) == 0)
    stop("No region labels matched --region_label_regex = ", opt$region_label_regex)
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

# The original plot used the complete region/filter/version description in the
# legend. Those strings are informative in a table but too long for a figure.
# Assign stable short IDs after sorting by the actual settings. The IDs are used
# only for display; the complete fields remain in every tabular output.
plot_key <- all_summaries %>%
  dplyr::distinct(combo_id, region_label, adjustment_version, subsample_tag) %>%
  dplyr::arrange(region_label, adjustment_version, subsample_tag) %>%
  dplyr::mutate(
    plot_id = sprintf("C%02d", dplyr::row_number()),
    plot_description = paste(region_label, adjustment_version, subsample_tag,
                             sep = " | ")
  )

all_summaries <- all_summaries %>%
  dplyr::left_join(
    dplyr::select(plot_key, combo_id, plot_id),
    by = "combo_id"
  )

# Slope-aware benchmark files are required for valid ranking. Older 08b output
# can still be read for plotting, but it cannot be declared topology-valid.
if (!"Mean_Slope" %in% colnames(all_summaries)) {
  warning("Mean_Slope is absent. These are old 08b results; rerun revised 08b ",
          "before choosing a final filter.")
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
  # New joint-dynamic structure from Script 05b.
  joint_ids <- file.path(
    pipeline_root, "05b_shared_complete_regions", opt$cpg_label, region_label,
    "joint_eligible_regionIDs.txt"
  )
  if (file.exists(joint_ids)) {
    return(tryCatch(length(readLines(joint_ids)),
                    error = function(e) NA_integer_))
  }

  # Legacy reference-selected structure from Script 02.
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
  dplyr::group_by(combo_id, plot_id, region_label, adjustment_version,
                  subsample_tag, Power) %>%
  dplyr::summarise(
    n_datasets    = dplyr::n_distinct(Dataset),
    n_meet_cutoff = dplyr::n_distinct(Dataset[which(
      Mean_SFT_R2 >= opt$scale_free_cutoff & !is.na(Mean_Slope) & Mean_Slope < 0
    )]),
    min_R2        = min(Mean_SFT_R2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(combo_id, Power)

best_per_combo <- rank_table %>%
  dplyr::group_by(combo_id, plot_id, region_label, adjustment_version,
                  subsample_tag) %>%
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
  # Successful candidates are ranked by the lowest valid power. If none pass,
  # the most useful ordering is the strongest worst-dataset R2, not the lowest
  # power or smallest region count.
  dplyr::arrange(
    dplyr::desc(cutoff_reached),
    dplyr::if_else(cutoff_reached, power_at_cutoff, Inf),
    dplyr::desc(min_R2_at_cutoff),
    n_regions_full
  )

# ----------------------------------------------------------------
# 7. Output directory and files
# ----------------------------------------------------------------
out_dir <- file.path(pipeline_root, "08c_softpower_benchmark_compare", opt$cpg_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_tsv(best_per_combo, file.path(out_dir, "region_filter_comparison.tsv"))
message("Saved: ", file.path(out_dir, "region_filter_comparison.tsv"))

plot_key_out <- plot_key %>%
  dplyr::left_join(
    dplyr::select(best_per_combo, plot_id, cutoff_reached, power_at_cutoff,
                  min_R2_at_cutoff, n_regions_full),
    by = "plot_id"
  ) %>%
  dplyr::arrange(plot_id)

readr::write_tsv(plot_key_out, file.path(out_dir, "region_filter_plot_key.tsv"))
message("Saved: ", file.path(out_dir, "region_filter_plot_key.tsv"))
message("\nRanked combos (best first):")
print(as.data.frame(best_per_combo))

# ----------------------------------------------------------------
# 8. Comparison PDF
#
# Page 1 uses short IDs so the legend remains readable even when region labels
# are long. Page 2 is a complete lookup table. The same lookup is also written
# to region_filter_plot_key.tsv for sorting/filtering outside the PDF.
# ----------------------------------------------------------------
pdf(file.path(out_dir, "Region_Filter_Comparison.pdf"), width = 13, height = 8.5)

p <- ggplot(all_summaries,
            aes(x = Power, y = Mean_SFT_R2, color = plot_id,
                linetype = adjustment_version,
                group = interaction(plot_id, Dataset))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1) +
  geom_hline(yintercept = opt$scale_free_cutoff, linetype = "dashed") +
  facet_wrap(~ Dataset) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = grid::unit(1.4, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  ) +
  guides(
    color = guide_legend(title = "Candidate ID", nrow = 2, byrow = TRUE),
    linetype = guide_legend(title = "Adjustment")
  ) +
  labs(x = "Soft Threshold (Power)", y = "Scale-Free R^2 (seed-averaged)",
       title = paste0("Region filter comparison: ", opt$cpg_label),
       subtitle = "Use the candidate key on page 2 (or region_filter_plot_key.tsv) for complete settings")

print(p)

# Create a dependency-free table page with grid graphics. Long descriptions
# wrap inside the page instead of expanding the plot legend horizontally.
key_display <- plot_key_out %>%
  dplyr::transmute(
    ID = plot_id,
    Settings = plot_description,
    `Pass cutoff` = ifelse(cutoff_reached, "YES", "NO"),
    `Selected power` = power_at_cutoff,
    `Worst R2` = sprintf("%.3f", min_R2_at_cutoff),
    `Full regions` = ifelse(is.na(n_regions_full), "NA",
                            format(n_regions_full, big.mark = ","))
  )

grid::grid.newpage()
grid::grid.text(
  paste0("Candidate key: ", opt$cpg_label),
  x = 0.03, y = 0.965, just = c("left", "top"),
  gp = grid::gpar(fontsize = 16, fontface = "bold")
)
grid::grid.text(
  "Full filter definitions corresponding to the short IDs on page 1",
  x = 0.03, y = 0.925, just = c("left", "top"),
  gp = grid::gpar(fontsize = 10, col = "grey30")
)

n_key <- nrow(key_display)
if (n_key > 0) {
  # Divide the available vertical space evenly. A minimum font size prevents a
  # very large comparison grid from producing unreadably small text; the TSV is
  # the authoritative key if more candidates are compared than fit comfortably.
  row_height <- min(0.075, 0.82 / n_key)
  key_font <- max(6.5, min(9.5, 7.5 / sqrt(max(1, n_key / 10))))
  y_start <- 0.875

  headers <- c("ID", "Settings (region | adjustment | subsample)",
               "Pass", "Power", "Worst R2", "Regions")
  x_pos <- c(0.03, 0.10, 0.69, 0.79, 0.87, 0.96)
  justs <- c("left", "left", "center", "center", "center", "right")
  for (j in seq_along(headers)) {
    grid::grid.text(headers[j], x = x_pos[j], y = y_start,
                    just = c(justs[j], "center"),
                    gp = grid::gpar(fontsize = 9, fontface = "bold"))
  }
  grid::grid.lines(x = c(0.03, 0.96), y = rep(y_start - 0.022, 2),
                   gp = grid::gpar(col = "grey45"))

  for (i in seq_len(n_key)) {
    y <- y_start - i * row_height
    vals <- c(key_display$ID[i], key_display$Settings[i],
              key_display$`Pass cutoff`[i], key_display$`Selected power`[i],
              key_display$`Worst R2`[i], key_display$`Full regions`[i])
    if (i %% 2 == 0) {
      grid::grid.rect(x = 0.495, y = y, width = 0.93,
                      height = row_height * 0.9,
                      gp = grid::gpar(fill = "grey96", col = NA))
    }
    for (j in seq_along(vals)) {
      grid::grid.text(vals[j], x = x_pos[j], y = y,
                      just = c(justs[j], "center"),
                      gp = grid::gpar(fontsize = key_font))
    }
  }
}

dev.off()

message("Saved: ", file.path(out_dir, "Region_Filter_Comparison.pdf"))

# ----------------------------------------------------------------
# 9. Run parameters
# ----------------------------------------------------------------
writeLines(c(
  paste("project_root:",            opt$project_root),
  paste("cpg_label:",                opt$cpg_label),
  paste("adjustment_version_filter:", ifelse(is.null(opt$adjustment_version), "ALL", opt$adjustment_version)),
  paste("subsample_tag_filter:",       ifelse(is.null(opt$subsample_tag), "ALL", opt$subsample_tag)),
  paste("region_label_regex:",         ifelse(is.null(opt$region_label_regex), "ALL", opt$region_label_regex)),
  paste("scale_free_cutoff:",        opt$scale_free_cutoff),
  paste("reference_dataset_label:",  reference_dataset_label),
  paste("n_combos_compared:",        length(combos)),
  paste("out_dir:",                  out_dir),
  paste("date:",                     as.character(Sys.time()))
), con = file.path(out_dir, "run_parameters.txt"))

message("\nScript 08c complete: region filter comparison finished")
message("Output: ", out_dir)
