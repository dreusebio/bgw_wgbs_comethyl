#!/usr/bin/env Rscript
# ============================================================
# 13B_regulatory_state_presentation.R
#
# Purpose
#   Presentation/reporting layer for regulatory-state overlap results
#   produced by 13A_core_regulatory_state_overlap.R.
#
# Reads 13A outputs and creates:
#     1) heatmaps (single module or combined modules)
#     2) per-source annotation summaries
#     3) stacked barplots of annotation composition
#     4) top-region / top-gene slide-friendly outputs
#
# Assumes 13A output structure:
#   <out_parent>/<run_tag>/module_<module>/
#       dominant_state_long_all_sources.csv
#       dominant_state_matrix_all_sources.csv
#       <source>/dominant_state_long.csv
#       <source>/dominant_state_matrix.csv
#
# Design goals
#   - no bedtools here
#   - no downloads here
#   - source-aware
#   - reproducible
# ============================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ============================================================
# 1) CLI helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

msg <- function(...) cat(sprintf(...), "\n")

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 2) source / annotation colors
# ============================================================
ROADMAP15_DESC_COLORS <- c(
  "Active TSS"                    = "#FF0000",
  "Flanking Active TSS"           = "#FF4500",
  "Transcribed at 5' and 3'"      = "#32CD32",
  "Strong transcription"          = "#008000",
  "Weak transcription"            = "#006400",
  "Genic enhancer"                = "#C2E105",
  "Enhancer"                      = "#FFFF00",
  "ZNF genes and repeats"         = "#66CDAA",
  "Heterochromatin"               = "#8A91D0",
  "Bivalent TSS"                  = "#CD5C5C",
  "Flanking bivalent TSS/enhancer"= "#E9967A",
  "Bivalent enhancer"             = "#BDB76B",
  "Repressed PolyComb"            = "#808080",
  "Weak repressed PolyComb"       = "#C0C0C0",
  "Quiescent/Low"                 = "#F7F7F7",
  "Missing"                       = "#EDEDED"
)

DEFAULT_OTHER_COL <- "#D0D0D0"

# ============================================================
# 3) loaders
# ============================================================
read_module_long <- function(module_dir) {
  f <- file.path(module_dir, "dominant_state_long_all_sources.csv")
  if (!file.exists(f)) {
    stop("Missing 13A long file: ", f)
  }

  df <- read_csv(f, show_col_types = FALSE)

  need <- c("module", "source", "dataset", "track_id", "track_label",
            "region_id", "annotation", "annotation_desc")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("Missing required columns in ", f, ": ", paste(missing, collapse = ", "))
  }

  if (!("gene_symbol" %in% names(df))) df$gene_symbol <- NA_character_
  if (!("membership" %in% names(df)))  df$membership  <- NA_real_
  if (!("overlap_bp" %in% names(df)))  df$overlap_bp  <- NA_real_

  df %>%
    mutate(
      module = as.character(module),
      source = as.character(source),
      dataset = as.character(dataset),
      track_id = as.character(track_id),
      track_label = as.character(track_label),
      region_id = as.character(region_id),
      gene_symbol = as.character(gene_symbol),
      annotation = as.character(annotation),
      annotation_desc = as.character(annotation_desc),
      membership = suppressWarnings(as.numeric(membership)),
      overlap_bp = suppressWarnings(as.numeric(overlap_bp))
    )
}

# ============================================================
# 4) ranking helpers
# ============================================================
rank_regions <- function(df) {
  df %>%
    group_by(module, region_id) %>%
    summarise(
      gene_symbol = {
        gs <- unique(gene_symbol[!is.na(gene_symbol) & gene_symbol != ""])
        if (length(gs) == 0) NA_character_ else gs[1]
      },
      membership = if (all(is.na(membership))) NA_real_ else max(membership, na.rm = TRUE),
      overlap_bp = if (all(is.na(overlap_bp))) NA_real_ else max(overlap_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      membership = ifelse(is.infinite(membership), NA_real_, membership),
      overlap_bp = ifelse(is.infinite(overlap_bp), NA_real_, overlap_bp)
    )
}

rank_genes <- function(df) {
  df %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    group_by(module, gene_symbol) %>%
    summarise(
      max_membership = if (all(is.na(membership))) NA_real_ else max(membership, na.rm = TRUE),
      max_overlap_bp = if (all(is.na(overlap_bp))) NA_real_ else max(overlap_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      max_membership = ifelse(is.infinite(max_membership), NA_real_, max_membership),
      max_overlap_bp = ifelse(is.infinite(max_overlap_bp), NA_real_, max_overlap_bp)
    )
}

limit_to_top_genes <- function(df, max_genes = NA_integer_) {
  if (is.na(max_genes) || !is.finite(max_genes) || max_genes <= 0) return(df)

  g <- rank_genes(df) %>%
    arrange(module, desc(max_membership), desc(max_overlap_bp), gene_symbol) %>%
    group_by(module) %>%
    slice_head(n = max_genes) %>%
    ungroup()

  df %>% semi_join(g, by = c("module", "gene_symbol"))
}

limit_to_top_regions <- function(df, max_regions = 50) {
  r <- rank_regions(df) %>%
    arrange(module, desc(membership), desc(overlap_bp), region_id) %>%
    group_by(module) %>%
    slice_head(n = max_regions) %>%
    ungroup()

  df %>% semi_join(r, by = c("module", "region_id"))
}

# ============================================================
# 5) matrix builder
# ============================================================
build_heatmap_matrix <- function(df, combine_modules = FALSE) {
  x <- df %>%
    mutate(
      row_label = case_when(
        !is.na(gene_symbol) & gene_symbol != "" ~ paste0(region_id, " | ", gene_symbol),
        TRUE ~ region_id
      ),
      col_label = case_when(
        combine_modules ~ paste(module, source, track_label, sep = " | "),
        TRUE ~ paste(source, track_label, sep = " | ")
      )
    )

  if (combine_modules) {
    x <- x %>%
      mutate(row_label = paste(module, row_label, sep = " :: "))
  }

  x %>%
    select(row_label, col_label, annotation_desc) %>%
    distinct() %>%
    pivot_wider(names_from = col_label, values_from = annotation_desc)
}

# ============================================================
# 6) ComplexHeatmap plot
# ============================================================
make_complex_heatmap <- function(mat_df, out_pdf, title = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("ComplexHeatmap not installed.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize not installed.")
  }

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
  })

  if (nrow(mat_df) == 0) {
    grDevices::pdf(out_pdf, width = 10, height = 7)
    grid::grid.newpage()
    grid::grid.text("No rows to plot", gp = grid::gpar(fontsize = 14, fontface = "bold"))
    grDevices::dev.off()
    return(invisible(out_pdf))
  }

  row_labels <- mat_df[[1]]
  m <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(m) <- row_labels

  m2 <- m
  m2[is.na(m2) | m2 == ""] <- "Missing"

  present <- sort(unique(as.vector(m2)))
  colors <- ROADMAP15_DESC_COLORS
  unknown <- setdiff(present, names(colors))
  if (length(unknown) > 0) {
    colors <- c(colors, setNames(rep(DEFAULT_OTHER_COL, length(unknown)), unknown))
  }

  grDevices::pdf(out_pdf, width = 14, height = 9)

  ht <- ComplexHeatmap::Heatmap(
    m2,
    name = "Annotation",
    col = colors,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 9),
    rect_gp = grid::gpar(col = "grey70", lwd = 0.5),
    column_title = title,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
  )

  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  grDevices::dev.off()

  invisible(out_pdf)
}

# ============================================================
# 7) summary tables
# ============================================================
make_annotation_summary <- function(df) {
  df %>%
    group_by(module, source, annotation_desc) %>%
    summarise(
      n_regions = n_distinct(region_id),
      n_tracks = n_distinct(track_id),
      .groups = "drop"
    ) %>%
    group_by(module, source) %>%
    mutate(frac_regions = n_regions / sum(n_regions)) %>%
    ungroup() %>%
    arrange(module, source, desc(n_regions), annotation_desc)
}

make_top_annotation_summary <- function(df) {
  df %>%
    group_by(module, source, region_id) %>%
    slice_max(order_by = overlap_bp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    group_by(module, source, annotation_desc) %>%
    summarise(
      n_regions = n_distinct(region_id),
      .groups = "drop"
    ) %>%
    group_by(module, source) %>%
    mutate(frac_regions = n_regions / sum(n_regions)) %>%
    ungroup() %>%
    arrange(module, source, desc(n_regions), annotation_desc)
}

# ============================================================
# 8) stacked barplot
# ============================================================
plot_annotation_composition <- function(sum_df, out_pdf, title = "Annotation composition") {
  if (nrow(sum_df) == 0) return(invisible(NULL))

  color_map <- ROADMAP15_DESC_COLORS
  unknown <- setdiff(unique(sum_df$annotation_desc), names(color_map))
  if (length(unknown) > 0) {
    color_map <- c(color_map, setNames(rep(DEFAULT_OTHER_COL, length(unknown)), unknown))
  }

  p <- ggplot(sum_df, aes(x = module, y = frac_regions, fill = annotation_desc)) +
    geom_col(position = "stack") +
    facet_wrap(~source, scales = "free_x") +
    scale_fill_manual(values = color_map) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      title = title,
      x = "Module",
      y = "Fraction of regions",
      fill = "Annotation"
    )

  ggsave(out_pdf, plot = p, width = 12, height = 8, units = "in", dpi = 300)
  invisible(out_pdf)
}

# ============================================================
# 9) source filter
# ============================================================
filter_sources <- function(df, source_arg = "") {
  if (is.null(source_arg) || is.na(source_arg) || !nzchar(source_arg)) return(df)
  keep <- split_csv(source_arg)
  df %>% filter(source %in% keep)
}

# ============================================================
# 10) main args
# ============================================================
project_root <- get_arg("--project_root", "")
region_version <- get_arg("--region_version", "")
if (!nzchar(project_root)) stop("--project_root is required")
if (!nzchar(region_version)) stop("--region_version is required")

base_adj <- file.path(project_root, "comethyl_output", "filter_regions", region_version, "methylation_adjustment")
default_state_root <- file.path(base_adj, "regulatory_state_overlap")

state_root <- get_arg("--state_root", default_state_root)
run_tag <- get_arg("--run_tag", "")
modules_in <- unique(split_csv(get_arg("--modules", "")))
if (length(modules_in) == 0) stop("Provide --modules")

source_arg <- get_arg("--source", "")
max_regions <- as.integer(get_arg("--max_regions", "50"))
max_genes   <- as.integer(get_arg("--max_genes", NA_character_))
combine_modules <- as_bool(get_arg("--combine_modules", "true"), TRUE)
make_per_module_heatmaps <- as_bool(get_arg("--make_per_module_heatmaps", "true"), TRUE)
make_combined_heatmap <- as_bool(get_arg("--make_combined_heatmap", "true"), TRUE)
make_barplots <- as_bool(get_arg("--make_barplots", "true"), TRUE)

if (!nzchar(run_tag)) {
  stop("--run_tag is required, e.g. v1_all_pcs")
}

run_root <- file.path(state_root, run_tag)
if (!dir.exists(run_root)) stop("run_tag folder not found: ", run_root)

out_dir <- file.path(run_root, "presentation")
safe_dir_create(out_dir)

msg("state_root: %s", state_root)
msg("run_tag: %s", run_tag)
msg("modules: %s", paste(modules_in, collapse = ", "))
msg("source filter: %s", ifelse(nzchar(source_arg), source_arg, "(all)"))

# ============================================================
# 11) load modules
# ============================================================
module_data <- list()

for (mod in modules_in) {
  module_dir <- file.path(run_root, paste0("module_", mod))
  if (!dir.exists(module_dir)) {
    msg("[WARN] Missing module dir: %s", module_dir)
    next
  }

  df <- read_module_long(module_dir)
  df <- filter_sources(df, source_arg)

  if (nrow(df) == 0) {
    msg("[WARN] No rows left after source filtering for module=%s", mod)
    next
  }

  df <- limit_to_top_genes(df, max_genes = max_genes)
  df <- limit_to_top_regions(df, max_regions = max_regions)

  module_data[[mod]] <- df
}

if (length(module_data) == 0) {
  stop("No module data available after loading/filtering.")
}

combined_df <- bind_rows(module_data)

# ============================================================
# 12) save summaries
# ============================================================
sum_all <- make_annotation_summary(combined_df)
write_csv(sum_all, file.path(out_dir, "annotation_summary_all.csv"))

sum_top <- make_top_annotation_summary(combined_df)
write_csv(sum_top, file.path(out_dir, "annotation_summary_top_regions.csv"))

# ============================================================
# 13) barplots
# ============================================================
if (isTRUE(make_barplots)) {
  plot_annotation_composition(
    sum_df = sum_all,
    out_pdf = file.path(out_dir, "annotation_composition_all.pdf"),
    title = "Annotation composition across selected modules"
  )

  plot_annotation_composition(
    sum_df = sum_top,
    out_pdf = file.path(out_dir, "annotation_composition_top_regions.pdf"),
    title = "Annotation composition in top selected regions"
  )
}

# ============================================================
# 14) heatmaps
# ============================================================
if (isTRUE(make_per_module_heatmaps)) {
  for (mod in names(module_data)) {
    mat_df <- build_heatmap_matrix(module_data[[mod]], combine_modules = FALSE)
    write_csv(mat_df, file.path(out_dir, paste0("heatmap_matrix_", mod, ".csv")))

    make_complex_heatmap(
      mat_df = mat_df,
      out_pdf = file.path(out_dir, paste0("heatmap_", mod, ".pdf")),
      title = paste0("Regulatory annotations | ", run_tag, " | ", mod)
    )
  }
}

if (isTRUE(make_combined_heatmap)) {
  mat_df <- build_heatmap_matrix(combined_df, combine_modules = combine_modules)
  write_csv(mat_df, file.path(out_dir, "heatmap_matrix_combined.csv"))

  make_complex_heatmap(
    mat_df = mat_df,
    out_pdf = file.path(out_dir, "heatmap_combined.pdf"),
    title = paste0("Regulatory annotations | ", run_tag, " | combined modules")
  )
}

msg("[✓] 13B complete.")