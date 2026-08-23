#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 18: Consensus Regulatory State Presentation
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Presentation and reporting layer for regulatory-state overlap
#   results produced by script 17:
#     1) Load dominant-state long tables from script 17 per module
#     2) Optionally filter to specific regulatory sources
#     3) Limit to top regions by module membership or overlap size
#     4) Generate per-module regulatory state heatmaps
#     5) Generate a combined multi-module heatmap
#     6) Generate stacked barplots of annotation composition
#     7) Save annotation summary tables
#
# NOTE
#   Module assignments are shared — presentation runs once per
#   adjustment version. No per-dataset loop needed.
#
# REQUIRED INPUTS
#   --project_root           : root directory of the analysis project
#   --regulatory_overlap_dir : path to script 17 output directory
#                              for this adjustment_version
#   --modules                : comma-separated module names to visualize
#                              (e.g. turquoise,blue,brown)
#
# OPTIONAL INPUTS
#   --adjustment_version       : label for output directory [default = unadjusted]
#   --source                   : comma-separated source filter (roadmap,roadmap18,encode)
#                                [default = all sources]
#   --max_regions              : maximum regions per module to display [default = 50]
#   --max_genes                : maximum genes per module to display [default = all]
#   --combine_modules          : include module label in combined heatmap TRUE/FALSE [default = TRUE]
#   --make_per_module_heatmaps : generate individual module heatmaps TRUE/FALSE [default = TRUE]
#   --make_combined_heatmap    : generate combined heatmap TRUE/FALSE [default = TRUE]
#   --make_barplots            : generate composition barplots TRUE/FALSE [default = TRUE]
#
# OUTPUTS
#   comethyl_output/consensus/18_regulatory_presentation/<adjustment_version>/
#       annotation_summary_all.csv
#       annotation_summary_top_regions.csv
#       annotation_composition_all.pdf
#       annotation_composition_top_regions.pdf
#       heatmap_matrix_<module>.csv
#       heatmap_<module>.pdf
#       heatmap_matrix_combined.csv
#       heatmap_combined.pdf
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - Regions are ranked by module membership score first, then
#     overlap base pairs, then region ID
#   - Roadmap 15-state colors follow the standard ChromHMM palette
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/18_plot_chromatin_heatmaps.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --regulatory_overlap_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/17_regulatory_overlap/v1_all_pcs \
#     --modules turquoise,blue,brown \
#     --adjustment_version v1_all_pcs \
#     --source roadmap \
#     --max_regions 50 \
#     --make_per_module_heatmaps TRUE \
#     --make_combined_heatmap TRUE \
#     --make_barplots TRUE
#
# PURPOSE
#   Presentation/reporting layer for regulatory-state overlap
#   results produced by script 17.
#
#   Reads script 17 outputs and creates:
#     1) Per-module heatmaps
#     2) Combined multi-module heatmap
#     3) Per-source annotation summaries
#     4) Stacked barplots of annotation composition
#
# NOTE
#   Module assignments are shared — presentation runs once per
#   adjustment version. No per-dataset loop needed.
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/18_regulatory_presentation/<adjustment_version>/
#       annotation_summary_all.csv
#       annotation_summary_top_regions.csv
#       annotation_composition_all.pdf
#       annotation_composition_top_regions.pdf
#       heatmap_matrix_<module>.csv
#       heatmap_<module>.pdf
#       heatmap_matrix_combined.csv
#       heatmap_combined.pdf
#       run_parameters.txt
#       sessionInfo.txt
# ================================================================
message("Starting Script 18 ")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ============================================================
# Helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- match(flag, args)
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

msg <- function(...) message(sprintf(...))

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

# ============================================================
# Annotation colors
# ============================================================
ROADMAP15_DESC_COLORS <- c(
  "Active TSS"                     = "#FF0000",
  "Flanking Active TSS"            = "#FF4500",
  "Transcribed at 5' and 3'"       = "#32CD32",
  "Strong transcription"           = "#008000",
  "Weak transcription"             = "#006400",
  "Genic enhancer"                 = "#C2E105",
  "Enhancer"                       = "#FFFF00",
  "ZNF genes and repeats"          = "#66CDAA",
  "Heterochromatin"                = "#8A91D0",
  "Bivalent TSS"                   = "#CD5C5C",
  "Flanking bivalent TSS/enhancer" = "#E9967A",
  "Bivalent enhancer"              = "#BDB76B",
  "Repressed PolyComb"             = "#808080",
  "Weak repressed PolyComb"        = "#C0C0C0",
  "Quiescent/Low"                  = "#F7F7F7",
  "Missing"                        = "#EDEDED"
)

DEFAULT_OTHER_COL <- "#D0D0D0"

# ============================================================
# Loaders
# ============================================================
read_module_long <- function(module_dir, mod) {
  # FIX: script 17 writes dominant_state_long_all_sources.csv
  f <- file.path(module_dir, "dominant_state_long_all_sources.csv")
  if (!file.exists(f))
    stop("Missing script 17 output: ", f)

  df <- read_csv(f, show_col_types = FALSE)

  need    <- c("module", "source", "dataset", "track_id", "track_label",
               "region_id", "annotation", "annotation_desc")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0)
    stop("Missing columns in ", basename(f), ": ",
         paste(missing, collapse = ", "))

  if (!"gene_symbol" %in% names(df)) df$gene_symbol <- NA_character_
  if (!"membership"  %in% names(df)) df$membership  <- NA_real_
  if (!"overlap_bp"  %in% names(df)) df$overlap_bp  <- NA_real_

  df %>%
    mutate(
      module          = as.character(module),
      source          = as.character(source),
      dataset         = as.character(dataset),
      track_id        = as.character(track_id),
      track_label     = as.character(track_label),
      region_id       = as.character(region_id),
      gene_symbol     = as.character(gene_symbol),
      annotation      = as.character(annotation),
      annotation_desc = as.character(annotation_desc),
      membership      = suppressWarnings(as.numeric(membership)),
      overlap_bp      = suppressWarnings(as.numeric(overlap_bp))
    )
}

# ============================================================
# Ranking helpers
# ============================================================
rank_regions <- function(df) {
  df %>%
    group_by(module, region_id) %>%
    summarise(
      gene_symbol = {
        gs <- unique(gene_symbol[!is.na(gene_symbol) & gene_symbol != ""])
        if (length(gs) == 0) NA_character_ else gs[1]
      },
      membership = if (all(is.na(membership))) NA_real_
                   else max(membership, na.rm = TRUE),
      overlap_bp = if (all(is.na(overlap_bp))) NA_real_
                   else max(overlap_bp, na.rm = TRUE),
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
      max_membership = if (all(is.na(membership))) NA_real_
                       else max(membership, na.rm = TRUE),
      max_overlap_bp = if (all(is.na(overlap_bp))) NA_real_
                       else max(overlap_bp, na.rm = TRUE),
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
# Matrix builder
# ============================================================
build_heatmap_matrix <- function(df, combine_modules = FALSE) {
  x <- df %>%
    mutate(
      row_label = case_when(
        !is.na(gene_symbol) & gene_symbol != "" ~
          paste0(region_id, " | ", gene_symbol),
        TRUE ~ region_id
      ),
      col_label = case_when(
        combine_modules ~ paste(module, source, track_label, sep = " | "),
        TRUE            ~ paste(source, track_label, sep = " | ")
      )
    )

  if (combine_modules)
    x <- x %>% mutate(row_label = paste(module, row_label, sep = " :: "))

  x %>%
    select(row_label, col_label, annotation_desc) %>%
    distinct() %>%
    pivot_wider(names_from = col_label, values_from = annotation_desc)
}

# ============================================================
# ComplexHeatmap
# ============================================================
make_complex_heatmap <- function(mat_df, out_pdf, title = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("ComplexHeatmap not installed.")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("circlize not installed.")

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
  })

  if (nrow(mat_df) == 0) {
    grDevices::pdf(out_pdf, width = 10, height = 7)
    grid::grid.newpage()
    grid::grid.text("No rows to plot",
                    gp = grid::gpar(fontsize = 14, fontface = "bold"))
    grDevices::dev.off()
    return(invisible(out_pdf))
  }

  row_labels <- mat_df[[1]]
  m          <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(m) <- row_labels

  m2 <- m
  m2[is.na(m2) | m2 == ""] <- "Missing"

  present <- sort(unique(as.vector(m2)))
  colors  <- ROADMAP15_DESC_COLORS
  unknown <- setdiff(present, names(colors))
  if (length(unknown) > 0)
    colors <- c(colors,
                setNames(rep(DEFAULT_OTHER_COL, length(unknown)), unknown))

  grDevices::pdf(out_pdf, width = 14, height = 9)
  ht <- ComplexHeatmap::Heatmap(
    m2,
    name               = "Annotation",
    col                = colors,
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    column_names_rot   = 45,
    row_names_gp       = grid::gpar(fontsize = 8),
    column_names_gp    = grid::gpar(fontsize = 9),
    rect_gp            = grid::gpar(col = "grey70", lwd = 0.5),
    column_title       = title,
    column_title_gp    = grid::gpar(fontsize = 12, fontface = "bold")
  )
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  grDevices::dev.off()
  invisible(out_pdf)
}

# ============================================================
# Summary tables
# ============================================================
make_annotation_summary <- function(df) {
  df %>%
    group_by(module, source, annotation_desc) %>%
    summarise(
      n_regions = n_distinct(region_id),
      n_tracks  = n_distinct(track_id),
      .groups   = "drop"
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
    summarise(n_regions = n_distinct(region_id), .groups = "drop") %>%
    group_by(module, source) %>%
    mutate(frac_regions = n_regions / sum(n_regions)) %>%
    ungroup() %>%
    arrange(module, source, desc(n_regions), annotation_desc)
}

# ============================================================
# Stacked barplot
# ============================================================
plot_annotation_composition <- function(sum_df, out_pdf,
                                         title = "Annotation composition") {
  if (nrow(sum_df) == 0) return(invisible(NULL))

  color_map <- ROADMAP15_DESC_COLORS
  unknown   <- setdiff(unique(sum_df$annotation_desc), names(color_map))
  if (length(unknown) > 0)
    color_map <- c(color_map,
                   setNames(rep(DEFAULT_OTHER_COL, length(unknown)), unknown))

  p <- ggplot(sum_df,
              aes(x = module, y = frac_regions, fill = annotation_desc)) +
    geom_col(position = "stack") +
    facet_wrap(~source, scales = "free_x") +
    scale_fill_manual(values = color_map) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    ) +
    labs(title = title, x = "Module",
         y = "Fraction of regions", fill = "Annotation")

  tryCatch(
    ggsave(out_pdf, plot = p, width = 12, height = 8,
           units = "in", dpi = 300),
    error = function(e)
      msg("[WARN] Barplot failed: %s", conditionMessage(e))
  )
  invisible(out_pdf)
}

# ============================================================
# Source filter
# ============================================================
filter_sources <- function(df, source_arg = "") {
  if (is.null(source_arg) || is.na(source_arg) || !nzchar(source_arg))
    return(df)
  df %>% filter(source %in% split_csv(source_arg))
}

# ============================================================
# Read arguments
# ============================================================
project_root          <- get_arg("--project_root",    "")
regulatory_overlap_dir <- get_arg("--regulatory_overlap_dir", "")
adjustment_version    <- get_arg("--adjustment_version", "unadjusted")
modules_in            <- unique(split_csv(get_arg("--modules", "")))

source_arg            <- get_arg("--source",               "")
max_regions           <- suppressWarnings(
                           as.integer(get_arg("--max_regions", "50")))
max_genes             <- suppressWarnings(
                           as.integer(get_arg("--max_genes",   NA_character_)))
combine_modules         <- as_bool(get_arg("--combine_modules",          "true"), TRUE)
make_per_module_heatmaps <- as_bool(get_arg("--make_per_module_heatmaps", "true"), TRUE)
make_combined_heatmap   <- as_bool(get_arg("--make_combined_heatmap",    "true"), TRUE)
make_barplots           <- as_bool(get_arg("--make_barplots",            "true"), TRUE)

# ============================================================
# Validate
# ============================================================
if (!nzchar(project_root))
  stop("--project_root is required")
if (!nzchar(regulatory_overlap_dir))
  stop("--regulatory_overlap_dir is required — path to script 17 output for this adjustment_version")
if (length(modules_in) == 0)
  stop("--modules is required (comma-separated module names, e.g. turquoise,blue)")
if (!dir.exists(project_root))
  stop("project_root not found: ", project_root)
if (!dir.exists(regulatory_overlap_dir))
  stop("regulatory_overlap_dir not found: ", regulatory_overlap_dir)

# ============================================================
# Output directory
# ============================================================
pipeline_root <- file.path(project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "18_regulatory_presentation",
                            adjustment_version)
safe_dir_create(step_dir)
message("Output directory: ", step_dir)

msg("regulatory_overlap_dir: %s", regulatory_overlap_dir)
msg("adjustment_version:     %s", adjustment_version)
msg("modules:                %s", paste(modules_in, collapse = ", "))
msg("source filter:          %s", ifelse(nzchar(source_arg), source_arg, "(all)"))
msg("max_regions:            %d", max_regions)
msg("max_genes:              %s", ifelse(is.na(max_genes), "all", max_genes))

# ============================================================
# Load module data from script 17 outputs
# ============================================================
module_data <- list()

for (mod in modules_in) {
  # FIX: script 17 writes to step_dir/module_<mod>/
  module_dir <- file.path(regulatory_overlap_dir, paste0("module_", mod))

  if (!dir.exists(module_dir)) {
    msg("[WARN] Module dir not found, skipping: %s", module_dir)
    next
  }

  df <- tryCatch(
    read_module_long(module_dir, mod),
    error = function(e) {
      msg("[WARN] Could not load module %s: %s", mod, conditionMessage(e))
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) {
    msg("[WARN] No data for module %s after loading", mod)
    next
  }

  df <- filter_sources(df, source_arg)
  if (nrow(df) == 0) {
    msg("[WARN] No rows after source filtering for module %s", mod)
    next
  }

  df <- limit_to_top_genes(df,    max_genes   = max_genes)
  df <- limit_to_top_regions(df,  max_regions = max_regions)

  module_data[[mod]] <- df
  msg("[OK] Loaded module %s: %d rows", mod, nrow(df))
}

if (length(module_data) == 0)
  stop("No module data available after loading/filtering. ",
       "Check --modules and --regulatory_overlap_dir.")

combined_df <- bind_rows(module_data)

# ============================================================
# Summary tables
# ============================================================
sum_all <- make_annotation_summary(combined_df)
write_csv(sum_all,
          file.path(step_dir, "annotation_summary_all.csv"))
message("Saved: annotation_summary_all.csv")

sum_top <- make_top_annotation_summary(combined_df)
write_csv(sum_top,
          file.path(step_dir, "annotation_summary_top_regions.csv"))
message("Saved: annotation_summary_top_regions.csv")

# ============================================================
# Barplots
# ============================================================
if (isTRUE(make_barplots)) {
  plot_annotation_composition(
    sum_df  = sum_all,
    out_pdf = file.path(step_dir, "annotation_composition_all.pdf"),
    title   = paste0("Annotation composition | ", adjustment_version)
  )
  message("Saved: annotation_composition_all.pdf")

  plot_annotation_composition(
    sum_df  = sum_top,
    out_pdf = file.path(step_dir,
                        "annotation_composition_top_regions.pdf"),
    title   = paste0("Annotation composition (top regions) | ",
                     adjustment_version)
  )
  message("Saved: annotation_composition_top_regions.pdf")
}

# ============================================================
# Per-module heatmaps
# ============================================================
if (isTRUE(make_per_module_heatmaps)) {
  for (mod in names(module_data)) {
    mat_df <- build_heatmap_matrix(module_data[[mod]],
                                   combine_modules = FALSE)
    write_csv(mat_df,
              file.path(step_dir,
                        paste0("heatmap_matrix_", mod, ".csv")))

    tryCatch(
      make_complex_heatmap(
        mat_df  = mat_df,
        out_pdf = file.path(step_dir, paste0("heatmap_", mod, ".pdf")),
        title   = paste0("Regulatory annotations | ",
                         adjustment_version, " | ", mod)
      ),
      error = function(e)
        msg("[WARN] Heatmap failed for module %s: %s",
            mod, conditionMessage(e))
    )
    message("Saved: heatmap_", mod, ".pdf")
  }
}

# ============================================================
# Combined heatmap
# ============================================================
if (isTRUE(make_combined_heatmap)) {
  mat_df <- build_heatmap_matrix(combined_df,
                                  combine_modules = combine_modules)
  write_csv(mat_df,
            file.path(step_dir, "heatmap_matrix_combined.csv"))

  tryCatch(
    make_complex_heatmap(
      mat_df  = mat_df,
      out_pdf = file.path(step_dir, "heatmap_combined.pdf"),
      title   = paste0("Regulatory annotations | ",
                       adjustment_version, " | combined modules")
    ),
    error = function(e)
      msg("[WARN] Combined heatmap failed: %s", conditionMessage(e))
  )
  message("Saved: heatmap_combined.pdf")
}

# ============================================================
# Logs
# ============================================================
write_lines_safe(
  c(
    paste0("project_root\t",            project_root),
    paste0("regulatory_overlap_dir\t",  regulatory_overlap_dir),
    paste0("adjustment_version\t",      adjustment_version),
    paste0("modules\t",                 paste(modules_in, collapse = ", ")),
    paste0("modules_loaded\t",          paste(names(module_data), collapse = ", ")),
    paste0("source_filter\t",           ifelse(nzchar(source_arg), source_arg, "all")),
    paste0("max_regions\t",             max_regions),
    paste0("max_genes\t",               ifelse(is.na(max_genes), "all", max_genes)),
    paste0("combine_modules\t",         combine_modules),
    paste0("make_per_module_heatmaps\t",make_per_module_heatmaps),
    paste0("make_combined_heatmap\t",   make_combined_heatmap),
    paste0("make_barplots\t",           make_barplots),
    paste0("date\t",                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ),
  file.path(step_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()),
           con = file.path(step_dir, "sessionInfo.txt"))

msg("Script 18 complete: consensus regulatory state presentation finished")
message("Outputs saved under: ", step_dir)