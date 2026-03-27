#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 10: Circos Plots and Module Summaries
#
# PURPOSE
#   - Load one or more module objects
#   - Save module color counts
#   - Save module eigengenes
#   - Generate one circos plot with all non-grey modules
#   - Generate one circos plot per non-grey module
#   - Generate a module region count barplot
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --modules_v1   : path to v1 Modules.rds
#
# OPTIONAL INPUTS
#   --modules_v2           : optional path to v2 Modules.rds
#   --modules_v3           : optional path to v3 Modules.rds
#   --genome               : genome for circos ideogram [default = hg38]
#   --single_module_color  : color for per-module circos plots [default = red]
#   --all_modules_width    : width for all-modules circos pdf [default = 4.0]
#   --all_modules_height   : height for all-modules circos pdf [default = 4.0]
#   --single_module_width  : width for single-module circos pdf [default = 5.5]
#   --single_module_height : height for single-module circos pdf [default = 5.5]
#   --counts_width         : width for module count pdf [default = 7]
#   --counts_height        : height for module count pdf [default = 6]
#
# OUTPUTS
#   project_root/comethyl_output/10_circos_and_module_summaries/<cpg_label>/<region_label>/<variant>/
#       Module_Color_Counts.csv
#       Module_Eigengenes.csv
#       Circos_ALLModules.pdf
#       Circos_Module_<module>.pdf
#       Module_Region_Counts.pdf
#       run_parameters.txt
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(dplyr)
  library(ggplot2)
  library(circlize)
  library(cowplot)
  library(scales)
})

# ------------------------------------------------------------
# Load helper.R
# ------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_file_arg) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) {
  stop("helper.R not found next to this script: ", helper_file)
}
source(helper_file)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_pdf <- function(path, width = 5.5, height = 5.5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
}

validate_modules_for_circos <- function(x, label = "modules object") {
  x <- validate_modules_object(x, label)

  if (is.null(x$regions) || !(is.data.frame(x$regions) || is.matrix(x$regions))) {
    stop(label, " is missing a valid $regions object.")
  }

  regions <- as.data.frame(x$regions, stringsAsFactors = FALSE)
  required_cols <- c("chr", "start", "end", "module")
  missing_cols <- setdiff(required_cols, colnames(regions))
  if (length(missing_cols) > 0) {
    stop(label, " $regions is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (nrow(regions) == 0) stop(label, " $regions has zero rows.")
  if (anyNA(regions$chr)) stop(label, " $regions$chr contains NA values.")
  if (anyNA(regions$start)) stop(label, " $regions$start contains NA values.")
  if (anyNA(regions$end)) stop(label, " $regions$end contains NA values.")
  if (anyNA(regions$module)) stop(label, " $regions$module contains NA values.")

  x$regions <- regions
  x
}

make_module_color_map <- function(mod_names) {
  mod_names <- as.character(mod_names)
  uniq <- sort(unique(mod_names))
  ok <- uniq %in% grDevices::colors()

  if (all(ok)) {
    setNames(uniq, uniq)
  } else {
    pal <- grDevices::hcl.colors(length(uniq), palette = "Set 3")
    names(pal) <- uniq
    if ("grey" %in% uniq) pal["grey"] <- "grey"
    pal
  }
}

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--modules_v1", type = "character",
              help = "Path to v1 Modules.rds"),

  make_option("--modules_v2", type = "character", default = NULL,
              help = "Optional path to v2 Modules.rds"),

  make_option("--modules_v3", type = "character", default = NULL,
              help = "Optional path to v3 Modules.rds"),

  make_option("--genome", type = "character", default = "hg38",
              help = "Genome for circos ideogram [default = hg38]"),

  make_option("--single_module_color", type = "character", default = "red",
              help = "Color for per-module circos plots [default = red]"),

  make_option("--all_modules_width", type = "double", default = 4.0,
              help = "Width for all-modules circos pdf [default = 4.0]"),

  make_option("--all_modules_height", type = "double", default = 4.0,
              help = "Height for all-modules circos pdf [default = 4.0]"),

  make_option("--single_module_width", type = "double", default = 5.5,
              help = "Width for single-module circos pdf [default = 5.5]"),

  make_option("--single_module_height", type = "double", default = 5.5,
              help = "Height for single-module circos pdf [default = 5.5]"),

  make_option("--counts_width", type = "double", default = 7,
              help = "Width for module count pdf [default = 7]"),

  make_option("--counts_height", type = "double", default = 6,
              help = "Height for module count pdf [default = 6]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$modules_v1)) stop("--modules_v1 is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$modules_v1)) stop("modules_v1 not found: ", opt$modules_v1)
if (!is.null(opt$modules_v2) && !file.exists(opt$modules_v2)) stop("modules_v2 not found: ", opt$modules_v2)
if (!is.null(opt$modules_v3) && !file.exists(opt$modules_v3)) stop("modules_v3 not found: ", opt$modules_v3)

# ------------------------------------------------------------
# Configure cache
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)
WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Derive lineage from modules_v1
# Expected:
#   .../07_module_detection/<cpg_label>/<region_label>/v1_all_pcs/Modules.rds
# ------------------------------------------------------------
v1_variant_dir <- dirname(opt$modules_v1)
v1_region_dir <- dirname(v1_variant_dir)
region_label <- basename(v1_region_dir)
cpg_label <- basename(dirname(v1_region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "10_circos_and_module_summaries")
out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Build variant inputs
# ------------------------------------------------------------
variant_inputs <- list(
  v1_all_pcs = opt$modules_v1
)

if (!is.null(opt$modules_v2)) {
  variant_inputs[[basename(dirname(opt$modules_v2))]] <- opt$modules_v2
}

if (!is.null(opt$modules_v3)) {
  variant_inputs[[basename(dirname(opt$modules_v3))]] <- opt$modules_v3
}

# ------------------------------------------------------------
# Run per variant
# ------------------------------------------------------------
for (variant_name in names(variant_inputs)) {
  message("\n==============================")
  message("Running circos/module summaries for variant: ", variant_name)
  message("==============================\n")

  modules <- validate_modules_for_circos(
    readRDS(variant_inputs[[variant_name]]),
    paste0(variant_name, " modules object")
  )

  MEs <- modules$MEs
  regions <- modules$regions

  variant_out_dir <- file.path(out_dir, variant_name)
  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  # A) Module color counts
  module_counts_df <- as.data.frame(table(modules$colors), stringsAsFactors = FALSE)
  colnames(module_counts_df) <- c("ModuleColor", "Count")
  write.csv(
    module_counts_df,
    file = file.path(variant_out_dir, "Module_Color_Counts.csv"),
    row.names = FALSE
  )

  # B) Eigengenes
  write.csv(
    MEs,
    file = file.path(variant_out_dir, "Module_Eigengenes.csv"),
    row.names = TRUE
  )

  # C) Prepare non-grey regions
  regions_col <- regions %>%
    dplyr::filter(.data$module != "grey") %>%
    dplyr::select(.data$chr, .data$start, .data$end, .data$module)

  module_colors <- NULL
  n_non_grey_modules <- 0L

  if (nrow(regions_col) > 0) {
    module_colors <- make_module_color_map(regions_col$module)
    n_non_grey_modules <- length(unique(as.character(regions_col$module)))

    # D) All-modules circos
    all_pdf <- file.path(variant_out_dir, "Circos_ALLModules.pdf")
    safe_pdf(all_pdf, width = opt$all_modules_width, height = opt$all_modules_height)

    circos.clear()
    circos.par(
      gap.degree = 2,
      cell.padding = c(0.007, 0, 0.007, 0),
      circle.margin = 0.00001
    )
    circos.initializeWithIdeogram(
      species = opt$genome,
      plotType = c("ideogram", "labels")
    )

    circos.genomicTrack(
      regions_col,
      ylim = c(0, 1),
      panel.fun = function(region, value, ...) {
        col_here <- module_colors[as.character(value$module)]
        circos.genomicRect(
          region, value,
          ytop = 1, ybottom = 0,
          col = col_here, border = col_here
        )
      }
    )

    circos.clear()
    dev.off()

    # E) One circos per module
    modules_to_plot <- sort(unique(as.character(regions_col$module)))

    for (m in modules_to_plot) {
      sub <- regions_col %>% dplyr::filter(.data$module == m)
      if (nrow(sub) == 0) next

      chr_idx <- unique(sub$chr)

      per_pdf <- file.path(variant_out_dir, paste0("Circos_Module_", m, ".pdf"))
      safe_pdf(per_pdf, width = opt$single_module_width, height = opt$single_module_height)

      circos.clear()
      circos.par(
        gap.degree = 2,
        cell.padding = c(0.007, 0, 0.007, 0),
        circle.margin = 0.00001
      )

      circos.initializeWithIdeogram(
        species = opt$genome,
        chromosome.index = chr_idx,
        plotType = c("ideogram", "labels")
      )

      circos.genomicTrack(
        sub,
        ylim = c(0, 1),
        track.height = 0.08,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(
            region, value,
            ytop = 1, ybottom = 0,
            col = opt$single_module_color,
            border = opt$single_module_color
          )
        }
      )

      circos.clear()
      dev.off()
    }

    # F) Module region counts barplot
    region_counts <- as.data.frame(sort(table(regions_col$module), decreasing = TRUE))
    colnames(region_counts) <- c("Module", "Regions")
    region_counts$Module <- as.character(region_counts$Module)
    region_counts$Module <- factor(region_counts$Module, levels = rev(region_counts$Module))

    barplot_obj <- ggplot(region_counts, aes(x = Module, y = Regions)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 18) +
      theme(
        legend.position  = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y     = element_blank()
      ) +
      ylab("Regions")

    rowColors <- ggplot(data.frame(
      x = 0,
      y = seq_len(nrow(region_counts)),
      color = rev(as.character(region_counts$Module))
    )) +
      geom_tile(aes(x = x, y = y, fill = color)) +
      scale_fill_identity() +
      theme_void() +
      theme(plot.margin = unit(c(0.1, -1.5, 0.1, 0.5), "lines"))

    gg <- cowplot::plot_grid(rowColors, barplot_obj, ncol = 2, rel_widths = c(0.08, 1))

    counts_pdf <- file.path(variant_out_dir, "Module_Region_Counts.pdf")
    ggsave(
      counts_pdf,
      plot = gg,
      dpi = 600,
      width = opt$counts_width,
      height = opt$counts_height,
      units = "in"
    )
  }

  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("modules_file:", variant_inputs[[variant_name]]),
      paste("genome:", opt$genome),
      paste("single_module_color:", opt$single_module_color),
      paste("n_samples:", nrow(MEs)),
      paste("n_modules:", ncol(MEs)),
      paste("n_regions_total:", nrow(regions)),
      paste("n_regions_non_grey:", nrow(regions_col)),
      paste("n_non_grey_modules:", n_non_grey_modules),
      paste("all_modules_width:", opt$all_modules_width),
      paste("all_modules_height:", opt$all_modules_height),
      paste("single_module_width:", opt$single_module_width),
      paste("single_module_height:", opt$single_module_height),
      paste("counts_width:", opt$counts_width),
      paste("counts_height:", opt$counts_height)
    ),
    file.path(variant_out_dir, "run_parameters.txt")
  )

  message("✓ Finished variant: ", variant_name)
  message("  Outputs in: ", variant_out_dir)
}

# ------------------------------------------------------------
# Save top-level run parameters
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("modules_v1:", opt$modules_v1),
    paste("modules_v2:", ifelse(is.null(opt$modules_v2), "NULL", opt$modules_v2)),
    paste("modules_v3:", ifelse(is.null(opt$modules_v3), "NULL", opt$modules_v3)),
    paste("genome:", opt$genome),
    paste("single_module_color:", opt$single_module_color),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL CIRCOS PLOTS + SUMMARIES COMPLETE ✓")
message("Outputs saved under:\n  ", out_dir)