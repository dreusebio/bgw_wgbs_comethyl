#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 13: Consensus Circos Plots and Module Summaries
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Visualizes the genomic distribution of consensus modules:
#     1) Load consensus region assignments from script 09 shared/
#     2) Save module color counts table
#     3) Generate one circos plot showing all non-grey modules
#     4) Generate individual circos plots per non-grey module
#     5) Generate a module region count barplot
#
# NOTE
#   Module assignments are shared across datasets in consensus
#   analysis — circos plots run once per adjustment version only.
#   No per-dataset loop is needed.
#
# REQUIRED INPUTS
#   --project_root          : root directory of the analysis project
#   --consensus_modules_rds : path to Consensus_Modules.rds from script 09
#   --consensus_regions_tsv : path to Consensus_Region_Assignments.tsv
#                             from script 09 shared/
#   --regions_file          : path to Filtered_Regions.txt from script 02
#                             (used to derive cpg_label and region_label)
#
# OPTIONAL INPUTS
#   --adjustment_version    : label matching script 09 run [default = unadjusted]
#   --genome                : genome assembly (hg38) [default = hg38]
#   --single_module_color   : color for per-module circos tracks [default = red]
#   --all_modules_width/height   : all-modules circos dimensions [default = 4.0 x 4.0]
#   --single_module_width/height : per-module circos dimensions [default = 5.5 x 5.5]
#   --counts_width/height        : barplot dimensions [default = 7.0 x 6.0]
#
# OUTPUTS
#   comethyl_output/consensus/13_circos_and_module_summaries/<cpg_label>/<region_label>/<adjustment_version>/
#       Module_Color_Counts.csv
#       Consensus_Region_Assignments.tsv   (copy from script 09 shared/)
#       Circos_ALLModules.pdf
#       Circos_Module_<module>.pdf         (one per non-grey module)
#       Module_Region_Counts.pdf
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/13_circos_plots_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --consensus_regions_tsv /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Region_Assignments.tsv \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --adjustment_version v1_all_pcs \
#     --genome hg38
# ================================================================
# SCRIPT 13: Consensus Circos Plots and Module Summaries
#
# PURPOSE
#   - Load Consensus_Modules.rds from script 09
#   - Load consensus region assignments from shared TSV
#   - Save module color counts
#   - Generate one circos plot with all non-grey modules
#   - Generate one circos plot per non-grey module
#   - Generate a module region count barplot
#
# NOTE
#   Module assignments are shared across all datasets in a consensus
#   analysis — circos plots and region summaries do not vary by dataset.
#   Outputs are organised by adjustment_version only.
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/13_circos_and_module_summaries/<adjustment_version>/
#       Module_Color_Counts.csv
#       Consensus_Region_Assignments.tsv   (copy from script 09 shared/)
#       Circos_ALLModules.pdf
#       Circos_Module_<module>.pdf
#       Module_Region_Counts.pdf
#       run_parameters.txt
#       sessionInfo.txt
# ================================================================
message("Starting Script 13 ✓")

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
script_file_arg <- commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
]
if (length(script_file_arg) == 0) stop("Could not determine script path.")
script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found: ", helper_file)
source(helper_file)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_pdf <- function(path, width = 5.5, height = 5.5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
}

make_module_color_map <- function(mod_names) {
  mod_names <- as.character(mod_names)
  uniq <- sort(unique(mod_names))
  ok   <- uniq %in% grDevices::colors()
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
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",            type = "character"),
  make_option("--consensus_modules_rds",   type = "character"),
  make_option("--consensus_regions_tsv",   type = "character",
              help = "Consensus_Region_Assignments.tsv from script 09 shared/"),
  make_option("--regions_file",            type = "character"),
  make_option("--adjustment_version",      type = "character", default = "unadjusted"),
  make_option("--genome",                  type = "character", default = "hg38"),
  make_option("--single_module_color",     type = "character", default = "red"),
  make_option("--all_modules_width",       type = "double",    default = 4.0),
  make_option("--all_modules_height",      type = "double",    default = 4.0),
  make_option("--single_module_width",     type = "double",    default = 5.5),
  make_option("--single_module_height",    type = "double",    default = 5.5),
  make_option("--counts_width",            type = "double",    default = 7.0),
  make_option("--counts_height",           type = "double",    default = 6.0)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root))          stop("--project_root is required")
if (is.null(opt$consensus_modules_rds)) stop("--consensus_modules_rds is required")
if (is.null(opt$consensus_regions_tsv)) stop("--consensus_regions_tsv is required")
if (is.null(opt$regions_file))          stop("--regions_file is required")
if (!dir.exists(opt$project_root))      stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
if (!file.exists(opt$consensus_regions_tsv)) stop("consensus_regions_tsv not found: ", opt$consensus_regions_tsv)
if (!file.exists(opt$regions_file))     stop("regions_file not found: ", opt$regions_file)

# ------------------------------------------------------------
# Configure
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption("CACHE",
  file.path(opt$project_root, ".cache"))
WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Output dir
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
step_dir      <- file.path(pipeline_root, "13_circos_and_module_summaries",
                            cpg_label, region_label, opt$adjustment_version)
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", step_dir)

# ------------------------------------------------------------
# Load consensus modules and regions
# ------------------------------------------------------------
consensusMods <- readRDS(opt$consensus_modules_rds)

all_colors    <- consensusMods$colors
module_colors_real <- unique(all_colors[all_colors != "grey"])
n_modules     <- length(module_colors_real)
message("Non-grey consensus modules: ", n_modules,
        " (", paste(sort(module_colors_real), collapse = ", "), ")")

regions <- read.delim(opt$consensus_regions_tsv,
                      stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("RegionID", "chr", "start", "end", "module")
missing_cols  <- setdiff(required_cols, colnames(regions))
if (length(missing_cols) > 0)
  stop("consensus_regions_tsv missing columns: ",
       paste(missing_cols, collapse = ", "))

message("Regions loaded: ", nrow(regions),
        " total, ",
        sum(regions$module != "grey"), " non-grey")

# ------------------------------------------------------------
# A) Module color counts
# ------------------------------------------------------------
module_counts_df <- as.data.frame(
  sort(table(regions$module), decreasing = TRUE),
  stringsAsFactors = FALSE
)
colnames(module_counts_df) <- c("ModuleColor", "Count")

write.csv(module_counts_df,
          file      = file.path(step_dir, "Module_Color_Counts.csv"),
          row.names = FALSE)
message("Saved: Module_Color_Counts.csv")

# ------------------------------------------------------------
# B) Copy region assignments for convenience
# ------------------------------------------------------------
file.copy(opt$consensus_regions_tsv,
          file.path(step_dir, "Consensus_Region_Assignments.tsv"),
          overwrite = TRUE)
message("Copied: Consensus_Region_Assignments.tsv")

# ------------------------------------------------------------
# C) Prepare non-grey regions for circos
# ------------------------------------------------------------
regions_col <- regions %>%
  dplyr::filter(.data$module != "grey") %>%
  dplyr::select(.data$chr, .data$start, .data$end, .data$module)

if (nrow(regions_col) == 0) {
  message("WARNING: No non-grey regions found. Skipping all circos plots.")
} else {

  color_map <- make_module_color_map(regions_col$module)

  # ── D) All-modules circos ─────────────────────────────────
  all_pdf <- file.path(step_dir, "Circos_ALLModules.pdf")
  tryCatch({
    safe_pdf(all_pdf, width = opt$all_modules_width,
                      height = opt$all_modules_height)
    circos.clear()
    circos.par(
      gap.degree     = 2,
      cell.padding   = c(0.007, 0, 0.007, 0),
      circle.margin  = 0.00001
    )
    circos.initializeWithIdeogram(
      species  = opt$genome,
      plotType = c("ideogram", "labels")
    )
    circos.genomicTrack(
      regions_col,
      ylim = c(0, 1),
      panel.fun = function(region, value, ...) {
        col_here <- color_map[as.character(value$module)]
        circos.genomicRect(region, value,
                           ytop = 1, ybottom = 0,
                           col = col_here, border = col_here)
      }
    )
    circos.clear()
    dev.off()
    message("Saved: Circos_ALLModules.pdf")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message("Circos_ALLModules failed: ", conditionMessage(e))
  })

  # ── E) Per-module circos ──────────────────────────────────
  modules_to_plot <- sort(unique(as.character(regions_col$module)))

  for (m in modules_to_plot) {
    sub <- dplyr::filter(regions_col, .data$module == m)
    if (nrow(sub) == 0) next

    chr_idx  <- unique(sub$chr)
    per_pdf  <- file.path(step_dir, paste0("Circos_Module_", m, ".pdf"))

    tryCatch({
      safe_pdf(per_pdf, width  = opt$single_module_width,
                        height = opt$single_module_height)
      circos.clear()
      circos.par(
        gap.degree    = 2,
        cell.padding  = c(0.007, 0, 0.007, 0),
        circle.margin = 0.00001
      )
      circos.initializeWithIdeogram(
        species           = opt$genome,
        chromosome.index  = chr_idx,
        plotType          = c("ideogram", "labels")
      )
      circos.genomicTrack(
        sub,
        ylim         = c(0, 1),
        track.height = 0.08,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(region, value,
                             ytop = 1, ybottom = 0,
                             col    = opt$single_module_color,
                             border = opt$single_module_color)
        }
      )
      circos.clear()
      dev.off()
      message("Saved: Circos_Module_", m, ".pdf")
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message("Circos_Module_", m, " failed: ", conditionMessage(e))
    })
  }

  # ── F) Module region counts barplot ───────────────────────
  region_counts <- as.data.frame(
    sort(table(regions_col$module), decreasing = TRUE)
  )
  colnames(region_counts) <- c("Module", "Regions")
  region_counts$Module <- factor(as.character(region_counts$Module),
                                 levels = rev(as.character(region_counts$Module)))

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
    x     = 0,
    y     = seq_len(nrow(region_counts)),
    color = rev(as.character(region_counts$Module))
  )) +
    geom_tile(aes(x = x, y = y, fill = color)) +
    scale_fill_identity() +
    theme_void() +
    theme(plot.margin = unit(c(0.1, -1.5, 0.1, 0.5), "lines"))

  gg <- cowplot::plot_grid(rowColors, barplot_obj,
                           ncol = 2, rel_widths = c(0.08, 1))

  counts_pdf <- file.path(step_dir, "Module_Region_Counts.pdf")
  tryCatch({
    ggsave(counts_pdf, plot = gg, dpi = 600,
           width = opt$counts_width, height = opt$counts_height, units = "in")
    message("Saved: Module_Region_Counts.pdf")
  }, error = function(e)
    message("Module_Region_Counts barplot failed: ", conditionMessage(e)))
}

# ------------------------------------------------------------
# Logs
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:",           opt$project_root),
    paste("consensus_modules_rds:",  opt$consensus_modules_rds),
    paste("consensus_regions_tsv:",  opt$consensus_regions_tsv),
    paste("adjustment_version:",     opt$adjustment_version),
    paste("genome:",                 opt$genome),
    paste("single_module_color:",    opt$single_module_color),
    paste("n_real_modules:",         n_modules),
    paste("module_colors:",          paste(sort(module_colors_real), collapse = ", ")),
    paste("n_regions_total:",        nrow(regions)),
    paste("n_regions_non_grey:",     sum(regions$module != "grey")),
    paste("all_modules_width:",      opt$all_modules_width),
    paste("all_modules_height:",     opt$all_modules_height),
    paste("single_module_width:",    opt$single_module_width),
    paste("single_module_height:",   opt$single_module_height),
    paste("counts_width:",           opt$counts_width),
    paste("counts_height:",          opt$counts_height),
    paste("date:",                   as.character(Sys.time()))
  ),
  file.path(step_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()),
           con = file.path(step_dir, "sessionInfo.txt"))

message("Script 13 complete: consensus circos plots and module summaries finished")
message("Outputs saved under: ", step_dir)