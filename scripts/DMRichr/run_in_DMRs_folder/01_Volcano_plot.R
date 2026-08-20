#!/usr/bin/env Rscript
## ============================================================================
## 01_Volcano_plot.R
##
## Full volcano plot for DMRichr output, using ALL tested candidate regions
## (background_annotated.xlsx) plus the significance cutoff, so non-
## significant regions are visible behind the DMRs.
##
## Parameterized version of the original script -- same plot design (colors,
## shapes, labels), but input/output paths and thresholds are now CLI
## arguments so this can run unattended per-timepoint in a SLURM array,
## matching 02_enrichment_plot.R's conventions.
##
## USAGE
##   Rscript 01_Volcano_plot.R --background background_annotated.xlsx \
##     --outdir volcano_results --prefix egwg_36-38wk
##
## Required packages: optparse, readxl, ggplot2, ggrepel, dplyr
## ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# ---- CLI arguments ---------------------------------------------------------

option_list <- list(
  make_option(c("-b", "--background"), type = "character", default = NULL,
              help = "Path to background_annotated.xlsx (all tested candidate regions) [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default %default]"),
  make_option(c("--prefix"), type = "character", default = "DMR",
              help = "Prefix for output file names, e.g. 'egwg_36-38wk' -> egwg_36-38wk_volcano_full.pdf [default %default]"),
  make_option(c("--pval_cutoff"), type = "double", default = 0.05,
              help = "P-value cutoff for significance [default %default]"),
  make_option(c("--n_labels"), type = "integer", default = 20,
              help = "Number of top hits (by p-value) to label on the plot [default %default]"),
  make_option(c("--fig_width"), type = "double", default = 8.5,
              help = "Figure width in inches [default %default]"),
  make_option(c("--fig_height"), type = "double", default = 7,
              help = "Figure height in inches [default %default]"),
  make_option(c("--dpi"), type = "integer", default = 300,
              help = "PNG resolution [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$background)) stop("Please provide --background <path to background_annotated.xlsx>")
if (!file.exists(opt$background)) stop("File not found: ", opt$background)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

P_CUTOFF <- opt$pval_cutoff
N_LABELS <- opt$n_labels

# ---- 1. Load data -----------------------------------------------------------

message("Loading: ", opt$background)
background <- read_excel(opt$background)

required_cols <- c("p.value", "difference", "direction", "geneSymbol")
missing_cols <- setdiff(required_cols, names(background))
if (length(missing_cols) > 0) {
  stop("background_annotated.xlsx is missing required column(s): ",
       paste(missing_cols, collapse = ", "))
}

background <- background %>%
  mutate(
    negLog10P = -log10(p.value),
    category = case_when(
      p.value >= P_CUTOFF                                   ~ "NS",
      p.value <  P_CUTOFF & direction == "Hypermethylated"   ~ "UP",
      p.value <  P_CUTOFF & direction == "Hypomethylated"    ~ "DOWN"
    ),
    category = factor(category, levels = c("NS", "DOWN", "UP"))
  )

dir_counts <- table(background$category)

# ---- 2. Identify significant DMRs -------------------------------------------

sig_regions <- background %>% filter(category != "NS")
n_sig_regions <- nrow(sig_regions)
message(n_sig_regions, " significant regions out of ", nrow(background),
        " tested (p < ", P_CUTOFF, ")")

# ---- 3. Identify top hits to label ------------------------------------------

top_hits <- sig_regions %>%
  slice_min(order_by = p.value, n = N_LABELS)

# ---- 4. Build the volcano plot ----------------------------------------------

volcano <- ggplot(background, aes(x = difference, y = negLog10P, color = category)) +
  geom_point(size = 1.3, alpha = 0.85, shape = 16) +   # shape 16 = solid dot, no border
  geom_hline(yintercept = -log10(P_CUTOFF), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    name = "Direction",
    values = c("UP" = "#E63946", "DOWN" = "#1D3D8F", "NS" = "#CCCCCC"),
    labels = c(
      "UP"   = paste0("Hypermethylated (n=", dir_counts["UP"], ")"),
      "DOWN" = paste0("Hypomethylated (n=", dir_counts["DOWN"], ")"),
      "NS"   = paste0("Not significant (n=", dir_counts["NS"], ")")
    )
  ) +
  geom_text_repel(
    data = top_hits,
    aes(label = geneSymbol),
    size = 3, fontface = "italic",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    x = "Methylation difference (%)",
    y = expression(-log[10]~"(p-value)"),
    title = paste0("Differentially Methylated Regions (n = ", n_sig_regions, ")")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",              # Direction legend directly under the title
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15)
  )

# ---- 5. Save -----------------------------------------------------------------

pdf_path <- file.path(opt$outdir, paste0(opt$prefix, "_volcano_full.pdf"))
png_path <- file.path(opt$outdir, paste0(opt$prefix, "_volcano_full.png"))

ggsave(pdf_path, volcano, width = opt$fig_width, height = opt$fig_height)
ggsave(png_path, volcano, width = opt$fig_width, height = opt$fig_height, dpi = opt$dpi)

message("Saved: ", pdf_path)
message("Saved: ", png_path)
message("Done.")