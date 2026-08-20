## ---------------------------------------------------------------------------
## Full volcano plot for DMRichr output, using ALL tested candidate regions
## (background_annotated.xlsx) plus the significance cutoff, so non-
## significant regions are visible behind the DMRs.
## ---------------------------------------------------------------------------
## Drop this in after your background/candidate regions object is written
## out in the DMRichr script (e.g. right after `background_annotated.xlsx`
## / the pre-filter regions object is generated).
## ---------------------------------------------------------------------------

library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- 1. Load data -----------------------------------------------------------
# From xlsx:
background <- read_excel("background_annotated.xlsx")

# If reading directly from the object already in memory in your DMRichr
# script (before filtering to sigRegions), skip read_excel and just use
# that data frame/GRanges (as.data.frame(regions)) instead.

P_CUTOFF <- 0.05

background <- background %>%
  mutate(
    negLog10P = -log10(p.value),
    category = case_when(
      p.value >= P_CUTOFF                      ~ "NS",
      p.value <  P_CUTOFF & direction == "Hypermethylated" ~ "UP",
      p.value <  P_CUTOFF & direction == "Hypomethylated"  ~ "DOWN"
    ),
    category = factor(category, levels = c("NS", "DOWN", "UP"))
  )

dir_counts <- table(background$category)

# ---- 2. Identify significant DMRs -------------------------------------------
sig_regions <- background %>% filter(category != "NS")
n_sig_regions <- nrow(sig_regions)

# ---- 3. Identify top hits to label ------------------------------------------
n_labels <- 20
top_hits <- sig_regions %>%
  slice_min(order_by = p.value, n = n_labels)

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

print(volcano)

# ---- 4. Save ------------------------------------------------------------
ggsave("DMR_volcano_plot_full.pdf", volcano, width = 8.5, height = 7)
ggsave("DMR_volcano_plot_full.png", volcano, width = 8.5, height = 7, dpi = 300)