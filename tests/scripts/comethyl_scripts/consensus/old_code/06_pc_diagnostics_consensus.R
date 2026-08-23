#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 06: PC Derivation and PC-Trait Diagnostics for Consensus
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   For a single dataset (run separately per dataset):
#     1) Load region-level methylation matrix (from script 05b or 03)
#     2) Load sample information file
#     3) Align samples between methylation and sample info
#     4) Remove non-numeric, all-NA, and zero-variance trait columns
#     5) Optionally resolve protected and technical traits
#     6) Build SVA model matrix and estimate number of surrogate
#        variables via sva::num.sv()
#     7) Derive PCs via SVD and save selected PCs for adjustment
#     8) Compute PC-trait correlation statistics (Pearson and bicor)
#     9) Generate diagnostic heatmaps and trait dendrogram

# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --input_meth     : path to Region_Methylation.rds from script 03,
#                      05, or 05b (recommended: 05b complete matrix)
#   --sample_info    : path to sample info file (.xlsx, .csv, .tsv, .txt)
#
# OPTIONAL INPUTS
#   --sample_id_col          : column in sample_info with sample IDs
#                              [default = use existing rownames]
#   --protected_traits_file  : text file, one protected trait per line
#   --technical_traits_file  : text file, one technical trait per line
#   --threads                : WGCNA threads [default = 4]
#   --top_trait_n            : top N traits for overview heatmap [default = 50]
#   --top_cell_n             : top N PC-trait cells for heatmap [default = 50]
#   --primary_cor_method     : bicor or pearson [default = bicor]
#   --selected_trait_mode    : protected | protected_and_technical | none
#                              [default = auto-detected from protected traits]
#   --plot_width             : plot width in inches [default = 12]
#   --plot_height            : plot height in inches [default = 10]
#   --plot_dpi               : plot DPI [default = 600]
#   --base_size              : base font size [default = 11]
#   --axis_text_size         : axis label font size [default = 8]
#
# OUTPUTS
#   comethyl_output/consensus/06_pc_diagnostics/<dataset_label>/<cpg_label>/<region_label>/
#       PCs.rds
#       PC_Metadata.tsv
#       Variance_Explained_First20_PCs.pdf
#       PC_Trait_Correlation_Stats_Pearson.tsv
#       PC_Trait_Correlation_Stats_Bicor.tsv
#       PC_Trait_Bicor_Heatmap_PUBLICATION.pdf
#       PC_Trait_Bicor_Heatmap_Top<N>Cells.pdf
#       Trait_Dendrogram_From_PC_Associations.pdf
#       sample_info_non_numeric_columns_removed.txt
#       sample_info_all_na_numeric_columns_removed.txt
#       sample_info_zero_variance_numeric_columns_removed.txt
#       traits_requested/found/missing_protected.txt (if provided)
#       traits_requested/found/missing_technical.txt (if provided)
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - Run this script once per dataset (Baseline, 36_38wks, Postpartum)
#   - dataset_label, cpg_label, and region_label are derived automatically
#     from the --input_meth file path; no manual labelling required
#   - Use the saved PCs.rds as --input_pcs in script 07
#   - Use PC_Trait_Correlation_Stats_Bicor.tsv as --pc_trait_stats in script 07
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/06_pc_diagnostics_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin4_methSD0p08/Baseline_Methylation_complete.rds \
#     --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/Baseline_sample_info.xlsx \
#     --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/protected_traits.txt \
#     --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/technical_traits.txt \
#     --primary_cor_method bicor \
#     --threads 4
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/06_pc_diagnostics/
#     <dataset_label>/
#       <cpg_label>/
#         <region_label>/
#           PCs.rds
#           PC_Metadata.tsv
#           Variance_Explained_First20_PCs.pdf
#           PC_Trait_Correlation_Stats_Pearson.tsv
#           PC_Trait_Correlation_Stats_Bicor.tsv
#           PC_Trait_Bicor_Heatmap_PUBLICATION.pdf
#           PC_Trait_Bicor_Heatmap_Top<N>Cells.pdf
#           Trait_Dendrogram_From_PC_Associations.pdf
#           run_parameters.txt
#           sessionInfo.txt
#
# PATH DERIVATION
#   dataset_label, cpg_label, region_label are derived automatically
#   from --input_meth file path:
#
#   From 05b output (recommended for consensus):
#     .../05b_shared_complete_regions/<dataset_label>/<cpg_label>/<region_label>/
#         <dataset_label>_Methylation_complete.rds
#
#   From 07 adjusted output:
#     .../07_methylation_adjustment/<dataset_label>/<cpg_label>/<region_label>/
#         <variant>/Adjusted_Region_Methylation.rds
#
#   From 03 raw output:
#     .../03_region_methylation/<dataset_label>/<cpg_label>/<region_label>/
#         Region_Methylation.rds
# ================================================================

message("Starting Script 06")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(sva)
})

# ----------------------------------------------------------------
# Load helper.R
# ----------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[
  grep("^--file=", commandArgs(trailingOnly = FALSE))
]
if (length(script_file_arg) == 0) stop("Could not determine script path.")
script_dir  <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found: ", helper_file)
source(helper_file)

# ----------------------------------------------------------------
# Parse args
# ----------------------------------------------------------------
option_list <- list(
  make_option("--project_root",           type = "character"),
  make_option("--input_meth",             type = "character"),
  make_option("--sample_info",            type = "character"),
  make_option("--sample_id_col",          type = "character", default = NULL),
  make_option("--protected_traits_file",  type = "character", default = NULL),
  make_option("--technical_traits_file",  type = "character", default = NULL),
  make_option("--threads",                type = "integer",   default = 4),
  make_option("--top_trait_n",            type = "integer",   default = 50),
  make_option("--top_cell_n",             type = "integer",   default = 50),
  make_option("--primary_cor_method",     type = "character", default = "bicor"),
  make_option("--selected_trait_mode",    type = "character", default = NULL),
  make_option("--plot_width",             type = "double",    default = 12),
  make_option("--plot_height",            type = "double",    default = 10),
  make_option("--plot_dpi",               type = "integer",   default = 600),
  make_option("--base_size",              type = "double",    default = 11),
  make_option("--axis_text_size",         type = "double",    default = 8)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------------------------------------------
# Validate
# ----------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$input_meth))   stop("--input_meth is required")
if (is.null(opt$sample_info))  stop("--sample_info is required")

if (!dir.exists(opt$project_root))  stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$input_meth))   stop("input_meth not found: ", opt$input_meth)
if (!file.exists(opt$sample_info))  stop("sample_info not found: ", opt$sample_info)

if (!is.null(opt$protected_traits_file) && !file.exists(opt$protected_traits_file))
  stop("protected_traits_file not found: ", opt$protected_traits_file)
if (!is.null(opt$technical_traits_file) && !file.exists(opt$technical_traits_file))
  stop("technical_traits_file not found: ", opt$technical_traits_file)

if (opt$threads < 1)     stop("--threads must be >= 1")
if (opt$top_trait_n < 1) stop("--top_trait_n must be >= 1")
if (opt$top_cell_n < 1)  stop("--top_cell_n must be >= 1")

primary_cor_method <- tolower(opt$primary_cor_method)
if (!primary_cor_method %in% c("bicor", "pearson"))
  stop("--primary_cor_method must be bicor or pearson")

if (!is.null(opt$selected_trait_mode) &&
    !tolower(opt$selected_trait_mode) %in% c("protected", "protected_and_technical", "none"))
  stop("--selected_trait_mode must be: protected, protected_and_technical, or none")

WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ----------------------------------------------------------------
# Derive dataset_label / cpg_label / region_label from input path
#
# Supported path patterns (walking up from the .rds file):
#
#   05b: .../05b_shared_complete_regions/<ds>/<cpg>/<region>/<ds>_Methylation_complete.rds
#   07:  .../07_methylation_adjustment/<ds>/<cpg>/<region>/<variant>/Adjusted_Region_Methylation.rds
#   03:  .../03_region_methylation/<ds>/<cpg>/<region>/Region_Methylation.rds
#
# In all cases the structure one level above <region> is:
#   parent = <region_label>   basename(input_dir)
#   grandp = <cpg_label>      basename(dirname(input_dir))
#   greatg = <dataset_label>  basename(dirname(dirname(input_dir)))
#
# For 07 adjusted the variant folder sits between region and the file,
# so we walk up one extra level.
# ----------------------------------------------------------------
input_filename <- basename(opt$input_meth)
input_dir      <- dirname(opt$input_meth)

# Detect if we are inside a variant sub-folder (07 output)
# by checking whether the immediate parent looks like a variant name
# (v1_all_pcs / v2_exclude_protected / v3_technical_pcs_only / unadjusted)
variant_pattern <- "^(v[0-9]|unadjusted)"
if (grepl(variant_pattern, basename(input_dir))) {
  # Walk up one extra level past the variant folder
  region_label  <- basename(dirname(input_dir))
  cpg_label     <- basename(dirname(dirname(input_dir)))
  dataset_label <- basename(dirname(dirname(dirname(input_dir))))
} else if (grepl("05b_shared_complete_regions", input_dir)) {
  # 05b output: dataset label is encoded in the filename
  dataset_label <- sub("_Methylation_complete\\.rds$", "", input_filename)
  dataset_label <- sub("\\.rds$", "", dataset_label)
  region_label  <- basename(input_dir)
  cpg_label     <- basename(dirname(input_dir))
} else {
  # Standard 03/05 output
  region_label  <- basename(input_dir)
  cpg_label     <- basename(dirname(input_dir))
  dataset_label <- basename(dirname(dirname(input_dir)))
}

message("Derived labels:")
message("  dataset_label : ", dataset_label)
message("  cpg_label     : ", cpg_label)
message("  region_label  : ", region_label)

# ----------------------------------------------------------------
# Output directory
# ----------------------------------------------------------------
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "06_pc_diagnostics")
out_dir       <- file.path(step_dir, dataset_label, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", out_dir)

# ----------------------------------------------------------------
# File paths
# ----------------------------------------------------------------
pcs_rds                  <- file.path(out_dir, "PCs.rds")
pc_meta_tsv              <- file.path(out_dir, "PC_Metadata.tsv")
var_plot_pdf             <- file.path(out_dir, "Variance_Explained_First20_PCs.pdf")
pc_trait_pearson_tsv     <- file.path(out_dir, "PC_Trait_Correlation_Stats_Pearson.tsv")
pc_trait_bicor_tsv       <- file.path(out_dir, "PC_Trait_Correlation_Stats_Bicor.tsv")
heatmap_pub_pdf          <- file.path(out_dir, "PC_Trait_Bicor_Heatmap_PUBLICATION.pdf")
heatmap_topcells_pdf     <- file.path(out_dir, paste0("PC_Trait_Bicor_Heatmap_Top",
                                                       opt$top_cell_n, "Cells.pdf"))
trait_dendro_pdf         <- file.path(out_dir, "Trait_Dendrogram_From_PC_Associations.pdf")
removed_non_numeric_file <- file.path(out_dir, "sample_info_non_numeric_columns_removed.txt")
removed_all_na_file      <- file.path(out_dir, "sample_info_all_na_numeric_columns_removed.txt")
removed_zero_var_file    <- file.path(out_dir, "sample_info_zero_variance_numeric_columns_removed.txt")

# ----------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------
label_sv_in_stats <- function(stats_df, n_sv) {
  stats_df$module <- as.character(stats_df$module)
  sv_pcs <- paste0("PC_", seq_len(n_sv))
  stats_df$module <- ifelse(
    stats_df$module %in% sv_pcs,
    paste0(stats_df$module, " [SV]"),
    stats_df$module
  )
  all_pcs <- paste0("PC_", sort(unique(as.integer(
    gsub("PC_([0-9]+).*", "\\1", unique(stats_df$module))
  ))))
  new_levels <- ifelse(all_pcs %in% sv_pcs, paste0(all_pcs, " [SV]"), all_pcs)
  stats_df$module <- factor(stats_df$module, levels = new_levels)
  stats_df
}

getPCs_generic <- function(meth,
                            mod     = matrix(1, nrow = ncol(meth), ncol = 1),
                            save    = TRUE,
                            file    = "PCs.rds",
                            verbose = TRUE) {
  if (verbose) message("[getPCs_generic] Determining number of PCs via sva::num.sv")

  mod_noNA <- mod
  if (ncol(mod_noNA) > 1) {
    intercept_col <- apply(mod_noNA, 2, function(x) length(unique(x)) == 1)
    trait_cols    <- mod_noNA[, !intercept_col, drop = FALSE]
    col_var       <- apply(trait_cols, 2, function(x) var(x, na.rm = TRUE))
    trait_cols    <- trait_cols[, !is.na(col_var) & col_var > 0, drop = FALSE]
    if (ncol(trait_cols) > 0) {
      trait_cols <- scale(trait_cols, center = TRUE, scale = TRUE)
      trait_cols[is.na(trait_cols) | is.infinite(trait_cols)] <- 0
      mod_noNA <- cbind(mod_noNA[, intercept_col, drop = FALSE], trait_cols)
    } else {
      mod_noNA <- mod_noNA[, intercept_col, drop = FALSE]
    }
  }
  mod_noNA[is.na(mod_noNA) | is.infinite(mod_noNA)] <- 0
  if (any(!is.finite(mod_noNA))) {
    warning("[getPCs_generic] Non-finite values remain; falling back to intercept-only")
    mod_noNA <- matrix(1, nrow = ncol(meth), ncol = 1)
  }

  n_pc           <- sva::num.sv(meth, mod = mod_noNA, method = "be", seed = 5)
  meth_t         <- t(meth)
  meth_t_centered <- scale(meth_t, center = TRUE, scale = FALSE)
  ss             <- svd(meth_t_centered)

  PCs <- ss$u[, seq_len(n_pc), drop = FALSE]
  dimnames(PCs) <- list(rownames(meth_t), paste0("PC_", seq_len(n_pc)))

  var_explained <- (ss$d^2) / sum(ss$d^2)
  if (save) saveRDS(PCs, file = file)

  list(PCs = PCs, singular_values = ss$d,
       var_explained = var_explained, n_pc = n_pc)
}

# ----------------------------------------------------------------
# Load methylation
# ----------------------------------------------------------------
meth <- readRDS(opt$input_meth)
if (!(is.matrix(meth) || is.data.frame(meth))) stop("input_meth must be matrix-like.")
meth <- as.matrix(meth)
if (!is.numeric(meth))          stop("Methylation matrix must be numeric.")
if (nrow(meth) < 1)             stop("Methylation matrix has zero rows.")
if (ncol(meth) < 2)             stop("Methylation matrix must have >= 2 samples.")
if (is.null(colnames(meth)))    stop("Methylation matrix must have sample IDs in colnames.")
if (anyDuplicated(colnames(meth)))
  stop("Duplicate methylation sample IDs found.")
n_na <- sum(is.na(meth))
if (n_na > 0)
  stop("Methylation matrix contains ", n_na, " NA values. ",
       "Run Script 05b first to obtain an NA-free matrix.")
message("Loaded methylation: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ----------------------------------------------------------------
# Load sample info
# ----------------------------------------------------------------
sample_info <- readSampleInfo(
  file          = opt$sample_info,
  sample_id_col = opt$sample_id_col,
  verbose       = TRUE
)
if (nrow(sample_info) < 2)          stop("sample_info must have >= 2 samples.")
if (is.null(rownames(sample_info))) stop("sample_info must have rownames as sample IDs.")
if (anyDuplicated(rownames(sample_info)))
  stop("Duplicate sample IDs in sample_info.")

# ----------------------------------------------------------------
# Align samples
# ----------------------------------------------------------------
common_samples <- intersect(colnames(meth), rownames(sample_info))
if (length(common_samples) == 0)
  stop("No overlapping sample IDs between methylation and sample_info.")
if (length(common_samples) < 5)
  warning("Only ", length(common_samples), " overlapping samples found.")

meth        <- meth[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]
message("After alignment: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ----------------------------------------------------------------
# Numeric traits only
# ----------------------------------------------------------------
is_num          <- vapply(sample_info, function(x) is.numeric(x) && is.atomic(x), logical(1))
removed_non_num <- names(sample_info)[!is_num]
sample_info_num <- sample_info[, is_num, drop = FALSE]
write_vector_file(removed_non_num, removed_non_numeric_file)

if (ncol(sample_info_num) == 0)
  stop("No numeric trait columns available.")

all_na   <- vapply(sample_info_num, function(x) all(is.na(x)), logical(1))
zero_var <- vapply(sample_info_num, function(x) {
  v <- x[!is.na(x)]; if (length(v) <= 1) TRUE else length(unique(v)) == 1
}, logical(1))

removed_all_na   <- names(sample_info_num)[all_na]
removed_zero_var <- names(sample_info_num)[zero_var]
sample_info_num  <- sample_info_num[, !(all_na | zero_var), drop = FALSE]

write_vector_file(removed_all_na,   removed_all_na_file)
write_vector_file(removed_zero_var, removed_zero_var_file)

if (ncol(sample_info_num) == 0)
  stop("No usable numeric traits after removing all-NA and zero-variance columns.")
message("Numeric traits retained: ", ncol(sample_info_num))

# ----------------------------------------------------------------
# Protected / technical traits
# ----------------------------------------------------------------
protected_requested <- readTraitFile(opt$protected_traits_file, verbose = TRUE)
technical_requested <- readTraitFile(opt$technical_traits_file, verbose = TRUE)

protected_resolved <- resolveTraits(protected_requested, colnames(sample_info_num),
                                    label = "protected_traits", verbose = TRUE)
technical_resolved <- resolveTraits(technical_requested, colnames(sample_info_num),
                                    label = "technical_traits", verbose = TRUE)

if (!is.null(opt$protected_traits_file)) {
  write_vector_file(protected_resolved$requested, file.path(out_dir, "traits_requested_protected.txt"))
  write_vector_file(protected_resolved$found,     file.path(out_dir, "traits_found_protected.txt"))
  write_vector_file(protected_resolved$missing,   file.path(out_dir, "traits_missing_protected.txt"))
}
if (!is.null(opt$technical_traits_file)) {
  write_vector_file(technical_resolved$requested, file.path(out_dir, "traits_requested_technical.txt"))
  write_vector_file(technical_resolved$found,     file.path(out_dir, "traits_found_technical.txt"))
  write_vector_file(technical_resolved$missing,   file.path(out_dir, "traits_missing_technical.txt"))
}

# ----------------------------------------------------------------
# Selected trait mode
# ----------------------------------------------------------------
selected_trait_mode <- opt$selected_trait_mode
if (is.null(selected_trait_mode)) {
  selected_trait_mode <- if (length(protected_resolved$found) > 0) "protected" else "none"
}
selected_trait_mode <- tolower(selected_trait_mode)

selected_traits <- character(0)
if (selected_trait_mode == "protected") {
  selected_traits <- protected_resolved$found
} else if (selected_trait_mode == "protected_and_technical") {
  selected_traits <- unique(c(protected_resolved$found, technical_resolved$found))
}

# ----------------------------------------------------------------
# Model matrix for num.sv
# ----------------------------------------------------------------
old_na_action <- getOption("na.action")
options(na.action = "na.pass")
on.exit(options(na.action = old_na_action), add = TRUE)

if (length(protected_resolved$found) > 0) {
  protected_df <- sample_info_num[, protected_resolved$found, drop = FALSE]
  mod_formula  <- paste("~", paste(colnames(protected_df), collapse = " + "))
  mod          <- model.matrix(as.formula(mod_formula), data = protected_df)
  message("Protected traits in SVA model: ", paste(colnames(protected_df), collapse = ", "))
} else {
  mod_formula  <- "~ 1"
  mod          <- matrix(1, nrow = ncol(meth), ncol = 1)
  rownames(mod) <- colnames(meth)
  colnames(mod) <- "(Intercept)"
  message("No protected traits found; using intercept-only model for num.sv")
}

if (nrow(mod) != ncol(meth))
  stop("Model matrix rows (", nrow(mod), ") != methylation samples (", ncol(meth), ").")

# ----------------------------------------------------------------
# Derive PCs
# ----------------------------------------------------------------
pc_out         <- getPCs_generic(meth = meth, mod = mod, save = TRUE, file = pcs_rds)
PCs            <- pc_out$PCs
var_explained  <- pc_out$var_explained
sing_vals      <- pc_out$singular_values
n_pc           <- pc_out$n_pc
n_total_pcs    <- length(var_explained)

# ----------------------------------------------------------------
# PC metadata
# ----------------------------------------------------------------
pc_meta <- data.frame(
  PC                    = paste0("PC_", seq_along(sing_vals)),
  singular_value        = sing_vals,
  variance_explained    = var_explained,
  rank                  = seq_along(sing_vals),
  is_surrogate_variable = seq_along(sing_vals) <= n_pc,
  stringsAsFactors      = FALSE
)
readr::write_tsv(pc_meta, pc_meta_tsv)

# ----------------------------------------------------------------
# Variance explained plot
# ----------------------------------------------------------------
n_plot_pcs  <- min(20, n_total_pcs)
var_df_plot <- data.frame(
  PC                 = seq_len(n_plot_pcs),
  variance_explained = var_explained[seq_len(n_plot_pcs)],
  is_sv              = seq_len(n_plot_pcs) <= n_pc
)

p_var <- ggplot(var_df_plot, aes(x = PC, y = variance_explained)) +
  geom_line(color = "grey60", linewidth = 0.7) +
  geom_point(aes(color = is_sv), size = 2.5) +
  { if (n_pc > 0 && n_pc < n_plot_pcs)
      geom_vline(xintercept = n_pc + 0.5, linetype = "dashed",
                 linewidth = 0.6, color = "firebrick") } +
  scale_color_manual(
    values = c("TRUE"  = "firebrick", "FALSE" = "grey40"),
    labels = c("TRUE"  = paste0("Surrogate Variable (num.sv n=", n_pc, ")"),
               "FALSE" = "PC (not selected as SV)"),
    name   = NULL
  ) +
  scale_x_continuous(breaks = seq_len(n_plot_pcs)) +
  theme_bw(base_size = opt$base_size) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(
    title    = paste0("Variance Explained (Top ", n_plot_pcs, " PCs | SVs: ",
                      n_pc, ") - ", dataset_label),
    subtitle = "Red points = surrogate variables used for methylation adjustment",
    x        = "Principal Component",
    y        = "Proportion of Variance Explained"
  )

ggsave(var_plot_pdf, plot = p_var, width = opt$plot_width,
       height = opt$plot_height, dpi = opt$plot_dpi, units = "in", limitsize = FALSE)
message("Saved: Variance_Explained_First20_PCs.pdf")

# ----------------------------------------------------------------
# PC-trait correlations
# ----------------------------------------------------------------
PCtraitCor_pearson <- getMEtraitCor(PCs, colData = sample_info_num,
                                    corType = "pearson", file = pc_trait_pearson_tsv)
PCtraitCor_bicor   <- getMEtraitCor(PCs, colData = sample_info_num,
                                    corType = "bicor",   file = pc_trait_bicor_tsv)

primary_stats    <- if (primary_cor_method == "bicor") PCtraitCor_bicor else PCtraitCor_pearson
primary_stats_sv <- label_sv_in_stats(primary_stats, n_sv = n_pc)

# ----------------------------------------------------------------
# Publication heatmap (selected traits)
# ----------------------------------------------------------------
if (length(selected_traits) > 0 && selected_trait_mode != "none") {
  message("Saving publication heatmap: ", heatmap_pub_pdf)
  plotPCTrait(
    MEtraitCor     = primary_stats_sv,
    subsetTraits   = selected_traits,
    label_mode     = "star",
    autoSize       = TRUE,
    cell_in        = 0.22,
    max_width      = 14,
    max_height     = 8,
    force_x_labels = TRUE,
    axis.text.size = opt$axis_text_size,
    base_size      = opt$base_size,
    dpi            = opt$plot_dpi,
    file           = heatmap_pub_pdf,
    verbose        = TRUE
  )
} else {
  message("Skipping publication heatmap: no selected traits (mode = ", selected_trait_mode, ")")
}

# ----------------------------------------------------------------
# Top-cells heatmap
# ----------------------------------------------------------------
message("Saving top-cells heatmap: ", heatmap_topcells_pdf)
plotPCTrait(
  MEtraitCor     = primary_stats_sv,
  topCells       = opt$top_cell_n,
  label_mode     = "star",
  autoSize       = FALSE,
  width          = 10,
  height         = 16,
  force_x_labels = TRUE,
  force_y_labels = TRUE,
  axis.text.size = 9,
  wrapY          = 40,
  base_size      = opt$base_size,
  dpi            = opt$plot_dpi,
  file           = heatmap_topcells_pdf,
  verbose        = TRUE
)

# ----------------------------------------------------------------
# Trait dendrogram
# ----------------------------------------------------------------
plotTraitDendrogramFromPC(
  pc_trait_stats = primary_stats_sv,
  output_file    = trait_dendro_pdf,
  cor_column     = primary_cor_method,
  p_column       = "p",
  width          = 20,
  height         = 12,
  verbose        = TRUE
)

# ----------------------------------------------------------------
# Run parameters
# ----------------------------------------------------------------
writeLines(
  c(
    paste("project_root:",           opt$project_root),
    paste("input_meth:",             opt$input_meth),
    paste("sample_info:",            opt$sample_info),
    paste("sample_id_col:",          ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
    paste("protected_traits_file:",  ifelse(is.null(opt$protected_traits_file), "NULL", opt$protected_traits_file)),
    paste("technical_traits_file:",  ifelse(is.null(opt$technical_traits_file), "NULL", opt$technical_traits_file)),
    paste("dataset_label:",          dataset_label),
    paste("cpg_label:",              cpg_label),
    paste("region_label:",           region_label),
    paste("threads:",                opt$threads),
    paste("top_trait_n:",            opt$top_trait_n),
    paste("top_cell_n:",             opt$top_cell_n),
    paste("primary_cor_method:",     primary_cor_method),
    paste("selected_trait_mode:",    selected_trait_mode),
    paste("n_selected_traits:",      length(selected_traits)),
    paste("model_formula:",          mod_formula),
    paste("n_samples:",              ncol(meth)),
    paste("n_regions:",              nrow(meth)),
    paste("n_total_pcs:",            n_total_pcs),
    paste("n_svs_selected:",         n_pc),
    paste("n_numeric_traits:",       ncol(sample_info_num)),
    paste("n_protected_found:",      length(protected_resolved$found)),
    paste("n_technical_found:",      length(technical_resolved$found)),
    paste("output_dir:",             out_dir),
    paste("date:",                   as.character(Sys.time()))
  ),
  con = file.path(out_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

message("Script 06 complete: PC diagnostics finished")
message("Output: ", out_dir)
message("  dataset : ", dataset_label)
message("  cpg     : ", cpg_label)
message("  region  : ", region_label)
message("  SVs selected by num.sv: ", n_pc, " of ", n_total_pcs, " total PCs")