#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 14: Consensus ME-Trait Pair Plots
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Generates raw scatter/violin/dot plots for individual
#   module eigengene vs trait pairs per dataset:
#     1) Load Consensus_Modules.rds and extract per-dataset MEs
#     2) Load ME-trait stats from script 12a per dataset
#     3) Load sample info per dataset
#     4) For each trait-set file, select module-trait pairs
#        (optionally filtering to significant pairs only)
#     5) Plot each pair as a scatter (continuous) or violin/dot
#        (categorical/binary) plot
#
# REQUIRED INPUTS
#   --project_root          : root directory of the analysis project
#   --consensus_modules_rds : path to Consensus_Modules.rds from script 09
#   --regions_file          : path to Filtered_Regions.txt from script 02
#   --dataset1_label        : label for dataset 1
#   --dataset1_sample_info  : sample info file for dataset 1
#   --dataset1_stats_file   : ME-trait stats TSV for dataset 1 from script 12a
#
# OPTIONAL INPUTS
#   --dataset2_label/sample_info/stats_file  : dataset 2 inputs
#   --dataset3_label/sample_info/stats_file  : dataset 3 inputs
#   --sample_id_col          : sample ID column [default = use rownames]
#   --adjustment_version     : label matching script 12a run [default = unadjusted]
#   --set_dir                : directory of .txt trait-set files
#   --set_file to set_file5  : individual trait-set .txt files
#   --trait_exclude_file     : text file with traits to exclude
#   --p_thresh               : significance threshold [default = 0.05]
#   --only_significant       : plot only significant pairs TRUE/FALSE [default = TRUE]
#   --top_n_pairs            : maximum number of pairs to plot per set
#   --cat_plot_style         : dot | violin | both [default = violin]
#   --continuous_width/height  : continuous plot size [default = 8 x 8]
#   --categorical_width/height : categorical plot size [default = 8 x 8]
#
# OUTPUTS
#   comethyl_output/consensus/14_me_trait_pair_plots/
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           <set_name>/
#               <dataset_label>_<set_name>_pair_list.tsv
#               <dataset_label>_<module>_vs_<trait>.pdf
#               run_parameters.txt
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - Trait sets are plain text files with one trait name per line
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/14_me_trait_pair_plots_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --dataset1_label Baseline \
#     --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/Baseline_sample_info.xlsx \
#     --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --adjustment_version v1_all_pcs \
#     --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/trait_sets/ \
#     --cat_plot_style violin \
#     --only_significant TRUE
# ================================================================
# SCRIPT 14: Consensus ME-Trait Pair Plots
#
# PURPOSE
#   - Load Consensus_Modules.rds from script 09
#   - Load per-dataset sample info files
#   - Load one or more ME-trait stats files from script 12A
#   - Load one or more trait-set text files
#   - Generate raw pair plots for module-trait pairs per dataset
#
# DESIGN
#   - stats_file per dataset decides which module-trait pairs exist
#   - sample_info + consensus MEs provide the raw values for plotting
#   - trait sets are supplied externally via text files
#   - module assignments are shared; MEs are dataset-specific
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/14_me_trait_pair_plots/<adjustment_version>/
#       <dataset_label>/<set_name>/
#           <dataset_label>_traits_requested.txt
#           <dataset_label>_traits_found.txt
#           <dataset_label>_traits_missing.txt
#           <dataset_label>_pair_list.tsv
#           <dataset_label>_<module>_vs_<trait>.pdf
#           run_parameters.txt
#       sessionInfo.txt
#       run_parameters.txt
# ================================================================
message("Starting Script 14 ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(ggplot2)
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
validate_stats_table <- function(df, label) {
  required <- c("module", "trait", "p")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0)
    stop(label, " missing columns: ", paste(missing, collapse = ", "))
  if (!any(c("bicor", "cor") %in% colnames(df)))
    stop(label, " must contain 'bicor' or 'cor' column")
  df
}

sanitize_filename <- function(x) {
  x <- gsub("[/\\?<>\\:*|\"'` ]+", "_", x)
  gsub("_+", "_", x)
}

parse_optional_int <- function(x, arg_name) {
  if (is.null(x) || is.na(x) || !nzchar(as.character(x))) return(NULL)
  out <- suppressWarnings(as.integer(x))
  if (is.na(out) || out < 1) stop(arg_name, " must be a positive integer if provided.")
  out
}

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",          type = "character"),
  make_option("--consensus_modules_rds", type = "character"),
  make_option("--regions_file",          type = "character"),
  make_option("--adjustment_version",    type = "character", default = "unadjusted"),

  # Dataset 1
  make_option("--dataset1_label",        type = "character"),
  make_option("--dataset1_sample_info",  type = "character"),
  make_option("--dataset1_stats_file",   type = "character"),
  # Dataset 2
  make_option("--dataset2_label",        type = "character", default = NULL),
  make_option("--dataset2_sample_info",  type = "character", default = NULL),
  make_option("--dataset2_stats_file",   type = "character", default = NULL),
  # Dataset 3
  make_option("--dataset3_label",        type = "character", default = NULL),
  make_option("--dataset3_sample_info",  type = "character", default = NULL),
  make_option("--dataset3_stats_file",   type = "character", default = NULL),

  make_option("--sample_id_col",         type = "character", default = NULL),

  # Trait sets
  make_option("--set_dir",               type = "character", default = NULL),
  make_option("--set_file",              type = "character", default = NULL),
  make_option("--set_file2",             type = "character", default = NULL),
  make_option("--set_file3",             type = "character", default = NULL),
  make_option("--set_file4",             type = "character", default = NULL),
  make_option("--set_file5",             type = "character", default = NULL),

  make_option("--trait_exclude_file",    type = "character", default = NULL),
  make_option("--p_thresh",              type = "double",    default = 0.05),
  make_option("--only_significant",      type = "character", default = "TRUE"),
  make_option("--top_n_pairs",           type = "character", default = NULL),
  make_option("--cat_plot_style",        type = "character", default = "violin"),
  make_option("--continuous_width",      type = "double",    default = 8),
  make_option("--continuous_height",     type = "double",    default = 8),
  make_option("--categorical_width",     type = "double",    default = 8),
  make_option("--categorical_height",    type = "double",    default = 8)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate
# ------------------------------------------------------------
if (is.null(opt$project_root))          stop("--project_root is required")
if (is.null(opt$consensus_modules_rds)) stop("--consensus_modules_rds is required")
if (is.null(opt$regions_file))          stop("--regions_file is required")
if (!dir.exists(opt$project_root))      stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
if (!file.exists(opt$regions_file))     stop("regions_file not found: ", opt$regions_file)
if (is.null(opt$dataset1_label))        stop("--dataset1_label is required")
if (is.null(opt$dataset1_sample_info))  stop("--dataset1_sample_info is required")
if (is.null(opt$dataset1_stats_file))   stop("--dataset1_stats_file is required")
if (!dir.exists(opt$project_root))      stop("project_root not found")
if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found")
if (!file.exists(opt$dataset1_sample_info))  stop("dataset1_sample_info not found")
if (!file.exists(opt$dataset1_stats_file))   stop("dataset1_stats_file not found")

only_significant <- parse_bool(opt$only_significant, "--only_significant")
top_n_pairs      <- parse_optional_int(opt$top_n_pairs, "--top_n_pairs")

cat_plot_style <- tolower(opt$cat_plot_style)
if (!cat_plot_style %in% c("dot", "violin", "both"))
  stop("--cat_plot_style must be dot, violin, or both")
if (opt$p_thresh <= 0 || opt$p_thresh >= 1)
  stop("--p_thresh must be > 0 and < 1")
if (!is.null(opt$trait_exclude_file) && !file.exists(opt$trait_exclude_file))
  stop("trait_exclude_file not found: ", opt$trait_exclude_file)

# Build dataset list — each dataset needs label, sample_info, stats_file
dataset_inputs <- list(
  list(label       = opt$dataset1_label,
       sample_info = opt$dataset1_sample_info,
       stats_file  = opt$dataset1_stats_file)
)
for (i in 2:3) {
  lbl <- opt[[paste0("dataset", i, "_label")]]
  si  <- opt[[paste0("dataset", i, "_sample_info")]]
  sf  <- opt[[paste0("dataset", i, "_stats_file")]]
  if (!is.null(lbl) || !is.null(si) || !is.null(sf)) {
    if (is.null(lbl) || is.null(si) || is.null(sf))
      stop("Provide --dataset", i, "_label, --dataset", i,
           "_sample_info, and --dataset", i, "_stats_file together")
    if (!file.exists(si)) stop("dataset", i, "_sample_info not found: ", si)
    if (!file.exists(sf)) stop("dataset", i, "_stats_file not found: ", sf)
    dataset_inputs[[length(dataset_inputs) + 1]] <-
      list(label = lbl, sample_info = si, stats_file = sf)
  }
}

set_files <- collect_set_files(
  set_dir   = opt$set_dir,
  set_file  = opt$set_file,
  set_file2 = opt$set_file2,
  set_file3 = opt$set_file3,
  set_file4 = opt$set_file4,
  set_file5 = opt$set_file5
)

# ------------------------------------------------------------
# Configure
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption("CACHE",
  file.path(opt$project_root, ".cache"))
WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Output dirs
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
step_dir      <- file.path(pipeline_root, "14_me_trait_pair_plots")
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)
message("Output root: ", step_dir)

# ------------------------------------------------------------
# Load consensus modules
# ------------------------------------------------------------
consensusMods  <- readRDS(opt$consensus_modules_rds)
all_colors     <- consensusMods$colors
module_colors  <- unique(all_colors[all_colors != "grey"])
n_modules      <- length(module_colors)
message("Non-grey consensus modules: ", n_modules)

if (n_modules == 0)
  stop("No non-grey modules found. Cannot generate pair plots.")

# Optional shared trait exclusion
trait_exclude_requested <- readTraitFile(opt$trait_exclude_file, verbose = TRUE)

# ------------------------------------------------------------
# Per-dataset loop
# ------------------------------------------------------------
for (ds in dataset_inputs) {
  ds_label <- ds$label
  message("\n==============================")
  message("Dataset: ", ds_label)
  message("==============================")

  # Get MEs — remove MEgrey
  MEs_full <- consensusMods$multiMEs[[ds_label]]$data
  grey_col  <- grep("^MEgrey$", colnames(MEs_full), value = TRUE)
  MEs_full  <- MEs_full[, !colnames(MEs_full) %in% grey_col, drop = FALSE]

  if (ncol(MEs_full) == 0) {
    message(ds_label, ": No real MEs — skipping")
    next
  }
  message(ds_label, ": ", nrow(MEs_full), " samples x ",
          ncol(MEs_full), " modules")

  # Load stats
  stats_df <- readr::read_tsv(ds$stats_file, show_col_types = FALSE)
  stats_df <- validate_stats_table(stats_df, paste0(ds_label, " stats_file"))
  stats_df$module <- as.character(stats_df$module)
  stats_df$trait  <- as.character(stats_df$trait)

  # Load sample info
  sample_info <- readSampleInfo(
    file          = ds$sample_info,
    sample_id_col = opt$sample_id_col,
    verbose       = TRUE
  )
  if (nrow(sample_info) < 2)
    stop(ds_label, ": sample_info must have >= 2 samples")

  # Trait exclusion
  trait_exclude_resolved <- resolveTraits(
    requested_traits = trait_exclude_requested,
    available_traits = colnames(sample_info),
    label            = paste0(ds_label, "_trait_exclude"),
    verbose          = TRUE
  )
  if (length(trait_exclude_resolved$found) > 0)
    sample_info <- sample_info[
      , !colnames(sample_info) %in% trait_exclude_resolved$found,
      drop = FALSE
    ]

  # Align samples
  common_samples <- intersect(rownames(MEs_full), rownames(sample_info))
  if (length(common_samples) == 0)
    stop(ds_label, ": No overlapping samples between MEs and sample_info")
  message(ds_label, ": ", length(common_samples), " overlapping samples")

  MEs_use         <- MEs_full[common_samples, , drop = FALSE]
  sample_info_use <- sample_info[common_samples, , drop = FALSE]

  # MEs_use has columns like "MEturquoise", "MEred", etc.
  # stats_df$module has bare color names like "turquoise", "red"
  # Build a lookup: bare color name -> ME column name
  me_col_lookup <- setNames(
    colnames(MEs_use),
    sub("^ME", "", colnames(MEs_use))
  )

  # ── Per trait set ─────────────────────────────────────────
  for (set_file in set_files) {
    set_name    <- get_set_name(set_file)
    set_out_dir <- file.path(step_dir, ds_label, cpg_label, region_label, opt$adjustment_version, set_name)
    dir.create(set_out_dir, recursive = TRUE, showWarnings = FALSE)

    requested_traits <- read_trait_set_file(set_file)
    found_traits     <- intersect(requested_traits, colnames(sample_info_use))
    missing_traits   <- setdiff(requested_traits, colnames(sample_info_use))

    write_vector_file(requested_traits,
      file.path(set_out_dir,
        paste0(ds_label, "_", set_name, "_traits_requested.txt")))
    write_vector_file(found_traits,
      file.path(set_out_dir,
        paste0(ds_label, "_", set_name, "_traits_found.txt")))
    write_vector_file(missing_traits,
      file.path(set_out_dir,
        paste0(ds_label, "_", set_name, "_traits_missing.txt")))

    if (length(found_traits) == 0) {
      message(ds_label, "/", set_name, ": no traits found — skipping")
      write_log_lines(
        c(paste("dataset:", ds_label),
          paste("set_name:", set_name),
          "status: skipped",
          "reason: no requested traits found in sample_info"),
        file.path(set_out_dir, "run_parameters.txt")
      )
      next
    }

    # Filter pairs — stats module column has bare colors, MEs have ME prefix
    pair_df <- stats_df %>%
      dplyr::filter(.data$trait  %in% found_traits,
                    .data$module %in% names(me_col_lookup))

    if (only_significant)
      pair_df <- dplyr::filter(pair_df, !is.na(.data$p), .data$p < opt$p_thresh)

    pair_df <- dplyr::arrange(pair_df, .data$p)

    if (!is.null(top_n_pairs) && nrow(pair_df) > top_n_pairs)
      pair_df <- dplyr::slice_head(pair_df, n = top_n_pairs)

    readr::write_tsv(pair_df,
      file.path(set_out_dir,
        paste0(ds_label, "_", set_name, "_pair_list.tsv")))

    if (nrow(pair_df) == 0) {
      message(ds_label, "/", set_name,
              ": no pairs after filtering — skipping plots")
      write_log_lines(
        c(paste("dataset:", ds_label),
          paste("set_name:", set_name),
          paste("only_significant:", only_significant),
          paste("p_thresh:", opt$p_thresh),
          "status: skipped",
          "reason: no pairs left after filtering"),
        file.path(set_out_dir, "run_parameters.txt")
      )
      next
    }

    message(ds_label, "/", set_name, ": plotting ", nrow(pair_df), " pairs")

    # ── Plot each pair ──────────────────────────────────────
    for (i in seq_len(nrow(pair_df))) {
      mod <- as.character(pair_df$module[i])
      tr  <- as.character(pair_df$trait[i])

      me_col <- me_col_lookup[mod]  # "MEturquoise" from "turquoise"
      if (is.na(me_col) || !(me_col %in% colnames(MEs_use))) next
      if (!(tr  %in% colnames(sample_info_use))) next

      ME_vec    <- MEs_use[[me_col]]
      trait_vec <- sample_info_use[[tr]]
      yl        <- auto_ylim(ME_vec)

      subtitle_text <- make_me_trait_subtitle(stats_df, module = mod, trait = tr)

      make_cat <- is.factor(trait_vec) ||
                  is.character(trait_vec) ||
                  is_binary_like(trait_vec)

      if (make_cat) {
        if (!is.factor(trait_vec)) trait_vec <- as.factor(trait_vec)
        if (nlevels(trait_vec) == 2)
          trait_vec <- relabel_binary_to_yesno(trait_vec)
      }

      # Dataset label prefixed to every plot filename
      pdf_file <- file.path(set_out_dir,
        paste0(ds_label, "_",
               sanitize_filename(mod), "_vs_",
               sanitize_filename(tr), ".pdf"))

      if (!make_cat) {
        p <- plotMEtraitScatter(
          ME      = ME_vec,
          trait   = trait_vec,
          ylim    = yl,
          xlab    = tr,
          ylab    = paste0(mod, " Module Eigengene"),
          save    = FALSE,
          axis.text.size = 16,
          verbose = FALSE
        )
        p <- add_plot_subtitle(p, subtitle_text)
        ggplot2::ggsave(pdf_file, plot = p, dpi = 600,
                        width  = opt$continuous_width,
                        height = opt$continuous_height,
                        units  = "in")

      } else if (cat_plot_style == "dot") {
        p <- plotMEtraitDot(
          ME      = ME_vec,
          trait   = trait_vec,
          ylim    = yl,
          xlab    = tr,
          ylab    = paste0(mod, " Module Eigengene"),
          save    = FALSE,
          axis.text.size = 16,
          verbose = FALSE
        )
        p <- add_plot_subtitle(p, subtitle_text)
        ggplot2::ggsave(pdf_file, plot = p, dpi = 600,
                        width  = opt$categorical_width,
                        height = opt$categorical_height,
                        units  = "in")

      } else if (cat_plot_style == "violin") {
        p <- plotMEtraitViolin(
          ME_vec     = ME_vec,
          trait_fac  = trait_vec,
          module_name = mod,
          trait_name = tr,
          ylim       = yl
        )
        p <- add_plot_subtitle(p, subtitle_text)
        if (!is.null(p))
          ggplot2::ggsave(pdf_file, plot = p, dpi = 600,
                          width  = opt$categorical_width,
                          height = opt$categorical_height,
                          units  = "in")

      } else if (cat_plot_style == "both") {
        grDevices::pdf(pdf_file,
                       width  = opt$categorical_width,
                       height = opt$categorical_height)
        tryCatch({
          p_dot <- plotMEtraitDot(
            ME      = ME_vec,
            trait   = trait_vec,
            ylim    = yl,
            xlab    = tr,
            ylab    = paste0(mod, " Module Eigengene"),
            save    = FALSE,
            axis.text.size = 16,
            verbose = FALSE
          )
          print(add_plot_subtitle(p_dot, subtitle_text))

          p_v <- plotMEtraitViolin(
            ME_vec      = ME_vec,
            trait_fac   = trait_vec,
            module_name = mod,
            trait_name  = tr,
            ylim        = yl
          )
          if (!is.null(p_v)) print(add_plot_subtitle(p_v, subtitle_text))
        }, error = function(e)
          message(ds_label, ": pair plot failed for ", mod, " vs ", tr,
                  " — ", conditionMessage(e)))
        dev.off()
      }
    }

    write_log_lines(
      c(
        paste("dataset_label:",      ds_label),
        paste("set_name:",           set_name),
        paste("set_file:",           set_file),
        paste("consensus_modules:",  opt$consensus_modules_rds),
        paste("stats_file:",         ds$stats_file),
        paste("sample_info:",        ds$sample_info),
        paste("only_significant:",   only_significant),
        paste("p_thresh:",           opt$p_thresh),
        paste("top_n_pairs:",        ifelse(is.null(top_n_pairs), "NULL", top_n_pairs)),
        paste("cat_plot_style:",     cat_plot_style),
        paste("n_traits_requested:", length(requested_traits)),
        paste("n_traits_found:",     length(found_traits)),
        paste("n_traits_missing:",   length(missing_traits)),
        paste("n_pairs_plotted:",    nrow(pair_df)),
        paste("date:",               as.character(Sys.time()))
      ),
      file.path(set_out_dir, "run_parameters.txt")
    )

    message("✓ ", ds_label, "/", set_name, ": done")
  }
}

# ------------------------------------------------------------
# Top-level logs
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:",          opt$project_root),
    paste("consensus_modules_rds:", opt$consensus_modules_rds),
    paste("adjustment_version:",    opt$adjustment_version),
    paste("datasets:",              paste(sapply(dataset_inputs, `[[`, "label"),
                                         collapse = ", ")),
    paste("set_files:",             paste(set_files, collapse = ", ")),
    paste("only_significant:",      only_significant),
    paste("p_thresh:",              opt$p_thresh),
    paste("top_n_pairs:",           ifelse(is.null(top_n_pairs), "NULL", top_n_pairs)),
    paste("cat_plot_style:",        cat_plot_style),
    paste("n_real_modules:",        n_modules),
    paste("date:",                  as.character(Sys.time()))
  ),
  file.path(step_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()),
           con = file.path(step_dir, "sessionInfo.txt"))

message("\n✓ Script 14 complete: consensus ME-trait pair plots finished")
message("Outputs saved under: ", step_dir)