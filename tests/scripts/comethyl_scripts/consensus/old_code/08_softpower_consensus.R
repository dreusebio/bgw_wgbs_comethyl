#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08: Consensus Soft-Threshold Power Analysis
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Loads adjusted methylation matrices from two or three datasets
#   (all sharing the same adjustment version) and:
#     1) Runs soft-threshold power analysis separately per dataset
#     2) Generates per-dataset softpower tables and sample dendrograms
#     3) Selects a single shared consensus power (smallest power where
#        all datasets meet the scale-free topology cutoff)
#     4) Generates a combined multi-dataset softpower plot
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --dataset1_label : label for dataset 1 (e.g. Baseline)
#   --dataset1_meth  : path to adjusted Region_Methylation.rds for dataset 1
#
# OPTIONAL INPUTS
#   --dataset2_label     : label for dataset 2
#   --dataset2_meth      : path to adjusted methylation RDS for dataset 2
#   --dataset3_label     : label for dataset 3
#   --dataset3_meth      : path to adjusted methylation RDS for dataset 3
#   --adjustment_version : label for this run (e.g. v1_all_pcs) [default = unadjusted]
#   --softpower_cor      : bicor or pearson [default = pearson]
#   --power_min          : minimum soft power to test [default = 1]
#   --power_max          : maximum soft power to test [default = 20]
#   --block_size         : WGCNA block size [default = 20000]
#   --gc_interval        : garbage collection interval [default = 19999]
#   --threads            : WGCNA threads [default = 4]
#   --plot_sample_dendro : save sample dendrogram per dataset TRUE/FALSE [default = TRUE]
#   --scale_free_cutoff  : R² cutoff for scale-free topology [default = 0.8]
#
# OUTPUTS
#   comethyl_output/consensus/08_soft_power/
#       shared/<cpg_label>/<region_label>/<adjustment_version>/
#           chosen_power.txt
#           combined_softpower_summary.tsv
#           Consensus_SoftPower_Combined.pdf
#           run_parameters.txt
#           sessionInfo.txt
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           SoftPower_<cor>.rds
#           SoftPower_<cor>_table.tsv
#           Sample_Dendrogram_<cor>.pdf
#
# NOTES
#   - cpg_label and region_label are derived from --dataset1_meth path
#   - The chosen_power.txt output feeds into script 09 (--chosen_power_file)
#   - Run once per adjustment version, passing all datasets together
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/08_softpower_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --dataset1_label Baseline \
#     --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline_Adjusted_Region_Methylation_allPCs.rds \
#     --dataset2_label 36_38wks \
#     --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36_38wks/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/36_38wks_Adjusted_Region_Methylation_allPCs.rds \
#     --dataset3_label Postpartum \
#     --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Postpartum_Adjusted_Region_Methylation_allPCs.rds \
#     --adjustment_version v1_all_pcs \
#     --softpower_cor bicor \
#     --scale_free_cutoff 0.8 \
#     --threads 4
# ================================================================
# SCRIPT 08: Consensus Soft-Threshold Power Analysis
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/08_soft_power/
#     <cpg_label>/
#       <region_label>/
#         <adjustment_version>/
#           chosen_power.txt
#           combined_softpower_summary.tsv
#           Consensus_SoftPower_Combined.pdf
#           run_parameters.txt
#           sessionInfo.txt
#           <dataset_label>/
#             SoftPower_<cor>.rds
#             SoftPower_<cor>_table.tsv
#             Sample_Dendrogram_<cor>.pdf
#
# PATH DERIVATION
#   cpg_label and region_label are derived from --dataset1_meth path.
#   Expected path structure (from script 07 output):
#     .../07_methylation_adjustment/<ds>/<cpg_label>/<region_label>/<variant>/...rds
#   OR from script 05b:
#     .../05b_shared_complete_regions/<ds>/<cpg_label>/<region_label>/...rds
# ================================================================

message("Starting Script 8 ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(comethyl)
  library(gridExtra)
})

# ----------------------------------------------------------------
# Parse args
# ----------------------------------------------------------------
option_list <- list(
  make_option("--project_root",       type = "character"),
  make_option("--dataset1_label",     type = "character"),
  make_option("--dataset1_meth",      type = "character"),
  make_option("--dataset2_label",     type = "character", default = NULL),
  make_option("--dataset2_meth",      type = "character", default = NULL),
  make_option("--dataset3_label",     type = "character", default = NULL),
  make_option("--dataset3_meth",      type = "character", default = NULL),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--softpower_cor",      type = "character", default = "pearson"),
  make_option("--power_min",          type = "integer",   default = 1),
  make_option("--power_max",          type = "integer",   default = 20),
  make_option("--block_size",         type = "integer",   default = 20000),
  make_option("--gc_interval",        type = "integer",   default = 19999),
  make_option("--threads",            type = "integer",   default = 4),
  make_option("--plot_sample_dendro", type = "character", default = "TRUE"),
  make_option("--scale_free_cutoff",  type = "double",    default = 0.8)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  if (is.logical(x)) return(x)
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true","t","1","yes","y")) return(TRUE)
  if (x2 %in% c("false","f","0","no","n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) stop(label, " must be matrix-like.")
  x <- as.matrix(x)
  if (!is.numeric(x))          stop(label, " must be numeric.")
  if (nrow(x) < 1)             stop(label, " has zero rows.")
  if (ncol(x) < 2)             stop(label, " must have >= 2 samples.")
  if (is.null(rownames(x)))    stop(label, " must have region IDs in rownames.")
  if (is.null(colnames(x)))    stop(label, " must have sample IDs in colnames.")
  if (anyDuplicated(rownames(x))) stop(label, " has duplicated region IDs.")
  if (anyDuplicated(colnames(x))) stop(label, " has duplicated sample IDs.")
  x
}

extract_sft_table <- function(sft) {
  if (is.data.frame(sft))      return(sft)
  if (is.list(sft) && !is.null(sft$fitIndices)) return(as.data.frame(sft$fitIndices))
  if (is.list(sft) && length(sft) >= 1)         return(as.data.frame(sft[[1]]))
  stop("Cannot extract SFT table from getSoftPower() output.")
}

choose_consensus_power <- function(all_summary_list, cutoff = 0.8) {
  combined <- dplyr::bind_rows(all_summary_list)
  power_summary <- combined %>%
    dplyr::group_by(Power) %>%
    dplyr::summarise(
      n_datasets      = dplyr::n(),
      n_meet_cutoff   = sum(SFT_R2 >= cutoff, na.rm = TRUE),
      min_SFT_R2      = min(SFT_R2, na.rm = TRUE),
      median_SFT_R2   = median(SFT_R2, na.rm = TRUE),
      mean_connectivity = mean(Mean_Connectivity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Power)

  all_ok <- dplyr::filter(power_summary, n_meet_cutoff == n_datasets)
  if (nrow(all_ok) > 0)
    return(list(power = all_ok$Power[1],
                reason = "smallest power where all datasets meet scale-free cutoff",
                summary = power_summary))

  fallback <- dplyr::arrange(power_summary, desc(n_meet_cutoff), Power)
  list(power   = fallback$Power[1],
       reason  = "fallback: smallest power maximizing datasets meeting cutoff",
       summary = power_summary)
}

make_combined_softpower_plot <- function(summary_df, out_file, cutoff = 0.8,
                                          chosen_power = NULL) {
  make_p <- function(y_var, y_lab, add_cutoff = FALSE) {
    p <- ggplot(summary_df, aes_string(x = "Power", y = y_var,
                                        color = "Dataset", group = "Dataset")) +
      geom_line() + geom_point() +
      { if (add_cutoff) geom_hline(yintercept = cutoff, linetype = "dashed") } +
      { if (!is.null(chosen_power))
          geom_vline(xintercept = chosen_power, linetype = "dotted",
                     color = "black", linewidth = 0.8) } +
      scale_x_continuous(breaks = unique(summary_df$Power)) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank()) +
      labs(x = "Soft Threshold (Power)", y = y_lab)
    if (add_cutoff && !is.null(chosen_power))
      p <- p + annotate("text", x = chosen_power + 0.3, y = 0.05,
                         label = paste0("power = ", chosen_power), hjust = 0, size = 3.5)
    p
  }
  grDevices::pdf(out_file, width = 14, height = 10)
  gridExtra::grid.arrange(
    make_p("SFT_R2",            "Scale-Free R²",     add_cutoff = TRUE),
    make_p("Mean_Connectivity", "Mean Connectivity"),
    make_p("Median_Connectivity","Median Connectivity"),
    make_p("Max_Connectivity",  "Max Connectivity"),
    ncol = 2
  )
  dev.off()
}

# ----------------------------------------------------------------
# Validate
# ----------------------------------------------------------------
if (is.null(opt$project_root))   stop("--project_root is required")
if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth))  stop("--dataset1_meth is required")
if (!dir.exists(opt$project_root)) stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$dataset1_meth)) stop("dataset1_meth not found: ", opt$dataset1_meth)

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)
if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth))
    stop("Provide both --dataset2_label and --dataset2_meth.")
  if (!file.exists(opt$dataset2_meth)) stop("dataset2_meth not found: ", opt$dataset2_meth)
}
if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth))
    stop("Provide both --dataset3_label and --dataset3_meth.")
  if (!file.exists(opt$dataset3_meth)) stop("dataset3_meth not found: ", opt$dataset3_meth)
}

softpower_cor      <- tolower(opt$softpower_cor)
if (!softpower_cor %in% c("bicor","pearson")) stop("--softpower_cor must be bicor or pearson")
plot_sample_dendro <- parse_bool(opt$plot_sample_dendro, "--plot_sample_dendro")
powerVector        <- seq(opt$power_min, opt$power_max)

WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ----------------------------------------------------------------
# Derive cpg_label / region_label from dataset1_meth
# Walk up from the .rds file:
#   .../07_methylation_adjustment/<ds>/<cpg>/<region>/<variant>/<file>.rds
#   .../05b_shared_complete_regions/<ds>/<cpg>/<region>/<file>.rds
#   .../03_region_methylation/<ds>/<cpg>/<region>/<file>.rds
# ----------------------------------------------------------------
d1_dir  <- dirname(opt$dataset1_meth)
variant_pattern <- "^(v[0-9]|unadjusted)"
if (grepl(variant_pattern, basename(d1_dir))) {
  # inside a variant sub-folder — walk one extra level up
  region_label <- basename(dirname(d1_dir))
  cpg_label    <- basename(dirname(dirname(d1_dir)))
} else {
  region_label <- basename(d1_dir)
  cpg_label    <- basename(dirname(d1_dir))
}

message("Derived labels:")
message("  cpg_label    : ", cpg_label)
message("  region_label : ", region_label)

# ----------------------------------------------------------------
# Output directory
# ----------------------------------------------------------------
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "08_soft_power")
shared_dir    <- file.path(step_dir, "shared", cpg_label, region_label, opt$adjustment_version)
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
message("Shared output directory: ", shared_dir)

# ----------------------------------------------------------------
# Build dataset list
# ----------------------------------------------------------------
dataset_inputs <- list(list(label = opt$dataset1_label, file = opt$dataset1_meth))
if (dataset2_provided)
  dataset_inputs[[2]] <- list(label = opt$dataset2_label, file = opt$dataset2_meth)
if (dataset3_provided)
  dataset_inputs[[length(dataset_inputs)+1]] <- list(label = opt$dataset3_label, file = opt$dataset3_meth)

# ----------------------------------------------------------------
# Run soft power per dataset
# ----------------------------------------------------------------
all_summary <- list()

for (ds in dataset_inputs) {
  ds_label <- ds$label
  message("\n============================== Dataset: ", ds_label)

  meth <- validate_meth_matrix(readRDS(ds$file), ds_label)

  ds_out_dir <- file.path(step_dir, ds_label, cpg_label, region_label, opt$adjustment_version)
  dir.create(ds_out_dir, recursive = TRUE, showWarnings = FALSE)

  sft_rds       <- file.path(ds_out_dir, paste0("SoftPower_", softpower_cor, ".rds"))
  sft_table_tsv <- file.path(ds_out_dir, paste0("SoftPower_", softpower_cor, "_table.tsv"))
  dendro_pdf    <- file.path(ds_out_dir, paste0("Sample_Dendrogram_", softpower_cor, ".pdf"))

  sft <- getSoftPower(meth, powerVector = powerVector, corType = softpower_cor,
                       file = sft_rds, blockSize = opt$block_size,
                       gcInterval = opt$gc_interval)

  sft_table     <- extract_sft_table(sft)
  cn_lower      <- tolower(colnames(sft_table))

  get_col <- function(candidates) {
    for (c in candidates) {
      hit <- match(c, cn_lower)
      if (!is.na(hit)) return(colnames(sft_table)[hit])
    }
    NA_character_
  }

  power_col  <- get_col(c("power"))
  fit_col    <- get_col(c("sft.r.sq","sft.rsq"))
  mean_col   <- get_col(c("mean.k.","mean.k"))
  median_col <- get_col(c("median.k.","median.k"))
  max_col    <- get_col(c("max.k.","max.k"))

  if (any(is.na(c(power_col, fit_col, mean_col, median_col, max_col))))
    stop("Could not identify expected columns in SFT table for dataset: ", ds_label)

  summary_df <- tibble::tibble(
    Dataset             = ds_label,
    Power               = sft_table[[power_col]],
    SFT_R2              = abs(sft_table[[fit_col]]),
    Mean_Connectivity   = sft_table[[mean_col]],
    Median_Connectivity = sft_table[[median_col]],
    Max_Connectivity    = sft_table[[max_col]]
  )

  readr::write_tsv(summary_df, sft_table_tsv)
  all_summary[[ds_label]] <- summary_df

  if (isTRUE(plot_sample_dendro)) {
    tryCatch({
      d <- getDendro(meth, distance = "euclidean")
      plotDendro(d, file = dendro_pdf, expandY = c(0.25, 0.08))
      message("Saved sample dendrogram: ", ds_label)
    }, error = function(e)
      warning("Sample dendrogram failed for ", ds_label, ": ", conditionMessage(e)))
  }

  message("Soft power done: ", ds_label)
}

# ----------------------------------------------------------------
# Combined summary and chosen power
# ----------------------------------------------------------------
combined_summary <- dplyr::bind_rows(all_summary)
readr::write_tsv(combined_summary,
                 file.path(shared_dir, "combined_softpower_summary.tsv"))

power_result  <- choose_consensus_power(all_summary, cutoff = opt$scale_free_cutoff)
chosen_power  <- power_result$power
chosen_reason <- power_result$reason

writeLines(
  c(paste("chosen_power:", chosen_power),
    paste("reason:", chosen_reason),
    paste("scale_free_cutoff:", opt$scale_free_cutoff),
    paste("adjustment_version:", opt$adjustment_version),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("datasets:", paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
    paste("date:", as.character(Sys.time()))),
  con = file.path(shared_dir, "chosen_power.txt")
)
message("Chosen consensus power: ", chosen_power, " (", chosen_reason, ")")

make_combined_softpower_plot(
  summary_df   = combined_summary,
  out_file     = file.path(shared_dir, "Consensus_SoftPower_Combined.pdf"),
  cutoff       = opt$scale_free_cutoff,
  chosen_power = chosen_power
)
message("Saved: Consensus_SoftPower_Combined.pdf")

# ----------------------------------------------------------------
# Run parameters
# ----------------------------------------------------------------
writeLines(c(
  paste("project_root:",        opt$project_root),
  paste("adjustment_version:",  opt$adjustment_version),
  paste("cpg_label:",           cpg_label),
  paste("region_label:",        region_label),
  paste("datasets:",            paste(sapply(dataset_inputs, `[[`, "label"), collapse = ", ")),
  paste("softpower_cor:",       softpower_cor),
  paste("power_min:",           opt$power_min),
  paste("power_max:",           opt$power_max),
  paste("scale_free_cutoff:",   opt$scale_free_cutoff),
  paste("block_size:",          opt$block_size),
  paste("gc_interval:",         opt$gc_interval),
  paste("threads:",             opt$threads),
  paste("plot_sample_dendro:",  plot_sample_dendro),
  paste("chosen_power:",        chosen_power),
  paste("chosen_reason:",       chosen_reason),
  paste("shared_dir:",          shared_dir),
  paste("date:",                as.character(Sys.time()))
), con = file.path(shared_dir, "run_parameters.txt"))

writeLines(capture.output(sessionInfo()), con = file.path(shared_dir, "sessionInfo.txt"))

message(" Script 08 complete: soft power analysis finished")
message("Shared output: ", shared_dir)