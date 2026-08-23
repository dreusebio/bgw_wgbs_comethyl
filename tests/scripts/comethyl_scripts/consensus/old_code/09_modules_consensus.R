#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 09: Consensus Module Detection
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Loads adjusted methylation matrices from two or three datasets
#   sharing the same reference region set and:
#     1) Ensures all datasets share identical region ordering
#     2) Builds a WGCNA multiExpr object
#     3) Runs blockwiseConsensusModules() to detect shared modules
#     4) Saves the consensus module object (shared across datasets)
#     5) Saves dataset-specific module eigengenes and dendrograms
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --regions_file   : path to Filtered_Regions.txt from script 02
#                      (used to derive cpg_label and region_label)
#   --dataset1_label : label for dataset 1 (e.g. Baseline)
#   --dataset1_meth  : path to adjusted Region_Methylation.rds for dataset 1
#
# OPTIONAL INPUTS
#   --dataset2_label       : label for dataset 2
#   --dataset2_meth        : path to adjusted methylation RDS for dataset 2
#   --dataset3_label       : label for dataset 3
#   --dataset3_meth        : path to adjusted methylation RDS for dataset 3
#   --adjustment_version   : label for this run [default = unadjusted]
#   --power                : soft-threshold power (overrides chosen_power_file)
#   --chosen_power_file    : path to chosen_power.txt from script 08
#   --consensus_cor        : bicor or pearson [default = bicor]
#   --network_type         : signed or unsigned [default = signed]
#   --tom_type             : signed or unsigned [default = signed]
#   --deep_split           : module detection sensitivity 0-4 [default = 4]
#   --min_module_size      : minimum regions per module [default = 5]
#   --merge_cut_height     : ME correlation threshold for merging [default = 0.1]
#   --max_block_size       : maximum block size for blockwise analysis [default = 40000]
#   --network_calibration  : calibration method [default = full quantile]
#   --save_consensus_toms  : save consensus TOMs TRUE/FALSE [default = TRUE]
#   --threads              : WGCNA threads [default = 4]
#   --max_p_outliers       : bicor outlier fraction [default = 0.1]
#
# OUTPUTS
#   comethyl_output/consensus/09_consensus_modules/
#       shared/<cpg_label>/<region_label>/<adjustment_version>/
#           Consensus_Modules.rds
#           Consensus_Module_Distribution.csv
#           Consensus_Region_Assignments.tsv
#           Consensus_Dendrogram_with_Colors.pdf
#           run_parameters.txt
#           sessionInfo.txt
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           <dataset_label>_Consensus_MEs.csv
#           <dataset_label>_Consensus_Module_Eigengene_Dendrogram.pdf
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - Consensus_Modules.rds feeds into scripts 10, 11, 12a, 13, 14
#   - Consensus_Region_Assignments.tsv feeds into scripts 13, 15
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/09_modules_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --dataset1_label Baseline \
#     --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline_Adjusted_Region_Methylation_allPCs.rds \
#     --dataset2_label 36_38wks \
#     --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36_38wks/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/36_38wks_Adjusted_Region_Methylation_allPCs.rds \
#     --dataset3_label Postpartum \
#     --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Postpartum_Adjusted_Region_Methylation_allPCs.rds \
#     --adjustment_version v1_all_pcs \
#     --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/08_soft_power/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/chosen_power.txt \
#     --consensus_cor bicor \
#     --threads 8
# ================================================================
# SCRIPT 09: Consensus Module Detection
#
# PURPOSE
#   - Load 2-3 dataset methylation matrices for the same adjustment version
#   - Ensure region sets are shared and ordered identically
#   - Build multiExpr object for WGCNA consensus analysis
#   - Run blockwiseConsensusModules()
#   - Save shared consensus module object
#   - Save dataset-specific eigengenes and module tables
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/09_consensus_modules/<adjustment_version>/
#       shared/
#           Consensus_Modules.rds
#           Consensus_Module_Distribution.csv
#           Consensus_Region_Assignments.tsv
#           Consensus_Dendrogram_with_Colors.pdf
#           run_parameters.txt
#           sessionInfo.txt
#
#       <dataset_label>/
#           <dataset_label>_Consensus_MEs.csv
#           <dataset_label>_Consensus_Module_Eigengene_Dendrogram.pdf
# ================================================================
message("Starting Script 9 ✓")
suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--regions_file", type = "character"),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_meth", type = "character"),
  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_meth", type = "character", default = NULL),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_meth", type = "character", default = NULL),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--power", type = "integer", default = NULL),
  make_option("--chosen_power_file", type = "character", default = NULL),
  make_option("--consensus_cor", type = "character", default = "bicor"),
  make_option("--network_type", type = "character", default = "signed"),
  make_option("--tom_type", type = "character", default = "signed"),
  make_option("--deep_split", type = "integer", default = 4),
  make_option("--min_module_size", type = "integer", default = 5),
  make_option("--merge_cut_height", type = "double", default = 0.1),
  make_option("--max_block_size", type = "integer", default = 40000),
  make_option("--network_calibration", type = "character", default = "full quantile"),
  make_option("--save_consensus_toms", type = "character", default = "TRUE"),
  make_option("--threads", type = "integer", default = 4),
  make_option("--max_p_outliers", type = "double", default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE")
}

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) stop(label, " must be matrix-like.")
  x <- as.matrix(x)
  if (!is.numeric(x)) stop(label, " must be numeric.")
  if (nrow(x) < 1) stop(label, " has zero rows.")
  if (ncol(x) < 2) stop(label, " must have at least 2 samples.")
  if (is.null(rownames(x))) stop(label, " must have region IDs in rownames.")
  if (is.null(colnames(x))) stop(label, " must have sample IDs in colnames.")
  if (anyDuplicated(rownames(x))) stop(label, " has duplicated region IDs.")
  if (anyDuplicated(colnames(x))) stop(label, " has duplicated sample IDs.")
  x
}

read_chosen_power <- function(file) {
  if (is.null(file) || !nzchar(file) || !file.exists(file)) return(NULL)
  x <- readLines(file, warn = FALSE)
  hit <- grep("^chosen_power:", x, value = TRUE)
  if (length(hit) == 0) return(NULL)
  as.integer(trimws(sub("^chosen_power:\\s*", "", hit[1])))
}

safe_hclust_plot <- function(hc, file, main = "", cex = 0.6, h = NULL) {
  grDevices::pdf(file, width = 10, height = 5)
  par(mar = c(0, 5, 1, 1))
  plot(hc, main = main, xlab = "", sub = "", cex = cex)
  if (!is.null(h)) abline(h = h, col = "red")
  dev.off()
}

# ------------------------------------------------------------
# Validate args
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$regions_file)) stop("--regions_file is required")
if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth)) stop("--dataset1_meth is required")

if (!dir.exists(opt$project_root)) stop("project_root not found")
if (!file.exists(opt$regions_file)) stop("regions_file not found")
if (!file.exists(opt$dataset1_meth)) stop("dataset1_meth not found")

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth)) {
    stop("Provide both dataset2_label and dataset2_meth")
  }
  if (!file.exists(opt$dataset2_meth)) stop("dataset2_meth not found")
}
if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth)) {
    stop("Provide both dataset3_label and dataset3_meth")
  }
  if (!file.exists(opt$dataset3_meth)) stop("dataset3_meth not found")
}

consensus_cor <- tolower(opt$consensus_cor)
if (!consensus_cor %in% c("bicor", "pearson")) stop("--consensus_cor must be bicor or pearson")

save_consensus_toms <- parse_bool(opt$save_consensus_toms, "--save_consensus_toms")

chosen_power <- opt$power
if (is.null(chosen_power)) {
  chosen_power <- read_chosen_power(opt$chosen_power_file)
}
if (is.null(chosen_power) || is.na(chosen_power)) {
  stop("Provide --power or a valid --chosen_power_file")
}

# WGCNA::enableWGCNAThreads(nThreads = opt$threads)
# Force conservative threading behavior
if ("disableWGCNAThreads" %in% getNamespaceExports("WGCNA")) {
  WGCNA::disableWGCNAThreads()
}
if (opt$threads > 1) {
  WGCNA::enableWGCNAThreads(nThreads = opt$threads)
}
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
step_dir      <- file.path(pipeline_root, "09_consensus_modules")
shared_dir    <- file.path(step_dir, "shared", cpg_label, region_label, opt$adjustment_version)
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load regions
# ------------------------------------------------------------
regions <- read.delim(opt$regions_file, stringsAsFactors = FALSE, check.names = FALSE)
required_region_cols <- c("RegionID", "chr", "start", "end")
missing_region_cols <- setdiff(required_region_cols, colnames(regions))
if (length(missing_region_cols) > 0) {
  stop("regions_file missing required columns: ", paste(missing_region_cols, collapse = ", "))
}

# ------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------
dataset_inputs <- list(
  list(label = opt$dataset1_label, meth_file = opt$dataset1_meth)
)
if (dataset2_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset2_label, meth_file = opt$dataset2_meth)
}
if (dataset3_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset3_label, meth_file = opt$dataset3_meth)
}

meth_list <- list()
for (ds in dataset_inputs) {
  m <- validate_meth_matrix(readRDS(ds$meth_file), ds$label)
  meth_list[[ds$label]] <- m
}

# ------------------------------------------------------------
# Ensure shared region set and ordering
# ------------------------------------------------------------
ref_regions <- colnames(meth_list[[dataset_inputs[[1]]$label]])

dim(meth_list[["Baseline"]])
# Should print: 43 11243
# [1]  43 11243   <- rows=samples, cols=regions

for (nm in names(meth_list)) {
  if (!identical(colnames(meth_list[[nm]]), ref_regions)) {
    if (!all(ref_regions %in% colnames(meth_list[[nm]]))) {
      stop("Dataset ", nm, " does not contain all reference regions.")
    }
    meth_list[[nm]] <- meth_list[[nm]][ref_regions, , drop = FALSE]
  }
}
head(ref_regions)
dim(ref_regions)
head(regions)
dim(regions)
regions_use <- regions[match(ref_regions, regions$RegionID), , drop = FALSE]
if (any(is.na(regions_use$RegionID))) {
  stop("Some methylation rownames are not present in regions_file")
}

# ------------------------------------------------------------
# Build multiExpr
# WGCNA expects samples x features
# ------------------------------------------------------------
multiExpr <- lapply(meth_list, function(m) {
  list(data = as.data.frame(m))      # already samples x regions — correct
})

exprSize <- checkSets(multiExpr)
nSets <- exprSize$nSets

# ------------------------------------------------------------
# DIAGNOSTIC: Pairwise correlation check (before module detection)
# ------------------------------------------------------------
# diag_dir <- file.path(step_dir, "diagnostics")
# dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
# message("Diagnostic outputs will be saved to: ", diag_dir)

# for (nm in names(meth_list)) {
#   message("Running correlation diagnostic for: ", nm)
  
#   # Matrix is samples x regions, so cor() directly gives region x region correlations
#   # Subsample regions to avoid memory issues with 11,243 x 11,243 matrix
#   set.seed(42)
#   n_regions_total <- ncol(meth_list[[nm]])
#   sample_cols <- if (n_regions_total > 1000) {
#     sample(n_regions_total, 1000)
#   } else {
#     seq_len(n_regions_total)
#   }
  
#   test_cor <- cor(meth_list[[nm]][, sample_cols], use = "pairwise.complete.obs")
#   upper_vals <- test_cor[upper.tri(test_cor)]

#   grDevices::pdf(
#     file.path(diag_dir, paste0(nm, "_pairwise_region_correlations.pdf")),
#     width = 7, height = 5
#   )
#   hist(
#     upper_vals,
#     breaks = 50,
#     main = paste("Pairwise Region Correlations —", nm),
#     xlab = "Pearson r",
#     col = "steelblue",
#     border = "white"
#   )
#   abline(v = median(upper_vals), col = "red", lwd = 2, lty = 2)
#   legend("topright",
#     legend = paste0("Median r = ", round(median(upper_vals), 3)),
#     col = "red", lty = 2, lwd = 2, bty = "n"
#   )
#   dev.off()

#   # Summary with correct row/col interpretation
#   writeLines(
#     c(
#       paste("Dataset:", nm),
#       paste("N samples:", nrow(meth_list[[nm]])),   # rows = samples
#       paste("N regions:", ncol(meth_list[[nm]])),   # cols = regions
#       paste("N regions sampled for diagnostic:", length(sample_cols)),
#       paste("Median pairwise r:", round(median(upper_vals), 4)),
#       paste("Mean pairwise r:",   round(mean(upper_vals), 4)),
#       paste("SD pairwise r:",     round(sd(upper_vals), 4)),
#       paste("% r > 0.3:", round(mean(upper_vals > 0.3) * 100, 2)),
#       paste("% r > 0.5:", round(mean(upper_vals > 0.5) * 100, 2))
#     ),
#     con = file.path(diag_dir, paste0(nm, "_correlation_summary.txt"))
#   )
# }
# ------------------------------------------------------------
# Run consensus modules
# ------------------------------------------------------------
consensusMods <- blockwiseConsensusModules(
  multiExpr,
  checkMissingData = FALSE,
  maxBlockSize = opt$max_block_size,
  corType = consensus_cor,
  maxPOutliers = opt$max_p_outliers,
  power = chosen_power,
  networkType = opt$network_type,
  checkPower = FALSE,
  TOMType = opt$tom_type,
  networkCalibration = opt$network_calibration,
  saveConsensusTOMs = save_consensus_toms,
  deepSplit = opt$deep_split,
  minModuleSize = opt$min_module_size,
  mergeCutHeight = opt$merge_cut_height,
  nThreads = opt$threads,
  verbose = 5
)

saveRDS(consensusMods, file = file.path(shared_dir, "Consensus_Modules.rds"))

# ------------------------------------------------------------
# Shared outputs
# ------------------------------------------------------------
module_dist <- as.data.frame(sort(table(consensusMods$colors), decreasing = TRUE))
colnames(module_dist) <- c("Module", "Regions")
write.csv(module_dist, file.path(shared_dir, "Consensus_Module_Distribution.csv"), row.names = FALSE)

consensus_regions <- regions_use
consensus_regions$module <- consensusMods$colors
write.table(
  consensus_regions,
  file = file.path(shared_dir, "Consensus_Region_Assignments.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

grDevices::pdf(file.path(shared_dir, "Consensus_Dendrogram_with_Colors.pdf"), width = 10, height = 5)
for (b in seq_along(consensusMods$dendrograms)) {
  block_genes <- consensusMods$blockGenes[[b]]
  block_colors <- consensusMods$colors[block_genes]
  plotDendroAndColors(
    dendro = consensusMods$dendrograms[[b]],
    colors = block_colors,
    groupLabels = "Consensus Modules",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = c(1, 5, 1, 0),
    main = paste("Block", b),
    cex.colorLabels = 1.0
  )
}
dev.off()
# ------------------------------------------------------------
# Dataset-specific outputs
# ------------------------------------------------------------
for (ds in names(consensusMods$multiMEs)) {
  ds_dir <- file.path(step_dir, ds, cpg_label, region_label, opt$adjustment_version)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  me_df <- consensusMods$multiMEs[[ds]]$data
  write.csv(me_df, file = file.path(ds_dir, paste0(ds, "_Consensus_MEs.csv")), row.names = TRUE)

  hc <- if (consensus_cor == "bicor") {
    hclust(as.dist(1 - bicor(me_df, maxPOutliers = opt$max_p_outliers)), method = "average")
  } else {
    hclust(as.dist(1 - cor(me_df, use = "pairwise.complete.obs")), method = "average")
  }

  safe_hclust_plot(
    hc = hc,
    file = file.path(ds_dir, paste0(ds, "_Consensus_Module_Eigengene_Dendrogram.pdf")),
    h = opt$merge_cut_height
  )
}

# ------------------------------------------------------------
# Logs
# ------------------------------------------------------------
writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("regions_file:", opt$regions_file),
    paste("adjustment_version:", opt$adjustment_version),
    paste("datasets:", paste(names(meth_list), collapse = ", ")),
    paste("consensus_cor:", consensus_cor),
    paste("power:", chosen_power),
    paste("network_type:", opt$network_type),
    paste("tom_type:", opt$tom_type),
    paste("deep_split:", opt$deep_split),
    paste("min_module_size:", opt$min_module_size),
    paste("merge_cut_height:", opt$merge_cut_height),
    paste("max_block_size:", opt$max_block_size),
    paste("network_calibration:", opt$network_calibration),
    paste("save_consensus_toms:", save_consensus_toms),
    paste("threads:", opt$threads),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(shared_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()), con = file.path(shared_dir, "sessionInfo.txt"))

message(" Script 09 complete: consensus modules built")