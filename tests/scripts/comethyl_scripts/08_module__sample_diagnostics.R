#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08: Module and Sample Diagnostics
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load one or more module objects from Script 07
#   - Generate:
#       1) module-module structure diagnostics from module eigengenes
#       2) sample-sample structure diagnostics based on module eigengenes
#       3) sample-by-module eigengene heatmap
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --modules_v1   : path to v1 Modules.rds from Script 07
#
# OPTIONAL INPUTS
#   --modules_v2            : optional path to v2 Modules.rds
#   --modules_v3            : optional path to v3 Modules.rds
#   --module_cor            : correlation for module structure:
#                             bicor or pearson [default = pearson]
#   --sample_dendro_distance: distance for sample dendrogram:
#                             euclidean, pearson, or bicor [default = euclidean]
#
# OUTPUTS
#   project_root/comethyl_output/08_module_and_sample_diagnostics/<cpg_label>/<region_label>/<variant>/
#       Module_ME_Dendrogram.pdf
#       Module_Correlation_Heatmap.pdf
#       Module_Correlation_Stats.tsv
#       Sample_ME_Dendrogram.pdf
#       Sample_Correlation_Heatmap.pdf
#       Sample_ME_Heatmap.pdf
#       run_parameters.txt
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
})

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

validate_modules_object <- function(x, label) {
  if (!is.list(x)) stop(label, " must be a list-like module object.")
  if (is.null(x$MEs)) stop(label, " is missing $MEs.")
  if (!(is.matrix(x$MEs) || is.data.frame(x$MEs))) stop(label, "$MEs must be matrix-like.")

  MEs <- as.matrix(x$MEs)
  if (!is.numeric(MEs)) stop(label, "$MEs must be numeric.")
  if (nrow(MEs) < 2) stop(label, "$MEs must have at least 2 samples.")
  if (ncol(MEs) < 1) stop(label, "$MEs must have at least 1 module.")
  if (is.null(rownames(MEs))) stop(label, "$MEs must have sample IDs in rownames.")
  if (anyDuplicated(rownames(MEs))) {
    dup_ids <- unique(rownames(MEs)[duplicated(rownames(MEs))])
    stop(label, "$MEs has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }

  x$MEs <- MEs
  x
}

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--modules_v1", type = "character",
              help = "Path to v1 Modules.rds from Script 07"),

  make_option("--modules_v2", type = "character", default = NULL,
              help = "Optional path to v2 Modules.rds from Script 07"),

  make_option("--modules_v3", type = "character", default = NULL,
              help = "Optional path to v3 Modules.rds from Script 07"),

  make_option("--module_cor", type = "character", default = "pearson",
              help = "Correlation for module structure: bicor or pearson [default = pearson]"),

  make_option("--sample_dendro_distance", type = "character", default = "euclidean",
              help = "Distance for sample dendrogram: euclidean, pearson, or bicor [default = euclidean]")
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

module_cor <- tolower(opt$module_cor)
if (!module_cor %in% c("pearson", "bicor")) {
  stop("--module_cor must be 'pearson' or 'bicor'")
}

sample_dendro_distance <- tolower(opt$sample_dendro_distance)
if (!sample_dendro_distance %in% c("euclidean", "pearson", "bicor")) {
  stop("--sample_dendro_distance must be one of: euclidean, pearson, bicor")
}

# ------------------------------------------------------------
# Configure cache and threads
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Derive lineage from modules_v1
# Expected input:
#   .../07_module_detection/<cpg_label>/<region_label>/v1_all_pcs/Modules.rds
# ------------------------------------------------------------
v1_variant_dir <- dirname(opt$modules_v1)
v1_region_dir <- dirname(v1_variant_dir)
region_label <- basename(v1_region_dir)
cpg_label <- basename(dirname(v1_region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "08_module_and_sample_diagnostics")
out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Build variant input list
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
  message("Running module/sample diagnostics for variant: ", variant_name)
  message("==============================\n")

  modules <- validate_modules_object(
    readRDS(variant_inputs[[variant_name]]),
    paste0(variant_name, " modules object")
  )

  MEs <- modules$MEs
  message("[", variant_name, "] MEs dimensions: ", nrow(MEs), " samples x ", ncol(MEs), " modules")

  variant_out_dir <- file.path(out_dir, variant_name)
  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  f_module_me_dendro <- file.path(variant_out_dir, "Module_ME_Dendrogram.pdf")
  f_module_cor_hm <- file.path(variant_out_dir, "Module_Correlation_Heatmap.pdf")
  f_module_stats <- file.path(variant_out_dir, "Module_Correlation_Stats.tsv")

  f_sample_me_dendro <- file.path(variant_out_dir, "Sample_ME_Dendrogram.pdf")
  f_sample_cor_hm <- file.path(variant_out_dir, "Sample_Correlation_Heatmap.pdf")
  f_sample_me_hm <- file.path(variant_out_dir, "Sample_ME_Heatmap.pdf")

  # Module-module structure
  moduleDendro <- getDendro(MEs, distance = module_cor)
  plotDendro(
    moduleDendro,
    labelSize = 4,
    nBreaks = 5,
    file = f_module_me_dendro
  )

  moduleCor <- getCor(MEs, corType = module_cor)
  plotHeatmap(
    moduleCor,
    rowDendro = moduleDendro,
    colDendro = moduleDendro,
    file = f_module_cor_hm
  )

  moduleCorStats <- getMEtraitCor(
    MEs,
    colData = MEs,
    corType = module_cor,
    robustY = TRUE,
    file = f_module_stats
  )

  # Sample-sample structure
  sampleDendro <- getDendro(MEs, transpose = TRUE, distance = sample_dendro_distance)
  plotDendro(
    sampleDendro,
    labelSize = 3,
    nBreaks = 5,
    file = f_sample_me_dendro
  )

  sampleCor <- getCor(MEs, transpose = TRUE, corType = module_cor)
  plotHeatmap(
    sampleCor,
    rowDendro = sampleDendro,
    colDendro = sampleDendro,
    file = f_sample_cor_hm
  )

  plotHeatmap(
    MEs,
    rowDendro = sampleDendro,
    colDendro = moduleDendro,
    legend.title = "Module\nEigennode",
    legend.position = c(0.37, 0.89),
    file = f_sample_me_hm
  )

  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("modules_file:", variant_inputs[[variant_name]]),
      paste("module_cor:", module_cor),
      paste("sample_dendro_distance:", sample_dendro_distance),
      paste("n_samples:", nrow(MEs)),
      paste("n_modules:", ncol(MEs)),
      paste("module_me_dendrogram:", f_module_me_dendro),
      paste("module_correlation_heatmap:", f_module_cor_hm),
      paste("module_correlation_stats:", f_module_stats),
      paste("sample_me_dendrogram:", f_sample_me_dendro),
      paste("sample_correlation_heatmap:", f_sample_cor_hm),
      paste("sample_me_heatmap:", f_sample_me_hm)
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
    paste("module_cor:", module_cor),
    paste("sample_dendro_distance:", sample_dendro_distance),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL MODULE/SAMPLE DIAGNOSTIC RUNS COMPLETE ✓")
message("Outputs saved under:\n  ", out_dir)