#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 07: Module Detection
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load filtered genomic regions
#   - Load one or more adjusted region-methylation matrices from Script 05
#   - Load one or more soft-threshold power objects from Script 06
#   - Run comethyl::getModules() for each available variant
#   - Save module objects, dendrograms, and BED files
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --regions_file : path to filtered regions file
#   --input_v1     : path to v1 adjusted methylation RDS from Script 05
#   --softpower_v1 : path to v1 soft power RDS from Script 06
#
# OPTIONAL INPUTS
#   --input_v2       : path to v2 adjusted methylation RDS
#   --softpower_v2   : path to v2 soft power RDS
#   --input_v3       : path to v3 adjusted methylation RDS
#   --softpower_v3   : path to v3 soft power RDS
#   --module_cor     : correlation used for module detection: bicor or pearson
#                      [default = pearson]
#   --deep_split     : dynamic tree cut deepSplit [default = 4]
#   --min_module_size: minimum module size [default = 10]
#   --merge_cut_height : merge cut height [default = 0.1]
#   --max_block_size : max block size [default = 40000]
#   --threads        : number of WGCNA threads [default = 4]
#
# OUTPUTS
#   project_root/comethyl_output/07_module_detection/<cpg_label>/<region_label>/<variant>/
#       Modules.rds
#       Region_Dendrograms.pdf
#       Modules.bed
#       run_parameters.txt
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(data.table)
})

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--regions_file", type = "character",
              help = "Path to filtered regions file"),

  make_option("--input_v1", type = "character",
              help = "Path to v1 adjusted methylation RDS"),

  make_option("--softpower_v1", type = "character",
              help = "Path to v1 soft power RDS"),

  make_option("--input_v2", type = "character", default = NULL,
              help = "Optional path to v2 adjusted methylation RDS"),

  make_option("--softpower_v2", type = "character", default = NULL,
              help = "Optional path to v2 soft power RDS"),

  make_option("--input_v3", type = "character", default = NULL,
              help = "Optional path to v3 adjusted methylation RDS"),

  make_option("--softpower_v3", type = "character", default = NULL,
              help = "Optional path to v3 soft power RDS"),

  make_option("--module_cor", type = "character", default = "pearson",
              help = "Correlation for module detection: bicor or pearson [default = pearson]"),

  make_option("--deep_split", type = "integer", default = 4,
              help = "deepSplit for getModules() [default = 4]"),

  make_option("--min_module_size", type = "integer", default = 10,
              help = "Minimum module size [default = 10]"),

  make_option("--merge_cut_height", type = "double", default = 0.1,
              help = "Merge cut height [default = 0.1]"),

  make_option("--max_block_size", type = "integer", default = 40000,
              help = "Maximum block size [default = 40000]"),

  make_option("--threads", type = "integer", default = 4,
              help = "Number of WGCNA threads [default = 4]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(label, " must be a matrix-like object.")
  }

  x <- as.matrix(x)

  if (!is.numeric(x)) stop(label, " must be numeric.")
  if (nrow(x) < 1) stop(label, " has zero rows.")
  if (ncol(x) < 2) stop(label, " must have at least 2 samples.")
  if (is.null(colnames(x))) stop(label, " must have sample IDs in colnames.")
  if (anyDuplicated(colnames(x))) {
    dup_ids <- unique(colnames(x)[duplicated(colnames(x))])
    stop(label, " has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }

  x
}

validate_softpower <- function(x, label) {
  if (!is.list(x)) stop(label, " must be a list-like soft power object.")
  if (is.null(x$powerEstimate) || is.na(x$powerEstimate)) {
    stop(label, " has missing powerEstimate.")
  }
  x
}

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$regions_file)) stop("--regions_file is required")
if (is.null(opt$input_v1)) stop("--input_v1 is required")
if (is.null(opt$softpower_v1)) stop("--softpower_v1 is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$regions_file)) stop("regions_file not found: ", opt$regions_file)
if (!file.exists(opt$input_v1)) stop("input_v1 not found: ", opt$input_v1)
if (!file.exists(opt$softpower_v1)) stop("softpower_v1 not found: ", opt$softpower_v1)

if (!is.null(opt$input_v2) && !file.exists(opt$input_v2)) stop("input_v2 not found: ", opt$input_v2)
if (!is.null(opt$softpower_v2) && !file.exists(opt$softpower_v2)) stop("softpower_v2 not found: ", opt$softpower_v2)
if (!is.null(opt$input_v3) && !file.exists(opt$input_v3)) stop("input_v3 not found: ", opt$input_v3)
if (!is.null(opt$softpower_v3) && !file.exists(opt$softpower_v3)) stop("softpower_v3 not found: ", opt$softpower_v3)

if (tolower(opt$module_cor) %in% c("bicor", "pearson") == FALSE) {
  stop("--module_cor must be 'bicor' or 'pearson'")
}
if (opt$deep_split < 0) stop("--deep_split must be >= 0")
if (opt$min_module_size < 2) stop("--min_module_size must be >= 2")
if (opt$merge_cut_height <= 0 || opt$merge_cut_height >= 1) {
  stop("--merge_cut_height must be > 0 and < 1")
}
if (opt$max_block_size < 1) stop("--max_block_size must be >= 1")
if (opt$threads < 1) stop("--threads must be >= 1")

if (!is.null(opt$input_v2) && is.null(opt$softpower_v2)) {
  stop("If --input_v2 is provided, --softpower_v2 must also be provided.")
}
if (is.null(opt$input_v2) && !is.null(opt$softpower_v2)) {
  stop("If --softpower_v2 is provided, --input_v2 must also be provided.")
}
if (!is.null(opt$input_v3) && is.null(opt$softpower_v3)) {
  stop("If --input_v3 is provided, --softpower_v3 must also be provided.")
}
if (is.null(opt$input_v3) && !is.null(opt$softpower_v3)) {
  stop("If --softpower_v3 is provided, --input_v3 must also be provided.")
}

module_cor <- tolower(opt$module_cor)

# ------------------------------------------------------------
# Configure cache and threads
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ------------------------------------------------------------
# Derive lineage from input_v1
# Expected input:
#   .../05_methylation_adjustment/<cpg_label>/<region_label>/v1_all_pcs/<file>.rds
# ------------------------------------------------------------
v1_variant_dir <- dirname(opt$input_v1)
v1_region_dir <- dirname(v1_variant_dir)
region_label <- basename(v1_region_dir)
cpg_label <- basename(dirname(v1_region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir <- file.path(pipeline_root, "07_module_detection")
out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Load and validate regions
# ------------------------------------------------------------
regions <- data.table::fread(opt$regions_file, data.table = FALSE)

required_region_cols <- c("RegionID")
missing_region_cols <- setdiff(required_region_cols, colnames(regions))
if (length(missing_region_cols) > 0) {
  stop("Regions file is missing required columns: ",
       paste(missing_region_cols, collapse = ", "))
}

message("Loaded regions: ", nrow(regions), " rows")

# ------------------------------------------------------------
# Build variant input list
# ------------------------------------------------------------
variant_inputs <- list(
  v1_all_pcs = list(
    meth_file = opt$input_v1,
    softpower_file = opt$softpower_v1
  )
)

if (!is.null(opt$input_v2)) {
  variant_inputs[[basename(dirname(opt$input_v2))]] <- list(
    meth_file = opt$input_v2,
    softpower_file = opt$softpower_v2
  )
}

if (!is.null(opt$input_v3)) {
  variant_inputs[[basename(dirname(opt$input_v3))]] <- list(
    meth_file = opt$input_v3,
    softpower_file = opt$softpower_v3
  )
}

# ------------------------------------------------------------
# Run module detection for each variant
# ------------------------------------------------------------
for (variant_name in names(variant_inputs)) {
  message("\n==============================")
  message("Running module detection for variant: ", variant_name)
  message("==============================\n")

  methAdj <- validate_meth_matrix(
    readRDS(variant_inputs[[variant_name]]$meth_file),
    paste0(variant_name, " adjusted methylation")
  )

  sft <- validate_softpower(
    readRDS(variant_inputs[[variant_name]]$softpower_file),
    paste0(variant_name, " soft power object")
  )

  variant_out_dir <- file.path(out_dir, variant_name)
  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  modules_rds <- file.path(variant_out_dir, "Modules.rds")
  dendro_pdf <- file.path(variant_out_dir, "Region_Dendrograms.pdf")
  bed_file <- file.path(variant_out_dir, "Modules.bed")

  message("[", variant_name, "] methylation dimensions: ",
          nrow(methAdj), " regions x ", ncol(methAdj), " samples")
  message("[", variant_name, "] using powerEstimate = ", sft$powerEstimate)
  message("[", variant_name, "] module correlation = ", module_cor)

  modules <- getModules(
    meth = methAdj,
    power = sft$powerEstimate,
    regions = regions,
    corType = module_cor,
    deepSplit = opt$deep_split,
    minModuleSize = opt$min_module_size,
    mergeCutHeight = opt$merge_cut_height,
    maxBlockSize = opt$max_block_size,
    nThreads = opt$threads,
    file = modules_rds,
    verbose = TRUE
  )

  plotRegionDendro(modules, file = dendro_pdf)
  getModuleBED(modules$regions, file = bed_file)

  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("meth_file:", variant_inputs[[variant_name]]$meth_file),
      paste("softpower_file:", variant_inputs[[variant_name]]$softpower_file),
      paste("regions_file:", opt$regions_file),
      paste("module_cor:", module_cor),
      paste("power_estimate:", sft$powerEstimate),
      paste("deep_split:", opt$deep_split),
      paste("min_module_size:", opt$min_module_size),
      paste("merge_cut_height:", opt$merge_cut_height),
      paste("max_block_size:", opt$max_block_size),
      paste("threads:", opt$threads),
      paste("n_regions_input:", nrow(methAdj)),
      paste("n_samples:", ncol(methAdj)),
      paste("modules_rds:", modules_rds),
      paste("dendrogram_pdf:", dendro_pdf),
      paste("bed_file:", bed_file)
    ),
    file.path(variant_out_dir, "run_parameters.txt")
  )

  message("✓ Finished variant: ", variant_name)
  message("  Modules RDS: ", modules_rds)
  message("  Dendrogram:  ", dendro_pdf)
  message("  BED:         ", bed_file)
}

# ------------------------------------------------------------
# Save top-level run parameters
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("regions_file:", opt$regions_file),
    paste("input_v1:", opt$input_v1),
    paste("softpower_v1:", opt$softpower_v1),
    paste("input_v2:", ifelse(is.null(opt$input_v2), "NULL", opt$input_v2)),
    paste("softpower_v2:", ifelse(is.null(opt$softpower_v2), "NULL", opt$softpower_v2)),
    paste("input_v3:", ifelse(is.null(opt$input_v3), "NULL", opt$input_v3)),
    paste("softpower_v3:", ifelse(is.null(opt$softpower_v3), "NULL", opt$softpower_v3)),
    paste("module_cor:", module_cor),
    paste("deep_split:", opt$deep_split),
    paste("min_module_size:", opt$min_module_size),
    paste("merge_cut_height:", opt$merge_cut_height),
    paste("max_block_size:", opt$max_block_size),
    paste("threads:", opt$threads),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL MODULE DETECTION RUNS COMPLETE ✓")
message("Outputs saved under:\n  ", out_dir)