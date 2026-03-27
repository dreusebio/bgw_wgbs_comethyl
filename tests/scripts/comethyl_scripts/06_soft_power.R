#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 06: Soft-Threshold Power Analysis
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load one or more adjusted region-methylation matrices from Script 05
#   - Run comethyl::getSoftPower() for each available adjustment variant
#   - Generate soft-threshold power plots
#   - Optionally generate sample dendrograms
#   - Save outputs in a reproducible step-specific folder structure
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --input_v1     : path to v1 adjusted methylation RDS from Script 05
#
# OPTIONAL INPUTS
#   --input_v2            : path to v2 adjusted methylation RDS from Script 05
#   --input_v3            : path to v3 adjusted methylation RDS from Script 05
#   --softpower_cor       : correlation for getSoftPower(): bicor or pearson
#                           [default = pearson]
#   --power_min           : minimum soft-threshold power [default = 1]
#   --power_max           : maximum soft-threshold power [default = 20]
#   --block_size          : block size passed to getSoftPower() [default = 20000]
#   --gc_interval         : gcInterval passed to getSoftPower() [default = 19999]
#   --threads             : number of WGCNA threads [default = 4]
#   --plot_sample_dendro  : TRUE/FALSE, whether to generate sample dendrograms
#                           [default = TRUE]
#
# OUTPUTS
#   project_root/comethyl_output/06_soft_power/<cpg_label>/<region_label>/<variant>/
#       SoftPower_<softpower_cor>.rds
#       SoftPower_<softpower_cor>_Plots.pdf
#       Sample_Dendrogram_<softpower_cor>.pdf      (optional)
#       run_parameters.txt
#
# NOTES
#   - v1 is required because it is the default adjustment output in Script 05.
#   - v2 and v3 are optional and will only run if input files are provided.
#   - Variant names are derived from the input file paths.
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
})

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--input_v1", type = "character",
              help = "Path to v1 adjusted methylation RDS from Script 05"),

  make_option("--input_v2", type = "character", default = NULL,
              help = "Optional path to v2 adjusted methylation RDS from Script 05"),

  make_option("--input_v3", type = "character", default = NULL,
              help = "Optional path to v3 adjusted methylation RDS from Script 05"),

  make_option("--softpower_cor", type = "character", default = "pearson",
              help = "Correlation for getSoftPower(): bicor or pearson [default = pearson]"),

  make_option("--power_min", type = "integer", default = 1,
              help = "Minimum soft-threshold power [default = 1]"),

  make_option("--power_max", type = "integer", default = 20,
              help = "Maximum soft-threshold power [default = 20]"),

  make_option("--block_size", type = "integer", default = 20000,
              help = "blockSize for getSoftPower() [default = 20000]"),

  make_option("--gc_interval", type = "integer", default = 19999,
              help = "gcInterval for getSoftPower() [default = 19999]"),

  make_option("--threads", type = "integer", default = 4,
              help = "Number of WGCNA threads [default = 4]"),

  make_option("--plot_sample_dendro", type = "character", default = "TRUE",
              help = "Generate sample dendrograms: TRUE or FALSE [default = TRUE]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  if (is.logical(x)) return(x)
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(label, " must be a matrix-like object (matrix/data.frame).")
  }

  x <- as.matrix(x)

  if (!is.numeric(x)) stop(label, " must be numeric.")
  if (nrow(x) < 1) stop(label, " has zero rows.")
  if (ncol(x) < 2) stop(label, " must have at least 2 samples (columns).")
  if (is.null(colnames(x))) stop(label, " must have sample IDs in colnames.")
  if (anyDuplicated(colnames(x))) {
    dup_ids <- unique(colnames(x)[duplicated(colnames(x))])
    stop(label, " has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }

  return(x)
}

run_softpower <- function(methAdj,
                          variant_name,
                          variant_out_dir,
                          powerVector,
                          corType = "pearson",
                          blockSize = 20000,
                          gcInterval = 19999,
                          plotSampleDendro = TRUE,
                          input_file = NULL) {

  message("\n==============================")
  message("Running soft power for variant: ", variant_name)
  message("Correlation method: ", corType)
  message("==============================\n")

  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  rds_file <- file.path(variant_out_dir, paste0("SoftPower_", corType, ".rds"))
  pdf_file <- file.path(variant_out_dir, paste0("SoftPower_", corType, "_Plots.pdf"))
  dendro_file <- file.path(variant_out_dir, paste0("Sample_Dendrogram_", corType, ".pdf"))

  sft <- getSoftPower(
    methAdj,
    powerVector = powerVector,
    corType = corType,
    file = rds_file,
    blockSize = blockSize,
    gcInterval = gcInterval
  )

  plotSoftPower(sft, file = pdf_file)

  dendro_status <- "not_requested"

  if (isTRUE(plotSampleDendro)) {
    dendro_status <- "attempted"
    ok <- TRUE

    tryCatch({
      d <- getDendro(methAdj, distance = "euclidean")
      plotDendro(d, file = dendro_file, expandY = c(0.25, 0.08))
    }, error = function(e) {
      ok <<- FALSE
      dendro_status <<- paste0("failed: ", conditionMessage(e))
      warning("Sample dendrogram failed for ", variant_name, ". Continuing. Error: ", conditionMessage(e))
    })

    if (ok) {
      dendro_status <- "completed"
      message("✓ Sample dendrogram saved: ", dendro_file)
    }
  }

  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("input_file:", ifelse(is.null(input_file), "NULL", input_file)),
      paste("softpower_cor:", corType),
      paste("power_vector:", paste(powerVector, collapse = ",")),
      paste("block_size:", blockSize),
      paste("gc_interval:", gcInterval),
      paste("n_regions:", nrow(methAdj)),
      paste("n_samples:", ncol(methAdj)),
      paste("plot_sample_dendro:", plotSampleDendro),
      paste("sample_dendrogram_status:", dendro_status),
      paste("softpower_rds:", rds_file),
      paste("softpower_plot:", pdf_file),
      paste("sample_dendrogram_file:", ifelse(isTRUE(plotSampleDendro), dendro_file, "not_requested"))
    ),
    file.path(variant_out_dir, "run_parameters.txt")
  )

  message("✓ Finished soft-power run for: ", variant_name)
  message("  RDS: ", rds_file)
  message("  PDF: ", pdf_file)

  invisible(sft)
}

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$input_v1)) stop("--input_v1 is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$input_v1)) stop("input_v1 not found: ", opt$input_v1)
if (!is.null(opt$input_v2) && !file.exists(opt$input_v2)) stop("input_v2 not found: ", opt$input_v2)
if (!is.null(opt$input_v3) && !file.exists(opt$input_v3)) stop("input_v3 not found: ", opt$input_v3)

softpower_cor <- tolower(opt$softpower_cor)
if (!softpower_cor %in% c("bicor", "pearson")) {
  stop("--softpower_cor must be 'bicor' or 'pearson'")
}
if (opt$power_min < 1) stop("--power_min must be >= 1")
if (opt$power_max < opt$power_min) stop("--power_max must be >= --power_min")
if (opt$block_size < 1) stop("--block_size must be >= 1")
if (opt$gc_interval < 0) stop("--gc_interval must be >= 0")
if (opt$threads < 1) stop("--threads must be >= 1")

plot_sample_dendro <- parse_bool(opt$plot_sample_dendro, "--plot_sample_dendro")
powerVector <- seq(opt$power_min, opt$power_max)

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
step_dir <- file.path(pipeline_root, "06_soft_power")
out_dir <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Load and validate adjusted methylation inputs
# ------------------------------------------------------------
methAdj_v1 <- validate_meth_matrix(readRDS(opt$input_v1), "v1 adjusted methylation")

variant_inputs <- list(
  v1_all_pcs = list(
    meth = methAdj_v1,
    input_file = opt$input_v1
  )
)

if (!is.null(opt$input_v2)) {
  methAdj_v2 <- validate_meth_matrix(readRDS(opt$input_v2), "v2 adjusted methylation")
  variant_inputs[[basename(dirname(opt$input_v2))]] <- list(
    meth = methAdj_v2,
    input_file = opt$input_v2
  )
}

if (!is.null(opt$input_v3)) {
  methAdj_v3 <- validate_meth_matrix(readRDS(opt$input_v3), "v3 adjusted methylation")
  variant_inputs[[basename(dirname(opt$input_v3))]] <- list(
    meth = methAdj_v3,
    input_file = opt$input_v3
  )
}

# ------------------------------------------------------------
# Run soft-power for each available variant
# ------------------------------------------------------------
for (variant_name in names(variant_inputs)) {
  variant_out_dir <- file.path(out_dir, variant_name)

  run_softpower(
    methAdj = variant_inputs[[variant_name]]$meth,
    variant_name = variant_name,
    variant_out_dir = variant_out_dir,
    powerVector = powerVector,
    corType = softpower_cor,
    blockSize = opt$block_size,
    gcInterval = opt$gc_interval,
    plotSampleDendro = plot_sample_dendro,
    input_file = variant_inputs[[variant_name]]$input_file
  )
}

# ------------------------------------------------------------
# Save top-level run parameters
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("input_v1:", opt$input_v1),
    paste("input_v2:", ifelse(is.null(opt$input_v2), "NULL", opt$input_v2)),
    paste("input_v3:", ifelse(is.null(opt$input_v3), "NULL", opt$input_v3)),
    paste("softpower_cor:", softpower_cor),
    paste("power_min:", opt$power_min),
    paste("power_max:", opt$power_max),
    paste("power_vector:", paste(powerVector, collapse = ",")),
    paste("block_size:", opt$block_size),
    paste("gc_interval:", opt$gc_interval),
    paste("threads:", opt$threads),
    paste("plot_sample_dendro:", plot_sample_dendro),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("\nALL SOFT-POWER RUNS COMPLETE (", softpower_cor, ") ✓")
message("Outputs saved under:\n  ", out_dir)