#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08: Module and Sample Diagnostics
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load one or more module objects from Script 07
#   - Optionally align module eigengenes to sample metadata
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
#   --sample_info            : optional metadata xlsx with sample IDs as rownames
#   --module_cor            : correlation for module structure:
#                             bicor or pearson [default = pearson]
#   --sample_dendro_distance: distance for sample dendrogram:
#                             euclidean, pearson, or bicor [default = euclidean]
#   --max_p_outliers        : maxPOutliers for bicor calls [default = 0.1]
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
#
# NOTES
#   - Grey module (MEgrey) is removed before plotting.
#   - All plot calls are wrapped in tryCatch() so one failure does
#     not halt the rest of the script.
#   - Sample-level plots use transpose = TRUE when computing
#     dendrograms and correlations from the ME matrix because
#     MEs are stored as samples x modules.
# ================================================================
message("Starting Script 08")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(openxlsx)
  library(dplyr)
})

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

safe_close_pdf <- function() {
  cur_dev <- grDevices::dev.cur()
  if (!is.null(cur_dev) && cur_dev > 1) {
    try(grDevices::dev.off(), silent = TRUE)
  }
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
  if (is.null(colnames(MEs))) stop(label, "$MEs must have module names in colnames.")

  if (anyDuplicated(rownames(MEs))) {
    dup_ids <- unique(rownames(MEs)[duplicated(rownames(MEs))])
    stop(label, "$MEs has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }

  if (anyDuplicated(colnames(MEs))) {
    dup_mods <- unique(colnames(MEs)[duplicated(colnames(MEs))])
    stop(label, "$MEs has duplicated module names. Example: ",
         paste(head(dup_mods, 10), collapse = ", "))
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

  make_option("--sample_info", type = "character", default = NULL,
              help = "Optional metadata xlsx with sample IDs as rownames for sample alignment/logging"),

  make_option("--module_cor", type = "character", default = "pearson",
              help = "Correlation for module structure: bicor or pearson [default = pearson]"),

  make_option("--sample_dendro_distance", type = "character", default = "euclidean",
              help = "Distance for sample dendrogram: euclidean, pearson, or bicor [default = euclidean]"),

  make_option("--max_p_outliers", type = "double", default = 0.1,
              help = "maxPOutliers for bicor calls [default = 0.1]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$modules_v1))   stop("--modules_v1 is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$modules_v1))  stop("modules_v1 not found: ", opt$modules_v1)

if (!is.null(opt$modules_v2) && !file.exists(opt$modules_v2)) {
  stop("modules_v2 not found: ", opt$modules_v2)
}
if (!is.null(opt$modules_v3) && !file.exists(opt$modules_v3)) {
  stop("modules_v3 not found: ", opt$modules_v3)
}
if (!is.null(opt$sample_info) && !file.exists(opt$sample_info)) {
  stop("sample_info not found: ", opt$sample_info)
}

module_cor <- tolower(opt$module_cor)
if (!module_cor %in% c("pearson", "bicor")) {
  stop("--module_cor must be 'pearson' or 'bicor'")
}

sample_dendro_distance <- tolower(opt$sample_dendro_distance)
if (!sample_dendro_distance %in% c("euclidean", "pearson", "bicor")) {
  stop("--sample_dendro_distance must be one of: euclidean, pearson, bicor")
}

if (!is.numeric(opt$max_p_outliers) || opt$max_p_outliers < 0 || opt$max_p_outliers > 1) {
  stop("--max_p_outliers must be between 0 and 1")
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
# Optional sample metadata
# ------------------------------------------------------------
colData_num <- NULL

if (!is.null(opt$sample_info)) {
  colData <- openxlsx::read.xlsx(opt$sample_info, rowNames = TRUE)
  message("Loaded sample_info: ", opt$sample_info)
  message("colData dimensions: ", nrow(colData), " samples x ", ncol(colData), " traits")

  colData_num <- colData %>%
    dplyr::mutate(dplyr::across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::select(where(~ is.numeric(.) && is.atomic(.)))

  if (nrow(colData_num) == 0) {
    stop("sample_info was loaded but has zero rows after processing.")
  }

  message("colData_num dimensions: ", nrow(colData_num), " samples x ", ncol(colData_num), " numeric traits")
}

# ------------------------------------------------------------
# Derive lineage from modules_v1
# Expected input:
#   .../07_module_detection/<cpg_label>/<region_label>/v1_all_pcs/Modules.rds
# ------------------------------------------------------------
v1_variant_dir <- dirname(opt$modules_v1)
v1_region_dir  <- dirname(v1_variant_dir)
region_label   <- basename(v1_region_dir)
cpg_label      <- basename(dirname(v1_region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir      <- file.path(pipeline_root, "08_module_and_sample_diagnostics")
out_dir       <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Build variant input list
# ------------------------------------------------------------
variant_inputs <- list(v1_all_pcs = opt$modules_v1)

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

  MEs <- as.data.frame(modules$MEs)

  # ----------------------------------------------------------
  # Remove grey module before all plotting
  # ----------------------------------------------------------
  grey_col <- grep("^MEgrey$", colnames(MEs), value = TRUE)
  if (length(grey_col) > 0) {
    MEs <- MEs[, !colnames(MEs) %in% grey_col, drop = FALSE]
    message("[", variant_name, "] MEgrey removed (", ncol(MEs), " real MEs remaining)")
  } else {
    message("[", variant_name, "] No MEgrey column found (", ncol(MEs), " real MEs remaining)")
  }

  variant_out_dir <- file.path(out_dir, variant_name)
  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  if (ncol(MEs) < 2) {
    message("[", variant_name, "] Only ", ncol(MEs), " real ME(s) after grey removal — skipping all plots")
    writeLines(
      paste("Skipped: only", ncol(MEs), "real ME(s) after removing grey."),
      con = file.path(variant_out_dir, paste0(variant_name, "_plots_SKIPPED.txt"))
    )
    next
  }

  # ----------------------------------------------------------
  # Align samples to metadata if provided
  # ----------------------------------------------------------
  if (!is.null(colData_num)) {
    common_samples <- intersect(rownames(MEs), rownames(colData_num))

    if (length(common_samples) == 0) {
      stop("[", variant_name, "] No overlapping samples between module eigengenes and sample_info metadata.")
    }

    MEs_use <- MEs[common_samples, , drop = FALSE]
    colData_use <- colData_num[common_samples, , drop = FALSE]

    message("[", variant_name, "] Sample alignment to metadata: ",
            length(common_samples), " overlapping samples retained from ",
            nrow(MEs), " ME samples and ", nrow(colData_num), " metadata samples")
  } else {
    MEs_use <- MEs
    colData_use <- NULL
    common_samples <- rownames(MEs_use)
    message("[", variant_name, "] No sample_info provided; using all ", nrow(MEs_use), " ME samples")
  }

  if (nrow(MEs_use) < 2) {
    message("[", variant_name, "] Fewer than 2 samples after alignment — skipping plots")
    writeLines(
      "Skipped: fewer than 2 samples after metadata alignment.",
      con = file.path(variant_out_dir, paste0(variant_name, "_plots_SKIPPED.txt"))
    )
    next
  }

  message("[", variant_name, "] MEs used: ", nrow(MEs_use), " samples x ", ncol(MEs_use), " modules")

  f_module_me_dendro <- file.path(variant_out_dir, "Module_ME_Dendrogram.pdf")
  f_module_cor_hm    <- file.path(variant_out_dir, "Module_Correlation_Heatmap.pdf")
  f_module_stats     <- file.path(variant_out_dir, "Module_Correlation_Stats.tsv")
  f_sample_me_dendro <- file.path(variant_out_dir, "Sample_ME_Dendrogram.pdf")
  f_sample_cor_hm    <- file.path(variant_out_dir, "Sample_Correlation_Heatmap.pdf")
  f_sample_me_hm     <- file.path(variant_out_dir, "Sample_ME_Heatmap.pdf")

  # ==========================================================
  # MODULE-LEVEL diagnostics
  # ==========================================================
  moduleDendro <- tryCatch({
    getDendro(MEs_use, distance = module_cor)
  }, error = function(e) {
    message("[", variant_name, "] Module dendrogram failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(moduleDendro)) {
    tryCatch({
      plotDendro(moduleDendro, labelSize = 4, nBreaks = 5, file = f_module_me_dendro)
      message("[", variant_name, "] Saved module ME dendrogram")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed saving module dendrogram: ", conditionMessage(e))
    })
  }

  moduleCor <- tryCatch({
    getCor(MEs_use, corType = module_cor, maxPOutliers = opt$max_p_outliers)
  }, error = function(e) {
    message("[", variant_name, "] Module correlation failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(moduleCor) && !is.null(moduleDendro)) {
    tryCatch({
      plotHeatmap(
        moduleCor,
        rowDendro = moduleDendro,
        colDendro = moduleDendro,
        file      = f_module_cor_hm
      )
      message("[", variant_name, "] Saved module correlation heatmap")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed module correlation heatmap: ", conditionMessage(e))
    })
  }

  tryCatch({
    getMEtraitCor(
      MEs_use,
      colData  = MEs_use,
      corType  = module_cor,
      robustY  = TRUE,
      file     = f_module_stats
    )
    message("[", variant_name, "] Saved module correlation stats")
  }, error = function(e) {
    message("[", variant_name, "] Failed module correlation stats: ", conditionMessage(e))
    writeLines(
      paste("Failed:", conditionMessage(e)),
      con = file.path(variant_out_dir, "Module_Correlation_Stats_FAILED.txt")
    )
  })

  # ==========================================================
  # SAMPLE-LEVEL diagnostics
  # IMPORTANT: transpose = TRUE because MEs are samples x modules
  # ==========================================================
  sampleDendro <- tryCatch({
    getDendro(MEs_use, transpose = TRUE, distance = sample_dendro_distance)
  }, error = function(e) {
    message("[", variant_name, "] Sample dendrogram failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(sampleDendro)) {
    tryCatch({
      plotDendro(sampleDendro, labelSize = 3, nBreaks = 5, file = f_sample_me_dendro)
      message("[", variant_name, "] Saved sample ME dendrogram")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed saving sample dendrogram: ", conditionMessage(e))
    })
  }

  sampleCor <- tryCatch({
    getCor(MEs_use, transpose = TRUE, corType = module_cor, maxPOutliers = opt$max_p_outliers)
  }, error = function(e) {
    message("[", variant_name, "] Sample correlation failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(sampleCor) && !is.null(sampleDendro)) {
    tryCatch({
      plotHeatmap(
        sampleCor,
        rowDendro = sampleDendro,
        colDendro = sampleDendro,
        file      = f_sample_cor_hm
      )
      message("[", variant_name, "] Saved sample correlation heatmap")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed sample correlation heatmap: ", conditionMessage(e))
    })
  }

  tryCatch({
    plotHeatmap(
      MEs_use,
      rowDendro       = sampleDendro,
      colDendro       = moduleDendro,
      legend.title    = "Module\nEigengene",
      legend.position = c(0.37, 0.89),
      file            = f_sample_me_hm
    )
    message("[", variant_name, "] Saved sample x module eigengene heatmap")
  }, error = function(e) {
    safe_close_pdf()
    message("[", variant_name, "] Sample x module heatmap failed: ", conditionMessage(e))
    writeLines(
      paste("Failed:", conditionMessage(e)),
      con = file.path(variant_out_dir, "Sample_ME_Heatmap_FAILED.txt")
    )
  })

  # ----------------------------------------------------------
  # Per-variant run parameters
  # ----------------------------------------------------------
  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("modules_file:", variant_inputs[[variant_name]]),
      paste("sample_info:", ifelse(is.null(opt$sample_info), "NULL", opt$sample_info)),
      paste("module_cor:", module_cor),
      paste("sample_dendro_distance:", sample_dendro_distance),
      paste("max_p_outliers:", opt$max_p_outliers),
      paste("n_samples_before_metadata_alignment:", nrow(MEs)),
      paste("n_samples_after_metadata_alignment:", nrow(MEs_use)),
      paste("n_traits_numeric:", ifelse(is.null(colData_use), 0, ncol(colData_use))),
      paste("n_modules:", ncol(MEs_use)),
      paste("module_me_dendrogram:", f_module_me_dendro),
      paste("module_correlation_heatmap:", f_module_cor_hm),
      paste("module_correlation_stats:", f_module_stats),
      paste("sample_me_dendrogram:", f_sample_me_dendro),
      paste("sample_correlation_heatmap:", f_sample_cor_hm),
      paste("sample_me_heatmap:", f_sample_me_hm),
      paste("date:", as.character(Sys.time()))
    ),
    file.path(variant_out_dir, "run_parameters.txt")
  )

  message("Finished variant: ", variant_name)
  message("  Outputs in: ", variant_out_dir)
}

# ------------------------------------------------------------
# Top-level run parameters
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("modules_v1:", opt$modules_v1),
    paste("modules_v2:", ifelse(is.null(opt$modules_v2), "NULL", opt$modules_v2)),
    paste("modules_v3:", ifelse(is.null(opt$modules_v3), "NULL", opt$modules_v3)),
    paste("sample_info:", ifelse(is.null(opt$sample_info), "NULL", opt$sample_info)),
    paste("module_cor:", module_cor),
    paste("sample_dendro_distance:", sample_dendro_distance),
    paste("max_p_outliers:", opt$max_p_outliers),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("Script 08 complete: module and sample diagnostics finished")
message("Outputs saved under:\n  ", out_dir)