#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 10: Consensus Module Diagnostics
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Loads the consensus module object from script 09 and generates
#   a suite of diagnostic plots and tables:
#     1) Shared consensus eigengene network plot (across all datasets)
#     2) Per-dataset ME correlation heatmaps
#     3) Per-dataset sample dendrograms from consensus MEs
#     4) Per-dataset module dendrograms from consensus MEs
#     5) Per-dataset sample correlation heatmaps
#     6) Per-dataset sample × module eigengene heatmaps
#     7) Module-module correlation statistics tables
#
# REQUIRED INPUTS
#   --project_root          : root directory of the analysis project
#   --consensus_modules_rds : path to Consensus_Modules.rds from script 09
#   --regions_file          : path to Filtered_Regions.txt from script 02
#                             (used to derive cpg_label and region_label)
#
# OPTIONAL INPUTS
#   --adjustment_version : label matching script 09 run [default = unadjusted]
#   --consensus_cor      : bicor or pearson [default = bicor]
#   --max_p_outliers     : bicor outlier fraction [default = 0.1]
#
# OUTPUTS
#   comethyl_output/consensus/10_consensus_diagnostics/
#       shared/<cpg_label>/<region_label>/<adjustment_version>/
#           Consensus_Eigengene_Networks.pdf
#       <dataset_label>/<cpg_label>/<region_label>/<adjustment_version>/
#           ME_Correlation_Heatmap.pdf
#           Sample_Dendrogram_from_Consensus_MEs.pdf
#           Module_Dendrogram_from_Consensus_MEs.pdf
#           Sample_Correlation_Heatmap.pdf
#           Sample_ME_Heatmap.pdf
#           Module_Module_Correlation_Stats.tsv
#
# NOTES
#   - cpg_label and region_label are derived from --regions_file path
#   - These diagnostics help assess module quality and sample clustering
#     before proceeding to trait association analysis
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/10_diagnostics_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Modules.rds \
#     --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
#     --adjustment_version v1_all_pcs \
#     --consensus_cor bicor
# ================================================================
# SCRIPT 10: Consensus Diagnostics
#
# PURPOSE
#   - Load Consensus_Modules.rds
#   - Save consensus eigengene network plots
#   - Save per-dataset ME correlation heatmaps
#   - Save per-dataset sample dendrograms from consensus MEs
#   - Save per-dataset module dendrograms from consensus MEs
#   - Save per-dataset sample correlation heatmaps
#   - Save per-dataset sample x module eigengene heatmaps
#   - Save module-module correlation stats tables
# ================================================================
message("Starting Script 10 ✓")
suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
  library(comethyl)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_close_pdf <- function() {
  cur_dev <- grDevices::dev.cur()
  if (!is.null(cur_dev) && cur_dev > 1) {
    try(grDevices::dev.off(), silent = TRUE)
  }
}

write_tsv <- function(df, file) {
  utils::write.table(
    df,
    file = file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root",          type = "character"),
  make_option("--consensus_modules_rds", type = "character"),
  make_option("--regions_file",          type = "character"),
  make_option("--adjustment_version",    type = "character", default = "unadjusted"),
  make_option("--consensus_cor",         type = "character", default = "bicor"),
  make_option("--max_p_outliers",        type = "double",    default = 0.1),
  make_option("--helper_functions",      type = "character", default = NULL)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$project_root))
  stop("--project_root is required")
if (is.null(opt$consensus_modules_rds))
  stop("--consensus_modules_rds is required")
if (is.null(opt$regions_file))
  stop("--regions_file is required")
if (!dir.exists(opt$project_root))
  stop("project_root not found: ", opt$project_root)
if (!file.exists(opt$consensus_modules_rds))
  stop("consensus_modules_rds not found: ", opt$consensus_modules_rds)
if (!file.exists(opt$regions_file))
  stop("regions_file not found: ", opt$regions_file)

consensus_cor <- tolower(opt$consensus_cor)
if (!consensus_cor %in% c("bicor", "pearson"))
  stop("--consensus_cor must be 'bicor' or 'pearson'")

# ------------------------------------------------------------
# Load helper functions if provided
# ------------------------------------------------------------
if (!is.null(opt$helper_functions)) {
  if (!file.exists(opt$helper_functions)) {
    stop("--helper_functions file not found: ", opt$helper_functions)
  }
  source(opt$helper_functions)
}

required_helpers <- c("getDendro", "plotDendro", "getCor", "plotHeatmap")
missing_helpers <- required_helpers[!vapply(required_helpers, exists, logical(1), mode = "function")]

if (length(missing_helpers) > 0) {
  stop(
    "Missing required helper functions: ",
    paste(missing_helpers, collapse = ", "),
    ". Supply them via --helper_functions or load the package/script that defines them."
  )
}

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------
consensusMods <- readRDS(opt$consensus_modules_rds)

if (is.null(consensusMods$colors)) {
  stop("Consensus object does not contain $colors")
}
if (is.null(consensusMods$multiMEs)) {
  stop("Consensus object does not contain $multiMEs")
}


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
step_dir      <- file.path(pipeline_root, "10_consensus_diagnostics")
shared_dir    <- file.path(step_dir, "shared", cpg_label, region_label, opt$adjustment_version)

dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", step_dir)

# ------------------------------------------------------------
# Check how many real modules exist (excluding grey)
# ------------------------------------------------------------
all_colors <- consensusMods$colors
module_colors <- unique(all_colors[all_colors != "grey"])
n_modules <- length(module_colors)

message("Non-grey consensus modules detected: ", n_modules)

# ------------------------------------------------------------
# Shared eigengene networks
# ------------------------------------------------------------
if (n_modules >= 2) {
  tryCatch({
    consensusMEs <- consensusOrderMEs(consensusMods$multiMEs)

    grDevices::pdf(
      file.path(shared_dir, "Consensus_Eigengene_Networks.pdf"),
      width = 8,
      height = 7
    )
    par(cex = 0.8)

    plotEigengeneNetworks(
      consensusMEs,
      setLabels        = names(consensusMods$multiMEs),
      plotDendrograms  = FALSE,
      marHeatmap       = c(3, 3, 2, 1),
      zlimPreservation = c(0.5, 1),
      xLabelsAngle     = 90
    )

    safe_close_pdf()
    message("Saved: Consensus_Eigengene_Networks.pdf")
  }, error = function(e) {
    safe_close_pdf()
    message("WARNING: Failed to generate consensus eigengene networks: ", conditionMessage(e))
    writeLines(
      paste("Failed to generate Consensus_Eigengene_Networks.pdf:", conditionMessage(e)),
      con = file.path(shared_dir, "Consensus_Eigengene_Networks_FAILED.txt")
    )
  })
} else {
  message("WARNING: Fewer than 2 real modules found. Skipping eigengene network plot.")
  writeLines(
    paste(
      "Skipped: only", n_modules, "non-grey module(s) detected.",
      "plotEigengeneNetworks requires >= 2 modules."
    ),
    con = file.path(shared_dir, "Consensus_Eigengene_Networks_SKIPPED.txt")
  )
}

# ------------------------------------------------------------
# Per-dataset outputs
# ------------------------------------------------------------
for (ds in names(consensusMods$multiMEs)) {

  message("Processing dataset: ", ds)

  ds_dir <- file.path(step_dir, ds, cpg_label, region_label, opt$adjustment_version)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(consensusMods$multiMEs[[ds]]$data)) {
    message(ds, ": No $data found in multiMEs entry. Skipping.")
    writeLines(
      "Skipped: no eigengene data found in consensusMods$multiMEs[[ds]]$data",
      con = file.path(ds_dir, paste0(ds, "_plots_SKIPPED.txt"))
    )
    next
  }

  MEs_use <- consensusMods$multiMEs[[ds]]$data

  if (!is.data.frame(MEs_use) && !is.matrix(MEs_use)) {
    message(ds, ": MEs object is not a matrix/data.frame. Skipping.")
    writeLines(
      "Skipped: MEs is not a matrix/data.frame",
      con = file.path(ds_dir, paste0(ds, "_plots_SKIPPED.txt"))
    )
    next
  }

  MEs_use <- as.data.frame(MEs_use)

  grey_col <- grep("^MEgrey$", colnames(MEs_use), value = TRUE)
  if (length(grey_col) > 0) {
    MEs_use <- MEs_use[, !colnames(MEs_use) %in% grey_col, drop = FALSE]
    message(ds, ": MEgrey removed before plotting (", ncol(MEs_use), " real MEs remaining)")
  } else {
    message(ds, ": No MEgrey column found (", ncol(MEs_use), " real MEs remaining)")
  }

  if (ncol(MEs_use) < 2) {
    message(ds, ": Only ", ncol(MEs_use), " real ME(s) — skipping correlation/dendrogram plots")
    writeLines(
      paste("Skipped: only", ncol(MEs_use), "real ME(s) after removing grey."),
      con = file.path(ds_dir, paste0(ds, "_plots_SKIPPED.txt"))
    )
    next
  }

  corType_use <- consensus_cor

  f_sample_dendro   <- file.path(ds_dir, paste0(ds, "_Consensus_ME_Sample_Dendrogram.pdf"))
  f_module_dendro   <- file.path(ds_dir, paste0(ds, "_Consensus_ME_Module_Dendrogram.pdf"))
  f_me_cor_hm       <- file.path(ds_dir, paste0(ds, "_Consensus_ME_Correlation_Heatmap.pdf"))
  f_sample_cor_hm   <- file.path(ds_dir, paste0(ds, "_Consensus_Sample_Correlation_Heatmap.pdf"))
  f_sample_me_hm    <- file.path(ds_dir, paste0(ds, "_Consensus_Sample_ME_Heatmap.pdf"))
  f_module_cor_tsv  <- file.path(ds_dir, paste0(ds, "_Consensus_Module_Correlation_Stats.tsv"))

  # ==========================================================
  # 1) Sample dendrogram
  # ==========================================================
  sampleDendro <- tryCatch({
    getDendro(MEs_use, transpose = TRUE, distance = corType_use)
  }, error = function(e) {
    message(ds, ": Sample dendrogram failed — ", conditionMessage(e))
    NULL
  })

  if (!is.null(sampleDendro)) {
    tryCatch({
      plotDendro(
        sampleDendro,
        labelSize = 0.8,
        nBreaks   = 5,
        file      = f_sample_dendro
      )
      message(ds, ": Saved sample dendrogram")
    }, error = function(e) {
      message(ds, ": Failed saving sample dendrogram — ", conditionMessage(e))
    })
  }

  # ==========================================================
  # 2) Module dendrogram
  # ==========================================================
  moduleDendro <- tryCatch({
    getDendro(MEs_use, distance = corType_use)
  }, error = function(e) {
    message(ds, ": Module dendrogram failed — ", conditionMessage(e))
    NULL
  })

  if (!is.null(moduleDendro)) {
    tryCatch({
      plotDendro(
        moduleDendro,
        labelSize = 4,
        nBreaks   = 5,
        file      = f_module_dendro
      )
      message(ds, ": Saved module dendrogram")
    }, error = function(e) {
      message(ds, ": Failed saving module dendrogram — ", conditionMessage(e))
    })
  }

  # ==========================================================
  # 3) Module-module correlation heatmap
  # ==========================================================
  me_cor <- tryCatch({
    getCor(MEs_use, corType = corType_use, maxPOutliers = opt$max_p_outliers)
  }, error = function(e) {
    message(ds, ": ME correlation calculation failed — ", conditionMessage(e))
    NULL
  })

  if (is.null(me_cor)) {
    writeLines(
      "Failed: could not compute ME correlation matrix",
      con = file.path(ds_dir, paste0(ds, "_ME_correlation_FAILED.txt"))
    )
    next
  }

  tryCatch({
    plotHeatmap(
      me_cor,
      rowDendro = moduleDendro,
      colDendro = moduleDendro,
      file      = f_me_cor_hm
    )
    message(ds, ": Saved ME correlation heatmap")
  }, error = function(e) {
    message(ds, ": Failed ME correlation heatmap — ", conditionMessage(e))
  })

  # ==========================================================
  # 4) Sample correlation heatmap
  # ==========================================================
  sampleCor <- tryCatch({
    getCor(MEs_use, transpose = TRUE, corType = corType_use, maxPOutliers = opt$max_p_outliers)
  }, error = function(e) {
    message(ds, ": Sample correlation calculation failed — ", conditionMessage(e))
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
      message(ds, ": Saved sample correlation heatmap")
    }, error = function(e) {
      message(ds, ": Sample correlation heatmap failed — ", conditionMessage(e))
    })
  }

  # ==========================================================
  # 5) Sample x Module ME heatmap
  # ==========================================================
  tryCatch({
    plotHeatmap(
      MEs_use,
      rowDendro       = sampleDendro,
      colDendro       = moduleDendro,
      legend.title    = "Module\nEigengene",
      legend.position = c(0.37, 0.89),
      file            = f_sample_me_hm
    )
    message(ds, ": Saved sample x module eigengene heatmap")
  }, error = function(e) {
    message(ds, ": Sample x module heatmap failed — ", conditionMessage(e))
    writeLines(
      paste("Failed:", conditionMessage(e)),
      con = file.path(ds_dir, paste0(ds, "_Consensus_Sample_ME_Heatmap_FAILED.txt"))
    )
  })

  # ==========================================================
  # 6) Module-module correlation stats table
  # ==========================================================
  tryCatch({
    cor_long <- data.frame(
      module_1 = rep(colnames(me_cor), times = ncol(me_cor)),
      module_2 = rep(colnames(me_cor), each = nrow(me_cor)),
      correlation = as.vector(me_cor),
      stringsAsFactors = FALSE
    )

    cor_long <- cor_long[cor_long$module_1 != cor_long$module_2, , drop = FALSE]

    cor_long$pair_id <- apply(
      cor_long[, c("module_1", "module_2")],
      1,
      function(x) paste(sort(x), collapse = "__")
    )
    cor_long <- cor_long[!duplicated(cor_long$pair_id), , drop = FALSE]
    cor_long$pair_id <- NULL

    cor_long <- cor_long[order(abs(cor_long$correlation), decreasing = TRUE), , drop = FALSE]

    write_tsv(cor_long, f_module_cor_tsv)
    message(ds, ": Saved module correlation stats")
  }, error = function(e) {
    message(ds, ": Module correlation stats failed — ", conditionMessage(e))
    writeLines(
      paste("Failed:", conditionMessage(e)),
      con = file.path(ds_dir, paste0(ds, "_Consensus_Module_Correlation_Stats_FAILED.txt"))
    )
  })
}

# ------------------------------------------------------------
# Logs
# ------------------------------------------------------------
writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("consensus_modules_rds:", opt$consensus_modules_rds),
    paste("adjustment_version:", opt$adjustment_version),
    paste("consensus_cor:", consensus_cor),
    paste("max_p_outliers:", opt$max_p_outliers),
    paste("n_real_modules:", n_modules),
    paste("module_colors:", paste(module_colors, collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(shared_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(shared_dir, "sessionInfo.txt")
)

message("Script 10 complete: consensus diagnostics finished")