#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 16: Consensus Module Enrichment
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Runs pathway and gene-set enrichment analysis for selected
#   consensus modules using gene lists from script 15:
#     1) Load annotation files from script 15
#     2) Select modules via explicit list OR stats-based selection
#        (modules significant for any trait in any dataset)
#     3) Run per-module enrichment: KEGG, GO (BP/MF/CC), Reactome,
#        and/or Enrichr
#     4) Generate per-module result tables (.tsv/.xlsx) and dotplots
#     5) Generate summary dotplots across all selected modules
#
# NOTE
#   Module assignments are shared across datasets — enrichment runs
#   once per adjustment version. When using stats-based selection,
#   modules significant in ANY provided dataset are included.
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --annotation_dir : path to script 15 output directory
#                      (comethyl_output/consensus/15_module_annotation/
#                       <cpg_label>/<region_label>/<adjustment_version>/)
#
# SELECTION MODE (provide one):
#   --module_list_file : .txt file OR directory of .txt files with module
#                        names to enrich (one per line)
#   OR
#   --trait_list_file       : .txt file OR directory of .txt files with
#                             trait names for stats-based selection
#   --dataset1_stats_file   : ME-trait stats TSV from script 12a (dataset 1)
#   --dataset2_stats_file   : ME-trait stats TSV from script 12a (dataset 2)
#   --dataset3_stats_file   : ME-trait stats TSV from script 12a (dataset 3)
#   --dataset1_label        : label for dataset 1 [default = dataset1]
#   --dataset2_label        : label for dataset 2
#   --dataset3_label        : label for dataset 3
#   --p_thresh              : p-value threshold for module selection [default = 0.05]
#
# OPTIONAL INPUTS — ENRICHMENT METHODS
#   --do_kegg          : run KEGG TRUE/FALSE [default = TRUE]
#   --do_go            : run GO BP/MF/CC TRUE/FALSE [default = TRUE]
#   --do_reactome      : run Reactome TRUE/FALSE [default = TRUE]
#   --do_enrichr       : run Enrichr TRUE/FALSE [default = FALSE]
#   --enrichr_dbs      : comma-separated Enrichr DB names
#   --organism         : KEGG organism code [default = hsa]
#
# OPTIONAL INPUTS — ENRICHR SETTINGS
#   --use_enrichr_background : use background gene universe TRUE/FALSE [default = TRUE]
#   --min_bg_genes           : minimum background size [default = 200]
#   --enrichr_retries        : Enrichr retry count [default = 3]
#   --enrichr_sleep          : seconds between retries [default = 1]
#
# OPTIONAL INPUTS — PLOT SETTINGS
#   --plot_p_col      : p.adjust or pvalue for summary plots [default = p.adjust]
#   --plot_top_n      : top N terms per module in summary plots [default = 10]
#   --plot_sig_cutoff : significance cutoff for summary plots [default = 0.05]
#   --plot_width      : summary plot width in inches [default = 14]
#   --plot_height     : summary plot height in inches [default = 10]
#
# OUTPUTS
#   comethyl_output/consensus/16_enrichment/<cpg_label>/<region_label>/<adjustment_version>/[<list_stem>/]
#       selected_modules.txt
#       selected_modules_source.tsv
#       selected_modules_summary.tsv
#       <module>_KEGG.tsv / .xlsx / _dotplot.pdf
#       <module>_GO_BP/MF/CC.tsv / .xlsx / _dotplot.pdf
#       <module>_Reactome.tsv / .xlsx / _dotplot.pdf
#       <module>_enrichr.xlsx
#       Summary_KEGG_top<n>_<pcol><cutoff>.pdf
#       Summary_GO_BP/MF/CC_top<n>_<pcol><cutoff>.pdf
#       Summary_Reactome_top<n>_<pcol><cutoff>.pdf
#       Summary_Enrichr_<db>_top<n>_<pcol><cutoff>.pdf
#       run_log.txt
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - cpg_label and region_label are derived from --annotation_dir path
#   - If --module_list_file or --trait_list_file is a directory, each
#     .txt file is processed separately in its own output subfolder
#     named after the file stem
#   - Enrichr requires internet access; all other methods run offline
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/16_enrichment_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --annotation_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/15_module_annotation/cov3_75pct/covMin4_methSD0p08/v1_all_pcs \
#     --adjustment_version v1_all_pcs \
#     --trait_list_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/config/outcome_traits.txt \
#     --dataset1_label Baseline \
#     --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Baseline/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --dataset2_label 36_38wks \
#     --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/36_38wks/ME_Trait_Correlation_Stats_Bicor.tsv \
#     --p_thresh 0.05 \
#     --do_kegg TRUE \
#     --do_go TRUE \
#     --do_reactome TRUE \
#     --plot_top_n 10
# PURPOSE
#   - Load annotated regions from script 15
#   - Select modules via explicit list OR stats-based (trait + p-threshold)
#     across one or more datasets — union of significant modules
#   - Run KEGG, GO (BP/MF/CC), Reactome, and/or Enrichr per module
#   - Generate per-module result tables and dotplots
#   - Generate summary dotplots across all selected modules
#   - Write run logs
#
# NOTE
#   Module assignments are shared across datasets — enrichment runs
#   once per adjustment version. When using stats-based selection,
#   modules significant in ANY dataset are included.
#
# SELECTION MODES
#   Provide either:
#     --module_list_file  : explicit list of module names (file OR directory of .txt files)
#   OR:
#     --trait_list_file + --dataset1_stats_file (+ dataset2/3)
#       modules where trait is in list AND p <= p_thresh in ANY dataset
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/16_enrichment/
#     <cpg_label>/
#       <region_label>/
#         <adjustment_version>/
#           [<list_stem>/]
#             selected_modules.txt
#             selected_modules_source.tsv
#             selected_modules_summary.tsv
#             <module>_KEGG.tsv / .xlsx / _dotplot.pdf
#             <module>_GO_BP/MF/CC.tsv / .xlsx / _dotplot.pdf
#             <module>_Reactome.tsv / .xlsx / _dotplot.pdf
#             <module>_enrichr.xlsx
#             Summary_KEGG_top<n>_<pcol><cutoff>.pdf
#             Summary_GO_BP/MF/CC_top<n>_<pcol><cutoff>.pdf
#             Summary_Reactome_top<n>_<pcol><cutoff>.pdf
#             Summary_Enrichr_<db>_top<n>_<pcol><cutoff>.pdf
#             run_log.txt
#             run_parameters.txt
#             sessionInfo.txt
#
# PATH DERIVATION
#   cpg_label and region_label derived from --annotation_dir path:
#     .../15_module_annotation/<cpg_label>/<region_label>/<adjustment_version>/
# ================================================================

message("Starting Script 16 ")

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(scales)
})

# ================================================================
# Helpers
# ================================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x); if (!nzchar(x)) return(NULL); x
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("true","t","1","yes","y")
}

stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument: ", label, call. = FALSE)
}

validate_file_exists <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

validate_dir_exists <- function(path, label) {
  if (!dir.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  message(txt)
  cat(txt, "\n", file = logfile, append = TRUE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

qc_gene_list <- function(x, label = "genes") {
  x0 <- x
  x  <- trimws(as.character(x))
  x  <- x[!is.na(x) & x != ""]
  x  <- unique(x)
  ok <- grepl("^[A-Za-z0-9._-]+$", x)
  message(sprintf("[%s] input n=%d -> cleaned unique n=%d", label, length(x0), length(x)))
  message(sprintf("[%s] invalid-format count=%d", label, sum(!ok)))
  list(clean = x, ok = ok)
}

# ================================================================
# File-or-directory resolver (same pattern as regular 12b)
# ================================================================
resolve_list_input <- function(path) {
  if (is.null(path)) return(NULL)
  if (file.exists(path) && !dir.exists(path))
    return(setNames(normalizePath(path), ""))
  if (dir.exists(path)) {
    txts  <- sort(list.files(path, pattern = "\\.txt$", full.names = TRUE))
    if (length(txts) == 0)
      stop("Directory has no .txt files: ", path, call. = FALSE)
    stems <- tools::file_path_sans_ext(basename(txts))
    return(setNames(normalizePath(txts), stems))
  }
  stop("Path does not exist: ", path, call. = FALSE)
}

read_list_file <- function(path) {
  x <- readLines(path, warn = FALSE)
  unique(trimws(x[nzchar(trimws(x))]))
}

# ================================================================
# Annotation loader
# ================================================================
load_annotation_files <- function(annotation_dir) {
  af  <- file.path(annotation_dir, "Annotated_Regions.tsv")
  glf <- file.path(annotation_dir, "Module_Gene_List.tsv")
  smf <- file.path(annotation_dir, "Module_Gene_Summary.tsv")
  bgf <- file.path(annotation_dir, "Background_Genes.tsv")

  validate_file_exists(af,  "Annotated_Regions.tsv")
  validate_file_exists(glf, "Module_Gene_List.tsv")
  validate_file_exists(smf, "Module_Gene_Summary.tsv")

  annotated  <- read.delim(af,  stringsAsFactors = FALSE, check.names = FALSE)
  gene_list  <- read.delim(glf, stringsAsFactors = FALSE, check.names = FALSE)
  summary    <- read.delim(smf, stringsAsFactors = FALSE, check.names = FALSE)

  bg <- NULL
  if (file.exists(bgf)) {
    bg_df <- read.delim(bgf, stringsAsFactors = FALSE, check.names = FALSE)
    if ("gene_symbol" %in% names(bg_df)) bg <- bg_df$gene_symbol
  }

  list(annotated = annotated, gene_list = gene_list,
       summary = summary, background_genes = bg)
}

# ================================================================
# Stats-based module selection
# ================================================================
standardize_stats_columns <- function(df, label) {
  cn <- tolower(colnames(df))
  for (col in c("module","trait")) {
    if (!col %in% cn) stop(label, " must contain '", col, "' column.", call. = FALSE)
    names(df)[match(col, cn)] <- col; cn <- tolower(colnames(df))
  }
  p_cands <- c("p","pvalue","p.value","p_value")
  p_hit   <- p_cands[p_cands %in% cn][1]
  if (is.na(p_hit))
    stop(label, " must have a p-value column.", call. = FALSE)
  names(df)[match(p_hit, tolower(colnames(df)))] <- "p"
  df
}

read_stats_file <- function(path, label) {
  ext <- tolower(tools::file_ext(path))
  df  <- switch(ext,
    xlsx = openxlsx::read.xlsx(path),
    xls  = openxlsx::read.xlsx(path),
    tsv  = , txt = read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    csv  = read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    stop("Unsupported stats file format: ", path, call. = FALSE)
  )
  standardize_stats_columns(df, label)
}

select_modules_from_stats <- function(stats_df, trait_list, p_thresh) {
  stats_df %>%
    dplyr::mutate(p = as.numeric(p)) %>%
    dplyr::filter(!is.na(p), trait %in% trait_list, p <= p_thresh) %>%
    dplyr::pull(module) %>% unique() %>% as.character()
}

# ================================================================
# Gene list builder
# ================================================================
build_gene_vector <- function(gene_list_df, module_name) {
  gene_list_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene_symbol) %>%
    as.character() %>% trimws() %>% unique() %>%
    .[!is.na(.) & . != ""]
}

# ================================================================
# KEGG enrichment
# ================================================================
run_kegg_clusterprofiler <- function(gene_symbols, out_prefix, organism = "hsa") {
  gene_symbols <- qc_gene_list(gene_symbols, "KEGG")
  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    message("  [KEGG] No Entrez IDs mapped — skipping")
    return(invisible(NULL))
  }

  kegg_res <- suppressMessages(
    clusterProfiler::enrichKEGG(gene = unique(gene_df$ENTREZID),
                                  organism = organism,
                                  pvalueCutoff = 1, qvalueCutoff = 1)
  )
  df <- as.data.frame(kegg_res)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  write.table(df, paste0(out_prefix, "_KEGG.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(df, paste0(out_prefix, "_KEGG.xlsx"), rowNames = FALSE)

  top <- df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(25, nrow(df)))
  if (nrow(top) > 0) {
    top$Description <- factor(top$Description, levels = rev(top$Description))
    p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
      geom_point() + theme_bw(base_size = 12) +
      labs(title = "Top KEGG pathways", x = "-log10(adj. p)", y = NULL)
    ggsave(paste0(out_prefix, "_KEGG_dotplot.pdf"), plot = p, width = 9, height = 6)
  }
  message("  [KEGG] Done: ", nrow(df), " pathways")
  invisible(df)
}

# ================================================================
# GO enrichment (BP, MF, CC)
# ================================================================
run_go_clusterprofiler <- function(gene_symbols, out_prefix) {
  gene_symbols <- qc_gene_list(gene_symbols, "GO")
  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    message("  [GO] No Entrez IDs mapped — skipping")
    return(invisible(NULL))
  }
  entrez_ids <- unique(gene_df$ENTREZID)

  for (ont in c("BP","MF","CC")) {
    res <- tryCatch(
      suppressMessages(
        clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
                                   ont = ont, pvalueCutoff = 1, qvalueCutoff = 1,
                                   readable = TRUE)
      ),
      error = function(e) { message("  [GO ", ont, "] Failed: ", conditionMessage(e)); NULL }
    )
    if (is.null(res)) next
    df <- as.data.frame(res)
    if (nrow(df) == 0) next

    write.table(df, paste0(out_prefix, "_GO_", ont, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    openxlsx::write.xlsx(df, paste0(out_prefix, "_GO_", ont, ".xlsx"), rowNames = FALSE)

    top <- df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(25, nrow(df)))
    if (nrow(top) > 0) {
      top$Description <- factor(top$Description, levels = rev(top$Description))
      p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
        geom_point() + theme_bw(base_size = 12) +
        labs(title = paste0("Top GO ", ont, " terms"), x = "-log10(adj. p)", y = NULL)
      ggsave(paste0(out_prefix, "_GO_", ont, "_dotplot.pdf"), plot = p, width = 9, height = 6)
    }
    message("  [GO ", ont, "] Done: ", nrow(df), " terms")
  }
  invisible(NULL)
}

# ================================================================
# Reactome enrichment
# ================================================================
run_reactome_clusterprofiler <- function(gene_symbols, out_prefix) {
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    message("  [Reactome] ReactomePA not installed — skipping")
    return(invisible(NULL))
  }
  gene_symbols <- qc_gene_list(gene_symbols, "Reactome")
  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    message("  [Reactome] No Entrez IDs mapped — skipping")
    return(invisible(NULL))
  }

  res <- tryCatch(
    suppressMessages(
      ReactomePA::enrichPathway(gene = unique(gene_df$ENTREZID),
                                 organism = "human", pvalueCutoff = 1,
                                 qvalueCutoff = 1, readable = TRUE)
    ),
    error = function(e) { message("  [Reactome] Failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(res)) return(invisible(NULL))
  df <- as.data.frame(res)
  if (nrow(df) == 0) return(invisible(NULL))

  write.table(df, paste0(out_prefix, "_Reactome.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(df, paste0(out_prefix, "_Reactome.xlsx"), rowNames = FALSE)

  top <- df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = min(25, nrow(df)))
  if (nrow(top) > 0) {
    top$Description <- factor(top$Description, levels = rev(top$Description))
    p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
      geom_point() + theme_bw(base_size = 12) +
      labs(title = "Top Reactome Pathways", x = "-log10(adj. p)", y = NULL)
    ggsave(paste0(out_prefix, "_Reactome_dotplot.pdf"), plot = p, width = 9, height = 6)
  }
  message("  [Reactome] Done: ", nrow(df), " pathways")
  invisible(df)
}

# ================================================================
# Enrichr
# ================================================================
init_enrichr <- function() {
  options(
    enrichR.sites.base.address = "https://maayanlab.cloud/",
    enrichR.base.address       = "https://maayanlab.cloud/Enrichr/",
    speedrichr.base.address    = "https://maayanlab.cloud/speedrichr/api/",
    enrichR.live               = TRUE,
    enrichR.quiet              = TRUE,
    modEnrichR.use             = TRUE,
    enrichR.sites              = c("Enrichr", "FlyEnrichr", "WormEnrichr",
                                   "YeastEnrichr", "FishEnrichr", "OxEnrichr")
  )
}

validate_enrichr_dbs_safe <- function(dbs) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning("enrichR not installed; using requested list as-is.")
    return(list(ok = dbs, bad = character(0)))
  }
  tryCatch({
    init_enrichr(); Sys.sleep(0.5)
    avail <- as.character(enrichR::listEnrichrDbs()$libraryName)
    list(ok = intersect(dbs, avail), bad = setdiff(dbs, avail))
  }, error = function(e) {
    warning("Could not validate Enrichr DBs: ", conditionMessage(e))
    list(ok = dbs, bad = character(0))
  })
}
run_enrichr_simple <- function(gene_symbols, out_prefix, dbs,
                               background_genes = NULL,
                               retries          = 3,
                               sleep_time       = 1) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning("enrichR package not installed; skipping Enrichr.")
    return(invisible(NULL))
  }

  fg_qc        <- qc_gene_list(gene_symbols, "foreground")
  gene_symbols <- fg_qc$clean
  if (length(gene_symbols) == 0) return(invisible(NULL))
  if (length(dbs) == 0)          return(invisible(NULL))

  if (!is.null(background_genes)) {
    bg_qc            <- qc_gene_list(background_genes, "background")
    background_genes <- bg_qc$clean
    if (length(background_genes) == 0) background_genes <- NULL
  }

  init_enrichr()

  res <- NULL; last_err <- NULL

  for (attempt in seq_len(retries)) {
    message("[Enrichr] Attempt ", attempt, " of ", retries,
            " for: ", basename(out_prefix))

    res <- tryCatch({
      if (!is.null(background_genes) && length(background_genes) > 0) {
        enrichR::enrichr(genes = gene_symbols, databases = dbs,
                         background = background_genes)
      } else {
        enrichR::enrichr(genes = gene_symbols, databases = dbs)
      }
    }, error = function(e) { last_err <<- conditionMessage(e); NULL })

    if (!is.null(res) && length(res) > 0) break
    if (attempt < retries) {
      message("[Enrichr] Failed attempt ", attempt, ": ", last_err,
              " — retrying in ", sleep_time, "s")
      Sys.sleep(sleep_time)
    }
  }

  if (is.null(res) || length(res) == 0) {
    warning("Enrichr failed for ", out_prefix,
            if (!is.null(last_err)) paste0(": ", last_err) else ".")
    return(invisible(NULL))
  }

  # ---- Save combined Excel workbook ----
  wb <- openxlsx::createWorkbook(); wrote_any <- FALSE

  for (nm in names(res)) {
    df <- res[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    sheet <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df)
    wrote_any <- TRUE
  }

  if (wrote_any) {
    openxlsx::saveWorkbook(wb, paste0(out_prefix, "_enrichr.xlsx"), overwrite = TRUE)
    message("[Enrichr] Saved: ", basename(out_prefix), "_enrichr.xlsx")

    # ---- Per-database plotEnrich() PDFs ----
    # Requires enrichR::plotEnrich(); mirrors the old plot_and_save() logic.
    for (nm in names(res)) {
      df <- res[[nm]]
      if (is.null(df) || nrow(df) == 0) next

      db_safe  <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 60)
      pdf_file <- paste0(out_prefix, "_enrichr_", db_safe, ".pdf")
      mod_name <- basename(out_prefix)   # e.g. "skyblue2"

      tryCatch({
        grDevices::pdf(pdf_file, height = 7, width = 15)
        print(
          enrichR::plotEnrich(df, showTerms = 25, numChar = 75,
                              y = "Count", orderBy = "P.value") +
            ggplot2::ggtitle(paste(nm, "for", mod_name, "module"))
        )
        grDevices::dev.off()
        message("[Enrichr] Plot saved: ", basename(pdf_file))
      }, error = function(e) {
        # Ensure graphics device is always closed on error
        tryCatch(grDevices::dev.off(), error = function(e2) NULL)
        message("[Enrichr] Plot failed for DB '", nm, "': ", conditionMessage(e))
      })
    }
  }

  Sys.sleep(sleep_time)
  invisible(res)
}

# ================================================================
# Summary dotplot helpers (matching regular 12b)
# ================================================================
load_method_results <- function(out_dir, modules, file_suffix,
                                 p_col = "p.adjust",
                                 count_col = "Count",
                                 desc_col = "Description") {
  all_data <- data.frame()
  for (m in modules) {
    f  <- file.path(out_dir, paste0(m, file_suffix))
    if (!file.exists(f)) next
    df <- tryCatch(openxlsx::read.xlsx(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    actual_p <- if (p_col %in% colnames(df)) p_col else
                if ("p.adjust" %in% colnames(df)) "p.adjust" else
                if ("pvalue"   %in% colnames(df)) "pvalue" else { next }
    needed <- c(desc_col, actual_p, count_col)
    if (!all(needed %in% colnames(df))) next
    df <- df %>%
      dplyr::select(Description = !!desc_col,
                    p_col_val   = !!actual_p,
                    gene_count  = !!count_col) %>%
      dplyr::mutate(module = m)
    all_data <- dplyr::bind_rows(all_data, df)
  }
  all_data
}

make_enrichment_dotplot <- function(df, title,
                                     source_tag  = NULL,
                                     p_col_label = "adj. p",
                                     sig_cutoff  = 0.05,
                                     top_n       = 10,
                                     max_chars   = 55,
                                     plot_width  = 14,
                                     plot_height = 10,
                                     base_size   = 13,
                                     out_file    = NULL) {
  if (nrow(df) == 0) { message("[plot] No data for: ", title); return(invisible(NULL)) }

  df$Description <- ifelse(nchar(df$Description) > max_chars,
                            paste0(substr(df$Description, 1, max_chars - 3), "..."),
                            df$Description)

  df <- df %>%
    dplyr::filter(!is.na(p_col_val), p_col_val <= sig_cutoff) %>%
    dplyr::group_by(module) %>%
    dplyr::arrange(p_col_val, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  if (nrow(df) == 0) {
    message("[plot] No terms pass sig_cutoff=", sig_cutoff, " for: ", title)
    return(invisible(NULL))
  }

  df <- df %>%
    dplyr::arrange(module, p_col_val) %>%
    dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

  subtitle_text <- paste0("Top ", top_n, " terms per module | ",
                           p_col_label, " \u2264 ", sig_cutoff,
                           if (!is.null(source_tag)) paste0(" | ", source_tag) else "")

  p <- ggplot(df, aes(x = module, y = Description,
                       size = gene_count, fill = p_col_val)) +
    geom_point(shape = 23, color = "grey30") +
    scale_size(range = c(3,10), name = "Gene Count") +
    scale_fill_gradient(low = "#B22222", high = "#FFC0CB", name = p_col_label,
                        limits = c(0, sig_cutoff), oob = scales::squish) +
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(size = base_size - 1),
          axis.title  = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "right") +
    labs(title = title, subtitle = subtitle_text)

  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = plot_width, height = plot_height,
           dpi = 300, units = "in")
    message("[plot] Saved: ", basename(out_file))
  }
  invisible(p)
}

plot_kegg_summary <- function(out_dir, modules, p_col = "p.adjust",
                               sig_cutoff = 0.05, top_n = 10,
                               plot_width = 14, plot_height = 10) {
  p_col_label <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe  <- gsub("\\.", "", p_col)
  df <- load_method_results(out_dir, modules, "_KEGG.xlsx", p_col = p_col)
  if (nrow(df) == 0) return(invisible(NULL))
  make_enrichment_dotplot(df, title = "Top KEGG Pathways per Module",
    source_tag = "clusterProfiler", p_col_label = p_col_label,
    sig_cutoff = sig_cutoff, top_n = top_n,
    plot_width = plot_width, plot_height = plot_height,
    out_file = file.path(out_dir,
      sprintf("Summary_KEGG_top%d_%s%s.pdf", top_n, p_col_safe, sig_cutoff)))
}

plot_go_summary <- function(out_dir, modules, p_col = "p.adjust",
                             sig_cutoff = 0.05, top_n = 10,
                             plot_width = 14, plot_height = 10) {
  p_col_label <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe  <- gsub("\\.", "", p_col)
  for (ont in c("BP","MF","CC")) {
    df <- load_method_results(out_dir, modules,
                               paste0("_GO_", ont, ".xlsx"), p_col = p_col)
    if (nrow(df) == 0) next
    make_enrichment_dotplot(df, title = paste0("Top GO ", ont, " Terms per Module"),
      source_tag = "clusterProfiler", p_col_label = p_col_label,
      sig_cutoff = sig_cutoff, top_n = top_n,
      plot_width = plot_width, plot_height = plot_height,
      out_file = file.path(out_dir,
        sprintf("Summary_GO_%s_top%d_%s%s.pdf", ont, top_n, p_col_safe, sig_cutoff)))
  }
}

plot_reactome_summary <- function(out_dir, modules, p_col = "p.adjust",
                                   sig_cutoff = 0.05, top_n = 10,
                                   plot_width = 14, plot_height = 10) {
  p_col_label <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe  <- gsub("\\.", "", p_col)
  df <- load_method_results(out_dir, modules, "_Reactome.xlsx", p_col = p_col)
  if (nrow(df) == 0) return(invisible(NULL))
  make_enrichment_dotplot(df, title = "Top Reactome Pathways per Module",
    source_tag = "ReactomePA", p_col_label = p_col_label,
    sig_cutoff = sig_cutoff, top_n = top_n,
    plot_width = plot_width, plot_height = plot_height,
    out_file = file.path(out_dir,
      sprintf("Summary_Reactome_top%d_%s%s.pdf", top_n, p_col_safe, sig_cutoff)))
}

plot_enrichr_summary <- function(out_dir, modules, databases,
                                  p_col = "p.adjust", sig_cutoff = 0.05,
                                  top_n = 10, plot_width = 14, plot_height = 10) {
  enrichr_p_col <- if (p_col == "p.adjust") "Adjusted.P.value" else "P.value"
  p_col_label   <- if (p_col == "p.adjust") "adj. p" else "raw p"
  p_col_safe    <- gsub("\\.", "", p_col)

  for (db in databases) {
    all_data <- data.frame()
    for (m in modules) {
      f <- file.path(out_dir, paste0(m, "_enrichr.xlsx"))
      if (!file.exists(f)) next
      sheet <- substr(gsub("[^A-Za-z0-9]","_", db), 1, 31)
      df <- tryCatch(openxlsx::read.xlsx(f, sheet = sheet), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next
      actual_p <- if (enrichr_p_col %in% colnames(df)) enrichr_p_col else
                  if ("Adjusted.P.value" %in% colnames(df)) "Adjusted.P.value" else "P.value"
      if (!all(c("Term", actual_p, "Overlap") %in% colnames(df))) next
      df$gene_count  <- suppressWarnings(as.integer(sub("/.*$","", as.character(df$Overlap))))
      df$Description <- gsub(" \\(GO:\\d+\\)$","", df$Term)
      df <- df %>%
        dplyr::select(Description, p_col_val = !!actual_p, gene_count) %>%
        dplyr::mutate(module = m)
      all_data <- dplyr::bind_rows(all_data, df)
    }
    if (nrow(all_data) == 0) next
    db_safe  <- gsub("[^A-Za-z0-9]+","_", db)
    make_enrichment_dotplot(all_data,
      title      = paste0("Top Enrichr Terms: ", db),
      source_tag = "Enrichr / maayanlab.cloud",
      p_col_label = p_col_label,
      sig_cutoff = sig_cutoff, top_n = top_n,
      plot_width = plot_width, plot_height = plot_height,
      out_file = file.path(out_dir,
        sprintf("Summary_Enrichr_%s_top%d_%s%s.pdf", db_safe, top_n, p_col_safe, sig_cutoff)))
  }
}

# ================================================================
# Read arguments
# ================================================================
project_root      <- trim_or_null(get_arg("--project_root"))
annotation_dir    <- trim_or_null(get_arg("--annotation_dir"))
adjustment_version <- trim_or_null(get_arg("--adjustment_version", "unadjusted"))

module_list_input <- trim_or_null(get_arg("--module_list_file"))
trait_list_input  <- trim_or_null(get_arg("--trait_list_file"))

dataset1_label      <- trim_or_null(get_arg("--dataset1_label", "dataset1"))
dataset1_stats_file <- trim_or_null(get_arg("--dataset1_stats_file"))
dataset2_label      <- trim_or_null(get_arg("--dataset2_label"))
dataset2_stats_file <- trim_or_null(get_arg("--dataset2_stats_file"))
dataset3_label      <- trim_or_null(get_arg("--dataset3_label"))
dataset3_stats_file <- trim_or_null(get_arg("--dataset3_stats_file"))

p_thresh    <- as.numeric(get_arg("--p_thresh", "0.05"))
do_kegg     <- to_bool(get_arg("--do_kegg",     "true"),  default = TRUE)
do_go       <- to_bool(get_arg("--do_go",       "true"),  default = TRUE)
do_reactome <- to_bool(get_arg("--do_reactome", "true"),  default = TRUE)
do_enrichr  <- to_bool(get_arg("--do_enrichr",  "false"), default = FALSE)

enrichr_dbs <- split_csv(get_arg("--enrichr_dbs",
  "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025"))

organism               <- trim_or_null(get_arg("--organism", "hsa"))
use_enrichr_background <- to_bool(get_arg("--use_enrichr_background","true"), TRUE)
min_bg_genes           <- as.integer(get_arg("--min_bg_genes", "200"))
enrichr_retries        <- as.integer(get_arg("--enrichr_retries", "3"))
enrichr_sleep          <- as.numeric(get_arg("--enrichr_sleep", "1"))

plot_p_col      <- trim_or_null(get_arg("--plot_p_col",      "p.adjust"))
plot_top_n      <- as.integer(get_arg( "--plot_top_n",       "10"))
plot_sig_cutoff <- as.numeric(get_arg( "--plot_sig_cutoff",  "0.05"))
plot_width      <- as.numeric(get_arg( "--plot_width",       "14"))
plot_height     <- as.numeric(get_arg( "--plot_height",      "10"))

# ================================================================
# Validate
# ================================================================
stop_if_missing(project_root,   "--project_root")
stop_if_missing(annotation_dir, "--annotation_dir")

validate_dir_exists(project_root,   "project_root")
validate_dir_exists(annotation_dir, "annotation_dir")

if (!plot_p_col %in% c("p.adjust","pvalue"))
  stop("--plot_p_col must be 'p.adjust' or 'pvalue'", call. = FALSE)

# Resolve list inputs
module_list_files <- resolve_list_input(module_list_input)
trait_list_files  <- resolve_list_input(trait_list_input)

if (!is.null(module_list_files)) {
  if (!is.null(trait_list_files))
    message("Note: --module_list_file provided; --trait_list_file will be ignored.")
} else {
  if (is.null(trait_list_files))
    stop("Provide --module_list_file OR --trait_list_file with stats_file(s).", call. = FALSE)
  if (is.null(dataset1_stats_file))
    stop("--dataset1_stats_file required when using --trait_list_file", call. = FALSE)
  validate_file_exists(dataset1_stats_file, "--dataset1_stats_file")
  if (!is.null(dataset2_stats_file)) validate_file_exists(dataset2_stats_file, "--dataset2_stats_file")
  if (!is.null(dataset3_stats_file)) validate_file_exists(dataset3_stats_file, "--dataset3_stats_file")
}

if (do_enrichr) {
  db_check    <- validate_enrichr_dbs_safe(enrichr_dbs)
  if (length(db_check$bad) > 0)
    warning("Dropping invalid Enrichr DBs: ", paste(db_check$bad, collapse = ", "))
  enrichr_dbs <- db_check$ok
  if (length(enrichr_dbs) == 0)
    stop("No valid Enrichr DBs remain.", call. = FALSE)
}

setwd(project_root)

# ================================================================
# Derive cpg_label / region_label from annotation_dir path:
#   .../15_module_annotation/<cpg_label>/<region_label>/<adjustment_version>/
# ================================================================
{
  adj_dir      <- annotation_dir
  region_label <- basename(dirname(adj_dir))
  cpg_label    <- basename(dirname(dirname(adj_dir)))
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}

# ================================================================
# Load annotation files
# ================================================================
anno_files       <- load_annotation_files(annotation_dir)
gene_list_df     <- anno_files$gene_list
module_summary   <- anno_files$summary
background_genes <- anno_files$background_genes
available_modules <- unique(as.character(gene_list_df$module))

# ================================================================
# Build run configs (stem × list approach matching regular 12b)
# ================================================================
run_configs <- list()

if (!is.null(module_list_files)) {
  for (i in seq_along(module_list_files)) {
    run_configs[[length(run_configs)+1]] <- list(
      stem             = names(module_list_files)[i],
      module_list_file = unname(module_list_files[i]),
      trait_list       = NULL
    )
  }
} else {
  for (i in seq_along(trait_list_files)) {
    run_configs[[length(run_configs)+1]] <- list(
      stem             = names(trait_list_files)[i],
      module_list_file = NULL,
      trait_list       = read_list_file(unname(trait_list_files[i]))
    )
  }
}

# ================================================================
# Enrichr background prep (shared across all configs)
# ================================================================
bg_to_use <- NULL; bg_ok <- FALSE

if (do_enrichr && isTRUE(use_enrichr_background) && !is.null(background_genes)) {
  bg_clean  <- qc_gene_list(background_genes, "background")
  mapped    <- AnnotationDbi::select(org.Hs.eg.db, keys = bg_clean,
                                      keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID"))
  mapped_ok <- unique(mapped$SYMBOL[!is.na(mapped$ENTREZID)])
  if (length(bg_clean) >= min_bg_genes) {
    bg_to_use <- bg_clean; bg_ok <- TRUE
    message("Enrichr background: ", length(bg_clean), " genes (",
            round(100*length(mapped_ok)/length(bg_clean), 1), "% map to Entrez)")
  } else {
    message("Background too small (", length(bg_clean), " < ", min_bg_genes,
            ") — disabling Enrichr background mode.")
  }
}

# ================================================================
# Main loop: for each config
# ================================================================
for (cfg in run_configs) {

  stem_label <- if (nzchar(cfg$stem)) paste0(" [list: ", cfg$stem, "]") else ""
  message("\n======================================================")
  message("Running config", stem_label)
  message("======================================================")

  # ---- Output dir ----
  pipeline_root <- file.path(project_root, "comethyl_output", "consensus")
  base_out_dir  <- file.path(pipeline_root, "16_enrichment",
                              cpg_label, region_label, adjustment_version)
  out_dir <- if (nzchar(cfg$stem))
    file.path(base_out_dir, cfg$stem) else base_out_dir
  safe_dir_create(out_dir)

  log_file    <- file.path(out_dir, "run_log.txt")
  params_file <- file.path(out_dir, "run_parameters.txt")

  append_log(log_file, "Starting enrichment", stem_label)
  append_log(log_file, "annotation_dir: ", annotation_dir)
  append_log(log_file, "out_dir: ", out_dir)

  # ---- Module selection ----
  selected_modules <- character(0)
  selection_mode   <- NULL
  module_source_rows <- list()

  if (!is.null(cfg$module_list_file)) {
    selected_modules <- read_list_file(cfg$module_list_file)
    selection_mode   <- "module_list_file"
    append_log(log_file, "Selection: explicit list (", length(selected_modules), " modules)")

  } else {
    selection_mode <- "stats_based_union_across_datasets"
    trait_list     <- cfg$trait_list
    append_log(log_file, "Trait list: ", paste(trait_list, collapse = ", "))

    dataset_stats <- list()
    if (!is.null(dataset1_stats_file)) dataset_stats[[dataset1_label]] <- dataset1_stats_file
    if (!is.null(dataset2_stats_file)) dataset_stats[[dataset2_label]] <- dataset2_stats_file
    if (!is.null(dataset3_stats_file)) dataset_stats[[dataset3_label]] <- dataset3_stats_file

    for (ds_label in names(dataset_stats)) {
      sf       <- dataset_stats[[ds_label]]
      stats_df <- read_stats_file(sf, paste0(ds_label, " stats_file"))
      mods     <- select_modules_from_stats(stats_df, trait_list = trait_list, p_thresh = p_thresh)
      append_log(log_file, ds_label, ": ", length(mods), " modules selected (p <= ", p_thresh, ")")
      for (m in mods)
        module_source_rows[[length(module_source_rows)+1]] <-
          data.frame(module = m, dataset = ds_label, stringsAsFactors = FALSE)
      selected_modules <- union(selected_modules, mods)
    }
  }

  # Filter to annotated, non-grey
  selected_modules <- unique(as.character(selected_modules[nzchar(selected_modules)]))
  not_found        <- setdiff(selected_modules, available_modules)
  selected_modules <- intersect(selected_modules, available_modules)
  selected_modules <- selected_modules[!grepl("^grey$", selected_modules, ignore.case = TRUE)]

  if (length(not_found) > 0)
    append_log(log_file, "WARNING: modules not in annotation (skipped): ",
               paste(not_found, collapse = ", "))
  append_log(log_file, "Final selected modules (", length(selected_modules), "): ",
             paste(sort(selected_modules), collapse = ", "))

  # Save selection outputs
  write_lines_safe(selected_modules, file.path(out_dir, "selected_modules.txt"))

  if (length(module_source_rows) > 0) {
    source_df <- dplyr::bind_rows(module_source_rows) %>%
      dplyr::filter(module %in% selected_modules) %>%
      dplyr::arrange(module, dataset)
    write.table(source_df, file.path(out_dir, "selected_modules_source.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  sel_summary <- module_summary %>% dplyr::filter(module %in% selected_modules)
  write.table(sel_summary, file.path(out_dir, "selected_modules_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (length(selected_modules) == 0) {
    append_log(log_file, "No modules selected — skipping enrichment.")
    next
  }

  # ---- Per-module enrichment ----
  for (m in selected_modules) {
    gene_vec   <- build_gene_vector(gene_list_df, m)
    out_prefix <- file.path(out_dir, m)

    append_log(log_file, "Module ", m, ": ", length(gene_vec), " genes")
    message("\n--- Module: ", m, " (", length(gene_vec), " genes) ---")

    if (do_kegg)
      tryCatch(run_kegg_clusterprofiler(gene_vec, out_prefix, organism),
               error = function(e) append_log(log_file, "KEGG failed [", m, "]: ", conditionMessage(e)))

    if (do_go)
      tryCatch(run_go_clusterprofiler(gene_vec, out_prefix),
               error = function(e) append_log(log_file, "GO failed [", m, "]: ", conditionMessage(e)))

    if (do_reactome)
      tryCatch(run_reactome_clusterprofiler(gene_vec, out_prefix),
               error = function(e) append_log(log_file, "Reactome failed [", m, "]: ", conditionMessage(e)))

    if (do_enrichr)
      tryCatch(run_enrichr_simple(gene_vec, out_prefix, enrichr_dbs,
                                   background_genes = if (bg_ok) bg_to_use else NULL,
                                   retries = enrichr_retries, sleep_time = enrichr_sleep),
               error = function(e) append_log(log_file, "Enrichr failed [", m, "]: ", conditionMessage(e)))
  }

  # ---- Summary dotplots ----
  message("\n--- Summary plots ---")
  if (do_kegg)     plot_kegg_summary(out_dir, selected_modules, p_col = plot_p_col,
                                      sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                                      plot_width = plot_width, plot_height = plot_height)
  if (do_go)       plot_go_summary(out_dir, selected_modules, p_col = plot_p_col,
                                    sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                                    plot_width = plot_width, plot_height = plot_height)
  if (do_reactome) plot_reactome_summary(out_dir, selected_modules, p_col = plot_p_col,
                                          sig_cutoff = plot_sig_cutoff, top_n = plot_top_n,
                                          plot_width = plot_width, plot_height = plot_height)
  if (do_enrichr && length(enrichr_dbs) > 0)
                   plot_enrichr_summary(out_dir, selected_modules, enrichr_dbs,
                                         p_col = plot_p_col, sig_cutoff = plot_sig_cutoff,
                                         top_n = plot_top_n, plot_width = plot_width,
                                         plot_height = plot_height)

  # ---- Run parameters ----
  write_lines_safe(c(
    paste0("timestamp\t",              timestamp_now()),
    paste0("project_root\t",           project_root),
    paste0("annotation_dir\t",         annotation_dir),
    paste0("adjustment_version\t",     adjustment_version),
    paste0("cpg_label\t",              cpg_label),
    paste0("region_label\t",           region_label),
    paste0("list_stem\t",              cfg$stem),
    paste0("selection_mode\t",         selection_mode),
    paste0("module_list_file\t",       ifelse(is.null(cfg$module_list_file),"", cfg$module_list_file)),
    paste0("trait_list_file\t",        ifelse(is.null(trait_list_input),"", trait_list_input)),
    paste0("dataset1_label\t",         ifelse(is.null(dataset1_label),"", dataset1_label)),
    paste0("dataset1_stats_file\t",    ifelse(is.null(dataset1_stats_file),"", dataset1_stats_file)),
    paste0("dataset2_label\t",         ifelse(is.null(dataset2_label),"", dataset2_label)),
    paste0("dataset2_stats_file\t",    ifelse(is.null(dataset2_stats_file),"", dataset2_stats_file)),
    paste0("dataset3_label\t",         ifelse(is.null(dataset3_label),"", dataset3_label)),
    paste0("dataset3_stats_file\t",    ifelse(is.null(dataset3_stats_file),"", dataset3_stats_file)),
    paste0("p_thresh\t",               p_thresh),
    paste0("do_kegg\t",                do_kegg),
    paste0("do_go\t",                  do_go),
    paste0("do_reactome\t",            do_reactome),
    paste0("do_enrichr\t",             do_enrichr),
    paste0("organism\t",               organism),
    paste0("enrichr_dbs\t",            paste(enrichr_dbs, collapse = ",")),
    paste0("use_enrichr_background\t", use_enrichr_background),
    paste0("background_ok\t",          bg_ok),
    paste0("min_bg_genes\t",           min_bg_genes),
    paste0("enrichr_retries\t",        enrichr_retries),
    paste0("enrichr_sleep\t",          enrichr_sleep),
    paste0("plot_p_col\t",             plot_p_col),
    paste0("plot_top_n\t",             plot_top_n),
    paste0("plot_sig_cutoff\t",        plot_sig_cutoff),
    paste0("plot_width\t",             plot_width),
    paste0("plot_height\t",            plot_height),
    paste0("n_selected_modules\t",     length(selected_modules)),
    paste0("selected_modules\t",       paste(sort(selected_modules), collapse = ", ")),
    paste0("out_dir\t",                out_dir)
  ), params_file)

  append_log(log_file, "Finished config", stem_label)
}

writeLines(capture.output(sessionInfo()),
           con = file.path(
             file.path(project_root, "comethyl_output", "consensus",
                       "16_enrichment", cpg_label, region_label, adjustment_version),
             "sessionInfo.txt"))

message("Script 16 complete: consensus module enrichment finished")