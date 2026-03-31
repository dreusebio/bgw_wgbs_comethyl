#!/usr/bin/env Rscript

# ============================================================
# 12B_enrich_selected_modules.R
#
# Purpose:
#   Reproducible enrichment step for selected comethyl modules.
#
# What this script does:
#   - loads annotated module-region output from 12A
#   - determines which modules to enrich
#   - runs KEGG and/or Enrichr
#   - writes per-module enrichment tables
#   - writes summary tables
#   - writes run logs
#
# What this script does NOT do:
#   - region annotation
#   - module-trait analysis
#   - presentation-specific filtering
#
# Ways to select modules:
#   1) --module_list_file
#      Plain text file, one module per line
#
#   2) --stats_file + --trait_list_file
#      Uses module-trait stats and selects modules with:
#        trait in trait_list_file
#        p <= p_thresh
#
# Inputs:
#   --project_root            required
#   --annotation_dir_v1       required
#   --annotation_dir_v2       optional
#   --annotation_dir_v3       optional
#
#   --module_list_file        optional
#   --stats_file_v1           optional
#   --stats_file_v2           optional
#   --stats_file_v3           optional
#   --trait_list_file         optional
#   --p_thresh                optional, default 0.05
#
#   --do_kegg                 optional, true/false
#   --do_enrichr              optional, true/false
#   --enrichr_dbs             optional, comma-separated
#   --organism                optional, default hsa
#
# Expected 12A files in each annotation dir:
#   Annotated_Regions.tsv
#   Module_Gene_List.tsv
#   Module_Gene_Summary.tsv
#
# Output structure per variant:
#   <annotation_dir>/enrichment/
#     selected_modules.txt
#     selected_modules_summary.tsv
#     run_parameters.txt
#     run_log.txt
#     <MODULE>_KEGG.tsv
#     <MODULE>_KEGG.xlsx
#     <MODULE>_KEGG_dotplot.pdf
#     <MODULE>_enrichr.xlsx
# ============================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# ============================================================
# 1) Basic helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  x
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
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

timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = logfile, append = TRUE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

# ============================================================
# 2) Input readers
# ============================================================
read_module_list_file <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

read_trait_list_file <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

load_annotation_files <- function(annotation_dir) {
  annotated_regions_file <- file.path(annotation_dir, "Annotated_Regions.tsv")
  gene_list_file         <- file.path(annotation_dir, "Module_Gene_List.tsv")
  module_summary_file    <- file.path(annotation_dir, "Module_Gene_Summary.tsv")

  validate_file_exists(annotated_regions_file, "Annotated_Regions.tsv")
  validate_file_exists(gene_list_file, "Module_Gene_List.tsv")
  validate_file_exists(module_summary_file, "Module_Gene_Summary.tsv")

  annotated_regions <- read.delim(
    annotated_regions_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  gene_list <- read.delim(
    gene_list_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  module_summary <- read.delim(
    module_summary_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  req_regions <- c("RegionID", "chr", "start", "end", "module")
  miss_regions <- setdiff(req_regions, colnames(annotated_regions))
  if (length(miss_regions) > 0) {
    stop("Annotated_Regions.tsv missing required columns: ",
         paste(miss_regions, collapse = ", "), call. = FALSE)
  }

  req_gene <- c("module", "gene_symbol")
  miss_gene <- setdiff(req_gene, colnames(gene_list))
  if (length(miss_gene) > 0) {
    stop("Module_Gene_List.tsv missing required columns: ",
         paste(miss_gene, collapse = ", "), call. = FALSE)
  }

  list(
    annotated_regions = annotated_regions,
    gene_list = gene_list,
    module_summary = module_summary
  )
}

# ============================================================
# 3) Stats-based module selection
# ============================================================
standardize_stats_columns <- function(df) {
  cn <- tolower(colnames(df))

  if ("module" %in% cn) {
    names(df)[match("module", cn)] <- "module"
  } else {
    stop("stats file must contain module column.", call. = FALSE)
  }

  if ("trait" %in% cn) {
    names(df)[match("trait", cn)] <- "trait"
  } else {
    stop("stats file must contain trait column.", call. = FALSE)
  }

  if ("p" %in% cn) {
    names(df)[match("p", cn)] <- "p"
  } else if ("pvalue" %in% cn) {
    names(df)[match("pvalue", cn)] <- "p"
  } else if ("p.value" %in% cn) {
    names(df)[match("p.value", cn)] <- "p"
  } else if ("p_value" %in% cn) {
    names(df)[match("p_value", cn)] <- "p"
  } else {
    stop("stats file must contain p-value column named p, pvalue, p.value, or p_value.", call. = FALSE)
  }

  df
}

read_stats_file <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    df <- openxlsx::read.xlsx(path)
  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("csv")) {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported stats file format: ", path, call. = FALSE)
  }

  standardize_stats_columns(df)
}

select_modules_from_stats <- function(stats_df, trait_list, p_thresh) {
  stats_df %>%
    dplyr::mutate(p = as.numeric(p)) %>%
    dplyr::filter(!is.na(p), trait %in% trait_list, p <= p_thresh) %>%
    dplyr::pull(module) %>%
    unique() %>%
    as.character()
}

# ============================================================
# 4) Gene list builder
# ============================================================
build_gene_vector_for_module <- function(gene_list_df, module_name) {
  gene_list_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene_symbol) %>%
    as.character() %>%
    trimws() %>%
    unique() %>%
    .[!is.na(.) & . != ""]
}

# ============================================================
# 5) KEGG enrichment
# ============================================================
run_kegg_clusterprofiler <- function(gene_symbols, out_prefix, organism = "hsa") {
  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]

  if (length(gene_symbols) == 0) return(invisible(NULL))

  gene_df <- suppressMessages(
    clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  if (is.null(gene_df) || nrow(gene_df) == 0) return(invisible(NULL))

  kegg_res <- suppressMessages(
    clusterProfiler::enrichKEGG(
      gene = unique(gene_df$ENTREZID),
      organism = organism,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  )

  kegg_df <- as.data.frame(kegg_res)
  if (is.null(kegg_df) || nrow(kegg_df) == 0) return(invisible(NULL))

  write.table(
    kegg_df,
    file = paste0(out_prefix, "_KEGG.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  openxlsx::write.xlsx(kegg_df, file = paste0(out_prefix, "_KEGG.xlsx"), rowNames = FALSE)

  if (all(c("Description", "p.adjust") %in% colnames(kegg_df))) {
    top <- kegg_df %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = min(25, nrow(kegg_df)))

    if (nrow(top) > 0) {
      top$Description <- factor(top$Description, levels = rev(top$Description))

      p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
        geom_point() +
        theme_bw(base_size = 12) +
        labs(title = "Top KEGG pathways", x = "-log10(adj p)", y = NULL)

      ggsave(
        filename = paste0(out_prefix, "_KEGG_dotplot.pdf"),
        plot = p,
        width = 9,
        height = 6
      )
    }
  }

  invisible(kegg_df)
}

# ============================================================
# 6) Enrichr
# ============================================================
run_enrichr_simple <- function(gene_symbols, out_prefix, dbs, site = "Enrichr") {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning("enrichR package not installed; skipping Enrichr.")
    return(invisible(NULL))
  }

  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]

  if (length(gene_symbols) == 0) return(invisible(NULL))
  if (length(dbs) == 0) return(invisible(NULL))

  enrichR::setEnrichrSite(site)

  res <- tryCatch(
    {
      enrichR::enrichr(gene_symbols, dbs)
    },
    error = function(e) {
      warning("Enrichr failed for ", out_prefix, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res) || length(res) == 0) return(invisible(NULL))

  wb <- openxlsx::createWorkbook()

  wrote_any <- FALSE
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
  }

  invisible(res)
}

# ============================================================
# 7) Per-variant enrichment runner
# ============================================================
run_enrichment_for_variant <- function(annotation_dir,
                                       variant_label,
                                       module_list_file = NULL,
                                       stats_file = NULL,
                                       trait_list = NULL,
                                       p_thresh = 0.05,
                                       do_kegg = TRUE,
                                       do_enrichr = FALSE,
                                       enrichr_dbs = character(0),
                                       organism = "hsa",
                                       project_root = NULL) {

  validate_dir_exists(annotation_dir, paste0("annotation_dir_", variant_label))

  files <- load_annotation_files(annotation_dir)
  gene_list_df <- files$gene_list
  annotated_regions <- files$annotated_regions
  module_summary <- files$module_summary

  out_dir <- file.path(annotation_dir, "enrichment")
  safe_dir_create(out_dir)

  log_file <- file.path(out_dir, "run_log.txt")
  params_file <- file.path(out_dir, "run_parameters.txt")
  selected_modules_file <- file.path(out_dir, "selected_modules.txt")
  selected_modules_summary_file <- file.path(out_dir, "selected_modules_summary.tsv")

  append_log(log_file, "Starting enrichment for variant: ", variant_label)
  append_log(log_file, "annotation_dir: ", annotation_dir)

  selected_modules <- character(0)
  selection_mode <- NULL

  if (!is.null(module_list_file)) {
    validate_file_exists(module_list_file, "module_list_file")
    selected_modules <- read_module_list_file(module_list_file)
    selection_mode <- "module_list_file"
  } else if (!is.null(stats_file) && !is.null(trait_list) && length(trait_list) > 0) {
    validate_file_exists(stats_file, paste0("stats_file_", variant_label))
    stats_df <- read_stats_file(stats_file)
    selected_modules <- select_modules_from_stats(stats_df, trait_list = trait_list, p_thresh = p_thresh)
    selection_mode <- "stats_file_plus_trait_list"
  } else {
    stop(
      "For ", variant_label,
      ", provide either --module_list_file OR both stats_file and trait_list_file.",
      call. = FALSE
    )
  }

  selected_modules <- unique(as.character(selected_modules))
  selected_modules <- selected_modules[nzchar(selected_modules)]
  selected_modules <- selected_modules[selected_modules %in% unique(gene_list_df$module)]

  write_lines_safe(selected_modules, selected_modules_file)

  sel_summary <- module_summary %>%
    dplyr::filter(module %in% selected_modules)

  write.table(
    sel_summary,
    selected_modules_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  append_log(log_file, "selection_mode: ", selection_mode)
  append_log(log_file, "selected modules count: ", length(selected_modules))
  append_log(log_file, "selected modules: ", paste(selected_modules, collapse = ", "))

  if (length(selected_modules) == 0) {
    append_log(log_file, "No selected modules remained after filtering to annotated modules.")
    return(invisible(NULL))
  }

  for (m in selected_modules) {
    gene_vec <- build_gene_vector_for_module(gene_list_df, m)
    out_prefix <- file.path(out_dir, m)

    append_log(log_file, "Module ", m, ": genes=", length(gene_vec))

    if (do_kegg) {
      run_kegg_clusterprofiler(
        gene_symbols = gene_vec,
        out_prefix = out_prefix,
        organism = organism
      )
    }

    if (do_enrichr) {
      run_enrichr_simple(
        gene_symbols = gene_vec,
        out_prefix = out_prefix,
        dbs = enrichr_dbs,
        site = "Enrichr"
      )
    }
  }

  params <- c(
    paste0("timestamp\t", timestamp_now()),
    paste0("project_root\t", ifelse(is.null(project_root), "", project_root)),
    paste0("variant_label\t", variant_label),
    paste0("annotation_dir\t", annotation_dir),
    paste0("selection_mode\t", selection_mode),
    paste0("module_list_file\t", ifelse(is.null(module_list_file), "", module_list_file)),
    paste0("stats_file\t", ifelse(is.null(stats_file), "", stats_file)),
    paste0("p_thresh\t", p_thresh),
    paste0("do_kegg\t", do_kegg),
    paste0("do_enrichr\t", do_enrichr),
    paste0("organism\t", organism),
    paste0("enrichr_dbs\t", paste(enrichr_dbs, collapse = ",")),
    paste0("n_selected_modules\t", length(selected_modules))
  )
  write_lines_safe(params, params_file)

  append_log(log_file, "Finished enrichment for variant: ", variant_label)

  invisible(list(
    selected_modules = selected_modules,
    out_dir = out_dir
  ))
}

# ============================================================
# 8) Read arguments
# ============================================================
project_root      <- trim_or_null(get_arg("--project_root"))

annotation_dir_v1 <- trim_or_null(get_arg("--annotation_dir_v1"))
annotation_dir_v2 <- trim_or_null(get_arg("--annotation_dir_v2"))
annotation_dir_v3 <- trim_or_null(get_arg("--annotation_dir_v3"))

module_list_file  <- trim_or_null(get_arg("--module_list_file"))

stats_file_v1     <- trim_or_null(get_arg("--stats_file_v1"))
stats_file_v2     <- trim_or_null(get_arg("--stats_file_v2"))
stats_file_v3     <- trim_or_null(get_arg("--stats_file_v3"))

trait_list_file   <- trim_or_null(get_arg("--trait_list_file"))
p_thresh          <- as.numeric(get_arg("--p_thresh", "0.05"))

do_kegg           <- to_bool(get_arg("--do_kegg", "true"), default = TRUE)
do_enrichr        <- to_bool(get_arg("--do_enrichr", "false"), default = FALSE)

enrichr_dbs       <- split_csv(get_arg(
  "--enrichr_dbs",
  "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025,ClinVar_2025,KEGG_2021_Human"
))

organism          <- trim_or_null(get_arg("--organism", "hsa"))

# ============================================================
# 9) Validate top-level arguments
# ============================================================
stop_if_missing(project_root, "--project_root")
stop_if_missing(annotation_dir_v1, "--annotation_dir_v1")

validate_dir_exists(project_root, "project_root")
validate_dir_exists(annotation_dir_v1, "annotation_dir_v1")
if (!is.null(annotation_dir_v2)) validate_dir_exists(annotation_dir_v2, "annotation_dir_v2")
if (!is.null(annotation_dir_v3)) validate_dir_exists(annotation_dir_v3, "annotation_dir_v3")

if (!is.null(module_list_file)) {
  validate_file_exists(module_list_file, "module_list_file")
} else {
  if (is.null(trait_list_file)) {
    stop("Provide either --module_list_file OR --trait_list_file with stats_file(s).", call. = FALSE)
  }
  validate_file_exists(trait_list_file, "trait_list_file")

  validate_file_exists(stats_file_v1, "stats_file_v1")
  if (!is.null(annotation_dir_v2) && is.null(stats_file_v2)) {
    stop("annotation_dir_v2 provided but stats_file_v2 missing.", call. = FALSE)
  }
  if (!is.null(annotation_dir_v3) && is.null(stats_file_v3)) {
    stop("annotation_dir_v3 provided but stats_file_v3 missing.", call. = FALSE)
  }
}

if (do_enrichr && length(enrichr_dbs) == 0) {
  stop("do_enrichr=TRUE but no --enrichr_dbs provided.", call. = FALSE)
}

setwd(project_root)

trait_list <- NULL
if (!is.null(trait_list_file)) {
  trait_list <- read_trait_list_file(trait_list_file)
}

cat("project_root: ", project_root, "\n", sep = "")
cat("annotation_dir_v1: ", annotation_dir_v1, "\n", sep = "")
if (!is.null(annotation_dir_v2)) cat("annotation_dir_v2: ", annotation_dir_v2, "\n", sep = "")
if (!is.null(annotation_dir_v3)) cat("annotation_dir_v3: ", annotation_dir_v3, "\n", sep = "")
cat("do_kegg: ", do_kegg, "\n", sep = "")
cat("do_enrichr: ", do_enrichr, "\n", sep = "")
cat("organism: ", organism, "\n", sep = "")

# ============================================================
# 10) Run variants
# ============================================================
run_enrichment_for_variant(
  annotation_dir = annotation_dir_v1,
  variant_label = "v1",
  module_list_file = module_list_file,
  stats_file = stats_file_v1,
  trait_list = trait_list,
  p_thresh = p_thresh,
  do_kegg = do_kegg,
  do_enrichr = do_enrichr,
  enrichr_dbs = enrichr_dbs,
  organism = organism,
  project_root = project_root
)

if (!is.null(annotation_dir_v2)) {
  run_enrichment_for_variant(
    annotation_dir = annotation_dir_v2,
    variant_label = "v2",
    module_list_file = module_list_file,
    stats_file = stats_file_v2,
    trait_list = trait_list,
    p_thresh = p_thresh,
    do_kegg = do_kegg,
    do_enrichr = do_enrichr,
    enrichr_dbs = enrichr_dbs,
    organism = organism,
    project_root = project_root
  )
}

if (!is.null(annotation_dir_v3)) {
  run_enrichment_for_variant(
    annotation_dir = annotation_dir_v3,
    variant_label = "v3",
    module_list_file = module_list_file,
    stats_file = stats_file_v3,
    trait_list = trait_list,
    p_thresh = p_thresh,
    do_kegg = do_kegg,
    do_enrichr = do_enrichr,
    enrichr_dbs = enrichr_dbs,
    organism = organism,
    project_root = project_root
  )
}

cat("\n12B enrichment complete.\n")