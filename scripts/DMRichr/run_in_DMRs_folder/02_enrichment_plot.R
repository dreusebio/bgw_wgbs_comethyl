#!/usr/bin/env Rscript
# ============================================================================
# 02_enrichment_plot.R
#
# Flexible enrichment + plotting pipeline for a gene list (e.g. DMR-annotated
# genes from DMRichR). GO and KEGG default to fully OFFLINE methods
# (clusterProfiler::enrichGO / enrichKEGG against locally installed OrgDb),
# matching how compute nodes without outbound internet access need to run.
# Enrichr (GO/KEGG alt methods, GWAS Catalog, --extra_dbs) requires internet.
# The Enrichr server addresses are initialized explicitly for compatibility
# with enrichR versions that otherwise fail with:
#   "Must specify at least one of url or handle"
# API calls are retried and degrade gracefully rather than halting the script.
#
#   --mode presentation : one combined overview figure, few terms, large fonts
#   --mode publication   : separate panel per database, more terms, print fonts
#
# In addition to the plots, a single multi-sheet Excel workbook
# (<prefix>_full_results.xlsx) is written with the FULL, unfiltered result
# table for every database tested (GO BP/CC/MF, KEGG, GWAS Catalog, any
# --extra_dbs) -- every tested term, not just what clears --pval_cutoff or
# fits in --top_n, each with the specific genes behind it in a "Genes"
# column. Disable with --skip_full_results.
#
# IMPORTANT: enrichR is NEVER loaded with library(enrichR) here. Some enrichR
# versions run a connectivity check in .onAttach() (attachNamespace()) that
# hard-fails the whole script if the compute node can't resolve
# maayanlab.cloud -- which is normal on HPC compute nodes with no outbound
# DNS. Calling enrichR only via requireNamespace() + enrichR::enrichr(...)
# skips .onAttach entirely, so the package loads fine even offline; only the
# actual API call fails, and that failure is caught per-database.
#
# Required packages (CRAN): optparse, dplyr, stringr, forcats, ggplot2, ggsci,
#   glue, purrr, Hmisc, openxlsx
# Required packages (Bioconductor): org.Hs.eg.db (or your OrgDb of choice),
#   clusterProfiler (default GO/KEGG methods), rrvgo (GO slimming)
# Optional (internet required): enrichR -- only needed for --go_method enrichr,
#   --kegg_method enrichr, GWAS Catalog (on by default; use --skip_gwas to
#   avoid it on offline nodes), or --extra_dbs
#
#   install.packages(c("optparse","dplyr","stringr","forcats","ggplot2",
#                       "ggsci","glue","purrr","Hmisc","openxlsx","enrichR"))
#   BiocManager::install(c("rrvgo","org.Hs.eg.db","clusterProfiler"))
#
# USAGE EXAMPLES
# ---------------------------------------------------------------------------
# Fully offline: GO + KEGG via clusterProfiler, no GWAS Catalog, no internet:
#   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --skip_gwas
#
# Publication figure, GO + KEGG (clusterProfiler) + GWAS Catalog (needs internet):
#   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --mode publication
#
# Presentation overview slide, KEGG via Enrichr instead, top 4 terms:
#   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --mode presentation \
#     --kegg_method enrichr --top_n 4
#
# Add a database beyond KEGG/GWAS Catalog:
#   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx \
#     --extra_dbs "Human_Phenotype_Ontology,Reactome_2022"
#
# Skip the full-results workbook (plots only):
#   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --skip_full_results
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggsci)
  library(glue)
  library(purrr)
  library(Hmisc)
  library(openxlsx)
})
# enrichR is intentionally NOT library()'d here -- see header note.

# ---- CLI arguments ---------------------------------------------------------

option_list <- list(
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "Path to gene list: .txt (one symbol per line) or .xlsx (needs a geneSymbol column) [required]"),
  make_option(c("-o", "--outdir"), type = "character", default = "enrichment_results",
              help = "Output directory [default %default]"),
  make_option(c("--prefix"), type = "character", default = "enrichment",
              help = "Prefix for output file names [default %default]"),

  make_option(c("--skip_go"), action = "store_true", default = FALSE,
              help = "Skip GO enrichment entirely"),
  make_option(c("--go_method"), type = "character", default = "clusterProfiler",
              help = "GO method: 'clusterProfiler' (offline, default) or 'enrichr' (needs internet) [default %default]"),
  make_option(c("--go_dbs"), type = "character",
              default = "GO_Biological_Process_2018,GO_Cellular_Component_2018,GO_Molecular_Function_2018",
              help = "Comma-separated Enrichr GO database names, used only if --go_method enrichr [default %default]"),
  make_option(c("--slim_threshold"), type = "double", default = 0.7,
              help = "rrvgo similarity threshold for slimming GO terms: 0.9 large, 0.7 medium, 0.5 small, 0.4 tiny [default %default]"),

  make_option(c("--skip_kegg"), action = "store_true", default = FALSE,
              help = "Skip KEGG enrichment"),
  make_option(c("--kegg_method"), type = "character", default = "clusterProfiler",
              help = "KEGG method: 'clusterProfiler' (offline, default) or 'enrichr' (needs internet) [default %default]"),
  make_option(c("--kegg_db"), type = "character", default = "KEGG_2021_Human",
              help = "Enrichr KEGG database name, used only if --kegg_method enrichr [default %default]"),
  make_option(c("--organism"), type = "character", default = "hsa",
              help = "KEGG organism code, used only if --kegg_method clusterProfiler [default %default]"),

  make_option(c("--skip_gwas"), action = "store_true", default = FALSE,
              help = "Skip GWAS Catalog enrichment (recommended on compute nodes with no internet access)"),
  make_option(c("--gwas_db"), type = "character", default = "GWAS_Catalog_2023",
              help = "Enrichr GWAS Catalog database name -- confirm current name via enrichR::listEnrichrDbs() before running [default %default]"),

  make_option(c("--extra_dbs"), type = "character", default = NULL,
              help = "Comma-separated additional Enrichr database names to test and plot (needs internet; one panel each in publication mode)"),

  make_option(c("--mode"), type = "character", default = "publication",
              help = "'publication' (separate panels, conservative top-N) or 'presentation' (one combined overview, fewer terms, larger fonts) [default %default]"),
  make_option(c("--top_n"), type = "integer", default = NULL,
              help = "Terms shown per database/category in the PLOTS. Defaults to 7 for publication, 3 for presentation. Does not limit the full-results workbook."),
  make_option(c("--pval_cutoff"), type = "double", default = 0.05,
              help = "P-value cutoff for retaining terms in the PLOTS [default %default]. Does not limit the full-results workbook."),
  make_option(c("--use_adjusted"), action = "store_true", default = FALSE,
              help = "Filter and rank plotted terms on adjusted p-value instead of raw p-value"),

  make_option(c("--skip_full_results"), action = "store_true", default = FALSE,
              help = "Skip writing <prefix>_full_results.xlsx (the complete, unfiltered per-database results with gene lists)"),

  make_option(c("--annoDb"), type = "character", default = "org.Hs.eg.db",
              help = "OrgDb annotation package name, used for GO/KEGG offline enrichment and GO slimming [default %default]"),
  make_option(c("--fig_format"), type = "character", default = NULL,
              help = "Output figure format: pdf or png [default pdf]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$genes)) stop("Please provide --genes <path>")
opt$mode <- match.arg(opt$mode, c("publication", "presentation"))
opt$kegg_method <- match.arg(opt$kegg_method, c("clusterProfiler", "enrichr"))
opt$go_method   <- match.arg(opt$go_method,   c("clusterProfiler", "enrichr"))
if (is.null(opt$top_n)) opt$top_n <- if (opt$mode == "presentation") 3 else 7
if (is.null(opt$fig_format)) opt$fig_format <- "pdf"
opt$fig_format <- match.arg(tolower(opt$fig_format), c("pdf", "png"))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

pval_col <- if (opt$use_adjusted) "Adjusted.P.value" else "P.value"

font_cfg <- if (opt$mode == "presentation") {
  list(axis_text = 20, axis_title = 20, legend_text = 18, legend_title = 18)
} else {
  list(axis_text = 12, axis_title = 14, legend_text = 12, legend_title = 14)
}

# ---- Safe Enrichr wrapper: explicit server setup + retry handling -----------

init_enrichr <- function() {
  options(
    enrichR.sites.base.address = "https://maayanlab.cloud/",
    enrichR.base.address = "https://maayanlab.cloud/Enrichr/",
    speedrichr.base.address = "https://maayanlab.cloud/speedrichr/api/",
    enrichR.live = TRUE,
    enrichR.quiet = TRUE,
    modEnrichR.use = TRUE,
    enrichR.sites = c(
      "Enrichr",
      "FlyEnrichr",
      "WormEnrichr",
      "YeastEnrichr",
      "FishEnrichr",
      "OxEnrichr"
    )
  )
}

safe_enrichr <- function(genes, dbs, context = "Enrichr",
                         retries = 3L, sleep_time = 2) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning(glue::glue("enrichR package not installed -- skipping {context}"))
    return(NULL)
  }

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  dbs <- unique(trimws(as.character(dbs)))
  dbs <- dbs[!is.na(dbs) & nzchar(dbs)]

  if (length(genes) == 0L) {
    warning(glue::glue("No valid genes supplied for {context} -- skipping"))
    return(NULL)
  }
  if (length(dbs) == 0L) {
    warning(glue::glue("No databases supplied for {context} -- skipping"))
    return(NULL)
  }

  init_enrichr()
  last_error <- "unknown error"

  for (attempt in seq_len(retries)) {
    message(glue::glue("{context}: Enrichr attempt {attempt} of {retries}"))

    result <- tryCatch(
      enrichR::enrichr(genes = genes, databases = dbs),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(result) && length(result) > 0L) return(result)

    if (attempt < retries) {
      message(glue::glue(
        "{context} failed: {last_error}. Retrying in {sleep_time} seconds."
      ))
      Sys.sleep(sleep_time)
    }
  }

  warning(glue::glue(
    "{context} failed after {retries} attempts ({last_error}) -- skipping. ",
    "Confirm that this compute node can reach maayanlab.cloud."
  ))
  NULL
}

# ---- Full-results workbook: every tested term + the genes behind it --------

full_wb <- if (!opt$skip_full_results) openxlsx::createWorkbook() else NULL

# Adds one sheet to the full-results workbook. Standardizes gene-list columns
# across sources: clusterProfiler's "geneID" (Entrez, "/"-separated, or gene
# symbols if setReadable()/readable=TRUE was used) is renamed to "Genes" and
# re-delimited with ";" to match Enrichr's native "Genes" column.
add_full_sheet <- function(wb, df, sheet_name) {
  if (is.null(wb) || is.null(df) || nrow(df) == 0) return(invisible(NULL))

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if ("geneID" %in% names(df)) {
    df$geneID <- gsub("/", ";", df$geneID)
    names(df)[names(df) == "geneID"] <- "Genes"
  }

  safe_name <- substr(gsub("[^A-Za-z0-9]", "_", sheet_name), 1, 31)
  existing  <- tryCatch(openxlsx::sheets(wb), error = function(e) character(0))
  candidate <- safe_name
  i <- 1
  while (candidate %in% existing) {
    i <- i + 1
    candidate <- substr(paste0(safe_name, "_", i), 1, 31)
  }

  openxlsx::addWorksheet(wb, candidate)
  openxlsx::writeData(wb, candidate, df)
  message(glue::glue("Added full-results sheet: {candidate} ({nrow(df)} rows)"))
}

# ---- Read gene list ---------------------------------------------------------

read_gene_list <- function(path) {
  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    df <- openxlsx::read.xlsx(path)
    if (!"geneSymbol" %in% names(df)) stop("Expected a 'geneSymbol' column in ", path)
    unique(na.omit(df$geneSymbol))
  } else {
    unique(na.omit(readLines(path)))
  }
}

genes <- read_gene_list(opt$genes)
message(glue::glue("Loaded {length(genes)} unique gene symbols from {opt$genes}"))
if (length(genes) == 0) stop("Gene list is empty")

# ---- Shared helper: flatten any Enrichr/clusterProfiler-style result table -
# NOTE: this is only used to build the PLOTS (top_n, pval_cutoff applied
# here). The full-results workbook is written separately from the raw,
# unfiltered result tables via add_full_sheet(), before this trimming.

format_flat <- function(df, label, n, pcol = pval_col, cutoff = opt$pval_cutoff) {
  if (is.null(df) || nrow(df) == 0) {
    warning(glue::glue("No results to format for {label}"))
    return(NULL)
  }
  if (!pcol %in% names(df)) {
    warning(glue::glue("Column '{pcol}' not found for {label}; skipping"))
    return(NULL)
  }
  df %>%
    dplyr::filter(.data[[pcol]] < cutoff) %>%
    dplyr::arrange(.data[[pcol]]) %>%
    dplyr::slice(1:min(n, dplyr::n())) %>%
    dplyr::transmute(
      `Gene Ontology` = label,
      Term = Term,
      `-log10.p-value` = -log10(.data[[pcol]])
    )
}

# ---- GO enrichment, two sources feeding the same rrvgo slimming step -------

run_go_offline_flat <- function(genes, annoDb, wb = NULL) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Install clusterProfiler for --go_method clusterProfiler")
  message("Running GO enrichment via clusterProfiler::enrichGO() (offline)...")
  gene_ids <- suppressMessages(clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = annoDb))
  if (is.null(gene_ids) || nrow(gene_ids) == 0) {
    warning("No Entrez IDs mapped for GO enrichment")
    return(NULL)
  }
  purrr::map_dfr(c("BP", "CC", "MF"), function(ont) {
    res <- tryCatch(
      clusterProfiler::enrichGO(gene = unique(gene_ids$ENTREZID), OrgDb = annoDb,
                                 ont = ont, pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
      error = function(e) {
        warning(glue::glue("GO {ont} enrichment failed: {conditionMessage(e)}"))
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (nrow(df) == 0) return(NULL)
    # Full, unfiltered per-ontology results (all tested terms + genes) --
    # this is the "readable" table with gene symbols in geneID, since
    # readable = TRUE above already converted Entrez IDs.
    add_full_sheet(wb, df, glue::glue("GO_{ont}_full"))
    df %>% dplyr::transmute(`Gene Ontology` = ont, go = ID, p = pvalue)
  }) %>%
    dplyr::filter(p < opt$pval_cutoff)
}

run_go_enrichr_flat <- function(genes, go_dbs, wb = NULL) {
  raw <- safe_enrichr(genes, go_dbs, context = "GO enrichment (Enrichr)")
  if (is.null(raw)) return(NULL)

  purrr::iwalk(raw, function(df, db_name) {
    add_full_sheet(wb, df, glue::glue("GO_{db_name}_full"))
  })

  purrr::imap_dfr(raw, function(df, db_name) {
    dplyr::as_tibble(df) %>% dplyr::mutate(db = db_name)
  }) %>%
    dplyr::mutate(
      `Gene Ontology` = dplyr::case_when(
        stringr::str_detect(db, "Biological_Process") ~ "BP",
        stringr::str_detect(db, "Cellular_Component")  ~ "CC",
        stringr::str_detect(db, "Molecular_Function")   ~ "MF",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(`Gene Ontology`)) %>%
    dplyr::mutate(
      Term = stringr::str_extract(Term, "\\(GO.*\\)") %>% stringr::str_remove_all("[()]"),
      p = .data[[pval_col]]
    ) %>%
    dplyr::select(p, go = Term, `Gene Ontology`) %>%
    dplyr::filter(p < opt$pval_cutoff)
}

slim_go <- function(genes, method, go_dbs, annoDb, threshold, wb = NULL) {

  go_flat <- if (method == "clusterProfiler") {
    run_go_offline_flat(genes, annoDb, wb = wb)
  } else {
    run_go_enrichr_flat(genes, go_dbs, wb = wb)
  }

  if (is.null(go_flat) || nrow(go_flat) == 0) {
    if (method == "enrichr") {
      warning(paste(
        "GO enrichment returned no usable significant terms -- skipping GO slimming.",
        "Review the preceding Enrichr messages to distinguish an API failure",
        "from a successful analysis with no terms passing the cutoff."
      ))
    } else {
      warning("No significant GO terms found -- skipping GO slimming")
    }
    return(NULL)
  }

  slim_one_ontology <- function(ont_label) {
    sub <- go_flat %>% dplyr::filter(`Gene Ontology` == ont_label)
    if (nrow(sub) == 0) return(NULL)

    message(glue::glue("Slimming {nrow(sub)} {ont_label} terms with rrvgo (threshold = {threshold})..."))

    simMatrix <- rrvgo::calculateSimMatrix(sub$go, orgdb = annoDb, ont = ont_label, method = "Rel")
    reduced   <- rrvgo::reduceSimMatrix(simMatrix, setNames(-log10(sub$p), sub$go),
                                         threshold = threshold, orgdb = annoDb)

    reduced %>%
      dplyr::as_tibble() %>%
      dplyr::inner_join(sub, by = "go") %>%
      dplyr::filter(term == as.character(parentTerm)) %>%
      dplyr::mutate(`-log10.p-value` = -log10(p), `Gene Ontology` = ont_label) %>%
      dplyr::select(`Gene Ontology`, Term = term, `-log10.p-value`)
  }

  dplyr::bind_rows(slim_one_ontology("BP"), slim_one_ontology("CC"), slim_one_ontology("MF")) %>%
    dplyr::arrange(dplyr::desc(`-log10.p-value`))
}

# ---- KEGG: clusterProfiler (offline, default) or Enrichr -------------------

run_kegg <- function(genes, method, organism, kegg_db, annoDb, wb = NULL) {

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]

  if (method == "clusterProfiler") {

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Install clusterProfiler for --kegg_method clusterProfiler")
    }

    # annoDb must be an OrgDb object, not only the package name.
    if (is.character(annoDb)) {
      if (!requireNamespace(annoDb, quietly = TRUE)) {
        stop("Annotation package is not installed: ", annoDb)
      }
      annoDb_object <- getExportedValue(annoDb, annoDb)
    } else {
      annoDb_object <- annoDb
    }

    message("Converting ", length(genes), " gene symbols to Entrez IDs...")

    gene_ids <- suppressMessages(
      clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = annoDb_object)
    )

    if (is.null(gene_ids) || nrow(gene_ids) == 0) {
      warning("No input gene symbols mapped to Entrez IDs for KEGG.")
      return(NULL)
    }

    gene_ids <- gene_ids[!is.na(gene_ids$ENTREZID) & gene_ids$ENTREZID != "", , drop = FALSE]
    entrez_ids <- unique(as.character(gene_ids$ENTREZID))

    message(
      "Mapped ", length(unique(gene_ids$SYMBOL)), " of ", length(genes),
      " gene symbols to ", length(entrez_ids), " unique Entrez IDs."
    )

    if (length(entrez_ids) == 0) {
      warning("No valid Entrez IDs remain after filtering.")
      return(NULL)
    }

    message("Example Entrez IDs: ", paste(utils::head(entrez_ids, 10), collapse = ", "))
    message("Running KEGG enrichment via clusterProfiler::enrichKEGG()...")

    res <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = entrez_ids, organism = organism, keyType = "kegg",
        pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"
      ),
      error = function(e) {
        warning(glue::glue(
          "KEGG enrichment failed: {conditionMessage(e)}. ",
          "enrichKEGG() accesses rest.kegg.jp and therefore requires internet access."
        ))
        NULL
      }
    )

    if (is.null(res)) return(NULL)

    # Convert Entrez IDs -> gene symbols so the exported table (and
    # downstream "Genes" column) is readable, matching enrichGO's
    # readable = TRUE behavior.
    res <- tryCatch(
      clusterProfiler::setReadable(res, OrgDb = annoDb_object, keyType = "ENTREZID"),
      error = function(e) {
        message(glue::glue("Could not convert KEGG gene IDs to symbols: {conditionMessage(e)} -- keeping Entrez IDs"))
        res
      }
    )

    res_df <- as.data.frame(res)

    if (nrow(res_df) == 0) {
      warning("The genes mapped to Entrez IDs, but no KEGG pathways passed the enrichment criteria.")
      return(NULL)
    }

    add_full_sheet(wb, res_df, "KEGG_full")

    res_df %>%
      dplyr::rename(Term = Description, P.value = pvalue, Adjusted.P.value = p.adjust) %>%
      dplyr::as_tibble()

  } else if (method == "enrichr") {

    message(glue::glue("Running KEGG enrichment via Enrichr ({kegg_db})..."))

    raw <- safe_enrichr(genes, kegg_db, context = "KEGG enrichment (Enrichr)")
    if (is.null(raw)) return(NULL)

    result <- raw[[kegg_db]]

    if (is.null(result) || nrow(result) == 0) {
      warning(glue::glue("No KEGG enrichment results returned from {kegg_db}."))
      return(NULL)
    }

    add_full_sheet(wb, result, "KEGG_full")

    dplyr::as_tibble(result)

  } else {
    stop("Unknown KEGG method: ", method, ". Use 'clusterProfiler' or 'enrichr'.")
  }
}

# ---- Run all requested databases -------------------------------------------

go_slimmed <- if (!opt$skip_go) {
  slim_go(genes, opt$go_method, trimws(stringr::str_split(opt$go_dbs, ",")[[1]]),
          opt$annoDb, opt$slim_threshold, wb = full_wb)
} else NULL

go_top <- if (!is.null(go_slimmed) && nrow(go_slimmed) > 0) {
  go_slimmed %>%
    dplyr::group_by(`Gene Ontology`) %>%
    dplyr::slice_max(order_by = `-log10.p-value`, n = opt$top_n) %>%
    dplyr::ungroup()
} else NULL

kegg_flat <- if (!opt$skip_kegg) {
  kegg_raw <- run_kegg(genes, opt$kegg_method, opt$organism, opt$kegg_db, opt$annoDb, wb = full_wb)
  format_flat(kegg_raw, "KEGG", opt$top_n)
} else NULL

gwas_flat <- if (!opt$skip_gwas) {
  gwas_raw_list <- safe_enrichr(genes, opt$gwas_db, context = "GWAS Catalog enrichment")
  gwas_raw <- if (!is.null(gwas_raw_list)) gwas_raw_list[[opt$gwas_db]] else NULL
  if (!is.null(gwas_raw)) add_full_sheet(full_wb, gwas_raw, "GWAS_Catalog_full")
  if (is.null(gwas_raw)) NULL else format_flat(gwas_raw, "GWAS Catalog", opt$top_n)
} else NULL

extra_flat <- NULL
if (!is.null(opt$extra_dbs)) {
  extra_dbs <- trimws(stringr::str_split(opt$extra_dbs, ",")[[1]])
  extra_raw <- safe_enrichr(genes, extra_dbs, context = "extra database enrichment")
  if (!is.null(extra_raw)) {
    purrr::iwalk(extra_raw, function(df, db_name) {
      add_full_sheet(full_wb, df, glue::glue("{db_name}_full"))
    })
    extra_flat <- purrr::imap_dfr(extra_raw, function(df, db_name) {
      res <- format_flat(df, db_name, opt$top_n)
      if (is.null(res)) return(dplyr::tibble())
      res
    })
  }
}

# ---- Save the full-results workbook (independent of plotting) -------------

if (!opt$skip_full_results && !is.null(full_wb)) {
  n_sheets <- length(tryCatch(openxlsx::sheets(full_wb), error = function(e) character(0)))
  if (n_sheets > 0) {
    full_results_path <- file.path(opt$outdir, glue::glue("{opt$prefix}_full_results.xlsx"))
    openxlsx::saveWorkbook(full_wb, full_results_path, overwrite = TRUE)
    message(glue::glue("Saved full results workbook ({n_sheets} sheets): {full_results_path}"))
  } else {
    message("No full-results sheets to save (no databases returned data).")
  }
}

# ---- Plotting ---------------------------------------------------------------

make_barplot <- function(df, title = NULL, fonts = font_cfg) {
  df <- df %>%
    dplyr::mutate(
      Term = stringr::str_trim(Term),
      Term = Hmisc::capitalize(Term),
      Term = stringr::str_wrap(Term, 50),
      Term = factor(Term, levels = unique(Term[order(forcats::fct_rev(`Gene Ontology`), `-log10.p-value`)]))
    )

  ggplot2::ggplot(df, ggplot2::aes(x = Term, y = `-log10.p-value`,
                                    fill = `Gene Ontology`, group = `Gene Ontology`)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "black") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
    ggsci::scale_fill_d3() +
    ggplot2::labs(y = expression("-log"[10](p)), title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = fonts$axis_text, color = "black"),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = fonts$axis_title),
      legend.text = ggplot2::element_text(size = fonts$legend_text),
      legend.title = ggplot2::element_text(size = fonts$legend_title),
      plot.title = ggplot2::element_text(size = fonts$axis_title + 2, face = "bold", hjust = 0.5)
    )
}

save_plot <- function(p, filename, width, height) {
  path <- file.path(opt$outdir, filename)
  ggplot2::ggsave(
    filename = path,
    plot = p,
    device = opt$fig_format,
    width = width,
    height = height,
    units = "in",
    dpi = if (opt$mode == "presentation") 150 else 300
  )
  message(glue::glue("Saved {path}"))
}

ext <- opt$fig_format

all_empty <- all(sapply(list(go_top, kegg_flat, gwas_flat, extra_flat), function(x) is.null(x) || nrow(x) == 0))
if (all_empty) stop("No significant terms found across any requested database -- nothing to plot")

if (opt$mode == "presentation") {

  combined <- dplyr::bind_rows(go_top, kegg_flat, gwas_flat, extra_flat)
  p <- make_barplot(combined)
  save_plot(p, glue::glue("{opt$prefix}_overview.{ext}"), width = 11, height = max(4, nrow(combined) * 0.5))

} else {

  if (!is.null(go_top) && nrow(go_top) > 0) {
    p_go <- make_barplot(go_top, title = "Gene ontology")
    save_plot(p_go, glue::glue("{opt$prefix}_GO.{ext}"), width = 9, height = max(4, nrow(go_top) * 0.35))
  }

  if (!is.null(kegg_flat) && nrow(kegg_flat) > 0) {
    p_kegg <- make_barplot(kegg_flat, title = "KEGG pathways")
    save_plot(p_kegg, glue::glue("{opt$prefix}_KEGG.{ext}"), width = 9, height = max(3, nrow(kegg_flat) * 0.35))
  }

  if (!is.null(gwas_flat) && nrow(gwas_flat) > 0) {
    p_gwas <- make_barplot(gwas_flat, title = "GWAS Catalog")
    save_plot(p_gwas, glue::glue("{opt$prefix}_GWAS.{ext}"), width = 9, height = max(3, nrow(gwas_flat) * 0.35))
  }

  if (!is.null(extra_flat) && nrow(extra_flat) > 0) {
    for (db in unique(extra_flat$`Gene Ontology`)) {
      sub <- extra_flat %>% dplyr::filter(`Gene Ontology` == db)
      p_extra <- make_barplot(sub, title = db)
      save_plot(p_extra, glue::glue("{opt$prefix}_{db}.{ext}"), width = 9, height = max(3, nrow(sub) * 0.35))
    }
  }
}

message("Done.")
# Working. I wanted to add excel outputs
# #!/usr/bin/env Rscript
# # ============================================================================
# # 02_enrichment_plot.R
# #
# # Flexible enrichment + plotting pipeline for a gene list (e.g. DMR-annotated
# # genes from DMRichR). GO and KEGG default to fully OFFLINE methods
# # (clusterProfiler::enrichGO / enrichKEGG against locally installed OrgDb),
# # matching how compute nodes without outbound internet access need to run.
# # Enrichr (GO/KEGG alt methods, GWAS Catalog, --extra_dbs) requires internet.
# # The Enrichr server addresses are initialized explicitly for compatibility
# # with enrichR versions that otherwise fail with:
# #   "Must specify at least one of url or handle"
# # API calls are retried and degrade gracefully rather than halting the script.
# #
# #   --mode presentation : one combined overview figure, few terms, large fonts
# #   --mode publication   : separate panel per database, more terms, print fonts
# #
# # IMPORTANT: enrichR is NEVER loaded with library(enrichR) here. Some enrichR
# # versions run a connectivity check in .onAttach() (attachNamespace()) that
# # hard-fails the whole script if the compute node can't resolve
# # maayanlab.cloud -- which is normal on HPC compute nodes with no outbound
# # DNS. Calling enrichR only via requireNamespace() + enrichR::enrichr(...)
# # skips .onAttach entirely, so the package loads fine even offline; only the
# # actual API call fails, and that failure is caught per-database.
# #
# # Required packages (CRAN): optparse, dplyr, stringr, forcats, ggplot2, ggsci,
# #   glue, purrr, Hmisc, openxlsx (only if reading .xlsx gene lists)
# # Required packages (Bioconductor): org.Hs.eg.db (or your OrgDb of choice),
# #   clusterProfiler (default GO/KEGG methods), rrvgo (GO slimming)
# # Optional (internet required): enrichR -- only needed for --go_method enrichr,
# #   --kegg_method enrichr, GWAS Catalog (on by default; use --skip_gwas to
# #   avoid it on offline nodes), or --extra_dbs
# #
# #   install.packages(c("optparse","dplyr","stringr","forcats","ggplot2",
# #                       "ggsci","glue","purrr","Hmisc","openxlsx","enrichR"))
# #   BiocManager::install(c("rrvgo","org.Hs.eg.db","clusterProfiler"))
# #
# # USAGE EXAMPLES
# # ---------------------------------------------------------------------------
# # Fully offline: GO + KEGG via clusterProfiler, no GWAS Catalog, no internet:
# #   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --skip_gwas
# #
# # Publication figure, GO + KEGG (clusterProfiler) + GWAS Catalog (needs internet):
# #   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --mode publication
# #
# # Presentation overview slide, KEGG via Enrichr instead, top 4 terms:
# #   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx --mode presentation \
# #     --kegg_method enrichr --top_n 4
# #
# # Add a database beyond KEGG/GWAS Catalog:
# #   Rscript 02_enrichment_plot.R --genes DMRs_annotated.xlsx \
# #     --extra_dbs "Human_Phenotype_Ontology,Reactome_2022"
# # ============================================================================

# suppressPackageStartupMessages({
#   library(optparse)
#   library(dplyr)
#   library(stringr)
#   library(forcats)
#   library(ggplot2)
#   library(ggsci)
#   library(glue)
#   library(purrr)
#   library(Hmisc)
# })
# # enrichR is intentionally NOT library()'d here -- see header note.

# # ---- CLI arguments ---------------------------------------------------------

# option_list <- list(
#   make_option(c("-g", "--genes"), type = "character", default = NULL,
#               help = "Path to gene list: .txt (one symbol per line) or .xlsx (needs a geneSymbol column) [required]"),
#   make_option(c("-o", "--outdir"), type = "character", default = "enrichment_results",
#               help = "Output directory [default %default]"),
#   make_option(c("--prefix"), type = "character", default = "enrichment",
#               help = "Prefix for output file names [default %default]"),

#   make_option(c("--skip_go"), action = "store_true", default = FALSE,
#               help = "Skip GO enrichment entirely"),
#   make_option(c("--go_method"), type = "character", default = "clusterProfiler",
#               help = "GO method: 'clusterProfiler' (offline, default) or 'enrichr' (needs internet) [default %default]"),
#   make_option(c("--go_dbs"), type = "character",
#               default = "GO_Biological_Process_2018,GO_Cellular_Component_2018,GO_Molecular_Function_2018",
#               help = "Comma-separated Enrichr GO database names, used only if --go_method enrichr [default %default]"),
#   make_option(c("--slim_threshold"), type = "double", default = 0.7,
#               help = "rrvgo similarity threshold for slimming GO terms: 0.9 large, 0.7 medium, 0.5 small, 0.4 tiny [default %default]"),

#   make_option(c("--skip_kegg"), action = "store_true", default = FALSE,
#               help = "Skip KEGG enrichment"),
#   make_option(c("--kegg_method"), type = "character", default = "clusterProfiler",
#               help = "KEGG method: 'clusterProfiler' (offline, default) or 'enrichr' (needs internet) [default %default]"),
#   make_option(c("--kegg_db"), type = "character", default = "KEGG_2021_Human",
#               help = "Enrichr KEGG database name, used only if --kegg_method enrichr [default %default]"),
#   make_option(c("--organism"), type = "character", default = "hsa",
#               help = "KEGG organism code, used only if --kegg_method clusterProfiler [default %default]"),

#   make_option(c("--skip_gwas"), action = "store_true", default = FALSE,
#               help = "Skip GWAS Catalog enrichment (recommended on compute nodes with no internet access)"),
#   make_option(c("--gwas_db"), type = "character", default = "GWAS_Catalog_2023",
#               help = "Enrichr GWAS Catalog database name -- confirm current name via enrichR::listEnrichrDbs() before running [default %default]"),

#   make_option(c("--extra_dbs"), type = "character", default = NULL,
#               help = "Comma-separated additional Enrichr database names to test and plot (needs internet; one panel each in publication mode)"),

#   make_option(c("--mode"), type = "character", default = "publication",
#               help = "'publication' (separate panels, conservative top-N) or 'presentation' (one combined overview, fewer terms, larger fonts) [default %default]"),
#   make_option(c("--top_n"), type = "integer", default = NULL,
#               help = "Terms shown per database/category. Defaults to 7 for publication, 3 for presentation"),
#   make_option(c("--pval_cutoff"), type = "double", default = 0.05,
#               help = "P-value cutoff for retaining terms [default %default]"),
#   make_option(c("--use_adjusted"), action = "store_true", default = FALSE,
#               help = "Filter and rank terms on adjusted p-value instead of raw p-value"),

#   make_option(c("--annoDb"), type = "character", default = "org.Hs.eg.db",
#               help = "OrgDb annotation package name, used for GO/KEGG offline enrichment and GO slimming [default %default]"),
#   make_option(c("--fig_format"), type = "character", default = NULL,
#               help = "Output figure format: pdf or png [default pdf]")
# )

# opt <- parse_args(OptionParser(option_list = option_list))

# if (is.null(opt$genes)) stop("Please provide --genes <path>")
# opt$mode <- match.arg(opt$mode, c("publication", "presentation"))
# opt$kegg_method <- match.arg(opt$kegg_method, c("clusterProfiler", "enrichr"))
# opt$go_method   <- match.arg(opt$go_method,   c("clusterProfiler", "enrichr"))
# if (is.null(opt$top_n)) opt$top_n <- if (opt$mode == "presentation") 3 else 7
# if (is.null(opt$fig_format)) opt$fig_format <- "pdf"
# opt$fig_format <- match.arg(tolower(opt$fig_format), c("pdf", "png"))

# dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# pval_col <- if (opt$use_adjusted) "Adjusted.P.value" else "P.value"

# font_cfg <- if (opt$mode == "presentation") {
#   list(axis_text = 20, axis_title = 20, legend_text = 18, legend_title = 18)
# } else {
#   list(axis_text = 12, axis_title = 14, legend_text = 12, legend_title = 14)
# }

# # ---- Safe Enrichr wrapper: explicit server setup + retry handling -----------

# init_enrichr <- function() {
#   options(
#     enrichR.sites.base.address = "https://maayanlab.cloud/",
#     enrichR.base.address = "https://maayanlab.cloud/Enrichr/",
#     speedrichr.base.address = "https://maayanlab.cloud/speedrichr/api/",
#     enrichR.live = TRUE,
#     enrichR.quiet = TRUE,
#     modEnrichR.use = TRUE,
#     enrichR.sites = c(
#       "Enrichr",
#       "FlyEnrichr",
#       "WormEnrichr",
#       "YeastEnrichr",
#       "FishEnrichr",
#       "OxEnrichr"
#     )
#   )
# }

# safe_enrichr <- function(genes, dbs, context = "Enrichr",
#                          retries = 3L, sleep_time = 2) {
#   if (!requireNamespace("enrichR", quietly = TRUE)) {
#     warning(glue::glue("enrichR package not installed -- skipping {context}"))
#     return(NULL)
#   }

#   genes <- unique(trimws(as.character(genes)))
#   genes <- genes[!is.na(genes) & nzchar(genes)]
#   dbs <- unique(trimws(as.character(dbs)))
#   dbs <- dbs[!is.na(dbs) & nzchar(dbs)]

#   if (length(genes) == 0L) {
#     warning(glue::glue("No valid genes supplied for {context} -- skipping"))
#     return(NULL)
#   }
#   if (length(dbs) == 0L) {
#     warning(glue::glue("No databases supplied for {context} -- skipping"))
#     return(NULL)
#   }

#   init_enrichr()
#   last_error <- "unknown error"

#   for (attempt in seq_len(retries)) {
#     message(glue::glue("{context}: Enrichr attempt {attempt} of {retries}"))

#     result <- tryCatch(
#       enrichR::enrichr(genes = genes, databases = dbs),
#       error = function(e) {
#         last_error <<- conditionMessage(e)
#         NULL
#       }
#     )

#     if (!is.null(result) && length(result) > 0L) return(result)

#     if (attempt < retries) {
#       message(glue::glue(
#         "{context} failed: {last_error}. Retrying in {sleep_time} seconds."
#       ))
#       Sys.sleep(sleep_time)
#     }
#   }

#   warning(glue::glue(
#     "{context} failed after {retries} attempts ({last_error}) -- skipping. ",
#     "Confirm that this compute node can reach maayanlab.cloud."
#   ))
#   NULL
# }

# # ---- Read gene list ---------------------------------------------------------

# read_gene_list <- function(path) {
#   if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
#     if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Install openxlsx to read .xlsx gene lists")
#     df <- openxlsx::read.xlsx(path)
#     if (!"geneSymbol" %in% names(df)) stop("Expected a 'geneSymbol' column in ", path)
#     unique(na.omit(df$geneSymbol))
#   } else {
#     unique(na.omit(readLines(path)))
#   }
# }

# genes <- read_gene_list(opt$genes)
# message(glue::glue("Loaded {length(genes)} unique gene symbols from {opt$genes}"))
# if (length(genes) == 0) stop("Gene list is empty")

# # ---- Shared helper: flatten any Enrichr/clusterProfiler-style result table -

# format_flat <- function(df, label, n, pcol = pval_col, cutoff = opt$pval_cutoff) {
#   if (is.null(df) || nrow(df) == 0) {
#     warning(glue::glue("No results to format for {label}"))
#     return(NULL)
#   }
#   if (!pcol %in% names(df)) {
#     warning(glue::glue("Column '{pcol}' not found for {label}; skipping"))
#     return(NULL)
#   }
#   df %>%
#     dplyr::filter(.data[[pcol]] < cutoff) %>%
#     dplyr::arrange(.data[[pcol]]) %>%
#     dplyr::slice(1:min(n, dplyr::n())) %>%
#     dplyr::transmute(
#       `Gene Ontology` = label,
#       Term = Term,
#       `-log10.p-value` = -log10(.data[[pcol]])
#     )
# }

# # ---- GO enrichment, two sources feeding the same rrvgo slimming step -------

# run_go_offline_flat <- function(genes, annoDb) {
#   if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Install clusterProfiler for --go_method clusterProfiler")
#   message("Running GO enrichment via clusterProfiler::enrichGO() (offline)...")
#   gene_ids <- suppressMessages(clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = annoDb))
#   if (is.null(gene_ids) || nrow(gene_ids) == 0) {
#     warning("No Entrez IDs mapped for GO enrichment")
#     return(NULL)
#   }
#   purrr::map_dfr(c("BP", "CC", "MF"), function(ont) {
#     res <- tryCatch(
#       clusterProfiler::enrichGO(gene = unique(gene_ids$ENTREZID), OrgDb = annoDb,
#                                  ont = ont, pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
#       error = function(e) {
#         warning(glue::glue("GO {ont} enrichment failed: {conditionMessage(e)}"))
#         NULL
#       }
#     )
#     if (is.null(res)) return(NULL)
#     df <- as.data.frame(res)
#     if (nrow(df) == 0) return(NULL)
#     df %>% dplyr::transmute(`Gene Ontology` = ont, go = ID, p = pvalue)
#   }) %>%
#     dplyr::filter(p < opt$pval_cutoff)
# }

# run_go_enrichr_flat <- function(genes, go_dbs) {
#   raw <- safe_enrichr(genes, go_dbs, context = "GO enrichment (Enrichr)")
#   if (is.null(raw)) return(NULL)

#   purrr::imap_dfr(raw, function(df, db_name) {
#     dplyr::as_tibble(df) %>% dplyr::mutate(db = db_name)
#   }) %>%
#     dplyr::mutate(
#       `Gene Ontology` = dplyr::case_when(
#         stringr::str_detect(db, "Biological_Process") ~ "BP",
#         stringr::str_detect(db, "Cellular_Component")  ~ "CC",
#         stringr::str_detect(db, "Molecular_Function")   ~ "MF",
#         TRUE ~ NA_character_
#       )
#     ) %>%
#     dplyr::filter(!is.na(`Gene Ontology`)) %>%
#     dplyr::mutate(
#       Term = stringr::str_extract(Term, "\\(GO.*\\)") %>% stringr::str_remove_all("[()]"),
#       p = .data[[pval_col]]
#     ) %>%
#     dplyr::select(p, go = Term, `Gene Ontology`) %>%
#     dplyr::filter(p < opt$pval_cutoff)
# }

# slim_go <- function(genes, method, go_dbs, annoDb, threshold) {

#   go_flat <- if (method == "clusterProfiler") {
#     run_go_offline_flat(genes, annoDb)
#   } else {
#     run_go_enrichr_flat(genes, go_dbs)
#   }

#   if (is.null(go_flat) || nrow(go_flat) == 0) {
#     if (method == "enrichr") {
#       warning(paste(
#         "GO enrichment returned no usable significant terms -- skipping GO slimming.",
#         "Review the preceding Enrichr messages to distinguish an API failure",
#         "from a successful analysis with no terms passing the cutoff."
#       ))
#     } else {
#       warning("No significant GO terms found -- skipping GO slimming")
#     }
#     return(NULL)
#   }

#   slim_one_ontology <- function(ont_label) {
#     sub <- go_flat %>% dplyr::filter(`Gene Ontology` == ont_label)
#     if (nrow(sub) == 0) return(NULL)

#     message(glue::glue("Slimming {nrow(sub)} {ont_label} terms with rrvgo (threshold = {threshold})..."))

#     simMatrix <- rrvgo::calculateSimMatrix(sub$go, orgdb = annoDb, ont = ont_label, method = "Rel")
#     reduced   <- rrvgo::reduceSimMatrix(simMatrix, setNames(-log10(sub$p), sub$go),
#                                          threshold = threshold, orgdb = annoDb)

#     reduced %>%
#       dplyr::as_tibble() %>%
#       dplyr::inner_join(sub, by = "go") %>%
#       dplyr::filter(term == as.character(parentTerm)) %>%
#       dplyr::mutate(`-log10.p-value` = -log10(p), `Gene Ontology` = ont_label) %>%
#       dplyr::select(`Gene Ontology`, Term = term, `-log10.p-value`)
#   }

#   dplyr::bind_rows(slim_one_ontology("BP"), slim_one_ontology("CC"), slim_one_ontology("MF")) %>%
#     dplyr::arrange(dplyr::desc(`-log10.p-value`))
# }

# # ---- KEGG: clusterProfiler (offline, default) or Enrichr -------------------

# # run_kegg <- function(genes, method, organism, kegg_db, annoDb) {
# #   if (method == "clusterProfiler") {
# #     if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Install clusterProfiler for --kegg_method clusterProfiler")
# #     message("Running KEGG enrichment via clusterProfiler::enrichKEGG()...")
# #     gene_ids <- suppressMessages(clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = annoDb))
# #     if (is.null(gene_ids) || nrow(gene_ids) == 0) {
# #       warning("No Entrez IDs mapped for KEGG enrichment")
# #       return(NULL)
# #     }
# #     res <- tryCatch(
# #       clusterProfiler::enrichKEGG(gene = gene_ids$ENTREZID, organism = organism, pvalueCutoff = 1),
# #       error = function(e) {
# #         warning(glue::glue("KEGG enrichment failed: {conditionMessage(e)}. ",
# #                             "Note: enrichKEGG() itself calls out to rest.kegg.jp and also needs internet access."))
# #         NULL
# #       }
# #     )
# #     if (is.null(res)) return(NULL)
# #     res_df <- as.data.frame(res)
# #     if (nrow(res_df) == 0) return(NULL)
# #     res_df %>%
# #       dplyr::rename(Term = Description, P.value = pvalue, Adjusted.P.value = p.adjust) %>%
# #       dplyr::as_tibble()
# #   } else {
# #     message(glue::glue("Running KEGG enrichment via Enrichr ({kegg_db})..."))
# #     raw <- safe_enrichr(genes, kegg_db, context = "KEGG enrichment (Enrichr)")
# #     if (is.null(raw)) return(NULL)
# #     result <- raw[[kegg_db]]
# #     if (is.null(result)) return(NULL)
# #     dplyr::as_tibble(result)
# #   }
# # }

# run_kegg <- function(genes, method, organism, kegg_db, annoDb) {

#   genes <- unique(trimws(as.character(genes)))
#   genes <- genes[!is.na(genes) & genes != ""]

#   if (method == "clusterProfiler") {

#     if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
#       stop(
#         "Install clusterProfiler for ",
#         "--kegg_method clusterProfiler"
#       )
#     }

#     # annoDb must be an OrgDb object, not only the package name.
#     if (is.character(annoDb)) {
#       if (!requireNamespace(annoDb, quietly = TRUE)) {
#         stop("Annotation package is not installed: ", annoDb)
#       }

#       annoDb_object <- getExportedValue(annoDb, annoDb)
#     } else {
#       annoDb_object <- annoDb
#     }

#     message(
#       "Converting ", length(genes),
#       " gene symbols to Entrez IDs..."
#     )

#     gene_ids <- suppressMessages(
#       clusterProfiler::bitr(
#         genes,
#         fromType = "SYMBOL",
#         toType   = "ENTREZID",
#         OrgDb    = annoDb_object
#       )
#     )

#     if (is.null(gene_ids) || nrow(gene_ids) == 0) {
#       warning("No input gene symbols mapped to Entrez IDs for KEGG.")
#       return(NULL)
#     }

#     # Remove duplicate and missing identifiers and force character format.
#     gene_ids <- gene_ids[
#       !is.na(gene_ids$ENTREZID) &
#         gene_ids$ENTREZID != "",
#       ,
#       drop = FALSE
#     ]

#     entrez_ids <- unique(as.character(gene_ids$ENTREZID))

#     message(
#       "Mapped ",
#       length(unique(gene_ids$SYMBOL)),
#       " of ", length(genes),
#       " gene symbols to ",
#       length(entrez_ids),
#       " unique Entrez IDs."
#     )

#     if (length(entrez_ids) == 0) {
#       warning("No valid Entrez IDs remain after filtering.")
#       return(NULL)
#     }

#     message(
#       "Example Entrez IDs: ",
#       paste(utils::head(entrez_ids, 10), collapse = ", ")
#     )

#     message(
#       "Running KEGG enrichment via ",
#       "clusterProfiler::enrichKEGG()..."
#     )

#     res <- tryCatch(
#       clusterProfiler::enrichKEGG(
#         gene          = entrez_ids,
#         organism      = organism,
#         keyType       = "kegg",
#         pvalueCutoff  = 1,
#         qvalueCutoff  = 1,
#         pAdjustMethod = "BH"
#       ),
#       error = function(e) {
#         warning(
#           glue::glue(
#             "KEGG enrichment failed: {conditionMessage(e)}. ",
#             "enrichKEGG() accesses rest.kegg.jp and therefore ",
#             "requires internet access."
#           )
#         )
#         NULL
#       }
#     )

#     if (is.null(res)) {
#       return(NULL)
#     }

#     res_df <- as.data.frame(res)

#     if (nrow(res_df) == 0) {
#       warning(
#         "The genes mapped to Entrez IDs, but no KEGG pathways ",
#         "passed the enrichment criteria."
#       )
#       return(NULL)
#     }

#     res_df %>%
#       dplyr::rename(
#         Term             = Description,
#         P.value          = pvalue,
#         Adjusted.P.value = p.adjust
#       ) %>%
#       dplyr::as_tibble()

#   } else if (method == "enrichr") {

#     message(
#       glue::glue(
#         "Running KEGG enrichment via Enrichr ({kegg_db})..."
#       )
#     )

#     raw <- safe_enrichr(
#       genes,
#       kegg_db,
#       context = "KEGG enrichment (Enrichr)"
#     )

#     if (is.null(raw)) {
#       return(NULL)
#     }

#     result <- raw[[kegg_db]]

#     if (is.null(result) || nrow(result) == 0) {
#       warning(
#         glue::glue(
#           "No KEGG enrichment results returned from {kegg_db}."
#         )
#       )
#       return(NULL)
#     }

#     dplyr::as_tibble(result)

#   } else {
#     stop(
#       "Unknown KEGG method: ", method,
#       ". Use 'clusterProfiler' or 'enrichr'."
#     )
#   }
# }

# # ---- Run all requested databases -------------------------------------------

# go_slimmed <- if (!opt$skip_go) {
#   slim_go(genes, opt$go_method, trimws(stringr::str_split(opt$go_dbs, ",")[[1]]), opt$annoDb, opt$slim_threshold)
# } else NULL

# go_top <- if (!is.null(go_slimmed) && nrow(go_slimmed) > 0) {
#   go_slimmed %>%
#     dplyr::group_by(`Gene Ontology`) %>%
#     dplyr::slice_max(order_by = `-log10.p-value`, n = opt$top_n) %>%
#     dplyr::ungroup()
# } else NULL

# kegg_flat <- if (!opt$skip_kegg) {
#   kegg_raw <- run_kegg(genes, opt$kegg_method, opt$organism, opt$kegg_db, opt$annoDb)
#   format_flat(kegg_raw, "KEGG", opt$top_n)
# } else NULL

# gwas_flat <- if (!opt$skip_gwas) {
#   gwas_raw_list <- safe_enrichr(genes, opt$gwas_db, context = "GWAS Catalog enrichment")
#   gwas_raw <- if (!is.null(gwas_raw_list)) gwas_raw_list[[opt$gwas_db]] else NULL
#   if (is.null(gwas_raw)) NULL else format_flat(gwas_raw, "GWAS Catalog", opt$top_n)
# } else NULL

# extra_flat <- NULL
# if (!is.null(opt$extra_dbs)) {
#   extra_dbs <- trimws(stringr::str_split(opt$extra_dbs, ",")[[1]])
#   extra_raw <- safe_enrichr(genes, extra_dbs, context = "extra database enrichment")
#   if (!is.null(extra_raw)) {
#     extra_flat <- purrr::imap_dfr(extra_raw, function(df, db_name) {
#       res <- format_flat(df, db_name, opt$top_n)
#       if (is.null(res)) return(dplyr::tibble())
#       res
#     })
#   }
# }

# # ---- Plotting ---------------------------------------------------------------

# make_barplot <- function(df, title = NULL, fonts = font_cfg) {
#   df <- df %>%
#     dplyr::mutate(
#       Term = stringr::str_trim(Term),
#       Term = Hmisc::capitalize(Term),
#       Term = stringr::str_wrap(Term, 50),
#       Term = factor(Term, levels = unique(Term[order(forcats::fct_rev(`Gene Ontology`), `-log10.p-value`)]))
#     )

#   ggplot2::ggplot(df, ggplot2::aes(x = Term, y = `-log10.p-value`,
#                                     fill = `Gene Ontology`, group = `Gene Ontology`)) +
#     ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "black") +
#     ggplot2::coord_flip() +
#     ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
#     ggsci::scale_fill_d3() +
#     ggplot2::labs(y = expression("-log"[10](p)), title = title) +
#     ggplot2::theme_classic() +
#     ggplot2::theme(
#       axis.text = ggplot2::element_text(size = fonts$axis_text, color = "black"),
#       axis.title.y = ggplot2::element_blank(),
#       axis.title.x = ggplot2::element_text(size = fonts$axis_title),
#       legend.text = ggplot2::element_text(size = fonts$legend_text),
#       legend.title = ggplot2::element_text(size = fonts$legend_title),
#       plot.title = ggplot2::element_text(size = fonts$axis_title + 2, face = "bold", hjust = 0.5)
#     )
# }

# save_plot <- function(p, filename, width, height) {
#   path <- file.path(opt$outdir, filename)
#   ggplot2::ggsave(
#     filename = path,
#     plot = p,
#     device = opt$fig_format,
#     width = width,
#     height = height,
#     units = "in",
#     dpi = if (opt$mode == "presentation") 150 else 300
#   )
#   message(glue::glue("Saved {path}"))
# }

# ext <- opt$fig_format

# all_empty <- all(sapply(list(go_top, kegg_flat, gwas_flat, extra_flat), function(x) is.null(x) || nrow(x) == 0))
# if (all_empty) stop("No significant terms found across any requested database -- nothing to plot")

# if (opt$mode == "presentation") {

#   combined <- dplyr::bind_rows(go_top, kegg_flat, gwas_flat, extra_flat)
#   p <- make_barplot(combined)
#   save_plot(p, glue::glue("{opt$prefix}_overview.{ext}"), width = 11, height = max(4, nrow(combined) * 0.5))

# } else {

#   if (!is.null(go_top) && nrow(go_top) > 0) {
#     p_go <- make_barplot(go_top, title = "Gene ontology")
#     save_plot(p_go, glue::glue("{opt$prefix}_GO.{ext}"), width = 9, height = max(4, nrow(go_top) * 0.35))
#   }

#   if (!is.null(kegg_flat) && nrow(kegg_flat) > 0) {
#     p_kegg <- make_barplot(kegg_flat, title = "KEGG pathways")
#     save_plot(p_kegg, glue::glue("{opt$prefix}_KEGG.{ext}"), width = 9, height = max(3, nrow(kegg_flat) * 0.35))
#   }

#   if (!is.null(gwas_flat) && nrow(gwas_flat) > 0) {
#     p_gwas <- make_barplot(gwas_flat, title = "GWAS Catalog")
#     save_plot(p_gwas, glue::glue("{opt$prefix}_GWAS.{ext}"), width = 9, height = max(3, nrow(gwas_flat) * 0.35))
#   }

#   if (!is.null(extra_flat) && nrow(extra_flat) > 0) {
#     for (db in unique(extra_flat$`Gene Ontology`)) {
#       sub <- extra_flat %>% dplyr::filter(`Gene Ontology` == db)
#       p_extra <- make_barplot(sub, title = db)
#       save_plot(p_extra, glue::glue("{opt$prefix}_{db}.{ext}"), width = 9, height = max(3, nrow(sub) * 0.35))
#     }
#   }
# }

# message("Done.")
