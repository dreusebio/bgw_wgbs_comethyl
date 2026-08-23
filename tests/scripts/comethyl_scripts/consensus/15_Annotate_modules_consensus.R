#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 15: Consensus Module Annotation
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   Annotates all consensus regions (including grey module) with
#   nearest gene information using GREAT (via comethyl) or offline
#   nearest-gene annotation (EnsDb.Hsapiens.v86 or TxDb fallback):
#     1) Load Consensus_Region_Assignments.tsv from script 09 shared/
#     2) Annotate all regions with nearest gene symbols and IDs
#     3) Build per-module gene lists
#     4) Build module-level region and gene count summaries
#     5) Build background gene universe for enrichment analysis
#
# NOTE
#   Module assignments are shared across datasets — annotation runs
#   once per adjustment version. No per-dataset loop needed.
#
# REQUIRED INPUTS
#   --project_root          : root directory of the analysis project
#   --consensus_regions_tsv : path to Consensus_Region_Assignments.tsv
#                             from script 09 shared/
#
# OPTIONAL INPUTS
#   --adjustment_version : label matching script 09 run [default = unadjusted]
#   --genome             : genome assembly [default = hg38]
#   --annotation_mode    : auto | great | offline [default = auto]
#                          auto tries GREAT first, falls back to offline
#
# OUTPUTS
#   comethyl_output/consensus/15_module_annotation/<cpg_label>/<region_label>/<adjustment_version>/
#       Annotated_Regions.tsv
#       Annotated_Regions.xlsx
#       Module_Gene_List.tsv
#       Module_Gene_Summary.tsv
#       Module_Gene_Summary.xlsx
#       Background_Genes.txt
#       Background_Genes.tsv
#       run_parameters.txt
#       run_log.txt
#       sessionInfo.txt
#
# NOTES
#   - cpg_label and region_label are derived automatically from the
#     --consensus_regions_tsv path (walking up through shared/ and
#     the adjustment_version/ directories)
#   - The annotation outputs feed directly into script 16 (--annotation_dir)
#   - offline mode requires EnsDb.Hsapiens.v86 or TxDb.Hsapiens.UCSC.hg38.knownGene
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/15_Annotate_modules_consensus.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --consensus_regions_tsv /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/shared/Consensus_Region_Assignments.tsv \
#     --adjustment_version v1_all_pcs \
#     --genome hg38 \
#     --annotation_mode auto
# ================================================================
# SCRIPT 15: Consensus Module Annotation
#
# PURPOSE
#   - Load consensus region assignments from script 09 shared/
#   - Annotate all regions (all modules including grey)
#   - Write annotated region tables
#   - Write per-module gene lists
#   - Write module-level summary tables
#   - Write background/testable gene universe
#   - Write run logs
#
# NOTE
#   Module assignments are shared across datasets in a consensus
#   analysis — annotation runs once per adjustment version only.
#   No per-dataset loop needed.
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/15_module_annotation/
#     <cpg_label>/
#       <region_label>/
#         <adjustment_version>/
#           Annotated_Regions.tsv
#           Annotated_Regions.xlsx
#           Module_Gene_List.tsv
#           Module_Gene_Summary.tsv
#           Module_Gene_Summary.xlsx
#           Background_Genes.txt
#           Background_Genes.tsv
#           run_parameters.txt
#           run_log.txt
#           sessionInfo.txt
#
# PATH DERIVATION
#   cpg_label and region_label are derived automatically from
#   --consensus_regions_tsv path:
#     .../09_consensus_modules/<cpg_label>/<region_label>/<adjustment_version>/shared/
#         Consensus_Region_Assignments.tsv
# ================================================================

message("Starting Script 15 ✓")

suppressPackageStartupMessages({
  library(comethyl)
  library(dplyr)
  library(openxlsx)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(AnnotationDbi)
  library(stringr)
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

stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument: ", label, call. = FALSE)
}

validate_file_exists <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

validate_regions_df <- function(regions, label = "regions") {
  req     <- c("RegionID", "chr", "start", "end", "module")
  missing <- setdiff(req, colnames(regions))
  if (length(missing) > 0)
    stop(label, " missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

# ================================================================
# Annotation helpers (matching regular 12a style)
# ================================================================
offline_nearest_gene <- function(gr, verbose = TRUE) {
  have_ensdb <- requireNamespace("EnsDb.Hsapiens.v86",                  quietly = TRUE)
  have_txdb  <- requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene",  quietly = TRUE)
  have_org   <- requireNamespace("org.Hs.eg.db",                        quietly = TRUE)

  if (!have_org)
    stop("org.Hs.eg.db is required for offline annotation.", call. = FALSE)

  if (have_ensdb) {
    if (verbose) message("[offline] Using EnsDb.Hsapiens.v86")
    edb  <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    seqs <- unique(as.character(GenomeInfoDb::seqnames(gr)))
    seqs <- seqs[seqs %in% GenomeInfoDb::seqlevels(edb)]
    genes <- ensembldb::genes(edb, filter = AnnotationFilter::SeqNameFilter(seqs))
    ggr   <- GenomicRanges::GRanges(genes)
    map   <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = genes$gene_id,
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID")
    )
    ggr$SYMBOL   <- map$SYMBOL[match(genes$gene_id,   map$ENSEMBL)]
    ggr$ENTREZID <- map$ENTREZID[match(genes$gene_id, map$ENSEMBL)]

  } else if (have_txdb) {
    if (verbose) message("[offline] Using TxDb.Hsapiens.UCSC.hg38.knownGene")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    ggr  <- GenomicFeatures::genes(txdb)
    map  <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = ggr$gene_id,
      keytype = "ENTREZID",
      columns = "SYMBOL"
    )
    ggr$ENTREZID <- ggr$gene_id
    ggr$SYMBOL   <- map$SYMBOL[match(ggr$ENTREZID, map$ENTREZID)]

  } else {
    stop("Install EnsDb.Hsapiens.v86 or TxDb.Hsapiens.UCSC.hg38.knownGene", call. = FALSE)
  }

  hit <- GenomicRanges::distanceToNearest(gr, ggr, ignore.strand = TRUE)
  ng  <- ggr[S4Vectors::subjectHits(hit)]

  data.frame(
    chr           = as.character(GenomeInfoDb::seqnames(gr))[S4Vectors::queryHits(hit)],
    start         = as.integer(S4Vectors::start(gr))[S4Vectors::queryHits(hit)],
    end           = as.integer(S4Vectors::end(gr))[S4Vectors::queryHits(hit)],
    gene_symbol   = as.character(ng$SYMBOL),
    gene_entrezID = as.character(ng$ENTREZID),
    stringsAsFactors = FALSE
  )
}

annotate_offline_only <- function(regions_df, genome = "hg38",
                                   file_txt = NULL, verbose = TRUE) {
  validate_regions_df(regions_df)
  gr <- GenomicRanges::GRanges(
    seqnames = regions_df$chr,
    ranges   = IRanges::IRanges(start = regions_df$start, end = regions_df$end),
    RegionID = regions_df$RegionID
  )
  ng  <- offline_nearest_gene(gr, verbose = verbose)
  out <- regions_df %>%
    dplyr::left_join(
      ng %>% dplyr::select(chr, start, end, gene_symbol, gene_entrezID),
      by = c("chr", "start", "end")
    ) %>%
    dplyr::mutate(
      gene_description = NA_character_,
      gene_ensemblID   = NA_character_
    )
  if (!is.null(file_txt) && nzchar(file_txt))
    write.table(out, file = file_txt, sep = "\t", quote = FALSE, row.names = FALSE)
  out
}

annotate_regions_safe <- function(regions_df,
                                   genome          = "hg38",
                                   annotation_mode = c("auto", "great", "offline"),
                                   file_txt        = NULL,
                                   verbose         = TRUE) {
  annotation_mode <- match.arg(annotation_mode)

  if (annotation_mode == "offline")
    return(annotate_offline_only(regions_df, genome = genome,
                                  file_txt = file_txt, verbose = verbose))

  if (annotation_mode == "great") {
    if (verbose) message("[annotate] Using comethyl::annotateModule()")
    return(comethyl::annotateModule(regions_df, genome = genome, file = file_txt))
  }

  # auto: try GREAT, fall back to offline
  tryCatch({
    if (verbose) message("[annotate] Trying comethyl::annotateModule()")
    comethyl::annotateModule(regions_df, genome = genome, file = file_txt)
  }, error = function(e) {
    message("[annotate] GREAT failed: ", conditionMessage(e))
    message("[annotate] Falling back to offline nearest-gene annotation.")
    annotate_offline_only(regions_df, genome = genome,
                           file_txt = file_txt, verbose = verbose)
  })
}

# ================================================================
# Summary table helpers (matching regular 12a)
# ================================================================
make_module_gene_list_table <- function(annotated_regions) {
  annotated_regions %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::select(module, gene_symbol) %>%
    dplyr::distinct() %>%
    dplyr::arrange(module, gene_symbol)
}

make_module_summary_table <- function(annotated_regions) {
  gene_col <- "gene_symbol" %in% colnames(annotated_regions)
  annotated_regions %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      n_regions      = dplyr::n(),
      n_unique_genes = if (gene_col)
        length(unique(na.omit(gene_symbol[gene_symbol != ""])))
        else NA_integer_,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_regions), module)
}

make_background_gene_table <- function(annotated_regions) {
  if (!"gene_symbol" %in% colnames(annotated_regions))
    stop("annotated_regions must contain gene_symbol column.", call. = FALSE)
  tibble::tibble(
    gene_symbol = sort(unique(na.omit(as.character(annotated_regions$gene_symbol))))
  ) %>%
    dplyr::filter(gene_symbol != "")
}

# ================================================================
# Read arguments
# ================================================================
project_root          <- trim_or_null(get_arg("--project_root"))
consensus_regions_tsv <- trim_or_null(get_arg("--consensus_regions_tsv"))
adjustment_version    <- trim_or_null(get_arg("--adjustment_version", "unadjusted"))
genome                <- trim_or_null(get_arg("--genome",             "hg38"))
annotation_mode       <- trim_or_null(get_arg("--annotation_mode",    "auto"))

stop_if_missing(project_root,          "--project_root")
stop_if_missing(consensus_regions_tsv, "--consensus_regions_tsv")

if (!dir.exists(project_root))
  stop("project_root not found: ", project_root, call. = FALSE)
validate_file_exists(consensus_regions_tsv, "--consensus_regions_tsv")

annotation_mode <- match.arg(annotation_mode, choices = c("auto", "great", "offline"))

# ================================================================
# Derive cpg_label / region_label from consensus_regions_tsv path
#
# Expected path (Option B structure):
#   .../09_consensus_modules/shared/<cpg_label>/<region_label>/<adjustment_version>/
#       Consensus_Region_Assignments.tsv
# Walking up from the file:
#   dirname        = <adjustment_version>/
#   dirname x2     = <region_label>/
#   dirname x3     = <cpg_label>/
# ================================================================
{
  adj_dir      <- dirname(consensus_regions_tsv)            # <adjustment_version>/
  region_label <- basename(dirname(adj_dir))                # <region_label>
  cpg_label    <- basename(dirname(dirname(adj_dir)))       # <cpg_label>
  message("Derived cpg_label    : ", cpg_label)
  message("Derived region_label : ", region_label)
}

# ================================================================
# Output directory
# ================================================================
pipeline_root <- file.path(project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "15_module_annotation",
                            cpg_label, region_label, adjustment_version)
safe_dir_create(step_dir)
message("Output directory: ", step_dir)

log_file    <- file.path(step_dir, "run_log.txt")
params_file <- file.path(step_dir, "run_parameters.txt")

# ================================================================
# Configure cache
# ================================================================
AnnotationHub::setAnnotationHubOption("CACHE",
  file.path(project_root, ".cache"))

# ================================================================
# Load consensus regions
# ================================================================
append_log(log_file, "Loading consensus region assignments: ", consensus_regions_tsv)

regions <- read.delim(consensus_regions_tsv,
                       stringsAsFactors = FALSE, check.names = FALSE)
validate_regions_df(regions, "--consensus_regions_tsv")

append_log(log_file, "Regions loaded: ", nrow(regions))
append_log(log_file, "Unique modules: ",
           length(unique(as.character(regions$module))))
append_log(log_file, "Non-grey regions: ", sum(regions$module != "grey"))

module_colors_real <- unique(regions$module[regions$module != "grey"])
message("Non-grey modules: ", length(module_colors_real),
        " (", paste(sort(module_colors_real), collapse = ", "), ")")

# ================================================================
# Annotate
# ================================================================
append_log(log_file, "Running annotation (mode = ", annotation_mode, ")")

annotated_tsv  <- file.path(step_dir, "Annotated_Regions.tsv")
annotated_xlsx <- file.path(step_dir, "Annotated_Regions.xlsx")

annotated_regions <- suppressWarnings(
  annotate_regions_safe(
    regions_df      = regions,
    genome          = genome,
    annotation_mode = annotation_mode,
    file_txt        = annotated_tsv,
    verbose         = TRUE
  )
)

append_log(log_file, "Annotated regions: ", nrow(annotated_regions))
openxlsx::write.xlsx(annotated_regions, annotated_xlsx, rowNames = FALSE)
message("Saved: Annotated_Regions.tsv / .xlsx")

# ================================================================
# Gene list and summary tables
# ================================================================
module_gene_list <- make_module_gene_list_table(annotated_regions)
write.table(module_gene_list,
            file.path(step_dir, "Module_Gene_List.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: Module_Gene_List.tsv")

module_summary <- make_module_summary_table(annotated_regions)
write.table(module_summary,
            file.path(step_dir, "Module_Gene_Summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
openxlsx::write.xlsx(module_summary,
                     file.path(step_dir, "Module_Gene_Summary.xlsx"),
                     rowNames = FALSE)
message("Saved: Module_Gene_Summary.tsv / .xlsx")

# ================================================================
# Background gene universe
# ================================================================
background_genes <- make_background_gene_table(annotated_regions)
write_lines_safe(background_genes$gene_symbol,
                 file.path(step_dir, "Background_Genes.txt"))
write.table(background_genes,
            file.path(step_dir, "Background_Genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: Background_Genes.txt / .tsv (", nrow(background_genes), " genes)")

append_log(log_file, "Module gene list rows: ", nrow(module_gene_list))
append_log(log_file, "Module summary rows: ",   nrow(module_summary))
append_log(log_file, "Background genes: ",      nrow(background_genes))

# ================================================================
# Logs
# ================================================================
write_lines_safe(
  c(
    paste0("timestamp\t",              timestamp_now()),
    paste0("project_root\t",          project_root),
    paste0("consensus_regions_tsv\t", consensus_regions_tsv),
    paste0("adjustment_version\t",    adjustment_version),
    paste0("cpg_label\t",             cpg_label),
    paste0("region_label\t",          region_label),
    paste0("genome\t",                genome),
    paste0("annotation_mode\t",       annotation_mode),
    paste0("n_regions_total\t",       nrow(regions)),
    paste0("n_regions_non_grey\t",    sum(regions$module != "grey")),
    paste0("n_modules_real\t",        length(module_colors_real)),
    paste0("n_annotated_rows\t",      nrow(annotated_regions)),
    paste0("n_background_genes\t",    nrow(background_genes)),
    paste0("output_dir\t",            step_dir)
  ),
  params_file
)

writeLines(capture.output(sessionInfo()),
           con = file.path(step_dir, "sessionInfo.txt"))

append_log(log_file, "Script 15 complete")
message("Script 15 complete: consensus module annotation finished")
message("Output: ", step_dir)