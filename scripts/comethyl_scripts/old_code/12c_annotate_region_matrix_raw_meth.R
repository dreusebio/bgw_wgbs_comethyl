#!/usr/bin/env Rscript

# ============================================================
# 12c_annotate_region_matrix.R
#
# Purpose:
#   Annotate ALL regions in Filtered_Regions.txt with nearest gene
#   symbol, genomic location, and CpG context -- independent of
#   module assignment. Annotation is computed ONCE (it depends only
#   on genomic coordinates, not methylation values) and then merged
#   against up to three separate methylation matrices:
#     - the raw, pre-adjustment matrix (step 03 Region_Methylation.rds)
#     - the v1_all_pcs PC-adjusted matrix (step 05)
#     - the v2_exclude_protected_pcs PC-adjusted matrix (step 05)
#   Any subset of these three can be provided; each gets its own
#   combined output file.
#
# Input:
#   --filtered_regions       Filtered_Regions.txt from step 02
#                             (columns: RegionID, chr, start, end, ...)
#   --methylation_matrix     (optional) raw Region_Methylation.rds (step 03)
#   --methylation_matrix_v1  (optional) Adjusted_Region_Methylation_allPCs.rds (step 05)
#   --methylation_matrix_v2  (optional) Adjusted_Region_Methylation_excluding_protected_PCs_*.rds (step 05)
#   --modules_v1             (optional) Modules.rds for v1_all_pcs (step 07) -- adds a "module" column
#   --modules_v2             (optional) Modules.rds for v2_exclude_protected_pcs (step 07) -- adds a "module" column
#   Each provided matrix is also used to sanity-check that its
#   RegionIDs match 1:1 with filtered_regions before annotating.
#
# Output, under <project_root>/comethyl_output/
#         12c_region_methylation_matrix_annotation/<cpg_label>/<region_label>/:
#   Region_Annotation/                          (computed once, shared across variants)
#     Region_Gene_Annotation.tsv                RegionID, chr, start, end, gene_symbol, ...
#     Region_Gene_Annotation.xlsx
#     run_parameters.txt
#     run_log.txt
#     session_info.txt
#   Region_Annotation_with_Methylation.tsv/.xlsx      (if --methylation_matrix given)
#   v1_all_pcs/
#     Region_Annotation_with_Methylation.tsv/.xlsx    (if --methylation_matrix_v1 given)
#   v2_exclude_protected_pcs/
#     Region_Annotation_with_Methylation.tsv/.xlsx    (if --methylation_matrix_v2 given)
#
# Usage:
#   Rscript 12c_annotate_region_matrix.R \
#     --project_root /path/to/project \
#     --filtered_regions /path/to/Filtered_Regions.txt \
#     --methylation_matrix_v1 /path/to/Adjusted_Region_Methylation_allPCs.rds \
#     --methylation_matrix_v2 /path/to/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#     --modules_v1 /path/to/v1/Modules.rds \
#     --modules_v2 /path/to/v2/Modules.rds \
#     --genome hg38 \
#     --annotation_mode offline \
#     --cpg_island_file /path/to/local/cpgIslandExt.txt.gz \
#     [--methylation_matrix /path/to/Region_Methylation.rds] \  # optional raw matrix too
#     [--add_genomic_location TRUE] \
#     [--txdb_pkg TxDb.Hsapiens.UCSC.hg38.knownGene] \  # override; default auto-resolves from --genome
#     [--write_methylation_xlsx FALSE] \                # TRUE to also write combined tables as .xlsx (large/slow)
#     [--helper_file /path/to/helper.R]   # optional; defaults to
#                                          # helper.R alongside this script
#
# Depends on: helper.R (sourced for annotate_regions_safe() and
# friends -- see the "ANNOTATION HELPERS" section at the bottom of
# helper.R). genomic_location requires ChIPseeker + a local TxDb
# (BiocManager::install('ChIPseeker')). cpg_context requires
# --cpg_island_file (a local UCSC cpgIslandExt table); omit the flag
# to skip that column with no error.
# ============================================================

message("Starting Script 12c: Annotate Region Matrix")

suppressPackageStartupMessages({
  library(comethyl)
  library(dplyr)
  library(openxlsx)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(AnnotationDbi)
})

# ------------------------------------------------------------
# Shared helpers
# ------------------------------------------------------------
# get_arg()/trim_or_null() live inside helper.R itself, so we need a
# minimal standalone version here just to read --helper_file before
# helper.R has been sourced.
.get_arg_bootstrap <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) NULL
)
if (is.null(script_dir) || !nzchar(script_dir)) script_dir <- getwd()

# Default: helper.R lives alongside this script. Override with
# --helper_file /full/path/to/helper.R if it lives elsewhere.
helper_file <- .get_arg_bootstrap("--helper_file", file.path(script_dir, "helper.R"))
if (!file.exists(helper_file)) {
  stop("helper.R not found at: ", helper_file,
       "\nPass --helper_file /path/to/helper.R if it lives elsewhere.", call. = FALSE)
}
source(helper_file)

# ------------------------------------------------------------
# Read arguments
# ------------------------------------------------------------
project_root       <- trim_or_null(get_arg("--project_root"))
filtered_regions_f <- trim_or_null(get_arg("--filtered_regions"))
meth_matrix_f      <- trim_or_null(get_arg("--methylation_matrix"))            # optional: raw (step 03) matrix
meth_matrix_v1_f   <- trim_or_null(get_arg("--methylation_matrix_v1"))         # optional: step 05 all_pcs
meth_matrix_v2_f   <- trim_or_null(get_arg("--methylation_matrix_v2"))         # optional: step 05 exclude_protected_pcs
modules_v1_f       <- trim_or_null(get_arg("--modules_v1"))                    # optional: step 07 Modules.rds (all_pcs)
modules_v2_f       <- trim_or_null(get_arg("--modules_v2"))                    # optional: step 07 Modules.rds (exclude_protected_pcs)
genome             <- trim_or_null(get_arg("--genome", "hg38"))
annotation_mode    <- trim_or_null(get_arg("--annotation_mode", "auto"))
cpg_island_file    <- trim_or_null(get_arg("--cpg_island_file"))
add_genomic_location <- parse_bool(get_arg("--add_genomic_location", "TRUE"), "--add_genomic_location")
txdb_pkg           <- trim_or_null(get_arg("--txdb_pkg"))  # override; default NULL = auto-resolve from --genome
write_methylation_xlsx <- parse_bool(get_arg("--write_methylation_xlsx", "FALSE"), "--write_methylation_xlsx")
exclude_grey_module <- parse_bool(get_arg("--exclude_grey_module", "TRUE"), "--exclude_grey_module")

stop_if_missing(project_root, "--project_root")
stop_if_missing(filtered_regions_f, "--filtered_regions")

if (!dir.exists(project_root)) stop("project_root not found: ", project_root, call. = FALSE)
validate_file_exists(filtered_regions_f, "filtered_regions")
if (!is.null(meth_matrix_f))    validate_file_exists(meth_matrix_f, "methylation_matrix")
if (!is.null(meth_matrix_v1_f)) validate_file_exists(meth_matrix_v1_f, "methylation_matrix_v1")
if (!is.null(meth_matrix_v2_f)) validate_file_exists(meth_matrix_v2_f, "methylation_matrix_v2")
if (!is.null(modules_v1_f))     validate_file_exists(modules_v1_f, "modules_v1")
if (!is.null(modules_v2_f))     validate_file_exists(modules_v2_f, "modules_v2")
if (!is.null(cpg_island_file))  validate_file_exists(cpg_island_file, "cpg_island_file")

annotation_mode <- match.arg(annotation_mode, choices = c("auto", "great", "offline"))

setwd(project_root)

AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(project_root, ".cache")
)

# ------------------------------------------------------------
# Derive output dir from filtered_regions path
#   .../02_filter_regions/<cpg_label>/<region_label>/Filtered_Regions.txt
#   -> comethyl_output/12c_region_methylation_matrix_annotation/<cpg_label>/<region_label>/
#        Region_Annotation/                    (annotation only -- shared, computed once)
#        v1_all_pcs/                           (annotation + v1 sample methylation, if provided)
#        v2_exclude_protected_pcs/             (annotation + v2 sample methylation, if provided)
#   We mirror whatever <cpg_label>/<region_label> we can find; if the
#   path doesn't match that shape, fall back to a flat output dir.
# ------------------------------------------------------------
derive_base_dir <- function(filtered_regions_path, project_root) {
  region_dir <- dirname(filtered_regions_path)      # .../covMin4_methSD0p08 (or similar)
  cpg_dir    <- dirname(region_dir)                 # .../cov3_75pct
  region_label <- basename(region_dir)
  cpg_label    <- basename(cpg_dir)

  base_dir <- file.path(
    project_root, "comethyl_output", "12c_region_methylation_matrix_annotation",
    cpg_label, region_label
  )
  safe_dir_create(base_dir)
  base_dir
}

base_dir <- derive_base_dir(filtered_regions_f, project_root)
out_dir  <- file.path(base_dir, "Region_Annotation")
safe_dir_create(out_dir)

log_file    <- file.path(out_dir, "run_log.txt")
params_file <- file.path(out_dir, "run_parameters.txt")


append_log(log_file, "project_root: ", project_root)
append_log(log_file, "filtered_regions: ", filtered_regions_f)
append_log(log_file, "methylation_matrix (raw): ", ifelse(is.null(meth_matrix_f), "(not provided)", meth_matrix_f))
append_log(log_file, "methylation_matrix_v1: ", ifelse(is.null(meth_matrix_v1_f), "(not provided)", meth_matrix_v1_f))
append_log(log_file, "methylation_matrix_v2: ", ifelse(is.null(meth_matrix_v2_f), "(not provided)", meth_matrix_v2_f))
append_log(log_file, "modules_v1: ", ifelse(is.null(modules_v1_f), "(not provided)", modules_v1_f))
append_log(log_file, "modules_v2: ", ifelse(is.null(modules_v2_f), "(not provided)", modules_v2_f))
append_log(log_file, "genome: ", genome)
append_log(log_file, "annotation_mode: ", annotation_mode)
append_log(log_file, "add_genomic_location: ", add_genomic_location)
append_log(log_file, "txdb_pkg: ", ifelse(is.null(txdb_pkg), paste0("(auto-resolved from genome=", genome, ")"), txdb_pkg))
append_log(log_file, "cpg_island_file: ", ifelse(is.null(cpg_island_file), "(not provided)", cpg_island_file))
append_log(log_file, "base_dir: ", base_dir)
append_log(log_file, "out_dir: ", out_dir)

# ------------------------------------------------------------
# Load filtered regions (has RegionID, chr, start, end, ...)
# ------------------------------------------------------------
filtered_regions <- read.delim(filtered_regions_f, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("RegionID", "chr", "start", "end")
missing_cols <- setdiff(required_cols, colnames(filtered_regions))
if (length(missing_cols) > 0) {
  stop("filtered_regions is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

append_log(log_file, "Loaded filtered_regions rows: ", nrow(filtered_regions))

# ------------------------------------------------------------
# Sanity check: RegionIDs match the methylation matrix rows.
# Reused for the raw matrix and each PC-adjusted variant.
#
# Matrices from different pipeline steps aren't guaranteed to share
# the same orientation -- the raw step-03 matrix has regions in rows,
# but PC-adjustment tools (step 05) commonly transpose to
# samples-in-rows for their own internal processing and leave it that
# way on output. Auto-detect and correct rather than assuming.
# ------------------------------------------------------------
orient_regions_in_rows <- function(mat, filtered_regions, label, log_file) {
  region_ids <- filtered_regions$RegionID
  rn_match <- sum(rownames(mat) %in% region_ids)
  cn_match <- sum(colnames(mat) %in% region_ids)

  if (rn_match >= cn_match) {
    return(mat)  # already regions-in-rows (or ambiguous -- leave as-is;
                 # the mismatch warnings below will still flag it)
  }

  append_log(log_file, "NOTE [", label, "]: matrix looks like samples-in-rows/regions-in-columns (",
             cn_match, " of ", ncol(mat), " column names match RegionIDs vs ",
             rn_match, " of ", nrow(mat), " row names). Transposing to regions-in-rows.")
  t(mat)
}

check_matrix_regionids <- function(mat_path, label, filtered_regions, log_file) {
  if (is.null(mat_path)) return(invisible(NULL))

  mat <- readRDS(mat_path)
  mat <- orient_regions_in_rows(mat, filtered_regions, label, log_file)
  mat_ids <- rownames(mat)

  n_mat <- length(mat_ids)
  n_filt <- nrow(filtered_regions)

  if (n_mat != n_filt) {
    append_log(log_file, "WARNING [", label, "]: row count mismatch -- matrix has ", n_mat,
               " regions, filtered_regions has ", n_filt)
  }

  missing_from_filtered <- setdiff(mat_ids, filtered_regions$RegionID)
  if (length(missing_from_filtered) > 0) {
    append_log(log_file, "WARNING [", label, "]: ", length(missing_from_filtered),
               " RegionIDs are absent from filtered_regions (e.g. ",
               paste(head(missing_from_filtered, 5), collapse = ", "), ")")
  } else {
    append_log(log_file, "OK [", label, "]: all RegionIDs found in filtered_regions")
  }

  mat
}

meth_mat    <- check_matrix_regionids(meth_matrix_f, "raw", filtered_regions, log_file)
meth_mat_v1 <- check_matrix_regionids(meth_matrix_v1_f, "v1_all_pcs", filtered_regions, log_file)
meth_mat_v2 <- check_matrix_regionids(meth_matrix_v2_f, "v2_exclude_protected_pcs", filtered_regions, log_file)

# ------------------------------------------------------------
# Build regions_df for annotate_regions_safe()
#   No "module" column exists yet at this stage, so we relax the
#   required-columns check via req_cols.
# ------------------------------------------------------------
regions_df <- filtered_regions %>%
  dplyr::select(RegionID, chr, start, end)

annotated_regions <- suppressWarnings(
  annotate_regions_safe(
    regions_df = regions_df,
    genome = genome,
    annotation_mode = annotation_mode,
    file_txt = NULL,
    verbose = TRUE,
    req_cols = c("RegionID", "chr", "start", "end"),
    add_genomic_location = add_genomic_location,
    cpg_island_file = cpg_island_file,
    txdb_pkg = txdb_pkg
  )
)

append_log(log_file, "Annotated regions rows: ", nrow(annotated_regions))
append_log(log_file, "Regions with a gene_symbol: ",
           sum(!is.na(annotated_regions$gene_symbol) & annotated_regions$gene_symbol != ""))
if ("genomic_location" %in% colnames(annotated_regions)) {
  append_log(log_file, "Regions with genomic_location: ",
             sum(!is.na(annotated_regions$genomic_location)))
}
if ("cpg_context" %in% colnames(annotated_regions)) {
  append_log(log_file, "Regions with cpg_context: ",
             sum(!is.na(annotated_regions$cpg_context)))
}

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------
out_tsv  <- file.path(out_dir, "Region_Gene_Annotation.tsv")
out_xlsx <- file.path(out_dir, "Region_Gene_Annotation.xlsx")

write.table(annotated_regions, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
openxlsx::write.xlsx(annotated_regions, out_xlsx, rowNames = FALSE)

params <- c(
  paste0("timestamp\t", timestamp_now()),
  paste0("project_root\t", project_root),
  paste0("filtered_regions\t", filtered_regions_f),
  paste0("methylation_matrix_raw\t", ifelse(is.null(meth_matrix_f), "", meth_matrix_f)),
  paste0("methylation_matrix_v1\t", ifelse(is.null(meth_matrix_v1_f), "", meth_matrix_v1_f)),
  paste0("methylation_matrix_v2\t", ifelse(is.null(meth_matrix_v2_f), "", meth_matrix_v2_f)),
  paste0("modules_v1\t", ifelse(is.null(modules_v1_f), "", modules_v1_f)),
  paste0("modules_v2\t", ifelse(is.null(modules_v2_f), "", modules_v2_f)),
  paste0("genome\t", genome),
  paste0("annotation_mode\t", annotation_mode),
  paste0("add_genomic_location\t", add_genomic_location),
  paste0("cpg_island_file\t", ifelse(is.null(cpg_island_file), "", cpg_island_file)),
  paste0("write_methylation_xlsx\t", write_methylation_xlsx),
  paste0("exclude_grey_module\t", exclude_grey_module),
  paste0("base_dir\t", base_dir),
  paste0("out_dir\t", out_dir),
  paste0("n_regions\t", nrow(filtered_regions)),
  paste0("n_annotated\t", nrow(annotated_regions))
)
write_lines_safe(params, params_file)

append_log(log_file, "Wrote: ", out_tsv)
append_log(log_file, "Wrote: ", out_xlsx)

# ------------------------------------------------------------
# Merge sample-level methylation values into the (shared) annotation
# table and write one combined file per matrix provided:
#   meth_mat    -> base_dir/Region_Annotation_with_Methylation.*        (raw, step 03)
#   meth_mat_v1 -> base_dir/v1_all_pcs/Region_Annotation_with_Methylation.*
#   meth_mat_v2 -> base_dir/v2_exclude_protected_pcs/Region_Annotation_with_Methylation.*
# ------------------------------------------------------------
# ------------------------------------------------------------
# Load RegionID -> module assignment from a Modules.rds (step 07),
# for merging into the per-variant output. Only makes sense per
# variant -- module assignment doesn't exist for the raw/pre-PC-
# adjustment matrix, since modules are detected downstream of that.
# ------------------------------------------------------------
load_module_assignment <- function(modules_path, label, filtered_regions, log_file) {
  if (is.null(modules_path)) return(invisible(NULL))

  obj <- readRDS(modules_path)
  regions <- extract_regions_from_module_object(obj, label = modules_path)

  mod_df <- regions %>% dplyr::select(RegionID, module) %>% dplyr::distinct(RegionID, .keep_all = TRUE)

  n_matched <- sum(filtered_regions$RegionID %in% mod_df$RegionID)
  append_log(log_file, label, ": ", nrow(mod_df), " regions have a module assignment; ",
             n_matched, " of ", nrow(filtered_regions),
             " filtered_regions rows have a matching module (rest = NA, e.g. grey/unassigned or excluded during module detection)")

  mod_df
}

modules_v1 <- load_module_assignment(modules_v1_f, "modules_v1", filtered_regions, log_file)
modules_v2 <- load_module_assignment(modules_v2_f, "modules_v2", filtered_regions, log_file)

# ------------------------------------------------------------
# Merge sample-level methylation values (and, for v1/v2, module
# assignment) into the (shared) annotation table and write one
# combined file per matrix provided:
#   meth_mat    -> base_dir/Region_Annotation_with_Methylation.*        (raw, step 03; no module)
#   meth_mat_v1 -> base_dir/v1_all_pcs/Region_Annotation_with_Methylation.*        (+ module)
#   meth_mat_v2 -> base_dir/v2_exclude_protected_pcs/Region_Annotation_with_Methylation.*  (+ module)
# ------------------------------------------------------------
write_annotated_methylation <- function(annotated_regions, mat, variant_dir, label,
                                        write_xlsx, log_file, mod_df = NULL,
                                        exclude_grey = TRUE) {
  if (is.null(mat)) {
    append_log(log_file, "Skipping ", label, " (matrix not provided)")
    return(invisible(NULL))
  }

  safe_dir_create(variant_dir)

  merged <- annotated_regions

  if (!is.null(mod_df)) {
    merged <- merged %>% dplyr::left_join(mod_df, by = "RegionID")
    # keep module right after the core annotation columns, before sample columns
    merged <- merged %>% dplyr::relocate(module, .after = dplyr::all_of(if ("cpg_context" %in% names(merged)) "cpg_context" else "gene_symbol"))

    if (isTRUE(exclude_grey)) {
      n_before <- nrow(merged)
      # exact match only -- "grey60" etc. are real named modules, not the
      # WGCNA unassigned bucket, and must NOT be caught by this filter
      merged <- merged %>% dplyr::filter(is.na(module) | module != "grey")
      n_dropped <- n_before - nrow(merged)
      append_log(log_file, label, ": excluded ", n_dropped, " grey-module regions (",
                 nrow(merged), " of ", n_before, " remain)")
    }
  }

  meth_df <- as.data.frame(mat, check.names = FALSE)
  meth_df$RegionID <- rownames(mat)
  meth_df <- meth_df[, c("RegionID", setdiff(colnames(meth_df), "RegionID"))]

  merged <- merged %>% dplyr::left_join(meth_df, by = "RegionID")

  first_sample_col <- colnames(mat)[1]
  n_matched <- sum(!is.na(merged[[first_sample_col]]))
  append_log(log_file, label, ": regions matched to sample methylation values: ",
             n_matched, " of ", nrow(merged))

  out_tsv <- file.path(variant_dir, "Region_Annotation_with_Methylation.tsv")
  write.table(merged, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  append_log(log_file, "Wrote: ", out_tsv)

  if (isTRUE(write_xlsx)) {
    out_xlsx_meth <- file.path(variant_dir, "Region_Annotation_with_Methylation.xlsx")
    openxlsx::write.xlsx(merged, out_xlsx_meth, rowNames = FALSE)
    append_log(log_file, "Wrote: ", out_xlsx_meth)
  } else {
    append_log(log_file, "Skipped ", label, " .xlsx (", ncol(merged), " columns x ", nrow(merged),
               " rows -- large and slow as .xlsx). Use the .tsv, or pass ",
               "--write_methylation_xlsx TRUE to force it.")
  }

  invisible(merged)
}

write_annotated_methylation(annotated_regions, meth_mat, base_dir,
                            "raw", write_methylation_xlsx, log_file)
write_annotated_methylation(annotated_regions, meth_mat_v1, file.path(base_dir, "v1_all_pcs"),
                            "v1_all_pcs", write_methylation_xlsx, log_file, mod_df = modules_v1,
                            exclude_grey = exclude_grey_module)
write_annotated_methylation(annotated_regions, meth_mat_v2, file.path(base_dir, "v2_exclude_protected_pcs"),
                            "v2_exclude_protected_pcs", write_methylation_xlsx, log_file, mod_df = modules_v2,
                            exclude_grey = exclude_grey_module)

session_info_file <- file.path(out_dir, "session_info.txt")
write_session_info(
  session_info_file,
  extra_files = list(
    filtered_regions = filtered_regions_f,
    methylation_matrix_raw = meth_matrix_f,
    methylation_matrix_v1 = meth_matrix_v1_f,
    methylation_matrix_v2 = meth_matrix_v2_f,
    modules_v1 = modules_v1_f,
    modules_v2 = modules_v2_f,
    cpg_island_file = cpg_island_file
  )
)
append_log(log_file, "Wrote: ", session_info_file)

message("Script 12c complete: Region matrix annotation finished")