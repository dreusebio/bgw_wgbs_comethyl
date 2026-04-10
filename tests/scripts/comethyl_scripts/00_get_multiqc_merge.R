#!/usr/bin/env Rscript

# ============================================================
# MultiQC (text outputs) -> merged QC table
# depending with files in the multiqc, you can adjust these or rename
#    - req_files <- c(
#   "mqc_bismark_alignment_1.txt",
#   "multiqc_bismark_dedup.txt",
#   "multiqc_bismark_methextract.txt"
# bam2nuc is Optional if it exists
# Optional:
#   - Merge to external trait XLSX via --trait_path
# Join key:
#   - Sample_name (derived from sample id / leading digits)
#
# Usage:
#   Rscript 00_get_multiqc_merge.R --multiqc_data /path/to/multiqc_data
#   Rscript 00_get_multiqc_merge.R \
#     --multiqc_data /path/to/multiqc_data \
#     --trait_path /path/to/FINAL_DNAm_Data.xlsx \
#     --out /path/to/merged_qc.xlsx
# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/scripts/comethyl_scripts/00_get_multiqc_merge.R \
#   --multiqc_data /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/data/processed/09_multiqc/merged_multiqc_data \
#   --trait_path /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/data/metadata/FINAL_DNAm_Data_04082026.xlsx \
#   --out /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/data/metadata/merged_qc.xlsx
# # ============================================================

message("Starting Script 00")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(openxlsx)
  library(tibble)
})

# ----------------------------
# CLI args
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste0("Missing value after ", flag))
  args[idx + 1]
}

multiqc_data <- get_arg("--multiqc_data", default = NA_character_)
trait_path   <- get_arg("--trait_path",   default = NA_character_)
out_path     <- get_arg("--out",          default = "merged_multiqc_selected.xlsx")

# ----------------------------
# Resolve multiqc data directory
# ----------------------------
if (!is.na(multiqc_data) && nzchar(multiqc_data)) {
  if (!dir.exists(multiqc_data))
    stop("--multiqc_data directory not found: ", multiqc_data)
  data_dir <- normalizePath(multiqc_data)
  message("multiqc_data: ", data_dir)
} else {
  # Fall back to current working directory (original behaviour)
  data_dir <- getwd()
  message("--multiqc_data not provided; using current directory: ", data_dir)
}

# Helper to build full path into data_dir
data_path <- function(filename) file.path(data_dir, filename)

# ----------------------------
# Helpers
# ----------------------------
read_multiqc_txt <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    as.data.frame()
}

# "201_S1_L004_..." -> "201"; "201_merged_name_sorted" -> "201"
get_core <- function(x) str_extract(x, "^[0-9]+[A-Za-z]*")

# ----------------------------
# Required files
# ----------------------------
req_files <- c(
  "mqc_bismark_alignment_1.txt",
  "multiqc_bismark_dedup.txt",
  "multiqc_bismark_methextract.txt"
)

missing_files <- req_files[!file.exists(data_path(req_files))]
if (length(missing_files) > 0) {
  stop(
    "Missing required file(s) in ", data_dir, ":\n  - ",
    paste(missing_files, collapse = "\n  - ")
  )
}

# ============================================================
# 1) Alignment (lane-level) -> collapse to Sample_name
# ============================================================
align_raw <- read_multiqc_txt(data_path("mqc_bismark_alignment_1.txt"))

needed_align <- c("Sample", "Aligned Uniquely", "Aligned Ambiguously", "Did Not Align", "No Genomic Sequence")
miss_align <- setdiff(needed_align, names(align_raw))
if (length(miss_align) > 0) {
  stop("mqc_bismark_alignment_1.txt missing columns: ", paste(miss_align, collapse = ", "),
       "\nFound: ", paste(names(align_raw), collapse = ", "))
}

align_collapsed <- align_raw %>%
  mutate(
    Sample_name = get_core(Sample),
    SequencingNumber = str_extract(Sample, "(?<=_S)\\d+"),
    SequencingNumber = suppressWarnings(as.numeric(SequencingNumber))
  ) %>%
  filter(!is.na(Sample_name)) %>%
  group_by(Sample_name) %>%
  summarise(
    SequencingNumber = first(SequencingNumber),
    Aligned_Uniquely    = sum(`Aligned Uniquely`, na.rm = TRUE),
    Aligned_Ambiguously = sum(`Aligned Ambiguously`, na.rm = TRUE),
    Did_Not_Align       = sum(`Did Not Align`, na.rm = TRUE),
    No_Genomic_Sequence = sum(`No Genomic Sequence`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Total_Reads         = Aligned_Uniquely + Aligned_Ambiguously + Did_Not_Align + No_Genomic_Sequence,
    Aligned_Reads_Total = Aligned_Uniquely + Aligned_Ambiguously,
    Mapping_Efficiency  = if_else(Total_Reads > 0, 100 * Aligned_Uniquely / Total_Reads, NA_real_),
    Percent_Aligned_Total = if_else(Total_Reads > 0, 100 * Aligned_Reads_Total / Total_Reads, NA_real_),
    sample_id = paste0(Sample_name, "_merged_name_sorted")
  ) %>%
  relocate(Sample_name, sample_id)

head(align_collapsed)
align_collapsed$Sample_name

# ============================================================
# 2) Deduplication
# ============================================================
dedup_raw <- read_multiqc_txt(data_path("multiqc_bismark_dedup.txt"))

needed_dedup <- c("Sample", "dedup_reads_percent", "dup_reads_percent")
miss_dedup <- setdiff(needed_dedup, names(dedup_raw))
if (length(miss_dedup) > 0) {
  stop("multiqc_bismark_dedup.txt missing columns: ", paste(miss_dedup, collapse = ", "),
       "\nFound: ", paste(names(dedup_raw), collapse = ", "))
}

dedup_clean <- dedup_raw %>%
  mutate(Sample_name = get_core(Sample)) %>%
  transmute(
    Sample_name,
    Percent_Deduplicated_reads = dedup_reads_percent,
    Percent_Duplicate_Reads    = dup_reads_percent,
    Percent_Duplicated_Reads   = if_else(
      (Percent_Deduplicated_reads + Percent_Duplicate_Reads) > 0,
      100 * Percent_Duplicate_Reads / (Percent_Deduplicated_reads + Percent_Duplicate_Reads),
      NA_real_
    )
  )

head(dedup_clean)

# ============================================================
# 3) Methylation percent CpG
# ============================================================
meth_raw <- read_multiqc_txt(data_path("multiqc_bismark_methextract.txt"))

needed_meth <- c("Sample", "percent_cpg_meth")
miss_meth <- setdiff(needed_meth, names(meth_raw))
if (length(miss_meth) > 0) {
  stop("multiqc_bismark_methextract.txt missing columns: ", paste(miss_meth, collapse = ", "),
       "\nFound: ", paste(names(meth_raw), collapse = ", "))
}

meth_clean <- meth_raw %>%
  mutate(Sample_name = get_core(Sample)) %>%
  transmute(
    Sample_name,
    Percent_CpG_Methylated = percent_cpg_meth
  )

head(meth_clean)

# ============================================================
# Merge all tables by Sample_name
# ============================================================
tables_to_merge <- list(
  align_collapsed,
  dedup_clean,
  meth_clean
)

merged_multiqc <- tables_to_merge %>%
  discard(is.null) %>%
  reduce(full_join, by = "Sample_name")

final_cols <- c(
  "Sample_name",
  "sample_id",
  "SequencingNumber",
  "Percent_CpG_Methylated",
  "Percent_Duplicated_Reads",
  "Percent_Deduplicated_reads",
  "Percent_Duplicate_Reads",
  "Aligned_Uniquely",
  "Aligned_Ambiguously",
  "Aligned_Reads_Total",
  "Mapping_Efficiency",
  "Percent_Aligned_Total",
  "Total_Reads"
)

merged_multiqc_selected <- merged_multiqc %>%
  select(any_of(final_cols))

# ============================================================
# Merge with trait data
# ============================================================
merged_final <- merged_multiqc_selected

if (!is.na(trait_path) && nzchar(trait_path)) {
  if (!file.exists(trait_path)) stop("trait_path does not exist: ", trait_path)

  trait_data <- openxlsx::read.xlsx(trait_path, rowNames = TRUE)

  trait_df <- as.data.frame(trait_data) %>%
    rownames_to_column(var = "Sample_name")

  merged_final <- merged_multiqc_selected %>%
    inner_join(trait_df, by = "Sample_name")
}

# ---- SET Sample_name AS ROWNAMES (FINAL STEP) ----
merged_final <- merged_final %>%
  tibble::column_to_rownames(var = "Sample_name")

# ============================================================
# Write outputs
# ============================================================
openxlsx::write.xlsx(
  merged_final,
  out_path,
  overwrite = TRUE,
  rowNames  = TRUE
)

out_multiqc_only <- sub("\\.xlsx$", "_multiqc_only.xlsx", out_path)

openxlsx::write.xlsx(
  merged_multiqc_selected,
  out_multiqc_only,
  overwrite = TRUE,
  rowNames  = TRUE
)

message("Done.")
message("Wrote: ", out_path,          "  (rows=", nrow(merged_final),            ", cols=", ncol(merged_final),            ")")
message("Wrote: ", out_multiqc_only,  "  (rows=", nrow(merged_multiqc_selected), ", cols=", ncol(merged_multiqc_selected), ")")
message("Final columns present: ", paste(names(merged_multiqc_selected), collapse = ", "))