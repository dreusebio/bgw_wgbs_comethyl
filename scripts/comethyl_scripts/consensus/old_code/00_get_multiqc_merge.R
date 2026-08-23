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
# Change directory to folder with your outputs from multiqc and then run the scripts 
# e.g cd /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/data/09_multiqc/multiqc_data
#
# Usage:
#   Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_get_multiqc_merge.R
#   Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_get_multiqc_merge.R --trait_path /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/data/metadata/Sample_names_GROWELL.xlsx
# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_get_multiqc_merge.R --trait_path /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/data/metadata/Sample_names_GROWELL.xlsx --out /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/data/metadata/merged_qc.xlsx
# ============================================================
#!/usr/bin/env Rscript

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

trait_path <- get_arg("--trait_path", default = NA_character_)
out_path   <- get_arg("--out", default = "merged_multiqc_selected.xlsx")

# ----------------------------
# Helpers
# ----------------------------
read_multiqc_txt <- function(path) {
  # MultiQC text tables are tab-delimited with header
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    as.data.frame()
}

# "201_S1_L004_..." -> "201"; "201_merged_name_sorted" -> "201"
get_core <- function(x) str_extract(x, "^[0-9]+[A-Za-z]*")

# ----------------------------
# Required files (based on your heads)
# ----------------------------
req_files <- c(
  "mqc_bismark_alignment_1.txt",
  "multiqc_bismark_dedup.txt",
  "multiqc_bismark_methextract.txt"
)   # adjust if your file name differs

missing_files <- req_files[!file.exists(req_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required file(s) in current directory:\n  - ",
    paste(missing_files, collapse = "\n  - ")
  )
}

# ============================================================
# 1) Alignment (lane-level) -> collapse to Sample_name
# Header (from you):
# Sample, Aligned Uniquely, Aligned Ambiguously, Did Not Align, No Genomic Sequence
# ============================================================
align_raw <- read_multiqc_txt("mqc_bismark_alignment_1.txt")   # adjust if your file name differs

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

    # Sum across lanes
    Aligned_Uniquely    = sum(`Aligned Uniquely`, na.rm = TRUE),
    Aligned_Ambiguously = sum(`Aligned Ambiguously`, na.rm = TRUE),
    Did_Not_Align       = sum(`Did Not Align`, na.rm = TRUE),
    No_Genomic_Sequence = sum(`No Genomic Sequence`, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    # Bismark-style totals
    Total_Reads = Aligned_Uniquely + Aligned_Ambiguously + Did_Not_Align + No_Genomic_Sequence,

    # Keep both metrics (do not conflate them)
    Aligned_Reads_Total = Aligned_Uniquely + Aligned_Ambiguously,

    # Bismark "Mapping efficiency" = uniquely aligned / total
    Mapping_Efficiency = if_else(Total_Reads > 0, 100 * Aligned_Uniquely / Total_Reads, NA_real_),

    # (Optional) informative but not Bismark mapping efficiency
    Percent_Aligned_Total = if_else(Total_Reads > 0, 100 * Aligned_Reads_Total / Total_Reads, NA_real_),

    # keep your sample_id convention if you still want it
    sample_id = paste0(Sample_name, "_merged_name_sorted")
  ) %>%
  relocate(Sample_name, sample_id)

head(align_collapsed)
align_collapsed$Sample_name

# ============================================================
# 2) Deduplication
# Header (from you):
# Sample, dedup_reads_percent, dup_reads_percent
# ============================================================
dedup_raw <- read_multiqc_txt("multiqc_bismark_dedup.txt")   # adjust if your file name differs

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
    Percent_Duplicate_Reads = dup_reads_percent,
    Percent_Duplicated_Reads = if_else(
      (Percent_Deduplicated_reads + Percent_Duplicate_Reads) > 0,
      100 * Percent_Duplicate_Reads / (Percent_Deduplicated_reads + Percent_Duplicate_Reads),
      NA_real_
    )
  )

head(dedup_clean)
# ============================================================
# 3) Methylation percent CpG
# Header (from you):
# Sample, Methylated CpG
# ============================================================
meth_raw <- read_multiqc_txt("multiqc_bismark_methextract.txt")

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
# Merge all tables by Sample_name (skip NULL tables safely)
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
# merge with trait data (exact load as requested)
# trait_data <- openxlsx::read.xlsx(path, rowNames=T)
# and then rownames -> Sample_name
# ============================================================

merged_final <- merged_multiqc_selected

if (!is.na(trait_path) && nzchar(trait_path)) {
  if (!file.exists(trait_path)) stop("trait_path does not exist: ", trait_path)

  trait_data <- openxlsx::read.xlsx(trait_path, rowNames = TRUE)

  trait_df <- as.data.frame(trait_data) %>%
    rownames_to_column(var = "Sample_name")

  # If trait Sample_name has prefixes/suffixes and you want numeric core only, uncomment:
  # trait_df$Sample_name <- str_extract(trait_df$Sample_name, "^[0-9]+")

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
  rowNames = TRUE
)

out_multiqc_only <- sub("\\.xlsx$", "_multiqc_only.xlsx", out_path)

openxlsx::write.xlsx(
  merged_multiqc_selected,
  out_multiqc_only,
  overwrite = TRUE,
  rowNames = TRUE
)

message("Done.")
message("Wrote: ", out_path, "  (rows=", nrow(merged_final), ", cols=", ncol(merged_final), ")")
message("Wrote: ", out_multiqc_only, "  (rows=", nrow(merged_multiqc_selected), ", cols=", ncol(merged_multiqc_selected), ")")
message("Final columns present: ", paste(names(merged_multiqc_selected), collapse = ", "))
