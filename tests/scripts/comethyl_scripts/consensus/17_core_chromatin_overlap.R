#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 17: Consensus Regulatory State Overlap
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   For each requested module, overlaps consensus regions against
#   regulatory annotation reference tracks using bedtools:
#     1) Load Annotated_Regions.tsv from script 15
#     2) Build a BED file for each module's regions
#     3) Resolve regulatory annotation sources:
#          roadmap   : Roadmap Epigenomics ChromHMM 15-state
#          roadmap18 : ChromHMM 18-state (local files required)
#          encode    : ENCODE cCRE BED (local files required)
#     4) Auto-download missing Roadmap 15-state files if possible
#     5) Run bedtools intersect for each module vs each track
#     6) Collapse overlaps to dominant annotation per region × tissue
#     7) Save long and wide overlap tables per module and source
#
# NOTE
#   Module assignments are shared across datasets — regulatory
#   overlap runs once per adjustment version. No per-dataset loop.
#
# REQUIRED INPUTS
#   --project_root   : root directory of the analysis project
#   --annotation_dir : path to script 15 output directory
#   --modules        : comma-separated module names to process
#                      (e.g. turquoise,blue,brown)
#   --reference_root : directory where reference BED files are stored
#                      or will be downloaded to
#
# OPTIONAL INPUTS
#   --adjustment_version  : label matching script 15 run [default = unadjusted]
#   --source              : comma-separated sources: roadmap,roadmap18,encode
#                           [default = all three]
#   --bedtools            : path to bedtools binary [default = auto-detect on PATH]
#   --download_missing    : auto-download missing Roadmap files TRUE/FALSE [default = TRUE]
#   --roadmap18_dir       : custom directory for 18-state Roadmap BEDs
#   --encode_dir          : custom directory for ENCODE cCRE BEDs
#
# OUTPUTS
#   comethyl_output/consensus/17_regulatory_overlap/<adjustment_version>/
#       resolved_reference_tracks.csv
#       module_<module>/
#           regions_<module>.bed
#           region_metadata.csv
#           dominant_state_long_all_sources.csv
#           dominant_state_matrix_all_sources.csv
#           <source>/
#               dominant_state_long.csv
#               dominant_state_matrix.csv
#               raw_intersections/
#                   <track_id>.overlap.tsv
#
# NOTES
#   - bedtools must be installed and accessible on PATH or via --bedtools
#   - Roadmap 15-state files are auto-downloaded if not present and
#     --download_missing TRUE; all other sources require local files
#   - Output of this script feeds into script 18 (--regulatory_overlap_dir)
#
# EXAMPLE
#   Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/consensus/17_core_chromatin_overlap.R \
#     --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#     --annotation_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/15_module_annotation/cov3_75pct/covMin4_methSD0p08/v1_all_pcs \
#     --modules turquoise,blue,brown \
#     --reference_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/reference/chromatin \
#     --adjustment_version v1_all_pcs \
#     --source roadmap \
#     --download_missing TRUE
#
# PURPOSE
#   For each requested module:
#     1) Load Annotated_Regions.tsv from script 15
#     2) Build BED of module regions
#     3) Resolve regulatory annotation source(s):
#          - roadmap   : ChromHMM 15-state
#          - roadmap18 : ChromHMM 18-state
#          - encode    : cCRE BED
#     4) Auto-download missing reference files if possible
#     5) Run bedtools intersect module BED vs annotation tracks
#     6) Collapse overlaps to dominant annotation per region x tissue
#     7) Save long and wide tables
#
# NOTE
#   Module assignments are shared across datasets — regulatory
#   overlap runs once per adjustment version only.
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/17_regulatory_overlap/<adjustment_version>/
#       resolved_reference_tracks.csv
#       module_<MODULE>/
#           regions_<MODULE>.bed
#           region_metadata.csv
#           dominant_state_long_all_sources.csv
#           dominant_state_matrix_all_sources.csv
#           <source>/
#               dominant_state_long.csv
#               dominant_state_matrix.csv
#               raw_intersections/
#                   <track_id>.overlap.tsv
# ================================================================
message("Starting Script 17 ")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

# ============================================================
# Helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

msg <- function(...) message(sprintf(...))

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# bedtools
# ============================================================
resolve_bedtools <- function(bedtools_path = "") {
  if (!is.na(bedtools_path) && nzchar(bedtools_path)) {
    if (!file.exists(bedtools_path))
      stop("--bedtools path not found: ", bedtools_path)
    return(normalizePath(bedtools_path))
  }
  bt <- Sys.which("bedtools")
  if (bt == "")
    stop("bedtools not found on PATH. Provide --bedtools /full/path/to/bedtools")
  bt
}

# ============================================================
# Load annotated regions from script 15
# ============================================================
read_annotated_regions <- function(annotation_dir) {
  f <- file.path(annotation_dir, "Annotated_Regions.tsv")
  if (!file.exists(f))
    stop("Annotated_Regions.tsv not found in annotation_dir: ", annotation_dir)
  msg("[ANNOT] Loading: %s", f)

  df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)

  need    <- c("RegionID", "module", "chr", "start", "end")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0)
    stop("Annotated_Regions.tsv missing columns: ",
         paste(missing, collapse = ", "))

  if (!"gene_symbol" %in% names(df)) df$gene_symbol <- NA_character_
  if (!"membership"  %in% names(df)) df$membership  <- NA_real_

  df %>%
    mutate(
      RegionID    = as.character(RegionID),
      module      = as.character(module),
      chr         = as.character(chr),
      start       = suppressWarnings(as.integer(start)),
      end         = suppressWarnings(as.integer(end)),
      gene_symbol = as.character(gene_symbol),
      membership  = suppressWarnings(as.numeric(membership))
    )
}

build_region_meta_map <- function(mod_df, membership_col = "membership") {
  if (!membership_col %in% names(mod_df)) mod_df[[membership_col]] <- NA_real_

  mod_df %>%
    transmute(
      RegionID    = as.character(RegionID),
      gene_symbol = as.character(gene_symbol),
      membership  = suppressWarnings(as.numeric(.data[[membership_col]]))
    ) %>%
    mutate(gene_symbol = ifelse(is.na(gene_symbol), "", str_trim(gene_symbol))) %>%
    group_by(RegionID) %>%
    summarise(
      gene_symbol = {
        gs <- unique(gene_symbol[gene_symbol != ""])
        if (length(gs) == 0) NA_character_ else paste(gs, collapse = ";")
      },
      membership = {
        mm <- suppressWarnings(as.numeric(membership))
        if (all(is.na(mm))) NA_real_ else max(mm, na.rm = TRUE)
      },
      .groups = "drop"
    )
}

# ============================================================
# Source registry
# ============================================================
SUPPORTED_SOURCES <- c("roadmap", "encode", "roadmap18")

get_requested_sources <- function(source_arg = "") {
  if (is.null(source_arg) || is.na(source_arg) || !nzchar(source_arg))
    return(SUPPORTED_SOURCES)
  src <- split_csv(source_arg)
  bad <- setdiff(src, SUPPORTED_SOURCES)
  if (length(bad) > 0)
    stop("Unsupported --source value(s): ", paste(bad, collapse = ", "),
         ". Supported: ", paste(SUPPORTED_SOURCES, collapse = ", "))
  unique(src)
}

# ============================================================
# Reference metadata
# ============================================================
ROADMAP15_EIDS <- c("E081","E082","E071","E073","E091",
                    "E062","E029","E003","E020")
ROADMAP15_LABELS <- c(
  E081 = "FetalBrain_M",
  E082 = "FetalBrain_F",
  E071 = "Brain_Hippocampus",
  E073 = "Brain_DLPFC",
  E091 = "Placenta_Fetal",
  E062 = "Blood_PBMC",
  E029 = "Blood_Monocytes",
  E003 = "ESC_H1",
  E020 = "iPSC"
)

roadmap15_filename <- function(eid)
  paste0(eid, "_15_coreMarks_hg38lift_mnemonics.bed.gz")

roadmap15_url <- function(eid)
  paste0(
    "https://egg2.wustl.edu/roadmap/data/byFileType/",
    "chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/",
    roadmap15_filename(eid)
  )

discover_local_bed_files <- function(dir_path) {
  if (!dir.exists(dir_path)) return(character(0))
  list.files(dir_path, pattern = "\\.(bed|bed\\.gz)$", full.names = TRUE)
}

download_file_safe <- function(url, destfile) {
  safe_dir_create(dirname(destfile))
  ok  <- FALSE
  err <- NULL
  tryCatch({
    if (nzchar(Sys.which("curl"))) {
      st <- system2("curl", c("-L", "-f", "-o", destfile, url),
                    stdout = FALSE, stderr = FALSE)
      ok <- identical(st, 0L) && file.exists(destfile) &&
            file.info(destfile)$size > 0
    } else if (nzchar(Sys.which("wget"))) {
      st <- system2("wget", c("-O", destfile, url),
                    stdout = FALSE, stderr = FALSE)
      ok <- identical(st, 0L) && file.exists(destfile) &&
            file.info(destfile)$size > 0
    } else {
      utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
      ok <- file.exists(destfile) && file.info(destfile)$size > 0
    }
  }, error = function(e) { err <<- conditionMessage(e) })
  list(ok = ok, error = err, destfile = destfile, url = url)
}

# ============================================================
# Resolve source files
# ============================================================
resolve_roadmap15_files <- function(reference_root, download_missing = TRUE) {
  src_dir <- file.path(reference_root, "roadmap", "chromhmm_15state")
  safe_dir_create(src_dir)

  out <- tibble(source = character(), dataset = character(),
                id = character(), label = character(),
                file = character(), ok = logical(), note = character())

  for (eid in ROADMAP15_EIDS) {
    f      <- file.path(src_dir, roadmap15_filename(eid))
    if (!file.exists(f) && isTRUE(download_missing)) {
      msg("[DL][roadmap] Missing %s -> attempting download", basename(f))
      dl <- download_file_safe(roadmap15_url(eid), f)
      if (!isTRUE(dl$ok)) msg("[WARN][roadmap] Download failed for %s", eid)
    }
    this_ok <- file.exists(f) && file.info(f)$size > 0
    out <- bind_rows(out, tibble(
      source  = "roadmap",
      dataset = "ChromHMM_15state",
      id      = eid,
      label   = unname(ROADMAP15_LABELS[eid]),
      file    = f,
      ok      = this_ok,
      note    = ifelse(this_ok, "", "missing_or_download_failed")
    ))
  }
  out
}

resolve_roadmap18_files <- function(reference_root, roadmap18_dir = NULL) {
  src_dir <- if (!is.null(roadmap18_dir) && nzchar(roadmap18_dir))
    roadmap18_dir else
    file.path(reference_root, "roadmap18", "chromhmm_18state")

  files <- discover_local_bed_files(src_dir)
  if (length(files) == 0) {
    msg("[WARN][roadmap18] No local BED files in: %s", src_dir)
    return(tibble(source = "roadmap18", dataset = "ChromHMM_18state",
                  id = character(), label = character(),
                  file = character(), ok = logical(), note = character()))
  }
  tibble(
    source  = "roadmap18",
    dataset = "ChromHMM_18state",
    id      = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    label   = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    file    = files,
    ok      = TRUE,
    note    = ""
  )
}

resolve_encode_ccre_files <- function(reference_root, encode_dir = NULL) {
  src_dir <- if (!is.null(encode_dir) && nzchar(encode_dir))
    encode_dir else
    file.path(reference_root, "encode", "ccre")

  files <- discover_local_bed_files(src_dir)
  if (length(files) == 0) {
    msg("[WARN][encode] No local BED files in: %s", src_dir)
    return(tibble(source = "encode", dataset = "cCRE",
                  id = character(), label = character(),
                  file = character(), ok = logical(), note = character()))
  }
  tibble(
    source  = "encode",
    dataset = "cCRE",
    id      = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    label   = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    file    = files,
    ok      = TRUE,
    note    = ""
  )
}

resolve_all_reference_files <- function(requested_sources, reference_root,
                                         download_missing = TRUE,
                                         roadmap18_dir = NULL,
                                         encode_dir    = NULL) {
  pieces <- list()
  if ("roadmap"   %in% requested_sources)
    pieces[["roadmap"]]   <- resolve_roadmap15_files(reference_root, download_missing)
  if ("roadmap18" %in% requested_sources)
    pieces[["roadmap18"]] <- resolve_roadmap18_files(reference_root, roadmap18_dir)
  if ("encode"    %in% requested_sources)
    pieces[["encode"]]    <- resolve_encode_ccre_files(reference_root, encode_dir)
  bind_rows(pieces)
}

# ============================================================
# State normalization
# ============================================================
ROADMAP15_DESC <- c(
  "1_TssA"      = "Active TSS",
  "2_TssAFlnk"  = "Flanking Active TSS",
  "3_TxFlnk"    = "Transcribed at 5' and 3'",
  "4_Tx"        = "Strong transcription",
  "5_TxWk"      = "Weak transcription",
  "6_EnhG"      = "Genic enhancer",
  "7_Enh"       = "Enhancer",
  "8_ZNF/Rpts"  = "ZNF genes and repeats",
  "9_Het"       = "Heterochromatin",
  "10_TssBiv"   = "Bivalent TSS",
  "11_BivFlnk"  = "Flanking bivalent TSS/enhancer",
  "12_EnhBiv"   = "Bivalent enhancer",
  "13_ReprPC"   = "Repressed PolyComb",
  "14_ReprPCWk" = "Weak repressed PolyComb",
  "15_Quies"    = "Quiescent/Low"
)

normalize_annotation_label <- function(source, raw_state) {
  if (is.na(raw_state) || !nzchar(raw_state)) return(NA_character_)
  raw_state <- as.character(raw_state)
  if (source == "roadmap") {
    if (raw_state %in% names(ROADMAP15_DESC)) return(ROADMAP15_DESC[[raw_state]])
    s2 <- str_replace(raw_state, "^E[0-9]{3}_", "")
    if (s2 %in% names(ROADMAP15_DESC)) ret