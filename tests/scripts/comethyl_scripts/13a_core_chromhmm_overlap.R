#!/usr/bin/env Rscript
# ============================================================
# 13A_core_regulatory_state_overlap.R
#
# Purpose
#   For each annotation directory from Script 12A and requested module:
#     1) load Annotated_Regions.tsv from 12_annotation
#     2) build BED of module regions
#     3) resolve regulatory annotation source(s):
#          - roadmap   : ChromHMM 15-state
#          - roadmap18 : ChromHMM 18-state
#          - encode    : cCRE BED
#     4) auto-download missing reference files if possible
#     5) if download fails, warn and continue
#     6) bedtools intersect module BED vs annotation tracks
#     7) collapse overlaps to dominant annotation per region x tissue/source
#     8) save long and wide tables
#
# Output structure
#   <project_root>/comethyl_output/13a_regulatory_state_overlap/<cpg_label>/<region_label>/<variant_name>/
#     run_log.txt
#     run_parameters.txt
#     resolved_reference_tracks.csv
#     module_<module>/
#       region_metadata.csv
#       regions_<module>.bed
#       dominant_state_long_all_sources.csv
#       dominant_state_matrix_all_sources.csv
#       roadmap/
#         dominant_state_long.csv
#         dominant_state_matrix.csv
#         raw_intersections/
#       encode/
#         ...
#       roadmap18/
#         ...
# ============================================================

message("Starting Script 13a")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

# ============================================================
# 1) CLI helpers
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

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

msg <- function(...) cat(sprintf(...), "\n")

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

validate_dir_exists <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !dir.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

validate_file_exists <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

# ============================================================
# 2) Output directory helper
# ============================================================
derive_pipeline_dirs_from_annotation_dir <- function(annotation_dir, project_root, step_name) {
  variant_name <- basename(annotation_dir)
  region_label <- basename(dirname(annotation_dir))
  cpg_label    <- basename(dirname(dirname(annotation_dir)))

  pipeline_root <- file.path(project_root, "comethyl_output")
  step_dir <- file.path(pipeline_root, step_name)
  out_dir <- file.path(step_dir, cpg_label, region_label, variant_name)

  safe_dir_create(out_dir)

  list(
    pipeline_root = pipeline_root,
    step_dir = step_dir,
    out_dir = out_dir,
    variant_name = variant_name,
    region_label = region_label,
    cpg_label = cpg_label
  )
}

# ============================================================
# 3) bedtools
# ============================================================
resolve_bedtools <- function(bedtools_path = "") {
  if (!is.na(bedtools_path) && nzchar(bedtools_path)) {
    if (!file.exists(bedtools_path)) stop("Provided --bedtools does not exist: ", bedtools_path, call. = FALSE)
    return(normalizePath(bedtools_path))
  }
  bt <- Sys.which("bedtools")
  if (bt == "") stop("bedtools not found on PATH. Install it in Pixi or provide --bedtools /full/path/to/bedtools", call. = FALSE)
  bt
}

# ============================================================
# 4) read annotated regions directly from 12_annotation
# ============================================================
read_annotated_regions_from_annotation_dir <- function(annotation_dir) {
  file_path <- file.path(annotation_dir, "Annotated_Regions.tsv")
  validate_file_exists(file_path, "Annotated_Regions.tsv")

  msg("[ANNOT] Using annotated regions file: %s", file_path)

  df <- read.delim(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  need <- c("RegionID", "module", "chr", "start", "end")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0) {
    stop("Annotated_Regions.tsv missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!("gene_symbol" %in% names(df))) df$gene_symbol <- NA_character_
  if (!("membership" %in% names(df)))  df$membership  <- NA_real_

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
      RegionID = as.character(RegionID),
      gene_symbol = as.character(gene_symbol),
      membership = suppressWarnings(as.numeric(.data[[membership_col]]))
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
# 5) source registry
# ============================================================
SUPPORTED_SOURCES <- c("roadmap", "encode", "roadmap18")

get_requested_sources <- function(source_arg = "") {
  if (is.null(source_arg) || is.na(source_arg) || !nzchar(source_arg)) {
    return(SUPPORTED_SOURCES)
  }
  src <- split_csv(source_arg)
  bad <- setdiff(src, SUPPORTED_SOURCES)
  if (length(bad) > 0) {
    stop(
      "Unsupported --source value(s): ", paste(bad, collapse = ", "),
      ". Supported: ", paste(SUPPORTED_SOURCES, collapse = ", "),
      call. = FALSE
    )
  }
  unique(src)
}

# ============================================================
# 6) reference metadata
# ============================================================
ROADMAP15_EIDS <- c("E081","E082","E071","E073","E091","E062","E029","E003","E020")
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

roadmap15_filename <- function(eid) paste0(eid, "_15_coreMarks_hg38lift_mnemonics.bed.gz")
roadmap15_url <- function(eid) {
  paste0(
    "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/",
    "ChmmModels/coreMarks/jointModel/final/",
    roadmap15_filename(eid)
  )
}

discover_local_bed_files <- function(dir_path, pattern = "\\.(bed|bed.gz)$") {
  if (!dir.exists(dir_path)) return(character(0))
  list.files(dir_path, pattern = pattern, full.names = TRUE)
}

download_file_safe <- function(url, destfile) {
  safe_dir_create(dirname(destfile))
  ok <- FALSE
  err <- NULL

  tryCatch({
    if (nzchar(Sys.which("curl"))) {
      status <- system2("curl", c("-L", "-f", "-o", destfile, url), stdout = FALSE, stderr = FALSE)
      ok <- identical(status, 0L) && file.exists(destfile) && file.info(destfile)$size > 0
    } else if (nzchar(Sys.which("wget"))) {
      status <- system2("wget", c("-O", destfile, url), stdout = FALSE, stderr = FALSE)
      ok <- identical(status, 0L) && file.exists(destfile) && file.info(destfile)$size > 0
    } else {
      utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
      ok <- file.exists(destfile) && file.info(destfile)$size > 0
    }
  }, error = function(e) {
    err <<- conditionMessage(e)
  })

  list(ok = ok, error = err, destfile = destfile, url = url)
}

# ============================================================
# 7) resolve source files
# ============================================================
resolve_roadmap15_files <- function(reference_root, download_missing = TRUE) {
  src_dir <- file.path(reference_root, "roadmap", "chromhmm_15state")
  safe_dir_create(src_dir)

  out <- tibble(
    source = character(),
    dataset = character(),
    id = character(),
    label = character(),
    file = character(),
    ok = logical(),
    note = character()
  )

  for (eid in ROADMAP15_EIDS) {
    f <- file.path(src_dir, roadmap15_filename(eid))

    if (!file.exists(f) && isTRUE(download_missing)) {
      msg("[DL][roadmap] Missing %s -> attempting download", basename(f))
      dl <- download_file_safe(roadmap15_url(eid), f)
      if (!isTRUE(dl$ok)) {
        msg("[WARN][roadmap] Download failed for %s", eid)
      }
    }

    this_ok <- file.exists(f) && file.info(f)$size > 0
    out <- bind_rows(out, tibble(
      source = "roadmap",
      dataset = "ChromHMM_15state",
      id = eid,
      label = unname(ROADMAP15_LABELS[eid]),
      file = f,
      ok = this_ok,
      note = ifelse(this_ok, "", "missing_or_download_failed")
    ))
  }

  out
}

resolve_roadmap18_files <- function(reference_root, roadmap18_dir = NULL) {
  src_dir <- if (!is.null(roadmap18_dir) && nzchar(roadmap18_dir)) roadmap18_dir else
    file.path(reference_root, "roadmap18", "chromhmm_18state")

  files <- discover_local_bed_files(src_dir)
  if (length(files) == 0) {
    msg("[WARN][roadmap18] No local BED files found in: %s", src_dir)
    return(tibble(
      source = "roadmap18",
      dataset = "ChromHMM_18state",
      id = character(),
      label = character(),
      file = character(),
      ok = logical(),
      note = character()
    ))
  }

  tibble(
    source = "roadmap18",
    dataset = "ChromHMM_18state",
    id = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    label = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    file = files,
    ok = TRUE,
    note = ""
  )
}

resolve_encode_ccre_files <- function(reference_root, encode_dir = NULL) {
  src_dir <- if (!is.null(encode_dir) && nzchar(encode_dir)) encode_dir else
    file.path(reference_root, "encode", "ccre")

  files <- discover_local_bed_files(src_dir)
  if (length(files) == 0) {
    msg("[WARN][encode] No local BED files found in: %s", src_dir)
    return(tibble(
      source = "encode",
      dataset = "cCRE",
      id = character(),
      label = character(),
      file = character(),
      ok = logical(),
      note = character()
    ))
  }

  tibble(
    source = "encode",
    dataset = "cCRE",
    id = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    label = tools::file_path_sans_ext(tools::file_path_sans_ext(basename(files))),
    file = files,
    ok = TRUE,
    note = ""
  )
}

resolve_requested_reference_files <- function(
    requested_sources,
    reference_root,
    download_missing = TRUE,
    roadmap18_dir = NULL,
    encode_dir = NULL) {

  pieces <- list()

  if ("roadmap" %in% requested_sources) {
    pieces[["roadmap"]] <- resolve_roadmap15_files(
      reference_root = reference_root,
      download_missing = download_missing
    )
  }

  if ("roadmap18" %in% requested_sources) {
    pieces[["roadmap18"]] <- resolve_roadmap18_files(
      reference_root = reference_root,
      roadmap18_dir = roadmap18_dir
    )
  }

  if ("encode" %in% requested_sources) {
    pieces[["encode"]] <- resolve_encode_ccre_files(
      reference_root = reference_root,
      encode_dir = encode_dir
    )
  }

  bind_rows(pieces)
}

# ============================================================
# 8) state normalization
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
    if (s2 %in% names(ROADMAP15_DESC)) return(ROADMAP15_DESC[[s2]])
    return(s2)
  }

  if (source %in% c("roadmap18", "encode")) {
    return(raw_state)
  }

  raw_state
}

# ============================================================
# 9) bedtools intersection
# ============================================================
run_intersections_one_source <- function(bed_file, ref_tbl, out_dir, bedtools_bin) {
  safe_dir_create(out_dir)
  out_files <- character(0)

  for (i in seq_len(nrow(ref_tbl))) {
    ref_file <- ref_tbl$file[i]
    ref_id   <- ref_tbl$id[i]

    if (!file.exists(ref_file)) {
      msg("[WARN] Reference missing, skipping: %s", ref_file)
      next
    }

    out_tsv <- file.path(out_dir, paste0(ref_id, ".overlap.tsv"))

    cmd <- sprintf(
      "%s intersect -a %s -b %s -wa -wb > %s",
      shQuote(bedtools_bin),
      shQuote(bed_file),
      shQuote(ref_file),
      shQuote(out_tsv)
    )

    msg("[BEDTOOLS] %s vs %s", basename(bed_file), basename(ref_file))
    st <- system(cmd)
    if (st != 0) {
      msg("[WARN] bedtools failed for %s", ref_id)
      next
    }
    out_files <- c(out_files, out_tsv)
  }

  out_files
}

# ============================================================
# 10) read + collapse overlaps
# ============================================================
guess_overlap_columns <- function(nc) {
  if (nc < 8) stop("Intersect output has too few columns: ", nc, call. = FALSE)

  a_cols <- c("chrA", "startA", "endA", "region_id")
  remaining <- nc - 4

  if (remaining == 4) {
    b_cols <- c("chrB", "startB", "endB", "state")
  } else {
    extra_n <- remaining - 4
    b_cols <- c("chrB", "startB", "endB", paste0("b_extra", seq_len(extra_n)), "state")
  }
  c(a_cols, b_cols)
}

read_and_collapse_overlaps <- function(overlap_files, source_name, ref_tbl) {
  if (length(overlap_files) == 0) {
    return(tibble(
      source = character(),
      dataset = character(),
      track_id = character(),
      track_label = character(),
      region_id = character(),
      annotation = character(),
      overlap_bp = integer()
    ))
  }

  all <- lapply(overlap_files, function(f) {
    ref_id <- str_extract(basename(f), "^[^.]+")
    meta <- ref_tbl %>% filter(id == ref_id) %>% slice(1)

    if (!file.exists(f) || file.info(f)$size == 0) {
      return(NULL)
    }

    x <- suppressWarnings(read_tsv(f, col_names = FALSE, show_col_types = FALSE))
    if (nrow(x) == 0) return(NULL)

    names(x) <- guess_overlap_columns(ncol(x))

    x %>%
      mutate(
        source = source_name,
        dataset = meta$dataset,
        track_id = ref_id,
        track_label = meta$label,
        overlap_bp = pmax(0L, pmin(endA, endB) - pmax(startA, startB)),
        annotation = as.character(state)
      ) %>%
      filter(overlap_bp > 0) %>%
      select(source, dataset, track_id, track_label, region_id, annotation, overlap_bp)
  })

  all <- bind_rows(all)
  if (nrow(all) == 0) return(all)

  all %>%
    group_by(source, dataset, track_id, track_label, region_id) %>%
    slice_max(order_by = overlap_bp, n = 1, with_ties = FALSE) %>%
    ungroup()
}

# ============================================================
# 11) arguments
# ============================================================
project_root <- trim_or_null(get_arg("--project_root"))

annotation_dir_v1 <- trim_or_null(get_arg("--annotation_dir_v1"))
annotation_dir_v2 <- trim_or_null(get_arg("--annotation_dir_v2"))
annotation_dir_v3 <- trim_or_null(get_arg("--annotation_dir_v3"))

modules_in <- unique(split_csv(get_arg("--modules", "")))
requested_sources <- get_requested_sources(get_arg("--source", ""))

reference_root <- trim_or_null(get_arg(
  "--reference_root",
  file.path(project_root, "reference_data", "regulatory_annotations")
))
roadmap18_dir <- trim_or_null(get_arg("--roadmap18_dir", ""))
encode_dir <- trim_or_null(get_arg("--encode_dir", ""))
download_missing <- as_bool(get_arg("--download_missing", "true"), TRUE)

out_parent <- trim_or_null(get_arg(
  "--out_parent",
  file.path(project_root, "comethyl_output", "13a_regulatory_state_overlap")
))

membership_col <- trim_or_null(get_arg("--membership_col", "membership"))

bedtools_bin <- resolve_bedtools(get_arg("--bedtools", ""))

if (is.null(project_root)) stop("--project_root is required", call. = FALSE)
if (length(modules_in) == 0) stop("You must provide --modules", call. = FALSE)

validate_dir_exists(project_root, "project_root")

annotation_dirs <- c(
  v1_all_pcs = annotation_dir_v1,
  v2_exclude_protected_pcs = annotation_dir_v2,
  v3_technical_pcs_only = annotation_dir_v3
)
annotation_dirs <- annotation_dirs[!is.na(annotation_dirs) & nzchar(annotation_dirs)]

if (length(annotation_dirs) == 0) {
  stop("Provide at least one of --annotation_dir_v1, --annotation_dir_v2, --annotation_dir_v3", call. = FALSE)
}

for (nm in names(annotation_dirs)) {
  validate_dir_exists(annotation_dirs[[nm]], nm)
  validate_file_exists(file.path(annotation_dirs[[nm]], "Annotated_Regions.tsv"), paste0(nm, "/Annotated_Regions.tsv"))
}

setwd(project_root)

msg("project_root: %s", project_root)
msg("annotation_dirs:")
for (nm in names(annotation_dirs)) msg("  - %s : %s", nm, annotation_dirs[[nm]])
msg("modules: %s", paste(modules_in, collapse = ", "))
msg("requested_sources: %s", paste(requested_sources, collapse = ", "))
msg("download_missing: %s", download_missing)
msg("reference_root: %s", reference_root)
msg("out_parent: %s", out_parent)
msg("bedtools: %s", bedtools_bin)

safe_dir_create(reference_root)
safe_dir_create(out_parent)

# ============================================================
# 12) resolve sources once
# ============================================================
ref_tbl_all <- resolve_requested_reference_files(
  requested_sources = requested_sources,
  reference_root = reference_root,
  download_missing = download_missing,
  roadmap18_dir = roadmap18_dir,
  encode_dir = encode_dir
)

if (nrow(ref_tbl_all) == 0) {
  stop("No reference tracks resolved for requested source(s).", call. = FALSE)
}

write_csv(ref_tbl_all, file.path(out_parent, "resolved_reference_tracks.csv"))

usable_sources <- ref_tbl_all %>%
  filter(ok) %>%
  count(source, name = "n_tracks")

if (nrow(usable_sources) == 0) {
  stop("None of the requested sources are usable. Check local files and/or internet access.", call. = FALSE)
}

msg("Usable sources:")
for (i in seq_len(nrow(usable_sources))) {
  msg("  - %s : %d track(s)", usable_sources$source[i], usable_sources$n_tracks[i])
}

# ============================================================
# 13) process annotation directories
# ============================================================
for (variant_label in names(annotation_dirs)) {
  annotation_dir <- annotation_dirs[[variant_label]]

  dir_info <- derive_pipeline_dirs_from_annotation_dir(
    annotation_dir = annotation_dir,
    project_root = project_root,
    step_name = "13a_regulatory_state_overlap"
  )

  variant_out_dir <- dir_info$out_dir
  safe_dir_create(variant_out_dir)

  variant_log <- file.path(variant_out_dir, "run_log.txt")
  variant_params <- file.path(variant_out_dir, "run_parameters.txt")

  append_log(variant_log, "Starting variant: ", variant_label)
  append_log(variant_log, "annotation_dir: ", annotation_dir)
  append_log(variant_log, "pipeline_root: ", dir_info$pipeline_root)
  append_log(variant_log, "step_dir: ", dir_info$step_dir)
  append_log(variant_log, "cpg_label: ", dir_info$cpg_label)
  append_log(variant_log, "region_label: ", dir_info$region_label)
  append_log(variant_log, "variant_name: ", dir_info$variant_name)
  append_log(variant_log, "variant_out_dir: ", variant_out_dir)
  append_log(variant_log, "modules: ", paste(modules_in, collapse = ", "))
  append_log(variant_log, "sources: ", paste(requested_sources, collapse = ", "))
  append_log(variant_log, "bedtools: ", bedtools_bin)

  anno <- read_annotated_regions_from_annotation_dir(annotation_dir)

  append_log(variant_log, "Annotated rows loaded: ", nrow(anno))
  append_log(variant_log, "Unique modules found: ", length(unique(anno$module)))

  for (mod in modules_in) {
    msg("\n==================================================")
    msg("[RUN] variant=%s | module=%s", variant_label, mod)
    append_log(variant_log, "Processing module: ", mod)

    mod_df <- anno %>% filter(module == mod)
    if (nrow(mod_df) == 0) {
      msg("[INFO] No regions found for module '%s' in %s", mod, variant_label)
      append_log(variant_log, "No regions found for module: ", mod)
      next
    }

    out_dir <- file.path(variant_out_dir, paste0("module_", mod))
    safe_dir_create(out_dir)

    region_meta_map <- build_region_meta_map(mod_df, membership_col = membership_col)
    write_csv(region_meta_map, file.path(out_dir, "region_metadata.csv"))

    bed <- mod_df %>%
      select(chr, start, end, RegionID) %>%
      distinct() %>%
      filter(!is.na(chr), !is.na(start), !is.na(end), !is.na(RegionID)) %>%
      mutate(chr = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr))) %>%
      arrange(chr, start, end)

    if (nrow(bed) == 0) {
      msg("[WARN] No valid BED rows for module %s", mod)
      append_log(variant_log, "No valid BED rows for module: ", mod)
      next
    }

    bed_file <- file.path(out_dir, paste0("regions_", mod, ".bed"))
    write_tsv(bed, bed_file, col_names = FALSE)

    source_results <- list()

    for (src in requested_sources) {
      ref_tbl <- ref_tbl_all %>% filter(source == src, ok)
      if (nrow(ref_tbl) == 0) {
        msg("[WARN] No usable tracks for source=%s; skipping", src)
        append_log(variant_log, "No usable tracks for source: ", src)
        next
      }

      src_out_dir <- file.path(out_dir, src)
      ov_dir <- file.path(src_out_dir, "raw_intersections")
      safe_dir_create(ov_dir)

      overlap_files <- run_intersections_one_source(
        bed_file = bed_file,
        ref_tbl = ref_tbl,
        out_dir = ov_dir,
        bedtools_bin = bedtools_bin
      )

      if (length(overlap_files) == 0) {
        msg("[WARN] No overlap files created for source=%s", src)
        append_log(variant_log, "No overlap files created for source: ", src, " module: ", mod)
        next
      }

      dom <- read_and_collapse_overlaps(
        overlap_files = overlap_files,
        source_name = src,
        ref_tbl = ref_tbl
      )

      if (nrow(dom) == 0) {
        msg("[WARN] No collapsed overlaps for source=%s", src)
        append_log(variant_log, "No collapsed overlaps for source: ", src, " module: ", mod)
        next
      }

      dom2 <- dom %>%
        left_join(region_meta_map, by = c("region_id" = "RegionID")) %>%
        mutate(
          module = mod,
          annotation_desc = vapply(annotation, function(x) normalize_annotation_label(src, x), character(1))
        ) %>%
        relocate(module, .before = source) %>%
        relocate(gene_symbol, membership, .after = region_id)

      write_csv(dom2, file.path(src_out_dir, "dominant_state_long.csv"))

      mat <- dom2 %>%
        mutate(
          column_label = paste(source, track_label, sep = " | "),
          row_label = ifelse(!is.na(gene_symbol) & gene_symbol != "",
                             paste0(region_id, " | ", gene_symbol),
                             region_id)
        ) %>%
        select(row_label, column_label, annotation_desc) %>%
        distinct() %>%
        pivot_wider(names_from = column_label, values_from = annotation_desc)

      write_csv(mat, file.path(src_out_dir, "dominant_state_matrix.csv"))

      source_results[[src]] <- dom2
      msg("[OK] Completed source=%s for module=%s", src, mod)
      append_log(variant_log, "Completed source: ", src, " for module: ", mod)
    }

    if (length(source_results) == 0) {
      msg("[WARN] No usable source results for variant=%s module=%s", variant_label, mod)
      append_log(variant_log, "No usable source results for module: ", mod)
      next
    }

    combined <- bind_rows(source_results)
    write_csv(combined, file.path(out_dir, "dominant_state_long_all_sources.csv"))

    combined_mat <- combined %>%
      mutate(
        column_label = paste(source, track_label, sep = " | "),
        row_label = ifelse(!is.na(gene_symbol) & gene_symbol != "",
                           paste0(region_id, " | ", gene_symbol),
                           region_id)
      ) %>%
      select(row_label, column_label, annotation_desc) %>%
      distinct() %>%
      pivot_wider(names_from = column_label, values_from = annotation_desc)

    write_csv(combined_mat, file.path(out_dir, "dominant_state_matrix_all_sources.csv"))

    append_log(variant_log, "Finished module: ", mod, " with ", nrow(combined), " overlap rows")
  }

  params <- c(
    paste0("timestamp\t", timestamp_now()),
    paste0("project_root\t", project_root),
    paste0("annotation_dir\t", annotation_dir),
    paste0("pipeline_root\t", dir_info$pipeline_root),
    paste0("step_dir\t", dir_info$step_dir),
    paste0("cpg_label\t", dir_info$cpg_label),
    paste0("region_label\t", dir_info$region_label),
    paste0("variant_name\t", dir_info$variant_name),
    paste0("variant_out_dir\t", variant_out_dir),
    paste0("modules\t", paste(modules_in, collapse = ",")),
    paste0("requested_sources\t", paste(requested_sources, collapse = ",")),
    paste0("reference_root\t", reference_root),
    paste0("download_missing\t", download_missing),
    paste0("bedtools\t", bedtools_bin),
    paste0("membership_col\t", membership_col)
  )
  write_lines_safe(params, variant_params)

  append_log(variant_log, "Finished variant: ", variant_label)
}

cat("\nScript 13a complete: Core Chromatin State finished\n")