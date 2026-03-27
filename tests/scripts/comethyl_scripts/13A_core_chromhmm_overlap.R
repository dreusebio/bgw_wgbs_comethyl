#!/usr/bin/env Rscript
# ============================================================
# 13A_core_regulatory_state_overlap.R
#
# Core purpose
#   For each run_root and requested module:
#     1) load Annotated_all_regions
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
# Design goals
#   - no project-specific paths
#   - reproducible
#   - default runs all supported sources
#   - if one source fails, script does not crash
#   - only fails if no requested source is usable
# ============================================================
message("Starting ✓")

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

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

msg <- function(...) cat(sprintf(...), "\n")

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

pick_file <- function(files, pick = c("newest", "oldest")) {
  pick <- match.arg(pick)
  if (length(files) == 0) return(NA_character_)
  info <- file.info(files)
  ord <- order(info$mtime, decreasing = (pick == "newest"))
  files[ord][1]
}

find_newest_file <- function(root, pattern, pick = "newest") {
  files <- list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
  pick_file(files, pick = pick)
}

# ============================================================
# 2) bedtools
# ============================================================
resolve_bedtools <- function(bedtools_path = "") {
  if (!is.na(bedtools_path) && nzchar(bedtools_path)) {
    if (!file.exists(bedtools_path)) stop("Provided --bedtools does not exist: ", bedtools_path)
    return(normalizePath(bedtools_path))
  }
  bt <- Sys.which("bedtools")
  if (bt == "") stop("bedtools not found on PATH. Provide --bedtools /full/path/to/bedtools")
  bt
}

# ============================================================
# 3) read annotated regions
# ============================================================
read_annotated_regions <- function(run_root) {
  candidates <- c(
    file.path(run_root, "enrichment",  "Annotated_all_regions.csv"),
    file.path(run_root, "enrichment",  "Annotated_all_regions.txt"),
    file.path(run_root, "annotations", "Annotated_all_regions.csv"),
    file.path(run_root, "annotations", "Annotated_all_regions.txt")
  )

  file_path <- candidates[file.exists(candidates)][1]
  if (is.na(file_path) || !file.exists(file_path)) {
    file_path <- find_newest_file(run_root, "Annotated_.*regions\\.(csv|txt)$", pick = "newest")
  }
  if (is.na(file_path) || !file.exists(file_path)) {
    stop("Could not find Annotated_all_regions under: ", run_root)
  }

  msg("[ANNOT] Using annotated regions file: %s", file_path)

  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE)
  }

  need <- c("RegionID", "module", "chr", "start", "end")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    stop("Annotated regions missing required columns: ", paste(missing, collapse = ", "))
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
# 4) source registry
# ============================================================
SUPPORTED_SOURCES <- c("roadmap", "encode", "roadmap18")

get_requested_sources <- function(source_arg = "") {
  if (is.null(source_arg) || is.na(source_arg) || !nzchar(source_arg)) {
    return(SUPPORTED_SOURCES)
  }
  src <- split_csv(source_arg)
  bad <- setdiff(src, SUPPORTED_SOURCES)
  if (length(bad) > 0) {
    stop("Unsupported --source value(s): ", paste(bad, collapse = ", "),
         ". Supported: ", paste(SUPPORTED_SOURCES, collapse = ", "))
  }
  unique(src)
}

# ============================================================
# 5) reference metadata
# ============================================================
# NOTE:
# - roadmap 15-state URLs are real pattern from Roadmap bundle
# - roadmap18 and encode are left as configurable placeholders
#   so users can plug in local mirrors / preferred releases without
#   baking project-specific assumptions into the script.

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

# For roadmap18 and encode, use manifest files or explicit local directory content.
# This keeps the script reproducible without hard-coding unstable URLs.
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
# 6) resolve source files
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
# 7) state normalization
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

# For encode cCRE, keep label as-is unless user later adds a mapping table.
normalize_annotation_label <- function(source, raw_state) {
  if (is.na(raw_state) || !nzchar(raw_state)) return(NA_character_)

  raw_state <- as.character(raw_state)

  if (source == "roadmap") {
    if (raw_state %in% names(ROADMAP15_DESC)) return(ROADMAP15_DESC[[raw_state]])
    s2 <- str_replace(raw_state, "^E[0-9]{3}_", "")
    if (s2 %in% names(ROADMAP15_DESC)) return(ROADMAP15_DESC[[s2]])
    return(s2)
  }

  if (source == "roadmap18") {
    return(raw_state)
  }

  if (source == "encode") {
    return(raw_state)
  }

  raw_state
}

# ============================================================
# 8) bedtools intersection
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
# 9) read + collapse overlaps
# ============================================================
guess_overlap_columns <- function(nc) {
  if (nc < 8) stop("Intersect output has too few columns: ", nc)

  a_cols <- c("chrA", "startA", "endA", "region_id")
  remaining <- nc - 4

  if (remaining == 4) {
    b_cols <- c("chrB", "startB", "endB", "state")
  } else {
    # keep last field as annotation label, preserve extras
    extra_n <- remaining - 4
    b_cols <- c("chrB", "startB", "endB",
                paste0("b_extra", seq_len(extra_n)),
                "state")
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
# 10) main args
# ============================================================
project_root <- get_arg("--project_root", "")
region_version <- get_arg("--region_version", "")
if (!nzchar(project_root)) stop("--project_root is required")
if (!nzchar(region_version)) stop("--region_version is required")

setwd(project_root)

base_adj <- file.path(project_root, "comethyl_output", "filter_regions", region_version, "methylation_adjustment")

default_run_roots <- paste(
  file.path(base_adj, "v1_all_pcs"),
  file.path(base_adj, "v2_exclude_outcome_exposure_pcs"),
  file.path(base_adj, "v3_technical_pcs_only"),
  sep = ","
)

run_roots <- split_csv(get_arg("--run_roots", default_run_roots))
modules_in <- unique(split_csv(get_arg("--modules", "")))
if (length(modules_in) == 0) stop("You must provide --modules")

requested_sources <- get_requested_sources(get_arg("--source", ""))

reference_root <- get_arg("--reference_root", file.path(project_root, "reference_data", "regulatory_annotations"))
roadmap18_dir  <- get_arg("--roadmap18_dir", "")
encode_dir     <- get_arg("--encode_dir", "")
download_missing <- as_bool(get_arg("--download_missing", "true"), TRUE)

out_parent <- get_arg("--out_parent", file.path(base_adj, "regulatory_state_overlap"))
membership_col <- get_arg("--membership_col", "membership")

default_bedtools <- "/quobyte/lasallegrp/programs/.conda/NGS_Tools/bin/bedtools"
bedtools_bin <- resolve_bedtools(get_arg("--bedtools", default_bedtools))

msg("project_root: %s", project_root)
msg("region_version: %s", region_version)
msg("requested_sources: %s", paste(requested_sources, collapse = ", "))
msg("download_missing: %s", download_missing)
msg("reference_root: %s", reference_root)
msg("out_parent: %s", out_parent)
msg("bedtools: %s", bedtools_bin)

safe_dir_create(reference_root)
safe_dir_create(out_parent)

# ============================================================
# 11) resolve sources once
# ============================================================
ref_tbl_all <- resolve_requested_reference_files(
  requested_sources = requested_sources,
  reference_root = reference_root,
  download_missing = download_missing,
  roadmap18_dir = roadmap18_dir,
  encode_dir = encode_dir
)

if (nrow(ref_tbl_all) == 0) {
  stop("No reference tracks resolved for requested source(s).")
}

write_csv(ref_tbl_all, file.path(out_parent, "resolved_reference_tracks.csv"))

usable_sources <- ref_tbl_all %>%
  filter(ok) %>%
  count(source, name = "n_tracks")

if (nrow(usable_sources) == 0) {
  stop("None of the requested sources are usable. Check local files and/or internet access.")
}

msg("Usable sources:")
for (i in seq_len(nrow(usable_sources))) {
  msg("  - %s : %d track(s)", usable_sources$source[i], usable_sources$n_tracks[i])
}

# ============================================================
# 12) process run_roots
# ============================================================
for (rr in run_roots) {
  rr <- trimws(rr)
  tag <- basename(rr)

  if (!dir.exists(rr)) {
    msg("[SKIP] run_root not found: %s", rr)
    next
  }

  anno <- read_annotated_regions(rr)

  for (mod in modules_in) {
    msg("\n==================================================")
    msg("[RUN] run_root=%s | module=%s", tag, mod)

    mod_df <- anno %>% filter(module == mod)
    if (nrow(mod_df) == 0) {
      msg("[INFO] No regions found for module '%s' in %s", mod, rr)
      next
    }

    out_dir <- file.path(out_parent, tag, paste0("module_", mod))
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
      next
    }

    bed_file <- file.path(out_dir, paste0("regions_", mod, ".bed"))
    write_tsv(bed, bed_file, col_names = FALSE)

    source_results <- list()

    for (src in requested_sources) {
      ref_tbl <- ref_tbl_all %>% filter(source == src, ok)
      if (nrow(ref_tbl) == 0) {
        msg("[WARN] No usable tracks for source=%s; skipping", src)
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
        next
      }

      dom <- read_and_collapse_overlaps(
        overlap_files = overlap_files,
        source_name = src,
        ref_tbl = ref_tbl
      )

      if (nrow(dom) == 0) {
        msg("[WARN] No collapsed overlaps for source=%s", src)
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
    }

    if (length(source_results) == 0) {
      msg("[WARN] No usable source results for run_root=%s module=%s", tag, mod)
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
  }
}

msg("\n[✓] 13A complete.")