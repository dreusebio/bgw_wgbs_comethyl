#=====================================================================================
#
#  Code chunk 1 - Venn Diagrams
#
#=====================================================================================
#load_data

# --- Venn for 3 lists ---
library(openxlsx)
library(ggvenn)

# --- Define output base directory ---
base_dir <- "/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/DMRichR_overlap_between_all_3"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read annotation files ---
batch1 <- read.xlsx("/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/DMRichR_environmental_exposure_index/DMRs/DMRs_annotated.xlsx")
batch2 <- read.xlsx("/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/rwg_wfl_0to6mo/cytosine_reports/cutoff_005/DMRs/DMRs_annotated.xlsx")
batch3 <- read.xlsx("/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/rwg_wfl_6to12mo/cytosine_reports/cutoff_005/DMRs/DMRs_annotated.xlsx")

# --- Sanity check for column name ---
stopifnot(
  "geneSymbol" %in% names(batch1),
  "geneSymbol" %in% names(batch2),
  "geneSymbol" %in% names(batch3)
)

# --- Prepare venn list ---
venn_list <- list(
  Environmental_exposure = unique(na.omit(batch1$geneSymbol)),
  Rapid_weight_0_6       = unique(na.omit(batch2$geneSymbol)),
  Rapid_weight_6_12      = unique(na.omit(batch3$geneSymbol))
)

# --- Plot ---
p <- ggvenn(
  venn_list,
  fill_color    = c("skyblue", "orchid", "palegreen"),
  stroke_size   = 0.5,
  text_size     = 5,
  set_name_size = 5
)

# --- Save Venn diagram properly inside base_dir ---
out_file <- file.path(base_dir, "Common_between_all_3.pdf")

ggsave(
  filename = out_file,
  plot = p,
  width = 6,
  height = 6
)

# --- Confirmation message ---
cat("✅ Venn diagram successfully saved as:", normalizePath(out_file), "\n")
#=====================================================================================
#
#  Code chunk 2 - #next save geneSymbols
#
#=====================================================================================
#load_data



# pairwise overlaps
pair_01 <- intersect(venn_list[[1]], venn_list[[2]])
pair_02 <- intersect(venn_list[[1]], venn_list[[3]])
pair_12 <- intersect(venn_list[[2]], venn_list[[3]])

# triple overlap
triple <- Reduce(intersect, venn_list)

# --- Convert overlap vectors to data frames with a "geneSymbol" column ---
pair_01_df <- data.frame(geneSymbol = sort(unique(pair_01)))
pair_02_df <- data.frame(geneSymbol = sort(unique(pair_02)))
pair_12_df <- data.frame(geneSymbol = sort(unique(pair_12)))
triple_df  <- data.frame(geneSymbol = sort(unique(triple)))

# --- Save to Excel workbook ---
openxlsx::write.xlsx(
  list(
    Pair_Env_vs_0_6      = pair_01_df,
    Pair_Env_vs_6_12     = pair_02_df,
    Pair_0_6_vs_6_12     = pair_12_df,
    Triple_Env_0_6_6_12  = triple_df
  ),
  file = "Venn_overlaps.xlsx",
  asTable = TRUE
)


# --- Prepare overlap data frames ---
pair_01_df <- data.frame(geneSymbol = sort(unique(pair_01)))
pair_02_df <- data.frame(geneSymbol = sort(unique(pair_02)))
pair_12_df <- data.frame(geneSymbol = sort(unique(pair_12)))
triple_df  <- data.frame(geneSymbol = sort(unique(triple)))

# --- Named list for looping ---
overlap_list <- list(
  Pair_Env_vs_0_6     = pair_01_df,
  Pair_Env_vs_6_12    = pair_02_df,
  Pair_0_6_vs_6_12    = pair_12_df,
  Triple_Env_0_6_6_12 = triple_df
)

# --- Loop through and save each as its own file in its own folder ---
for (name in names(overlap_list)) {
  subfolder <- file.path(base_dir, name)
  dir.create(subfolder, showWarnings = FALSE)
  
  file_path <- file.path(subfolder, paste0(name, ".xlsx"))
  
  openxlsx::write.xlsx(
    overlap_list[[name]],
    file = file_path,
    asTable = TRUE
  )
  
  message("✅ Saved: ", file_path)
}

cat("All overlap Excel files saved successfully!\n")
#=====================================================================================
#
#  Code chunk 3 - Extract information
#
#=====================================================================================
#load_data

#this last I added it so that it could extract all infomation from the original spreadsheets.
# ---- helper: bind data frames with different columns (no dplyr needed) ----
bind_rows_any <- function(dflist) {
  all_cols <- Reduce(union, lapply(dflist, names))
  dflist2 <- lapply(dflist, function(df) {
    missing <- setdiff(all_cols, names(df))
    if (length(missing)) df[missing] <- NA
    df[all_cols]
  })
  do.call(rbind, dflist2)
}


# ---- helper: build labeled table for an overlap set ----
build_overlap_table <- function(overlap_vec, batches, labels, overlap_name) {
  # batches: list of data.frames (e.g., list(batch1, batch2))
  # labels:  character vector same length as batches (e.g., c("Environmental_exposure","Rapid_weight_0_6"))
  stopifnot(length(batches) == length(labels))

  # extract rows by geneSymbol from each batch and label source
  extracted <- Map(function(df, lab) {
    sub <- df[df$geneSymbol %in% overlap_vec, , drop = FALSE]
    if (nrow(sub) == 0) return(sub)       # allow empty
    sub$Source <- lab
    sub$Overlap <- overlap_name
    sub
  }, batches, labels)

  # bind even if columns differ
  out <- bind_rows_any(extracted)
  # move key columns to front if present
  key_front <- c("Overlap", "Source", "geneSymbol")
  front <- intersect(key_front, names(out))
  rest  <- setdiff(names(out), front)
  out[, c(front, rest), drop = FALSE]
}

# =========================
# Build and save per-overlap
# =========================
# --- helper: save an annotated table into its matching subfolder ---
save_annotated_in_subfolder <- function(df, overlap_name, base_dir) {
  subfolder <- file.path(base_dir, overlap_name)
  dir.create(subfolder, showWarnings = FALSE, recursive = TRUE)
  out_path  <- file.path(subfolder, paste0("DMRs_annotated_", overlap_name, ".xlsx"))
  openxlsx::write.xlsx(df, out_path, asTable = TRUE)
  message("✅ Saved annotated table: ", out_path, " (n=", nrow(df), ")")
}

# =========================
# Build and save per-overlap
# =========================

label1 <- "Environmental_exposure"
label2 <- "Rapid_weight_0_6"
label3 <- "Rapid_weight_6_12"

# Pair: Env vs 0–6
tbl_01 <- build_overlap_table(
  overlap_vec = pair_01,
  batches     = list(batch1, batch2),
  labels      = c(label1, label2),
  overlap_name= "Pair_Env_vs_0_6"
)
save_annotated_in_subfolder(tbl_01, "Pair_Env_vs_0_6", base_dir)

# Pair: Env vs 6–12
tbl_02 <- build_overlap_table(
  overlap_vec = pair_02,
  batches     = list(batch1, batch3),
  labels      = c(label1, label3),
  overlap_name= "Pair_Env_vs_6_12"
)
save_annotated_in_subfolder(tbl_02, "Pair_Env_vs_6_12", base_dir)

# Pair: 0–6 vs 6–12
tbl_12 <- build_overlap_table(
  overlap_vec = pair_12,
  batches     = list(batch2, batch3),
  labels      = c(label2, label3),
  overlap_name= "Pair_0_6_vs_6_12"
)
save_annotated_in_subfolder(tbl_12, "Pair_0_6_vs_6_12", base_dir)

# Triple: Env ∩ 0–6 ∩ 6–12
tbl_tri <- build_overlap_table(
  overlap_vec = triple,
  batches     = list(batch1, batch2, batch3),
  labels      = c(label1, label2, label3),
  overlap_name= "Triple_Env_0_6_6_12"
)
save_annotated_in_subfolder(tbl_tri, "Triple_Env_0_6_6_12", base_dir)

cat("🎉 Finished writing annotated overlap files to per-overlap folders in:\n",
    normalizePath(base_dir), "\n")

# # Map your existing objects to readable labels
# label1 <- "Environmental_exposure"
# label2 <- "Rapid_weight_0_6"
# label3 <- "Rapid_weight_6_12"

# # Pair: Env vs 0–6
# tbl_01 <- build_overlap_table(
#   overlap_vec = pair_01,
#   batches     = list(batch1, batch2),
#   labels      = c(label1, label2),
#   overlap_name= "Pair_Env_vs_0_6"
# )
# file_01 <- file.path(base_dir, paste0("DMR_annotated_Pair_Env_vs_0_6.xlsx"))
# openxlsx::write.xlsx(tbl_01, file_01, asTable = TRUE)
# message("✅ Saved: ", file_01)

# # Pair: Env vs 6–12
# tbl_02 <- build_overlap_table(
#   overlap_vec = pair_02,
#   batches     = list(batch1, batch3),
#   labels      = c(label1, label3),
#   overlap_name= "Pair_Env_vs_6_12"
# )
# file_02 <- file.path(base_dir, paste0("DMR_annotated_Pair_Env_vs_6_12.xlsx"))
# openxlsx::write.xlsx(tbl_02, file_02, asTable = TRUE)
# message("✅ Saved: ", file_02)

# # Pair: 0–6 vs 6–12
# tbl_12 <- build_overlap_table(
#   overlap_vec = pair_12,
#   batches     = list(batch2, batch3),
#   labels      = c(label2, label3),
#   overlap_name= "Pair_0_6_vs_6_12"
# )
# file_12 <- file.path(base_dir, paste0("DMR_annotated_Pair_0_6_vs_6_12.xlsx"))
# openxlsx::write.xlsx(tbl_12, file_12, asTable = TRUE)
# message("✅ Saved: ", file_12)

# # Triple: Env ∩ 0–6 ∩ 6–12
# tbl_tri <- build_overlap_table(
#   overlap_vec = triple,
#   batches     = list(batch1, batch2, batch3),
#   labels      = c(label1, label2, label3),
#   overlap_name= "Triple_Env_0_6_6_12"
# )
# file_tri <- file.path(base_dir, paste0("DMR_annotated_Triple_Env_0_6_6_12.xlsx"))
# openxlsx::write.xlsx(tbl_tri, file_tri, asTable = TRUE)
# message("✅ Saved: ", file_tri)

# cat("🎉 Finished writing annotated overlap files to:\n", normalizePath(base_dir), "\n")
