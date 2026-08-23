#!/usr/bin/env Rscript

# ============================================================================
# SYNTHETIC TEST FOR REVISED SCRIPT 05b
# ============================================================================
# This test creates tiny female/male matrices with known SD patterns, invokes
# Script 05b as a real command-line program, and verifies that joint_sd_all:
#   - removes a region variable only in females;
#   - removes a region variable only in males;
#   - removes an incomplete region;
#   - retains regions variable and complete in both datasets;
#   - writes identical row order for both output matrices.
#
# It tests pipeline logic, file naming and argument parsing without WGBS files.
# It does not test comethyl/WGCNA numerical behavior.
# ============================================================================

args <- commandArgs(trailingOnly = FALSE)
this_arg <- grep("^--file=", args, value = TRUE)
if (length(this_arg) != 1) stop("Could not identify test script path")
this_file <- normalizePath(sub("^--file=", "", this_arg))
script_dir <- dirname(this_file)
script05b <- file.path(script_dir, "05b_filter_shared_complete_regions_consensus.R")

test_root <- tempfile("joint_dynamic_test_")
dir.create(test_root, recursive = TRUE)
on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

regions_dir <- file.path(test_root, "input", "cov3_75pct", "covMin10_methSD0")
dir.create(regions_dir, recursive = TRUE)

region_ids <- paste0("R", 1:6)
regions <- data.frame(
  RegionID = region_ids,
  chr = "chr1",
  start = seq(100, 600, by = 100),
  end = seq(150, 650, by = 100)
)
regions_file <- file.path(regions_dir, "Filtered_Regions.txt")
write.table(regions, regions_file, sep = "\t", quote = FALSE, row.names = FALSE)

# R1 and R6 vary in both groups and should survive.
# R2 varies only in females; R3 only in males; R4 is invariant; R5 has an NA.
female <- rbind(
  R1 = c(0.10, 0.20, 0.30, 0.40),
  R2 = c(0.10, 0.20, 0.30, 0.40),
  R3 = c(0.20, 0.20, 0.20, 0.20),
  R4 = c(0.50, 0.50, 0.50, 0.50),
  R5 = c(0.10, NA,   0.30, 0.40),
  R6 = c(0.60, 0.70, 0.80, 0.90)
)
male <- rbind(
  R1 = c(0.15, 0.25, 0.35, 0.45),
  R2 = c(0.20, 0.20, 0.20, 0.20),
  R3 = c(0.10, 0.20, 0.30, 0.40),
  R4 = c(0.50, 0.50, 0.50, 0.50),
  R5 = c(0.10, 0.20, 0.30, 0.40),
  R6 = c(0.55, 0.65, 0.75, 0.85)
)
colnames(female) <- paste0("F", 1:4)
colnames(male) <- paste0("M", 1:4)

female_file <- file.path(test_root, "female.rds")
male_file <- file.path(test_root, "male.rds")
saveRDS(female, female_file)
saveRDS(male, male_file)

cmd_args <- c(
  script05b,
  "--project_root", test_root,
  "--dataset1_label", "females",
  "--dataset1_meth", female_file,
  "--dataset2_label", "males",
  "--dataset2_meth", male_file,
  "--regions_file", regions_file,
  "--selection_mode", "joint_sd_all",
  "--joint_meth_sd", "0.07",
  "--min_regions", "1",
  "--write_sd_table", "TRUE"
)

status <- system2(file.path(R.home("bin"), "Rscript"), cmd_args)
if (status != 0) stop("Script 05b synthetic execution failed with status ", status)

out_dir <- file.path(
  test_root, "comethyl_output", "consensus", "05b_shared_complete_regions",
  "cov3_75pct", "covMin10_methSD0_jointSDall0p07"
)

f_out <- readRDS(file.path(out_dir, "females_Methylation_jointEligible.rds"))
m_out <- readRDS(file.path(out_dir, "males_Methylation_jointEligible.rds"))
ids <- readLines(file.path(out_dir, "joint_eligible_regionIDs.txt"))

stopifnot(identical(ids, c("R1", "R6")))
stopifnot(identical(rownames(f_out), ids))
stopifnot(identical(rownames(m_out), ids))
stopifnot(file.exists(file.path(out_dir, "joint_variability_diagnostics.tsv")))
stopifnot(file.exists(file.path(out_dir, "joint_eligible_regions_4col.tsv")))

message("PASS: joint_sd_all retained the expected regions and identical order.")
