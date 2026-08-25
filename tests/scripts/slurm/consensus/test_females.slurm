#!/bin/bash
#SBATCH --job-name=Test_Comethl_BGW_consensus_females
#SBATCH --mail-user=gekuodza@health.ucdavis.edu
#SBATCH --time=5-00:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --output=00_std_err_logs/comethyl_pearson_report_consensus_females_BGW.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report_consensus_females_BGW.err
#SBATCH --cpus-per-task=30
#SBATCH --account=lasallegrp #genome-center-grp
#SBATCH --mem=250G
#SBATCH --mail-type=ALL
#SBATCH --chdir=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests

###################
# Run Information #
###################

set -euo pipefail

###################
# Run Information #
###################
start=$(date +%s)
hostname
date

####################
# Load Environment #
####################
source ~/.profile

# Pixi settings for cluster stability
export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME

mkdir -p "$PIXI_CACHE_DIR"

# Keep R isolated from user libraries/configs
unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Optional threading hygiene
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Move into Pixi workspace
cd /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/reproducibility/comethyl/pixi

# Sanity check
pixi run Rscript --vanilla -e 'cat("R version:", R.version.string, "\n"); library(comethyl); cat("comethyl OK\n")'

##############
# Run Script #
##############

##############
# Comethyl Scripts#
##############

# ================================================================
# Convenience variables
# ================================================================
PROJECT=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests
SCRIPTS=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/scripts/comethyl_scripts/consensus
DATA=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/data/metadata
CONFIG=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/config
OUT=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/tests/comethyl_output/consensus

# ----------------------------------------------------------------
# REGION SELECTION -- joint-SD redesign (Candidate A)
#
# The decisive SD filter now happens in Script 05b, after region
# methylation exists for BOTH sexes -- not in the reference-only
# Script 02. Script 02 applies coverage only (--meth_sd 0); Script
# 05b applies the joint SD rule across females AND males together.
#
# To try Candidate B instead, change REF_COV_MIN to 14 and rerun
# from Script 02 onward. To sweep joint SD (0.06 / 0.08) without
# rerunning CpG import or region methylation, change JOINT_SD only
# and rerun from Script 05b onward -- Scripts 00-05 stay untouched.
# ----------------------------------------------------------------
CPG_LABEL=cov3_75pct                 # from Script 01 (cov=3, per_sample=0.75) -- unchanged
REF_COV_MIN=4                     # Script 02 coverage prefilter (Candidate A; use 14 for Candidate B)
REF_METH_SD=0                        # Script 02 SD is deliberately permissive now
REF_REGION_LABEL=covMin${REF_COV_MIN}_methSD${REF_METH_SD}   # Script 02's own output label, e.g. covMin10_methSD0
JOINT_SD=0.05                      # Script 05b's joint SD threshold, decided across both sexes
JOINT_SD_LABEL=0p05                  # must match Script 05b's number_label() formatting of JOINT_SD
EFFECTIVE_REGION_LABEL=${REF_REGION_LABEL}_jointSDall${JOINT_SD_LABEL}   # e.g. covMin10_methSD0_jointSDall0p05
                                      # -- this is the region_label every script from 06 onward will derive

# # # # # ================================================================
# # # # # SCRIPT 00: Import CpG reports (one per dataset)
# # # # # UNCHANGED
# # # # # ================================================================

# pixi run Rscript $SCRIPTS/00_import_cpg_reports_consensus.R \
#   --project_root $PROJECT \
#   --dataset_label females \
#   --meta_file $DATA/merged_qc_females.xlsx \
#   --cyto_dir $PROJECT/data/processed/cytosine_reports \
#   --workers 8

# pixi run Rscript $SCRIPTS/00_import_cpg_reports_consensus.R \
#   --project_root $PROJECT \
#   --dataset_label males \
#   --meta_file $DATA/merged_qc_males.xlsx \
#   --cyto_dir $PROJECT/data/processed/cytosine_reports \
#   --workers 8


# # # ================================================================
# # # SCRIPT 01: Filter CpGs and call reference regions (females only)
# # # UNCHANGED
# # # Output: $OUT/01_reference_filter_regions/females/$CPG_LABEL/
# # # ================================================================

# pixi run Rscript $SCRIPTS/01_filter_cpgs_call_regions_consensus.R \
#   --project_root $PROJECT \
#   --dataset_label females \
#   --input_rds $OUT/00_cpg_extraction/females/Unfiltered_BSseq.rds \
#   --cov 3 \
#   --per_sample 0.75

# # # ================================================================
# # # SCRIPT 02: Filter reference regions
# # #
# # # CHANGED: coverage prefilter only. --meth_sd is now 0 (permissive) --
# # # the decisive SD filter has moved to Script 05b so a female-only SD
# # # cut can't remove regions that turn out to be dynamic in males.
# # # Output: $OUT/02_reference_region_filter/females/$CPG_LABEL/$REF_REGION_LABEL/
# # # ================================================================

# pixi run Rscript $SCRIPTS/02_filter_reference_regions_consensus.R \
#   --project_root $PROJECT \
#   --dataset_label females \
#   --input_regions $OUT/01_reference_filter_regions/females/$CPG_LABEL/Regions.txt \
#   --cov_min $REF_COV_MIN \
#   --meth_sd $REF_METH_SD

# # # ================================================================
# # # SCRIPT 03: Build reference region methylation (females only)
# # #
# # # CHANGED: --input_regions now points at the permissive-SD region
# # # set from the updated Script 02 call above.
# # # Output: $OUT/03_region_methylation/females/$CPG_LABEL/$REF_REGION_LABEL/
# # # ================================================================

# pixi run Rscript $SCRIPTS/03_build_reference_region_methylation_consensus.R \
#   --project_root $PROJECT \
#   --dataset_label females \
#   --input_bs $OUT/01_reference_filter_regions/females/$CPG_LABEL/Filtered_BSseq.rds \
#   --input_regions $OUT/02_reference_region_filter/females/$CPG_LABEL/$REF_REGION_LABEL/Filtered_Regions.rds

# # # ================================================================
# # # SCRIPT 04: Align other datasets to reference CpGs
# # # UNCHANGED -- depends on Script 01's CpG-level filtering, not
# # # Script 02's region-level covMin/methSD, so it does not need to
# # # rerun when only the region-selection candidate changes.
# # # Output: $OUT/04_align_to_reference_cpgs/<dataset_label>/
# # # ================================================================

# pixi run Rscript $SCRIPTS/04_align_other_datasets_to_reference_cpgs_consensus.R \
#   --project_root $PROJECT \
#   --reference_bs $OUT/01_reference_filter_regions/females/$CPG_LABEL/Filtered_BSseq.rds \
#   --dataset1_label males \
#   --dataset1_meta $DATA/merged_qc_males.xlsx \
#   --dataset1_cyto_dir $PROJECT/data/processed/cytosine_reports \
#   --cov 3 \
#   --per_sample 0.75 \
#   --workers 8

# # # ================================================================
# # # SCRIPT 05: Build other dataset region methylation
# # #
# # # CHANGED: --reference_regions now points at the permissive-SD
# # # region set from the updated Script 02 call. --dataset1_bs is
# # # unchanged (Script 04 output, not candidate-specific).
# # # Output: $OUT/03_region_methylation/males/$CPG_LABEL/$REF_REGION_LABEL/
# # # ================================================================

# pixi run Rscript $SCRIPTS/05_build_other_dataset_region_methylation_consensus.R \
#   --project_root $PROJECT \
#   --reference_regions $OUT/02_reference_region_filter/females/$CPG_LABEL/$REF_REGION_LABEL/Filtered_Regions.rds \
#   --dataset1_label males \
#   --dataset1_bs $OUT/04_align_to_reference_cpgs/males/Filtered_BSseq.rds

# # # ================================================================
# # # SCRIPT 05b: Filter shared, complete, jointly-variable regions
# # #
# # # CHANGED: this is now where the decisive SD filter happens.
# # # --selection_mode joint_sd_all keeps only regions that are complete
# # # AND have SD >= $JOINT_SD in BOTH females and males, so a female-only
# # # cut can't discard regions that are actually dynamic in males (or
# # # vice versa). Output filenames also changed:
# # #   <dataset>_Methylation_complete.rds -> <dataset>_Methylation_jointEligible.rds
# # # Output: $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/
# # #   females_Methylation_jointEligible.rds
# # #   males_Methylation_jointEligible.rds
# # #   joint_eligible_regions_4col.tsv   <- new --regions_file for Scripts 09-14
# # # ================================================================

# pixi run Rscript $SCRIPTS/05b_filter_shared_complete_regions_consensus.R \
#   --project_root $PROJECT \
#   --dataset1_label females \
#   --dataset1_meth $OUT/03_region_methylation/females/$CPG_LABEL/$REF_REGION_LABEL/Region_Methylation.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/03_region_methylation/males/$CPG_LABEL/$REF_REGION_LABEL/Region_Methylation.rds \
#   --regions_file $OUT/02_reference_region_filter/females/$CPG_LABEL/$REF_REGION_LABEL/Filtered_Regions.txt \
#   --selection_mode joint_sd_all \
#   --joint_meth_sd $JOINT_SD \
#   --min_regions 1000 \
#   --write_sd_table TRUE

# # # ================================================================
# # # SCRIPT 06: PC diagnostics (run separately per dataset)
# # #
# # # CHANGED: --input_meth now points at Script 05b's joint-eligible
# # # matrices. Script 06 auto-detects the "05b_shared_complete_regions"
# # # path and the "_Methylation_jointEligible.rds" filename to derive
# # # dataset_label/cpg_label/region_label -- no other flags change.
# # # Output: $OUT/06_pc_diagnostics/<dataset>/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/
# # # ================================================================

# pixi run Rscript $SCRIPTS/06_pc_diagnostics_consensus.R \
#   --project_root $PROJECT \
#   --input_meth $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/females_Methylation_jointEligible.rds \
#   --sample_info $DATA/merged_qc_females.xlsx \
#   --protected_traits_file $CONFIG/protected_traits.txt \
#   --technical_traits_file $CONFIG/technical_traits.txt \
#   --selected_trait_mode protected_and_technical \
#   --threads 59

# pixi run Rscript $SCRIPTS/06_pc_diagnostics_consensus.R \
#   --project_root $PROJECT \
#   --input_meth $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/males_Methylation_jointEligible.rds \
#   --sample_info $DATA/merged_qc_males.xlsx \
#   --protected_traits_file $CONFIG/protected_traits.txt \
#   --technical_traits_file $CONFIG/technical_traits.txt \
#   --selected_trait_mode protected_and_technical \
#   --threads 59


# # # # ================================================================
# # # # SCRIPT 07: Methylation adjustment (run separately per dataset)
# # # #
# # # # CHANGED: --input_meth points at Script 05b's joint-eligible
# # # # matrices (same auto-detection as Script 06). --input_pcs and
# # # # --pc_trait_stats point at Script 06's new output directory.
# # # # Output: $OUT/07_methylation_adjustment/<dataset>/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<variant>/
# # # # ================================================================

# pixi run Rscript $SCRIPTS/07_methylation_adjustment_consensus.R \
#   --project_root $PROJECT \
#   --input_meth $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/females_Methylation_jointEligible.rds \
#   --input_pcs $OUT/06_pc_diagnostics/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/PCs.rds \
#   --pc_trait_stats $OUT/06_pc_diagnostics/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file $CONFIG/protected_traits.txt \
#   --technical_traits_file $CONFIG/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 FALSE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05

# pixi run Rscript $SCRIPTS/07_methylation_adjustment_consensus.R \
#   --project_root $PROJECT \
#   --input_meth $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/males_Methylation_jointEligible.rds \
#   --input_pcs $OUT/06_pc_diagnostics/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/PCs.rds \
#   --pc_trait_stats $OUT/06_pc_diagnostics/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file $CONFIG/protected_traits.txt \
#   --technical_traits_file $CONFIG/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 FALSE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05


# # ================================================================
# # SCRIPT 08: Soft-threshold power analysis (one run per adjustment version)
# #
# # CHANGED: --dataset*_meth now point at Script 07's new output
# # directories (which already include $EFFECTIVE_REGION_LABEL).
# # Output: $OUT/08_soft_power/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# # ================================================================

# # v1
# pixi run Rscript $SCRIPTS/08_softpower_consensus.R \
#   --project_root $PROJECT \
#   --dataset1_label females \
#   --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
#   --adjustment_version v1_all_pcs \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 59 \
#   --scale_free_cutoff 0.7

# # v2
# pixi run Rscript $SCRIPTS/08_softpower_consensus.R \
#   --project_root $PROJECT \
#   --dataset1_label females \
#   --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 59 \
#   --scale_free_cutoff 0.7


# # ================================================================
# # SCRIPT 08b: Soft-power benchmark (subsampled, fast screening)
# #
# # CHANGED: --dataset*_meth paths updated to the new $EFFECTIVE_REGION_LABEL.
# # Output: $OUT/08b_softpower_benchmark/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/subsample20000/
# # ================================================================

# # v1
# pixi run Rscript $SCRIPTS/08b_softpower_benchmark_consensus.R \
#   --project_root $PROJECT \
#   --dataset1_label females \
#   --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
#   --adjustment_version v1_all_pcs \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --subsample_size 600 \
#   --n_seeds 5 \
#   --threads 8 \
#   --scale_free_cutoff 0.7

# # v2
# pixi run Rscript $SCRIPTS/08b_softpower_benchmark_consensus.R \
#   --project_root $PROJECT \
#   --dataset1_label females \
#   --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --subsample_size 600 \
#   --n_seeds 5 \
#   --threads 8 \
#   --scale_free_cutoff 0.7


# # ================================================================
# # SCRIPT 08c: Compare soft-power benchmarks across region filters
# #
# # CHANGED: --reference_dataset_label region-count lookup now checks
# # Script 05b's joint_eligible_regionIDs.txt first (new joint-dynamic
# # support), falling back to Script 02's Filtered_Regions.txt.
# # Output: $OUT/08c_softpower_benchmark_compare/$CPG_LABEL/
# # ================================================================

# pixi run Rscript $SCRIPTS/08c_softpower_benchmark_compare_consensus.R \
#   --project_root $PROJECT \
#   --cpg_label $CPG_LABEL \
#   --reference_dataset_label females \
#   --scale_free_cutoff 0.7

# # ================================================================
# # SCRIPT 09: Consensus module detection (one run per adjustment version)
# #
# # CHANGED: --regions_file now points at Script 05b's
# # joint_eligible_regions_4col.tsv, NOT Script 02's Filtered_Regions.txt.
# # Script 02's table is deliberately permissive now (SD=0) and still
# # contains regions Script 05b's joint SD filter removed -- passing it
# # here would silently reintroduce those regions.
# # --dataset*_meth and --chosen_power_file updated to $EFFECTIVE_REGION_LABEL.
# # Output: $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# # ================================================================

# # v1
# pixi run Rscript $SCRIPTS/09_modules_consensus.R \
#   --project_root $PROJECT \
#   --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
#   --dataset1_label females \
#   --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
#   --adjustment_version v1_all_pcs \
#   --chosen_power_file $OUT/08_soft_power/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/chosen_power.txt \
#   --consensus_cor pearson \
#   --network_type signed \
#   --tom_type signed \
#   --deep_split 1 \
#   --min_module_size 5 \
#   --merge_cut_height 0.2 \
#   --max_block_size 9500 \
#   --network_calibration "single quantile" \
#   --save_consensus_toms TRUE \
#   --threads 20 \
#   --max_p_outliers 0.1

# # v2
# pixi run Rscript $SCRIPTS/09_modules_consensus.R \
#   --project_root $PROJECT \
#   --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
#   --dataset1_label females \
#   --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label males \
#   --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --chosen_power_file $OUT/08_soft_power/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/chosen_power.txt \
#   --consensus_cor pearson \
#   --network_type signed \
#   --tom_type signed \
#   --deep_split 1 \
#   --min_module_size 5 \
#   --merge_cut_height 0.2 \
#   --max_block_size 9500 \
#   --network_calibration "single quantile" \
#   --save_consensus_toms TRUE \
#   --threads 20 \
#   --max_p_outliers 0.1

# ================================================================
# SCRIPT 10: Consensus diagnostics
#
# CHANGED: --consensus_modules_rds and --regions_file updated to
# $EFFECTIVE_REGION_LABEL / Script 05b's joint region table.
# Output: $OUT/10_consensus_diagnostics/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript $SCRIPTS/10_diagnostics_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v1_all_pcs \
  --consensus_cor pearson \
  --max_p_outliers 0.1
# v2
pixi run Rscript $SCRIPTS/10_diagnostics_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --consensus_cor pearson \
  --max_p_outliers 0.1

# ================================================================
# SCRIPT 11: Consensus module membership
#
# CHANGED: --consensus_modules_rds, --regions_file, --dataset*_meth
# updated to $EFFECTIVE_REGION_LABEL / Script 05b's joint region table.
# Output: $OUT/11_consensus_membership/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript $SCRIPTS/11_membership_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --dataset1_label females \
  --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
  --dataset2_label males \
  --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
  --adjustment_version v1_all_pcs \
  --membership_cor bicor \
  --max_p_outliers 0.1
# v2
pixi run Rscript $SCRIPTS/11_membership_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --dataset1_label females \
  --dataset1_meth $OUT/07_methylation_adjustment/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label males \
  --dataset2_meth $OUT/07_methylation_adjustment/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --membership_cor bicor \
  --max_p_outliers 0.1

# ================================================================
# SCRIPT 12A: ME-trait analysis
#
# CHANGED (inferred -- 12a not re-uploaded this round, applying the
# same regions_file/path convention as 09-11 since the old slurm
# already passed --regions_file here. Verify against 12a's actual
# derivation code once you upload the updated script.)
# Output: $OUT/12a_me_trait_analysis/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript $SCRIPTS/12a_me_trait_analysis_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v1_all_pcs \
  --dataset1_label females \
  --dataset1_sample_info $DATA/merged_qc_females.xlsx \
  --dataset2_label males \
  --dataset2_sample_info $DATA/merged_qc_males.xlsx \
  --run_bicor TRUE \
  --run_pearson FALSE \
  --p_thresh 0.05 \
  --top_n 50 \
  --module_dendro_distance bicor \
  --max_p_outliers 0.1

# v2
pixi run Rscript $SCRIPTS/12a_me_trait_analysis_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --dataset1_label females \
  --dataset1_sample_info $DATA/merged_qc_females.xlsx \
  --dataset2_label males \
  --dataset2_sample_info $DATA/merged_qc_males.xlsx \
  --run_bicor TRUE \
  --run_pearson FALSE \
  --p_thresh 0.05 \
  --top_n 50 \
  --module_dendro_distance bicor \
  --max_p_outliers 0.1

# ================================================================
# SCRIPT 12B: ME-trait presentation heatmaps
#
# CHANGED (inferred, same caveat as 12a above)
# Output: $OUT/12b_me_trait_presentation/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript $SCRIPTS/12b_me_trait_presentation_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v1_all_pcs \
  --dataset1_label females \
  --dataset1_stats_file $OUT/12a_me_trait_analysis/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file $OUT/12a_me_trait_analysis/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir $CONFIG/presentation_sets \
  --module_dendro_distance bicor \
  --p_thresh 0.05 \
  --top_n 50 \
  --full_width 12 \
  --full_height 7 \
  --top_width 9 \
  --top_height 6

# v2
pixi run Rscript $SCRIPTS/12b_me_trait_presentation_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --dataset1_label females \
  --dataset1_stats_file $OUT/12a_me_trait_analysis/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file $OUT/12a_me_trait_analysis/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir $CONFIG/presentation_sets \
  --module_dendro_distance bicor \
  --p_thresh 0.05 \
  --top_n 50 \
  --full_width 12 \
  --full_height 7 \
  --top_width 9 \
  --top_height 6

# ================================================================
# SCRIPT 13: Circos plots and module summaries
#
# CHANGED (inferred, same caveat as 12a above)
# Output: $OUT/13_circos_and_module_summaries/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript $SCRIPTS/13_circos_plots_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Modules.rds \
  --consensus_regions_tsv $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Region_Assignments.tsv \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v1_all_pcs \
  --genome hg38

# v2
pixi run Rscript $SCRIPTS/13_circos_plots_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --consensus_regions_tsv $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Region_Assignments.tsv \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --genome hg38

# ================================================================
# SCRIPT 14: ME-trait pair plots
#
# CHANGED (inferred, same caveat as 12a above)
# Output: $OUT/14_me_trait_pair_plots/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript $SCRIPTS/14_me_trait_pair_plots_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v1_all_pcs \
  --dataset1_label females \
  --dataset1_sample_info $DATA/merged_qc_females.xlsx \
  --dataset1_stats_file $OUT/12a_me_trait_analysis/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_sample_info $DATA/merged_qc_males.xlsx \
  --dataset2_stats_file $OUT/12a_me_trait_analysis/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir $CONFIG/presentation_sets \
  --cat_plot_style violin \
  --only_significant TRUE \
  --p_thresh 0.05

# v2
pixi run Rscript $SCRIPTS/14_me_trait_pair_plots_consensus.R \
  --project_root $PROJECT \
  --consensus_modules_rds $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file $OUT/05b_shared_complete_regions/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/joint_eligible_regions_4col.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --dataset1_label females \
  --dataset1_sample_info $DATA/merged_qc_females.xlsx \
  --dataset1_stats_file $OUT/12a_me_trait_analysis/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_sample_info $DATA/merged_qc_males.xlsx \
  --dataset2_stats_file $OUT/12a_me_trait_analysis/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir $CONFIG/presentation_sets \
  --cat_plot_style violin \
  --only_significant TRUE \
  --p_thresh 0.05

# ================================================================
# SCRIPT 15: Module annotation
#
# CHANGED (inferred -- 15 not re-uploaded this round; this script
# never took --regions_file, so only --consensus_regions_tsv changes.
# Verify against 15's actual derivation code once you upload it.)
# Output: $OUT/15_module_annotation/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript $SCRIPTS/15_Annotate_modules_consensus.R \
  --project_root $PROJECT \
  --consensus_regions_tsv $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/Consensus_Region_Assignments.tsv \
  --adjustment_version v1_all_pcs \
  --genome hg38 \
  --annotation_mode auto

# v2
pixi run Rscript $SCRIPTS/15_Annotate_modules_consensus.R \
  --project_root $PROJECT \
  --consensus_regions_tsv $OUT/09_consensus_modules/shared/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/Consensus_Region_Assignments.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --genome hg38 \
  --annotation_mode auto

# ================================================================
# SCRIPT 16: Module enrichment
#
# CHANGED (inferred, same caveat as Script 15 above)
# Output: $OUT/16_enrichment/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript $SCRIPTS/16_enrichr_selected_modules.R \
  --project_root $PROJECT \
  --annotation_dir $OUT/15_module_annotation/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs \
  --adjustment_version v1_all_pcs \
  --trait_list_file $CONFIG/presentation_sets \
  --dataset1_label females \
  --dataset1_stats_file $OUT/12a_me_trait_analysis/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file $OUT/12a_me_trait_analysis/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --p_thresh 0.05 \
  --use_enrichr_background FALSE \
  --enrichr_retries 5 \
  --enrichr_sleep 2 \
  --do_kegg TRUE \
  --do_go FALSE \
  --do_reactome FALSE \
  --do_enrichr TRUE \
  --plot_top_n 10 \
  --plot_p_col pvalue \
  --plot_sig_cutoff 0.05

# v2
pixi run Rscript $SCRIPTS/16_enrichr_selected_modules.R \
  --project_root $PROJECT \
  --annotation_dir $OUT/15_module_annotation/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs \
  --adjustment_version v2_exclude_protected_pcs \
  --trait_list_file $CONFIG/presentation_sets \
  --dataset1_label females \
  --dataset1_stats_file $OUT/12a_me_trait_analysis/females/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file $OUT/12a_me_trait_analysis/males/$CPG_LABEL/$EFFECTIVE_REGION_LABEL/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --p_thresh 0.05 \
  --use_enrichr_background FALSE \
  --enrichr_retries 5 \
  --enrichr_sleep 2 \
  --do_kegg TRUE \
  --do_go FALSE \
  --do_reactome FALSE \
  --do_enrichr TRUE \
  --plot_top_n 10 \
  --plot_p_col pvalue \
  --plot_sig_cutoff 0.05

####################
# End              #
####################

end=`date +%s`
runtime=$((end-start))
echo "Runtime: $runtime seconds"
