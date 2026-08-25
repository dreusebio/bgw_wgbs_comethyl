#!/bin/bash
#SBATCH --job-name=Comethl_BGW_consensus_females
#SBATCH --mail-user=gekuodza@health.ucdavis.edu
#SBATCH --time=10-00:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --output=00_std_err_logs/comethyl_pearson_report_consensus_females_BGW.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report_consensus_females_BGW.err
#SBATCH --cpus-per-task=45
#SBATCH --account=lasallegrp #genome-center-grp
#SBATCH --mem=700G
#SBATCH --mail-type=ALL
#SBATCH --chdir=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria

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
PROJECT=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria
SCRIPTS=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus
DATA=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata
CONFIG=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config
OUT=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus

# Reference filter labels (derived from scripts 01 and 02)
CPG_LABEL=cov3_75pct
# REGION_LABEL=covMin18_methSD0p07 #covMin18_methSD0p07 covMin10_methSD0p09
JOINT_SD=0.07
# # # ================================================================
# # # SCRIPT 00: Import CpG reports (one per dataset)
# # # ================================================================

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset_label females \
  --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
  --cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/cytosine_reports \
  --workers 8

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset_label males \
  --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/cytosine_reports \
  --workers 8


# # ================================================================
# # SCRIPT 01: Filter CpGs and call reference regions (females only)
# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/females/cov3_75pct/
# # ================================================================

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/01_filter_cpgs_call_regions_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset_label females \
  --input_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/00_cpg_extraction/females/Unfiltered_BSseq.rds \
  --cov 3 \
  --per_sample 0.75

# # ================================================================
# # SCRIPT 02: Filter reference regions
# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/
# # ================================================================

# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/02_filter_reference_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label females \
#   --input_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/females/cov3_75pct/Regions.txt \
#   --cov_min 18 \
#   --meth_sd 0.07

# # ================================================================
# # Script 02 -- permissive reference SD (females only, as always)
# # ================================================================
pixi run Rscript $SCRIPTS/02_filter_reference_regions_consensus.R \
  --project_root $PROJECT \
  --dataset_label females \
  --input_regions $PROJECT/comethyl_output/consensus/01_reference_filter_regions/females/$CPG_LABEL/Regions.txt \
  --cov_min 10 \
  --meth_sd 0
# # ================================================================
# # SCRIPT 03: Build reference region methylation (females only)
# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/03_region_methylation/females/cov3_75pct/covMin18_methSD0p07/
# # ================================================================

# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/03_build_reference_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label females \
#   --input_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/females/cov3_75pct/Filtered_BSseq.rds \
#   --input_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.rds

# Script 03 -- females region methylation for this candidate's regions
pixi run Rscript $SCRIPTS/03_build_reference_region_methylation_consensus.R \
  --project_root $PROJECT \
  --dataset_label females \
  --input_bs $PROJECT/comethyl_output/consensus/01_reference_filter_regions/females/$CPG_LABEL/Filtered_BSseq.rds \
  --input_regions $PROJECT/comethyl_output/consensus/02_reference_region_filter/females/$CPG_LABEL/covMin10_methSD0/Filtered_Regions.rds

# # ================================================================
# # SCRIPT 04: Align other datasets to reference CpGs
# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/04_align_to_reference_cpgs/<dataset_label>/
# # ================================================================

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/04_align_other_datasets_to_reference_cpgs_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --reference_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/females/cov3_75pct/Filtered_BSseq.rds \
  --dataset1_label males \
  --dataset1_meta /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --dataset1_cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/cytosine_reports \
  --cov 3 \
  --per_sample 0.75 \
  --workers 8

# # ================================================================
# # SCRIPT 05: Build other dataset region methylation
# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/03_region_methylation/<dataset>/cov3_75pct/covMin18_methSD0p07/
# # ================================================================

# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/05_build_other_dataset_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --reference_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.rds \
#   --dataset1_label males \
#   --dataset1_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/04_align_to_reference_cpgs/males/Filtered_BSseq.rds \

# Script 05 -- males region methylation, same region definitions (04 unchanged, reused)
pixi run Rscript $SCRIPTS/05_build_other_dataset_region_methylation_consensus.R \
  --project_root $PROJECT \
  --reference_regions $PROJECT/comethyl_output/consensus/02_reference_region_filter/females/$CPG_LABEL/covMin10_methSD0/Filtered_Regions.rds \
  --dataset1_label males \
  --dataset1_bs $PROJECT/comethyl_output/consensus/04_align_to_reference_cpgs/males/Filtered_BSseq.rds

# # ================================================================
# # SCRIPT 05b: Filter shared complete regions
# #
# # FIX: --regions_file must point to Filtered_Regions.txt from script 02,
# #      NOT a file inside 05b output. This drives cpg_label/region_label
# #      derivation for the output path:
# #      /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin18_methSD0p07/
# # ================================================================

# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/05b_filter_shared_complete_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset1_label females \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/03_region_methylation/females/cov3_75pct/covMin18_methSD0p07/Region_Methylation.rds \
#   --dataset2_label males \
#   --dataset2_meth  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/03_region_methylation/males/cov3_75pct/covMin18_methSD0p07/Region_Methylation.rds \
#   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt

# Script 05b -- joint SD filter, decided across both sexes
pixi run Rscript $SCRIPTS/05b_filter_shared_complete_regions_consensus.R \
  --project_root $PROJECT \
  --dataset1_label females \
  --dataset1_meth $PROJECT/comethyl_output/consensus/03_region_methylation/females/$CPG_LABEL/covMin10_methSD0/Region_Methylation.rds \
  --dataset2_label males \
  --dataset2_meth $PROJECT/comethyl_output/consensus/03_region_methylation/males/$CPG_LABEL/covMin10_methSD0/Region_Methylation.rds \
  --regions_file $PROJECT/comethyl_output/consensus/02_reference_region_filter/females/$CPG_LABEL/covMin10_methSD0/Filtered_Regions.txt \
  --selection_mode joint_sd_all \
  --joint_meth_sd $JOINT_SD \
  --min_regions 1000 \
  --write_sd_table TRUE
  
# # ================================================================
# # SCRIPT 06: PC diagnostics (run separately per dataset)
# #
# # FIX: --input_meth now includes cov3_75pct/covMin18_methSD0p07 in path
# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/<dataset>/cov3_75pct/covMin18_methSD0p07/
# # ================================================================

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin18_methSD0p07/females_Methylation_complete.rds \
  --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
  --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
  --selected_trait_mode protected_and_technical \
  --threads 59

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin18_methSD0p07/males_Methylation_complete.rds \
  --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
  --selected_trait_mode protected_and_technical \
  --threads 59


# # # ================================================================
# # # SCRIPT 07: Methylation adjustment (run separately per dataset)
# # #
# # # FIX: --input_meth updated to new 05b path (includes cpg/region levels)
# # # FIX: --input_pcs and --pc_trait_stats now include cov3_75pct/covMin18_methSD0p07
# # # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/<dataset>/cov3_75pct/covMin18_methSD0p07/<variant>/
# # # ================================================================

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin18_methSD0p07/females_Methylation_complete.rds \
  --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/females/cov3_75pct/covMin18_methSD0p07/PCs.rds \
  --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/females/cov3_75pct/covMin18_methSD0p07/PC_Trait_Correlation_Stats_Bicor.tsv \
  --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
  --run_v1 TRUE \
  --run_v2 TRUE \
  --run_v3 FALSE \
  --v2_cor_method bicor \
  --v2_p_thresh 0.05

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/cov3_75pct/covMin18_methSD0p07/males_Methylation_complete.rds \
  --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/males/cov3_75pct/covMin18_methSD0p07/PCs.rds \
  --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/males/cov3_75pct/covMin18_methSD0p07/PC_Trait_Correlation_Stats_Bicor.tsv \
  --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
  --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
  --run_v1 TRUE \
  --run_v2 TRUE \
  --run_v3 FALSE \
  --v2_cor_method bicor \
  --v2_p_thresh 0.05


# ================================================================
# SCRIPT 08: Soft-threshold power analysis (one run per adjustment version)
#
# FIX: --dataset*_meth now includes cov3_75pct/covMin18_methSD0p07 levels
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
  --adjustment_version v1_all_pcs \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 20 \
  --threads 59 \
  --scale_free_cutoff 0.8

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 20 \
  --threads 59 \
  --scale_free_cutoff 0.8

# # v3
# # pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --dataset1_label females \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v3_technical_pcs_only/Males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --dataset2_label males \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v3_technical_pcs_only/males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --adjustment_version v3_technical_pcs_only \
# #   --softpower_cor pearson \
# #   --power_min 1 \
# #   --power_max 20 \
# #   --threads 59 \
# #   --scale_free_cutoff 0.8


# ================================================================
# SCRIPT 08b: Soft-power benchmark (subsampled, fast screening)
#
# Run this BEFORE committing to a full Script 08 run on a new
# covMin/methSD combo, or any time you want to sanity-check a fit
# without risking the OOM kills seen on the full 765k-region matrix.
# Fewer threads than Script 08 on purpose -- more workers means less
# memory per worker, which is what triggered the OOM kills before.
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08b_softpower_benchmark/shared/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/subsample20000/
# ================================================================

# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08b_softpower_benchmark_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
  --adjustment_version v1_all_pcs \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 20 \
  --subsample_size 20000 \
  --n_seeds 5 \
  --threads 8 \
  --scale_free_cutoff 0.8

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08b_softpower_benchmark_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 20 \
  --subsample_size 20000 \
  --n_seeds 5 \
  --threads 8 \
  --scale_free_cutoff 0.8


# ================================================================
# SCRIPT 08c: Compare soft-power benchmarks across region filters
#
# Only compares combos that have already been run through 08b above.
# With only covMin18_methSD0p07 run so far this will show one combo;
# add more --cov_min/--meth_sd sweeps through Scripts 02-08b (same
# pattern as above, new version_label each time) and rerun 08c to
# get a real ranked comparison across region filters.
#
# --reference_dataset_label is set to females explicitly because
# Script 02 (region filtering) was only run on the reference dataset
# (females) in this pipeline -- males has no 02_reference_region_filter
# output to look up region counts from.
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08c_softpower_benchmark_compare/cov3_75pct/
# ================================================================

pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08c_softpower_benchmark_compare_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --cpg_label cov3_75pct \
  --reference_dataset_label females \
  --scale_free_cutoff 0.8

# ================================================================
# SCRIPT 09: Consensus module detection (one run per adjustment version)
#
# FIX: --regions_file now points to Filtered_Regions.txt from script 02
#      (NOT shared_complete_regions_4col.tsv from 05b)
# FIX: --dataset*_meth now includes cov3_75pct/covMin18_methSD0p07 levels
# FIX: --chosen_power_file now includes cov3_75pct/covMin18_methSD0p07 levels
# FIX: file is named chosen_power.txt (not chosen_power_pearson.txt)
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
  --adjustment_version v1_all_pcs \
  --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/chosen_power.txt \
  --consensus_cor pearson \
  --network_type signed \
  --tom_type signed \
  --deep_split 4 \
  --min_module_size 50 \
  --merge_cut_height 0.2 \
  --max_block_size 40000 \
  --network_calibration "single quantile" \
  --save_consensus_toms TRUE \
  --threads 20 \
  --max_p_outliers 0.1

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/chosen_power.txt \
  --consensus_cor pearson \
  --network_type signed \
  --tom_type signed \
  --deep_split 4 \
  --min_module_size 50 \
  --merge_cut_height 0.2 \
  --max_block_size 40000 \
  --network_calibration "single quantile" \
  --save_consensus_toms TRUE \
  --threads 20 \
  --max_p_outliers 0.1

# # #v3
# # pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
# #   --dataset1_label females \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v3_technical_pcs_only/Males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --dataset2_label males \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v3_technical_pcs_only/males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --adjustment_version v3_technical_pcs_only \
# #   --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/shared/cov3_75pct/covMin18_methSD0p07/v3_technical_pcs_only/chosen_power.txt \
# #   --consensus_cor pearson \
# #   --network_type signed \
# #   --tom_type signed \
# #   --deep_split 4 \
# #   --min_module_size 50 \
# #   --merge_cut_height 0.1 \
# #   --max_block_size 50000 \
# #   --network_calibration "single quantile" \
# #   --save_consensus_toms TRUE \
# #   --threads 59 \
# #   --max_p_outliers 0.1

# ================================================================
# SCRIPT 10: Consensus diagnostics
#
# FIX: --consensus_modules_rds now includes cov3_75pct/covMin18_methSD0p07
# FIX: added --regions_file (now required for output path derivation)
# FIX: --adjustment_version was v2_exclude_outcome_exposure_pcs — corrected
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/10_consensus_diagnostics/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/10_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v1_all_pcs \
  --consensus_cor pearson \
  --max_p_outliers 0.1
# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/10_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v2_exclude_protected_pcs \
  --consensus_cor pearson \
  --max_p_outliers 0.1

# ================================================================
# SCRIPT 11: Consensus module membership
#
# FIX: --consensus_modules_rds now includes cov3_75pct/covMin18_methSD0p07
# FIX: added --regions_file (now required)
# FIX: --dataset*_meth now includes cov3_75pct/covMin18_methSD0p07 levels
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/11_consensus_membership/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================

 
# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/11_membership_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
  --adjustment_version v1_all_pcs \
  --membership_cor bicor \
  --max_p_outliers 0.1
# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/11_membership_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --dataset1_label females \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label males \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --membership_cor bicor \
  --max_p_outliers 0.1

# ================================================================
# SCRIPT 12A: ME-trait analysis
#
# FIX: --consensus_modules_rds now includes cov3_75pct/covMin18_methSD0p07
# FIX: added --regions_file (now required)
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/12a_me_trait_analysis_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v1_all_pcs \
  --dataset1_label females \
  --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
  --dataset2_label males \
  --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --run_bicor TRUE \
  --run_pearson FALSE \
  --p_thresh 0.05 \
  --top_n 50 \
  --module_dendro_distance bicor \
  --max_p_outliers 0.1

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/12a_me_trait_analysis_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v2_exclude_protected_pcs \
  --dataset1_label females \
  --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
  --dataset2_label males \
  --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --run_bicor TRUE \
  --run_pearson FALSE \
  --p_thresh 0.05 \
  --top_n 50 \
  --module_dendro_distance bicor \
  --max_p_outliers 0.1

# ================================================================
# SCRIPT 12B: ME-trait presentation heatmaps
#
# FIX: --consensus_modules_rds now includes cov3_75pct/covMin18_methSD0p07
# FIX: added --regions_file (now required)
# FIX: --dataset*_stats_file now includes cov3_75pct/covMin18_methSD0p07 levels
# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12b_me_trait_presentation/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/12b_me_trait_presentation_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v1_all_pcs \
  --dataset1_label females \
  --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
  --module_dendro_distance bicor \
  --p_thresh 0.05 \
  --top_n 50 \
  --full_width 12 \
  --full_height 7 \
  --top_width 9 \
  --top_height 6

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/12b_me_trait_presentation_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v2_exclude_protected_pcs \
  --dataset1_label females \
  --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
  --module_dendro_distance bicor \
  --p_thresh 0.05 \
  --top_n 50 \
  --full_width 12 \
  --full_height 7 \
  --top_width 9 \
  --top_height 6

# ================================================================
# SCRIPT 13: Circos plots and module summaries

# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/13_circos_and_module_summaries/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================
 
# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/13_circos_plots_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Modules.rds \
  --consensus_regions_tsv /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Region_Assignments.tsv \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v1_all_pcs \
  --genome hg38

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/13_circos_plots_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --consensus_regions_tsv /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Region_Assignments.tsv \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v2_exclude_protected_pcs \
  --genome hg38

# ================================================================
# SCRIPT 14: ME-trait pair plots

# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/14_me_trait_pair_plots/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/14_me_trait_pair_plots_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v1_all_pcs \
  --dataset1_label females \
  --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
  --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
  --cat_plot_style violin \
  --only_significant TRUE \
  --p_thresh 0.05

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/14_me_trait_pair_plots_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Modules.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/females/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt \
  --adjustment_version v2_exclude_protected_pcs \
  --dataset1_label females \
  --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
  --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
  --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
  --cat_plot_style violin \
  --only_significant TRUE \
  --p_thresh 0.05

# ================================================================
# SCRIPT 15: Module annotation

# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/15_module_annotation/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================

# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/15_Annotate_modules_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_regions_tsv /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/Consensus_Region_Assignments.tsv \
  --adjustment_version v1_all_pcs \
  --genome hg38 \
  --annotation_mode auto

# v2
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/15_Annotate_modules_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --consensus_regions_tsv /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/shared/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/Consensus_Region_Assignments.tsv \
  --adjustment_version v2_exclude_protected_pcs \
  --genome hg38 \
  --annotation_mode auto

# ================================================================
# SCRIPT 16: Module enrichment

# Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/16_enrichment/cov3_75pct/covMin18_methSD0p07/<adjustment_version>/
# ================================================================
# v1
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/16_enrichr_selected_modules.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --annotation_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/15_module_annotation/cov3_75pct/covMin18_methSD0p07/v1_all_pcs \
  --adjustment_version v1_all_pcs \
  --trait_list_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
  --dataset1_label females \
  --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/females/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/males/cov3_75pct/covMin18_methSD0p07/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
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
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/16_enrichr_selected_modules.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --annotation_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/15_module_annotation/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs \
  --adjustment_version v2_exclude_protected_pcs \
  --trait_list_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
  --dataset1_label females \
  --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/females/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --dataset2_label males \
  --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/males/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
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

# # ================================================================
# # SCRIPT 17: Regulatory state overlap (requires bedtools)

# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/17_regulatory_overlap/<adjustment_version>/
# # NOTE: update --modules to match actual module colors from script 09
# # ================================================================

# # # v2
# # pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/17_core_chromatin_overlap.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --annotation_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/15_module_annotation/cov3_75pct/covMin18_methSD0p07/v2_exclude_protected_pcs \
# #   --modules turquoise,blue,brown \
# #   --reference_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/reference/chromatin \
# #   --adjustment_version v2_exclude_protected_pcs \
# #   --source roadmap \
# #   --download_missing TRUE

# # ================================================================
# # SCRIPT 18: Regulatory state presentation

# # Output: /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/18_regulatory_presentation/<adjustment_version>/
# # NOTE: update --modules to match actual module colors from script 09
# # ================================================================

# # # v2
# # pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/18_plot_chromatin_heatmaps.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --regulatory_overlap_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/17_regulatory_overlap/v2_exclude_protected_pcs \
# #   --modules turquoise,blue,brown \
# #   --adjustment_version v2_exclude_protected_pcs \
# #   --source roadmap \
# #   --max_regions 50 \
# #   --make_per_module_heatmaps TRUE \
# #   --make_combined_heatmap TRUE \
# #   --make_barplots TRUE

# ####################
# # End              #
# ####################

# end=`date +%s`
# runtime=$((end-start))
# echo "Runtime: $runtime seconds"

# #!/bin/bash
# #SBATCH --job-name=Cometh_growell_test
# #SBATCH --mail-user=gekuodza@ucdavis.edu
# #SBATCH --time=14-00:00:00
# #SBATCH --partition=high
# #SBATCH --ntasks=1
# #SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
# #SBATCH --error=00_std_err_logs/comethyl_pearson_report.err                   	
# #SBATCH --cpus-per-task=59  #use 30 or 50 cores for certain things
# #SBATCH --account=lasallegrp #genome-center-grp
# #SBATCH --mem=950G #Need up to 500G for certain tasks
# #SBATCH --mail-type=ALL  
# #SBATCH --chdir=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria

# ###################
# # Run Information #
# ###################

# start=`date +%s`

# hostname

# cd /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/reproducibility/pixi
# pwd

# ####################
# # Load Environment #
# ####################

# export PIXI_HOME=/quobyte/lasallegrp/home/George/.pixi
# export PATH="$PIXI_HOME/bin:$PATH"
# export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
# mkdir -p "$PIXI_CACHE_DIR"


# pixi run Rscript -e 'cat("R home: ", R.home(), "\n"); sessionInfo()'

# ##############
# # Comethyl Scripts#
# ##############

# # Script 00
# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label males \
#   --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
#   --cyto_dir  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/males \
#   --workers 8

#  pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label females \
#   --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/females \
#   --workers 8

# pixi run Rscript  scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label Postpartum \
#   --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/sample_info_postpartum.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/Postpartum \
#   --workers 8

# # Script 01
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/01_filter_cpgs_call_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label males \
#   --input_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/00_cpg_extraction/males/Unfiltered_BSseq.rds \
#   --cov 3 \
#   --per_sample 0.75

# # Script 02
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/02_filter_reference_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label males \
#   --input_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/males/cov3_75pct/Regions.txt \
#   --cov_min 10 \
#   --meth_sd 0.09

# # Script 03
# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/03_build_reference_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset_label males \
#   --input_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/males/cov3_75pct/Filtered_BSseq.rds \
#   --input_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/males/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.rds

# # Script 04
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/04_align_other_datasets_to_reference_cpgs_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --reference_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/01_reference_filter_regions/males/cov3_75pct/Filtered_BSseq.rds \
#   --dataset1_label females \
#   --dataset1_meta /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
#   --dataset1_cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/females \
#   --dataset2_label Postpartum \
#   --dataset2_meta  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/sample_info_postpartum.xlsx \
#   --dataset2_cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/processed/Postpartum \
#   --cov 1 \
#   --per_sample 0.01 \
#   --workers 8  

# # Script 05
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/05_build_other_dataset_region_methylation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --reference_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/males/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.rds \
#   --dataset1_label females \
#   --dataset1_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/04_align_to_reference_cpgs/females/Filtered_BSseq.rds \
#   --dataset2_label Postpartum \
#   --dataset2_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/04_align_to_reference_cpgs/Postpartum/Filtered_BSseq.rds

# # Script 05b
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/05b_filter_shared_complete_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset1_label males \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/03_region_methylation/males/cov3_75pct/covMin18_methSD0p07/Region_Methylation.rds \
#   --dataset2_label females \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/03_region_methylation/females/females/covMin18_methSD0p07/Region_Methylation.rds \
#   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/02_reference_region_filter/males/cov3_75pct/covMin18_methSD0p07/Filtered_Regions.txt

# # Script 06

# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/males_Methylation_complete.rds \
#   --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
#   --selected_trait_mode   protected_and_technical

# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/females_Methylation_complete.rds \
#   --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
#   --selected_trait_mode   protected_and_technical
  
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/Postpartum_Methylation_complete.rds \
#   --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/sample_info_postpartum.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
#   --selected_trait_mode   protected_and_technical

# #Script 07
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/males_Methylation_complete.rds \
#   --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/males/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/males/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 TRUE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05

# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/females_Methylation_complete.rds \
#   --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/females/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/females/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 TRUE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05

# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/Postpartum_Methylation_complete.rds \
#   --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/Postpartum/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/06_pc_diagnostics/Postpartum/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 TRUE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05

# # # script 08
# # #v1
# # pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --dataset1_label males \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
# #   --dataset2_label females \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
# #   --adjustment_version v1_all_pcs \
# #   --softpower_cor pearson \
# #   --power_min 1 \
# #   --power_max 20 \
# #   --threads 59 \
# #   --scale_free_cutoff 0.8

# #v2
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --dataset1_label males \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label females \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 59 \
#   --scale_free_cutoff 0.8

# # #v3
# #pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --dataset1_label males \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v3_technical_pcs_only/males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --dataset2_label females \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/v3_technical_pcs_only/Males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --adjustment_version v3_technical_pcs_only \
# #   --softpower_cor pearson \
# #   --power_min 1 \
# #   --power_max 20 \
# #   --threads 59 \
# #   --scale_free_cutoff 0.8

# # #script 09
# # #v1
# # pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/shared_complete_regions_4col.tsv \
# #   --dataset1_label males \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v1_all_pcs/males_Adjusted_Region_Methylation_allPCs.rds \
# #   --dataset2_label females \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/v1_all_pcs/females_Adjusted_Region_Methylation_allPCs.rds \
# #   --adjustment_version v1_all_pcs \
# #   --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/v1_all_pcs/chosen_power_pearson.txt \
# #   --consensus_cor pearson \
# #   --network_type signed \
# #   --tom_type signed \
# #   --deep_split 1 \
# #   --min_module_size 10 \
# #   --merge_cut_height 0.2 \
# #   --max_block_size 1000 \
# #   --network_calibration "single quantile" \
# #   --save_consensus_toms TRUE \
# #   --threads 59 \
# #   --max_p_outliers 0.1

# #v2

# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/shared_complete_regions_4col.tsv \
#   --dataset1_label males \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label females \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/v2_exclude_protected_pcs/chosen_power_pearson.txt \
#   --consensus_cor pearson \
#   --network_type signed \
#   --tom_type signed \
#   --deep_split 1 \
#   --min_module_size 10 \
#   --merge_cut_height 0.2 \
#   --max_block_size 50000 \
#   --network_calibration "single quantile" \
#   --save_consensus_toms TRUE \
#   --threads 16 \
#   --max_p_outliers 0.1

# #v3

# # pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
# #   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/shared_complete_regions_4col.tsv \
# #   --dataset1_label males \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v3_technical_pcs_only/males_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --dataset2_label females \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/05b_shared_complete_regions/females_Methylation_complete.rds \
# #   --adjustment_version v3_technical_pcs_only \
# #   --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/v3_technical_pcs_only/chosen_power_pearson.txt \
# #   --consensus_cor pearson \
# #   --network_type signed \
# #   --tom_type signed \
# #   --deep_split 1 \
# #   --min_module_size 10 \
# #   --merge_cut_height 0.1 \
# #   --max_block_size 1000 \
# #   --network_calibration "single quantile" \
# #   --save_consensus_toms TRUE \
# #   --threads 59 \
# #   --max_p_outliers 0.1

# # # --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/08_soft_power/v2_exclude_protected_pcs/chosen_power_pearson.txt \
# # # --power 15

# # #Script 10
# #v2
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/10_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --adjustment_version v2_exclude_outcome_exposure_pcs \
#   --consensus_cor pearson \
#   --max_p_outliers 0.1

# # #Script 11
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/11_membership_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --dataset1_label males \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/males/v2_exclude_protected_pcs/males_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label females \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/07_methylation_adjustment/females/v2_exclude_protected_pcs/females_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --membership_cor bicor \
#   --max_p_outliers 0.1

# #script 12
# #v2
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/12a_me_trait_analysis_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --dataset1_label males \
#   --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_males.xlsx \
#   --dataset2_label females \
#   --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/data/metadata/merged_qc_females.xlsx \
#   --run_bicor TRUE \
#   --run_pearson FALSE \
#   --p_thresh 0.05 \
#   --top_n 50 \
#   --module_dendro_distance bicor \
#   --max_p_outliers 0.1

# # Script 12b
# pixi run Rscript   /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/comethyl_scripts/consensus/12b_me_trait_presentation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --dataset1_label males \
#   --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/v2_exclude_protected_pcs/males/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --dataset2_label females \
#   --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/comethyl_output/consensus/12a_me_trait_analysis/v2_exclude_protected_pcs/females/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --set_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/config/presentation_sets \
#   --module_dendro_distance bicor \
#   --p_thresh 0.05 \
#   --top_n 50 \
#   --full_width 12 \
#   --full_height 12 \
#   --top_width 9 \
#   --top_height 6