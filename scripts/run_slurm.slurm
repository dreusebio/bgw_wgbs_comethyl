#!/bin/bash
#SBATCH --job-name=Cometh_BGW
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --time=23:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err                   	
#SBATCH --cpus-per-task=59  #use 30 or 50 cores for certain things
#SBATCH --account=lasallegrp #genome-center-grp
#SBATCH --mem=700G #Need up to 500G for certain tasks
#SBATCH --mail-type=ALL  
#SBATCH --chdir=/quobyte/lasallegrp/projects/bgw_wgbs_comethyl

###################
# Run Information #
###################

start=`date +%s`

hostname

####################
# Load Environment #
####################


# Source your config file so conda can be used (likely .profile, could be .bashrc)

source ~/.profile
conda activate /quobyte/lasallegrp/George/envs/Comethyl_v1.3.0_BACKUP/
# conda activate /quobyte/lasallegrp/programs/.conda/envs/Comethyl_v1.3.0_UPDATED/
unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# # quick sanity check (writes to your .out log)
# Rscript --vanilla -e 'cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n"); print(.libPaths()); library(SparseArray); library(comethyl); cat("OK\n")'

##############
# Comethyl Scripts#
##############

# # Script 00
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/00_import_cpg_reports.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --meta_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/processed/cytosine_reports \
#   --workers 4

# # Script 01
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/01_filter_cpgs_call_regions.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --input_rds /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/00_cpg_extraction/Unfiltered_BSseq.rds \
#   --cov 3 \
#   --per_sample 0.75


# #Script 02
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/02_filter_regions.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --input_regions /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/01_filter_regions/cov3_75pct/Regions.txt \
#   --cov_min 4 \
#   --meth_sd 0.08


# # Script 03
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/03_region_methylation.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --input_bs /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/01_filter_regions/cov3_75pct/Filtered_BSseq.rds \
#   --input_regions /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/02_filter_regions/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt

# # Script 04
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/04_pcs_and_trait_diagnostics.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --input_meth /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/03_region_methylation/cov3_75pct/covMin4_methSD0p08/Region_Methylation.rds \
#   --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/technical_traits.txt \
#   --threads 59



# # Script 05
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/05_adjust_region_methylation.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --input_meth /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/03_region_methylation/cov3_75pct/covMin4_methSD0p08/Region_Methylation.rds \
#   --input_pcs /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin4_methSD0p08/PCs.rds \
#   --pc_trait_stats /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/04_pc_diagnostics/cov3_75pct/covMin4_methSD0p08/PC_Trait_Correlation_Stats_Bicor.tsv \
#   --protected_traits_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/technical_traits.txt \
#   --run_v1 TRUE \
#   --run_v2 TRUE \
#   --run_v3 TRUE \
#   --v2_cor_method bicor \
#   --v2_p_thresh 0.05 \
#   --v3_p_thresh 0.05


# # Script 06
# Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/06_soft_power.R \
#   --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
#   --input_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
#   --input_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --input_v3 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin4_methSD0p08/v3_technical_pcs_only/Adjusted_Region_Methylation_technicalPCs_only.rds \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 59

#Script 07
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/07_module_detection.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --regions_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/02_filter_regions/cov3_75pct/covMin4_methSD0p08/Filtered_Regions.txt \
  --input_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
  --softpower_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/06_soft_power/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/SoftPower_pearson.rds \
  --input_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/05_methylation_adjustment/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --softpower_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/06_soft_power/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/SoftPower_pearson.rds \
  --module_cor pearson \
  --deep_split 4 \
  --min_module_size 50 \
  --merge_cut_height 0.20 \
  --threads 59

# script 08

Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/08_module__sample_diagnostics.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --modules_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --modules_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --module_cor pearson \
  --sample_dendro_distance euclidean \
  --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx
#script 09a
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/09a_me_trait_analysis.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --modules_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --run_bicor TRUE \
  --p_thresh 0.05 \
  --top_n 25 \
  --module_dendro_distance bicor \

#script 09b
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/09b_me_trait_presentation_heatmaps.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --stats_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --modules_rds /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --stats_file_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --modules_rds_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --set_dir /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/presentation_sets \
  --p_thresh 0.05 \
  --top_n 25


# #Script 10
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/10_circos_and_module_summaries.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --modules_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --modules_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --genome hg38 \
  --single_module_color red \
  --all_modules_width 4.0 \
  --all_modules_height 4.0 \
  --single_module_width 5.5 \
  --single_module_height 5.5 \
  --counts_width 7 \
  --counts_height 6



# #Script 11
#v1
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/11_me_trait_pair_plots.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --modules_rds /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --stats_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
  --set_dir /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/presentation_sets \
  --p_thresh 0.05 \
  --only_significant TRUE \
  --cat_plot_style violin \
  --top_n_pairs 50

#v2
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/11_me_trait_pair_plots.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --modules_rds /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --stats_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
  --set_dir /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/presentation_sets \
  --p_thresh 0.05 \
  --only_significant TRUE \
  --cat_plot_style violin \
  --top_n_pairs 50

#script 12
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/12a_annotate_modules.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --modules_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --genome hg38 \
  --annotation_mode auto


#script 12a
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/12a_annotate_modules.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --sample_info /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/data/metadata/merged_qc.xlsx \
  --modules_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Modules.rds \
  --modules_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/07_module_detection/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/Modules.rds \
  --genome hg38 \
  --annotation_mode auto

#script 12b
Rscript /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/scripts/comethyl_scripts/12b_enrich_selected_modules.R \
  --project_root /quobyte/lasallegrp/projects/bgw_wgbs_comethyl \
  --annotation_dir_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/12_annotation/cov3_75pct/covMin4_methSD0p08/v1_all_pcs \
  --annotation_dir_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/12_annotation/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs \
  --trait_list_file /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/config/presentation_sets/heatmap_primary_traits.txt \
  --stats_file_v1 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --stats_file_v2 /quobyte/lasallegrp/projects/bgw_wgbs_comethyl/comethyl_output/09a_me_trait_analysis/cov3_75pct/covMin4_methSD0p08/v2_exclude_protected_pcs/ME_Trait_Correlation_Stats_Bicor.tsv \
  --p_thresh 0.05 \
  --do_kegg true \
  --do_go true \
  --do_reactome true \
  --do_enrichr true \
  --enrichr_dbs GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human \
  --use_enrichr_background false \
  --enrichr_retries 5 \
  --enrichr_sleep 2 \
  --plot_p_col pvalue \
  --plot_top_n 10 \
  --plot_sig_cutoff 0.05

# Script 12b
