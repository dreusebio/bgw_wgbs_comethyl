#!/bin/bash

# ============================================================================
# run_enrichment.sbatch
#
# Run enrichment for three DMR timepoints as a SLURM array.
#
# Submit all timepoints:
#   mkdir -p 00_std_err_logs
#   sbatch run_enrichment.sbatch
#
# Override plotting mode:
#   sbatch --export=ALL,MODE=publication run_enrichment.sbatch
# ============================================================================

#SBATCH --job-name=DMRichR_enrichment
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --output=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/logs/00_std_err_logs/DMRichR_enrichment.out
#SBATCH --error=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/logs/00_std_err_logs/DMRichR_enrichment.err
#SBATCH --mail-type=ALL

set -euo pipefail

# ============================================================================
# Run information
# ============================================================================

START_TIME=$(date +%s)

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

# ============================================================================
# Pixi environment
# ============================================================================

export PIXI_HOME="/quobyte/lasallegrp/home/George/.pixi"
export PATH="${PIXI_HOME}/bin:${PATH}"

export PIXI_CACHE_DIR="/tmp/${USER}/pixi-cache"
mkdir -p "${PIXI_CACHE_DIR}"

THREADS="${SLURM_CPUS_PER_TASK:-1}"

#the environment variable PIXI_PROJECT is set to the reproducibility/pixi folder of the comethyl project, because the dmrichr version is old.
# So KEGG with clusterProfiler will not work with the dmrichr version. 
# The comethyl version is up to date and works with clusterProfiler KEGG enrichment.
PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/reproducibility/comethyl/pixi"

PIXI_MANIFEST="${PIXI_PROJECT}/pixi.toml"

if [[ ! -f "${PIXI_MANIFEST}" ]]; then
    echo "ERROR: Pixi manifest not found: ${PIXI_MANIFEST}" >&2
    exit 1
fi

echo "Allocated threads: ${THREADS}"

# Optional environment check
pixi run \
    --manifest-path "${PIXI_MANIFEST}" \
    Rscript -e 'cat("R home:", R.home(), "\n"); sessionInfo()'

# ============================================================================
# Project paths
# ============================================================================

SCRIPT_DIR="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/DMRichr/run_in_DMRs_folder"
ENRICHMENT_SCRIPT="${SCRIPT_DIR}/02_enrichment_plot.R"

# Select one analysis
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/DMRichR_environmental_exposure_index"

# Other options:
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/rwg_wfl_0to6mo/cytosine_reports/cutoff_005"
BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/rwg_wfl_6to12mo/cytosine_reports/cutoff_005"


GENES="${BASE}/DMRs/DMRs_annotated.xlsx"
OUTDIR="${BASE}/enrichment"
PREFIX=$(basename "${BASE}")

# Choose either presentation or publication.
MODE="${MODE:-presentation}"

if [[ "${MODE}" != "presentation" && "${MODE}" != "publication" ]]; then
    echo "ERROR: MODE must be 'presentation' or 'publication'." >&2
    echo "Received: ${MODE}" >&2
    exit 1
fi

if [[ ! -f "${GENES}" ]]; then
    echo "ERROR: gene file not found: ${GENES}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# ============================================================================
# Run enrichment
# ============================================================================

echo
echo "Running enrichment"
echo "  Analysis: $(basename "${BASE}")"
echo "  Prefix:   ${PREFIX}"
echo "  Input:    ${GENES}"
echo "  Output:   ${OUTDIR}"
echo "  Mode:     ${MODE}"
echo

pixi run \
    --manifest-path "${PIXI_MANIFEST}" \
    Rscript --vanilla "${ENRICHMENT_SCRIPT}" \
    --genes "${GENES}" \
    --outdir "${OUTDIR}" \
    --prefix "${PREFIX}" \
    --mode "${MODE}" \
    --kegg_method "clusterProfiler" \
    --top_n 9 \
    --go_dbs "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025" \
    --gwas_db "GWAS_Catalog_2025"


    # --go_method enrichr \
    # --kegg_method "clusterProfiler" \
# ============================================================================
# Completion information
# ============================================================================

END_TIME=$(date +%s)
ELAPSED_SECONDS=$((END_TIME - START_TIME))

echo
echo "Enrichment analysis (${PREFIX}) completed successfully."
echo "End time: $(date)"
echo "Elapsed time: ${ELAPSED_SECONDS} seconds"
