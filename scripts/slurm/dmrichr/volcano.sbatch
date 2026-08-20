#!/bin/bash

# ============================================================================
# run_volcano.sbatch
#
# Run the full volcano plot (01_Volcano_plot.R) for three DMR timepoints as
# a SLURM array. Mirrors run_enrichment.sbatch's structure and naming
# convention (<group>_<timepoint> prefix) so outputs from both scripts are
# consistent and easy to cross-reference.
#
# Submit all timepoints:
#   mkdir -p 00_std_err_logs
#   sbatch run_volcano.sbatch
# ============================================================================

#SBATCH --job-name=DMRichR_volcano
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/logs/00_std_err_logs/DMRichR_volcano.out
#SBATCH --error=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/logs/00_std_err_logs/DMRichR_volcano.err
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

# Same pixi project used for enrichment (comethyl reproducibility env).
# EDIT if readxl/ggplot2/ggrepel/dplyr/optparse aren't all installed there --
# these are common enough that they likely already are, but worth confirming
# with the environment check below before trusting a real run.
# PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/reproducibility/pixi"
PIXI_PROJECT="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/reproducibility/comethyl/pixi"


PIXI_MANIFEST="${PIXI_PROJECT}/pixi.toml"

if [[ ! -f "${PIXI_MANIFEST}" ]]; then
    echo "ERROR: Pixi manifest not found: ${PIXI_MANIFEST}" >&2
    exit 1
fi

# Optional environment check
pixi run \
    --manifest-path "${PIXI_MANIFEST}" \
    Rscript -e 'cat("R home:", R.home(), "\n"); sessionInfo()'

# ============================================================================
# Project paths
# ============================================================================

SCRIPT_DIR="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/scripts/DMRichr/run_in_DMRs_folder"
VOLCANO_SCRIPT="${SCRIPT_DIR}/01_Volcano_plot.R"

# Select one analysis -- same convention as run_enrichment.sbatch
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/DMRichR_environmental_exposure_index"

# Other options:
# BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/rwg_wfl_0to6mo/cytosine_reports/cutoff_005"
BASE="/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria/results/rwg_wfl_6to12mo/cytosine_reports/cutoff_005"

if [[ ! -f "${VOLCANO_SCRIPT}" ]]; then
    echo "ERROR: volcano script not found: ${VOLCANO_SCRIPT}" >&2
    exit 1
fi

# : background_annotated.xlsx lives alongside DMRs_annotated.xlsx.
# Edit these paths if that's not where DMRichR wrote it for your runs.
BACKGROUND="${BASE}/DMRs/background_annotated.xlsx"
OUTDIR="${BASE}/volcano"
PREFIX=$(basename "${BASE}")

# Optional overrides at submit time, e.g.:
#   sbatch --export=ALL,PVAL_CUTOFF=0.01,N_LABELS=30 run_volcano.sbatch
PVAL_CUTOFF="${PVAL_CUTOFF:-0.05}"
N_LABELS="${N_LABELS:-20}"

if [[ ! -f "${BACKGROUND}" ]]; then
    echo "ERROR: background file not found: ${BACKGROUND}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# ============================================================================
# Run volcano plot
# ============================================================================

echo
echo "Running volcano plot"
echo "  Analysis:    $(basename "${BASE}")"
echo "  Prefix:      ${PREFIX}"
echo "  Input:       ${BACKGROUND}"
echo "  Output:      ${OUTDIR}"
echo "  P cutoff:    ${PVAL_CUTOFF}"
echo

cmd=(pixi run --manifest-path "${PIXI_MANIFEST}" Rscript --vanilla \
    "${VOLCANO_SCRIPT}" \
    --background "${BACKGROUND}" \
    --outdir "${OUTDIR}" \
    --prefix "${PREFIX}" \
    --n_labels "${N_LABELS}" \
    --pval_cutoff "${PVAL_CUTOFF}" )

printf '%q ' "${cmd[@]}"; echo
"${cmd[@]}"

# ============================================================================
# Completion information
# ============================================================================

END_TIME=$(date +%s)
ELAPSED_SECONDS=$((END_TIME - START_TIME))

echo
echo "Volcano plot (${PREFIX}) completed successfully."
echo "End time: $(date)"
echo "Elapsed time: ${ELAPSED_SECONDS} seconds"
