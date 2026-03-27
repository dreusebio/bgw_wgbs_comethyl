# bgw_wgbs_comethyl

Bioinformatics analysis of prenatal family and community exposures and their association with postnatal DNA methylation signatures related to child metabolic health.

---

Data type(s): Whole Genome Bisulfite Sequencing (WGBS)
Study / cohort: Babies GROWELL
Genome build: hg38 
Primary outputs: Co-methylation modules, ME–trait associations, figures

---

Project overview

This project investigates how prenatal family and community-level exposures influence early childhood DNA methylation (DNAm) patterns, particularly in genomic regions related to energy regulation and fat deposition.

Children (ages 1–4) from mothers with pre-pregnancy BMI 25–40 were profiled using WGBS. Prenatal exposures were quantified using the GRASP Place & Health burden percentile, derived from census tract data during pregnancy.

The central hypothesis is that higher prenatal exposure burden is associated with DNAm alterations in pathways relevant to metabolic health (e.g., obesity risk).

Study design
Design: Cross-sectional
Sample size: 49 child participants
Biological sample: Saliva
Age range: 1–4 years
Exposure variable: community exposures
Outcomes:
Genome-wide DNAm
Region-level DNAm (energy balance / adiposity pathways)
Statistical approach: Comethyl

---

This repository implements a WGBS → region-level methylation → co-methylation (Comethyl) pipeline:

Preprocessing
Bismark alignment and CpG extraction
Generation of CpG reports
Region definition
CpG clustering (coverage + variability thresholds)
Region-level methylation
Aggregation into regions
Filtering (coverage, SD thresholds)
Adjustment
Principal component (PC)-based correction
Network analysis
Co-methylation module detection (WGCNA framework)
Association testing
Module eigengenes (MEs) vs traits (GRASP, covariates)
Correlation (bicor/pearson) + regression models
Downstream interpretation
Functional enrichment
Pathway analysis (energy metabolism, adipogenesis)

---

Repository structure

- `data/` - raw + processed data + codebooks, data dictionaries, provenance notes (**not commited**)
- `scripts/` — entrypoints + SLURM submit scripts
- `analysis/` — downstream statistics/figures +  configuration + sample sheets + generated logs
- `docs/` — methods + workflow documentation
- `results/` — generated outputs
- `test/` —  Test data and vignettes

---

Quick start
Clone
git clone https://github.com/dreusebio/bgw_wgbs_comethyl.git
cd bgw_wgbs_comethyl

---

License

MIT License (see LICENSE)

---

Reproducibility standard

See `docs/lab_reproducibility_standard.md`.

---

