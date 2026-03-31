# bgw_wgbs_comethyl

Bioinformatics analysis of prenatal family and community exposures and their association with postnatal DNA methylation signatures related to child metabolic health.

---

## Project summary

**Data type:** Whole Genome Bisulfite Sequencing (WGBS)
**Study / cohort:** Babies GROWELL
**Genome build:** hg38
**Primary outputs:** Co-methylation modules, ME–trait associations, figures

---

## Project overview

This project investigates how prenatal family and community-level exposures influence early childhood DNA methylation (DNAm) patterns, particularly in genomic regions related to energy regulation and fat deposition.

Children (ages 1–4) from mothers with pre-pregnancy BMI 25–40 were profiled using WGBS. Prenatal exposures were quantified using the **GRASP Place & Health burden percentile**, derived from census tract data during pregnancy.

**Central hypothesis:**
Higher prenatal exposure burden is associated with DNAm alterations in pathways relevant to metabolic health (e.g., obesity risk).

---

## Study design

* **Design:** Cross-sectional
* **Sample size:** 49 child participants
* **Biological sample:** Saliva
* **Age range:** 1–4 years

**Exposure variable:**

* Community exposure burden (GRASP)

**Outcomes:**

* Genome-wide DNAm
* Region-level DNAm (energy balance / adiposity pathways)

**Statistical framework:** Comethyl (WGCNA-based)

---

## Analysis workflow

This repository implements a full WGBS → co-methylation pipeline:

### 1. Preprocessing

* Bismark alignment
* CpG extraction
* Generation of CpG reports

### 2. Region definition

* CpG clustering (coverage + variability thresholds)

### 3. Region-level methylation

* Aggregation into regions
* Filtering (coverage, SD thresholds)

### 4. Adjustment

* Principal component (PC)-based correction

### 5. Network analysis

* Co-methylation module detection (WGCNA framework)

### 6. Association testing

* Module eigengenes (MEs) vs traits
* Correlation (bicor / pearson)
* Regression models

### 7. Downstream interpretation

* Functional enrichment
* Pathway analysis (energy metabolism, adipogenesis)

---

## Repository structure

```text
bgw_wgbs_comethyl/
├── data/            # raw + processed data (**not committed**)
├── scripts/         # entrypoints + SLURM scripts
├── analysis/        # downstream analyses + configs
├── docs/            # methods + documentation
├── results/         # outputs
├── test/            # test data 
├── reproducibility/ # environment + setup
```

---

## Quick start (reproducible environment)

This project uses **Pixi** for reproducible environments.

### 1. Clone repository

```bash
git clone https://github.com/dreusebio/bgw_wgbs_comethyl.git
cd bgw_wgbs_comethyl
```

---

### 2. Configure environment 


```bash
export PIXI_HOME=/path/to/your/pixi_home
export PIXI_CACHE_DIR=/scratch/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME
```
#### (cluster users e.g UC Davis Hive)
```bash
export PIXI_HOME=/quobyte/lasallegrp/George/.pixi
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
unset RATTLER_CACHE_DIR
unset XDG_CACHE_HOME
```

---

### 3. Install environment

```bash
cd reproducibility/pixi
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

---

### 4. Install comethyl

```bash
pixi run install-comethyl
```

---

### 5. Test installation

```bash
pixi run Rscript -e 'library(comethyl); packageVersion("comethyl")'
```

Expected:

```text
[1] "1.3.0"
```

---

## Running analysis

All scripts should be executed through Pixi:

```bash
pixi run Rscript scripts/<your_script>.R
```

Example:

```bash
pixi run Rscript scripts/comethyl_scripts/00_get_multiqc_merge.R
```

---

## Reproducibility

This project is fully reproducible using:

```text
reproducibility/
  env/
    environment.portable.yml
    install_comethyl.R
  pixi/
    pixi.toml
    pixi.lock
```

### Key features

* `pixi.lock` → guarantees identical environments
* `install_comethyl.R` → handles Bioconductor + GitHub dependencies
* `comethyl` pinned to specific commit → reproducible results

---

## Platform support

| Platform        | Support                         |
| --------------- | ------------------------------- |
| Linux (cluster) | ✅ full                          |
| macOS (Intel)   | ✅ full                          |
| macOS (M1/M2)   | ⚠️ partial (Docker recommended) |
| Windows         | ⚠️ use WSL2 or Docker           |

---

## Troubleshooting

### Pixi install fails (network/cache issues)

```bash
export PIXI_CACHE_DIR=/tmp/$USER/pixi-cache
```

---

### Reinstall environment

```bash
rm -rf /tmp/$USER/pixi-cache
pixi install --concurrent-downloads 1 --concurrent-solves 1
```

---

### Missing R/Bioconductor packages

```bash
pixi run install-comethyl
```

---

## License

MIT License (see `LICENSE`)

---

## Reproducibility standard

See:

```text
docs/lab_reproducibility_standard.md
```

---

## Contact

George Eusebio Kuodza
Postdoctoral Scholar, UC Davis
[gekuodza@health.ucdavis.edu](mailto:gekuodza@health.ucdavis.edu)

---

## Future directions

* Docker container for full portability
* Snakemake / Nextflow integration
* Demo dataset for quick validation

---
