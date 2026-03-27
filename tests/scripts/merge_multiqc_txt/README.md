# MultiQC Data Merging

## Overview

In this project, multiple quality control (QC) runs were generated using the same Epigenerator pipeline version (v1.14).

Because:
- the same pipeline version (v1.14) was used across runs
- the same QC tools and parsing structure were expected
- intermediate QC files were not retained

the MultiQC summary tables in `multiqc_data/*.txt` were merged directly across runs.

---

## Rationale

Under normal circumstances, MultiQC should be rerun from the original QC outputs such as FastQC, Bismark, Cutadapt, and Picard logs.

In this case, only the `multiqc_data/` summary tables were available, so re-running MultiQC was not possible. A controlled merge of the tabular outputs was therefore used as a practical workaround.

---

## Method

A custom script (`merge_multiqc_txt.py`) was used to:

1. identify matching `.txt` files across multiple `multiqc_data/` folders
2. verify that columns match exactly before merging
3. append rows across runs for compatible files
4. remove duplicate samples when the same sample ID appears in more than one run
5. retain the first occurrence of each sample based on input order
6. add a `multiqc_run_label` column to track the run of origin
7. skip incompatible or non-tabular files safely

---

## Sample handling

If the same sample ID appears in multiple runs:

- the sample is assumed to represent the same biological sample
- only the first occurrence is retained
- duplicate entries from later runs are removed

This means that input order determines priority.

For example:

```bash
python3 tests/scripts/merge_multiqc_txt/merge_multiqc_txt.py \
  --inputs tests/comethyl_minidata/09_multiqc/multiqc_data tests/comethyl_minidata/09_multiqc/multiqc_data_1 \
  --labels run1 run2 \
  --output tests/comethyl_minidata/09_multiqc/merged_multiqc_data