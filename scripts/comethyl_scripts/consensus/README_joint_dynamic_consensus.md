# Revised consensus Comethyl workflow: jointly dynamic regions

## Design decision

The joint-dynamic behavior is implemented as an **optional general pipeline
mode**, not as a hidden project-specific rewrite.

- `--selection_mode complete_only` preserves the historical behavior.
- `--selection_mode joint_sd_all` is recommended when every retained feature
  must be variable in every consensus dataset, such as males versus females.

This makes the option reusable for sex, tissue, cohort or timepoint consensus
analyses without changing results for existing users.

## Most important workflow change

The decisive SD filter moves from the reference-only Script 02 to Script 05b,
after region methylation exists for every dataset.

1. Scripts 00-01 remain unchanged.
2. Script 02 applies region coverage but uses a permissive reference SD:

   ```bash
   --covMin 10 --methSD 0
   ```

3. Scripts 03-05 calculate the same region definitions in females and males.
4. Script 05b requires complete data and SD >= the chosen threshold in both:

   ```bash
   --selection_mode joint_sd_all --joint_meth_sd 0.07
   ```

5. Scripts 06 onward use only Script 05b outputs.

This prevents a female-only SD filter from removing regions that could be
dynamic in males before male methylation is evaluated.

## Recommended first candidates

Do not build a large coverage-by-SD grid immediately. Start with two coverage
prefilters and apply the same joint SD rule:

| Candidate | Script 02 coverage | Script 02 SD | Script 05b joint SD |
|---|---:|---:|---:|
| A | 10 | 0 | 0.07 |
| B | 14 | 0 | 0.07 |

If too many or too few regions remain, compare joint SD 0.06 and 0.08 using the
already-created region methylation matrices. Changing only Script 05b's joint
SD does not require rerunning CpG import or region methylation.

## Example Script 05b call for GROWell sexes

```bash
pixi run Rscript scripts/consensus/05b_filter_shared_complete_regions_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset1_label females \
  --dataset1_meth /path/to/females/Region_Methylation.rds \
  --dataset2_label males \
  --dataset2_meth /path/to/males/Region_Methylation.rds \
  --regions_file /path/to/Filtered_Regions.txt \
  --selection_mode joint_sd_all \
  --joint_meth_sd 0.07 \
  --min_regions 1000 \
  --write_sd_table TRUE
```

The output label records the selection rule, for example:

```text
covMin10_methSD0_jointSDall0p07
```

Therefore a historical run cannot be silently overwritten.

## How to test whether it is likely to work

### Stage 1: Inspect Script 05b QC

Review:

- `joint_variability_summary.tsv`
- `joint_variability_diagnostics.tsv`
- `run_parameters.txt`

Confirm that both datasets contribute a reasonable number of SD-passing
regions and that the final intersection is large enough.

### Stage 2: Run the low-memory repeated benchmark

Run Script 08b on the adjusted Script 07 matrices using 20,000-30,000 regions
and five seeds:

```bash
pixi run Rscript scripts/consensus/08b_softpower_benchmark_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_bgw_comethyl_Victoria \
  --dataset1_label females \
  --dataset1_meth /path/to/females_Adjusted_Region_Methylation_allPCs.rds \
  --dataset2_label males \
  --dataset2_meth /path/to/males_Adjusted_Region_Methylation_allPCs.rds \
  --adjustment_version v1_all_pcs \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 30 \
  --subsample_size 30000 \
  --n_seeds 5 \
  --threads 8
```

The revised benchmark requires both:

- seed-averaged scale-free R2 >= the requested cutoff; and
- a negative seed-averaged slope.

If no common power passes, it writes `chosen_power: NA`. It no longer chooses
power 1 merely because all powers tied at zero passing datasets.

### Stage 3: Confirm on the full region set

Only the best candidate proceeds to Script 08. Script 08b is a screening tool;
Script 08 remains the full-data confirmation.

### Stage 4: Construct modules

Run Script 09 only after Script 08 writes a numeric common power. Pass the
`joint_eligible_regions_4col.tsv` file from Script 05b as `--regions_file`.

## Adjustment-version guard

Only run a downstream adjustment variant when its files exist for **all**
datasets. In the GROWell log, male v2 was intentionally skipped because all
male PCs were classified as protected. A wrapper must therefore skip the v2
consensus call rather than assume the male v2 file exists.

## Thread and memory guidance

Do not request 45 Slurm CPUs while allowing WGCNA to start 59 workers. Set
`--threads` no higher than `SLURM_CPUS_PER_TASK`; 8-16 workers is a reasonable
starting point for the benchmark. A larger worker count can increase peak
memory through concurrent temporary objects.

## Scientific reporting

The joint SD filter establishes a common eligible feature set. It does not
guarantee identical correlations or a successful consensus topology. If the
male curve remains consistently unsuitable after joint filtering, report that
the common eligible regions did not support a stable shared topology rather
than repeatedly tuning filters until a desired R2 appears.

## Files changed materially

- `02_filter_reference_regions_consensus.R`: documents permissive reference SD
  for joint mode and warns about premature reference-only SD filtering.
- `05b_filter_shared_complete_regions_consensus.R`: adds joint-SD modes,
  diagnostics, collision-safe labels and identical final matrices.
- `08_softpower_consensus.R`: requires negative slope and returns `NA` instead
  of the invalid power-1 fallback.
- `08b_softpower_benchmark_consensus.R`: repeated subsampling now retains and
  evaluates slope.
- `08c_softpower_benchmark_compare_consensus.R`: ranks candidates using R2 and
  slope; old benchmark files are clearly flagged for rerunning.

All other supplied scripts are included in full so this folder is a complete
00-11 handoff.
