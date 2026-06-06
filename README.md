# snakemake_colocalisation

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Snakemake pipeline for Bayesian colocalisation of GWAS signals with eQTL and pQTL datasets. Originally created by Nikos Panousis and Kaur Alasoo; adapted by Laura Fachal, Tobi Alegbe and Bradley Harris.

This version of the pipeline was used in Fachal et al. (2026) *Resolving inflammatory bowel disease risk variants to genes and cell types* (medRxiv: https://www.medrxiv.org/content/10.64898/2026.05.13.26352926v2).

---

## Overview

This pipeline tests for colocalisation between GWAS summary statistics and molecular QTL datasets (eQTL, pQTL) using the [coloc](https://github.com/chr1swallace/coloc) R package. It is designed to run at scale on an HPC cluster via Snakemake.

---

## Repository structure

```
scripts/
├── run_coloc_eQTL_per_eqtl_study.smk       # Per-study Snakemake workflow
├── GWAS_run_coloc_eQTL_LF.R               # Core coloc execution script
├── functions_me_eQTL_bh.R                 # Helper functions
├── coloc_config_eQTL.yaml                  # Config template (standard)
├── coloc_config_eQTL_TensorQTL.yaml        # Config template (TensorQTL input)
├── coloc_config_eQTL_per_eqtl_study.yaml   # Config template (per-study)
├── cluster_config.yaml                     # HPC cluster resource config
├── join_files_per_eqtl_study.sh            # Merge results across conditions
├── eQTL_catalogue/                         # Scripts to prepare eQTL Catalogue data
├── hu_2021/                                # Scripts to prepare Hu 2021 scRNA eQTL data
├── macromap/                               # Scripts to prepare macromap data
├── pQTL_decode/                            # Scripts to prepare deCODE pQTL data
└── pQTL_sparc/                             # Scripts to prepare SPARC pQTL data

gwas_sumstats/                              # GWAS summary statistics (not included)
eQTL_datasets.csv                          # Registry of available eQTL datasets
qtl_datasets.xlsx                          # Registry of available QTL datasets
submit_snakemake_per_eqtl_study.sh         # Top-level job submission script
```

---

## Software environment

The Singularity container is available on Zenodo (https://zenodo.org/records/20560488). Download the `.sif` file and update the container path in `scripts/run_coloc_eQTL_per_eqtl_study.smk`:
```bash
singularity shell colocalization_singularity.sif
```

---

## Running the pipeline

### Step 1 — Prepare GWAS summary statistics


Format the full summary statistics file with the following columns:
```
RSid  Chr  Pos  Eff_allele  MAF  pval  beta  OR  log_OR  se  z.score  Disease  PubmedID  used_file
```

Sort, compress, and index:
```bash
bgzip file_sorted.txt
tabix -s2 -b3 -e3 -S1 file_sorted.txt.gz
```

Create a top-hits file (variants with p ≤ 1×10⁻⁵):
```bash
for i in *file_sorted.txt.gz; do
  g1=$(echo $i | cut -f1 -d ".")
  awk '$6<=1E-5' <(zcat ${i}) | bgzip -c > ${g1}.top_hits.txt.gz
done
```

Create a summary file (tab-separated: short name, file stem, disease category):
```
CD    cd_allchr_summary_stats_sorted    Autoimmune_Inflammatory_disease
IBD   ibd_allchr_summary_stats_sorted   Autoimmune_Inflammatory_disease
UC    uc_allchr_summary_stats_sorted     Autoimmune_Inflammatory_disease
```

Initialise a file to track all GWAS:
```bash
touch gwas_sumstats_final/gwas_files.txt
```

Place all files in `gwas_sumstats_final/`.


### Step 2 — Prepare QTL datasets

#### 2A. Publicly available datasets

Pre-processing scripts for the following datasets are included in `scripts/`:

| Dataset | Directory | Description |
|---------|-----------|-------------|
| eQTL Catalogue V7 | `eQTL_catalogue/` | Multi-tissue eQTLs including GTEx V8 |
| Hu 2021 | `hu_2021/` | Human Cell Atlas single-cell eQTLs |
| Macromap | `macromap/` | Macrophage single-cell eQTLs |
| deCODE pQTL | `pQTL_decode/` | Plasma protein QTLs |
| SPARC pQTL | `pQTL_sparc/` | Plasma protein QTLs |

Each directory contains scripts to generate: nominal summary statistics, permuted results, sample size files, and variant info files.

#### 2B. Custom (non-public) datasets

Files must be placed in `input/<study_name>/` with the following structure:

**1. Nominal summary statistics** — one gzipped, bgzipped and tabix-indexed file per condition:
```
variant  r2  pvalue  molecular_trait_object_id  molecular_trait_id  maf  gene_id  median_tpm  beta  se  an  ac  chromosome  position  ref  alt  type  rsid
```

**2. Permuted results** — one file per condition listing index variants for significant eGenes (FDR < 0.05):
```
Phenotype_ID  Chromosome_phe  TSS  TSS_end  Strand  Total_no_variants_cis  ...  Corrected_pvalue  std.err
```

**3. Sample sizes** — one tab-separated single-row file per condition:
```
<condition_name>    <n_samples>
```

**4. Variant info** — one bgzipped, tabix-indexed file per condition (no header):
```
<chr>  <pos>  <id>  <ref>  <alt>  <type>  <MAC>  <total_allele_number>
```

**5. Conditions file** — single-column file (no header) listing each condition name, one per row.

Scripts for generating these files from TensorQTL output are in `scripts/GWAS_run_coloc_eQTL_LF.R`.

---

### Step 3 — Configure and run the pipeline

Three config templates are provided in `scripts/`:

| Config file | Use case |
|-------------|----------|
| `coloc_config_eQTL.yaml` | Standard eQTL input |
| `coloc_config_eQTL_TensorQTL.yaml` | TensorQTL-formatted input |
| `coloc_config_eQTL_per_eqtl_study.yaml` | Per-study submission |

Edit the appropriate config file for your use case:

| Parameter | Description |
|-----------|-------------|
| `general_path` | Path to the cloned repository |
| `input_dir` | Path to QTL input files (from Step 2) |
| `coloc_output` | Path for results output |
| `coloc_phenotypes` | Set to `"featureCounts"` |
| `coloc_window` | Window around lead variant to test (e.g. `2e6`) |
| `gwas_dir` | Path to GWAS summary statistics (from Step 1) |
| `gwas_traits` | List of GWAS trait names to test |
| `gwas_list` | Corresponding GWAS summary statistic files |
| `chunks` | Job chunks (default: `["1"..."20"]`) |
| `conditions_file` | Path to conditions file (from Step 2B) |

> **Important:** Also update the config file path at the top of `scripts/run_coloc_eQTL_per_eqtl_study.smk`.

Submit the pipeline:
```bash
bsub -M 2000 -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
  -o sm_logs/snakemake_master-%J-output.log \
  -e sm_logs/snakemake_master-%J-error.log \
  -q oversubscribed -J "snakemake_master_COLOC" \
  < submit_snakemake_per_eqtl_study.sh
```

Once complete, merge results across conditions:
```bash
bash scripts/join_files_per_eqtl_study.sh
```

---

## Citation

> Fachal L, et al. (2026). Resolving inflammatory bowel disease risk variants to genes and cell types. medRxiv. https://www.medrxiv.org/content/10.64898/2026.05.13.26352926v2
