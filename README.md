# Bulk RNA-seq Pipeline (edgeR + voom/limma)

![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white) ![Bioconductor](https://img.shields.io/badge/Bioconductor-3B4BAA?logo=bioconductor&logoColor=white) ![Conda](https://img.shields.io/badge/conda-44A833?logo=anaconda&logoColor=white) ![Nextflow](https://img.shields.io/badge/Nextflow-0DC09D?logo=nextflow&logoColor=white)


## Overview
This script performs bulk RNA-seq differential expression analysis from a **raw counts matrix** and **sample metadata**:
- **edgeR (TMM)** normalization and log2-CPM calculation  
- Low-expression gene filtering  
- **limma-voom** linear modeling with empirical Bayes moderation  
- QC plots (boxplots + voom mean–variance trend) and DEG result tables  

## Inputs
### 1) Counts matrix (`--counts`)
- Formats: `.tsv / .csv / .txt / .xlsx / .xls`
- Rows = genes/features, columns = samples
- One column is gene IDs (auto-detected; override with `--feature-col`)
- Counts must be **non-negative** and sample column names must be **unique**

### 2) Metadata (`--meta`)
- Formats: `.tsv / .csv / .txt / .xlsx / .xls`
- Must contain:
  - `--sample-col`: sample IDs that match **counts column names**
  - `--group-col`: group/condition labels (must have ≥2 groups)
- `--ref-level` must exist in the `--group-col` values

## Install (Rscript usage)
This creates a fuller interactive environment for running the script directly:
```bash
conda env create -f environment-full.yml
conda activate rnaseq-r
```

## Run
### Text inputs (TSV/CSV)
```bash
Rscript bulk_rnaseq_pipeline.R \
  --counts counts.tsv \
  --meta metadata.tsv \
  --outdir results \
  --sample-col SampleID \
  --group-col PATHDX \
  --ref-level Healthy
```

### Excel metadata (XLSX)
```bash
Rscript bulk_rnaseq_pipeline.R \
  --counts counts.tsv \
  --meta metadata.xlsx \
  --meta-sheet 1 \
  --outdir results \
  --sample-col "Sample ID" \
  --group-col PATHDX \
  --ref-level Healthy
  ```

## Install (Nextflow wrapper) 
This repo includes a Nextflow (DSL2) wrapper (main.nf) that runs the same R pipeline reproducibly.

## Nextflow environment
The Nextflow pipeline is configured to use the lightweight conda environment at:
```text
    envs/environment.yml
```

You do not need environment-full.yml to run the Nextflow pipeline.

## Requirements

  -Java 17+
  -Nextflow
  -Conda (or Mamba)

## Run
```bash
nextflow run main.nf -profile conda \
  --id GSE203206 \
  --counts data/GSE203206_Subramaniam.ADRC_brain.counts.tsv \
  --meta   data/GSE203206_Subramaniam.ADRC_brain.metadata.tsv \
  --outdir results \
  --sample_col SampleID \
  --group_col PATHDX \
  --ref_level Healthy
```