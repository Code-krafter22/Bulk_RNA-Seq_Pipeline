# Bulk RNA-seq Pipeline (edgeR + voom/limma)

![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white) ![Bioconductor](https://img.shields.io/badge/Bioconductor-3B4BAA?logo=bioconductor&logoColor=white) ![Conda](https://img.shields.io/badge/conda-44A833?logo=anaconda&logoColor=white)

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

## Install
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
  --group-col Condition \
  --ref-level Control
```

### Excel metadata (XLSX)
```bash
Rscript bulk_rnaseq_pipeline.R \
  --counts counts.tsv \
  --meta metadata.xlsx \
  --meta-sheet 1 \
  --outdir results \
  --sample-col "Sample ID" \
  --group-col Condition \
  --ref-level HC
  ```