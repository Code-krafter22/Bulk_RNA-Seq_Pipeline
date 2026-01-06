# Bulk RNA-seq Pipeline (edgeR + voom/limma)

![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white) ![Bioconductor](https://img.shields.io/badge/Bioconductor-3B4BAA?logo=bioconductor&logoColor=white) ![Conda](https://img.shields.io/badge/conda-44A833?logo=anaconda&logoColor=white)

## What this script does
Given a **raw counts matrix** (genes × samples) and a **metadata table**, the script:
- Normalizes counts with **edgeR (TMM)** and computes log2-CPM
- Filters low-expression genes
- Runs differential expression with **limma-voom** (linear model + empirical Bayes)
- Writes QC plots (boxplots + voom trend) and DEG result tables

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

## Install (conda recommended)
```bash
conda create -n rnaseq-r -y -c conda-forge -c bioconda \
  r-base r-optparse r-ggplot2 r-reshape2 r-data.table r-readxl \
  bioconductor-edger bioconductor-limma
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