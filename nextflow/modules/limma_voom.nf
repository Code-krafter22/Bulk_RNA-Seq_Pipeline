process LIMMA_VOOM {
  tag "${id}"

  publishDir "${params.outdir}/${id}", mode: 'copy', overwrite: true

  // Choose ONE: conda or container (conda recommended to start)
  conda "${projectDir}/envs/environment.yml"

  input:
    tuple val(id), path(counts), path(meta)

  output:
    tuple val(id), path("out/**"), emit: results

  script:
    """
    set -euo pipefail

    mkdir -p out

    ARGS=()
    ARGS+=( --counts "${counts}" )
    ARGS+=( --meta "${meta}" )
    ARGS+=( --outdir out )

    ARGS+=( --sample-col "${params.sample_col}" )
    ARGS+=( --group-col "${params.group_col}" )
    ARGS+=( --ref-level "${params.ref_level}" )

    ARGS+=( --filter-logcpm "${params.filter_logcpm}" )
    ARGS+=( --filter-min-samples "${params.filter_min_samples}" )
    ARGS+=( --alpha "${params.alpha}" )
    ARGS+=( --lfc "${params.lfc}" )

    ARGS+=( --counts-sheet "${params.counts_sheet}" )
    ARGS+=( --meta-sheet "${params.meta_sheet}" )

    if [[ -n "${params.feature_col}" ]]; then
      ARGS+=( --feature-col "${params.feature_col}" )
    fi

    if [[ -n "${params.covariates}" ]]; then
      ARGS+=( --covariates "${params.covariates}" )
    fi

    Rscript ${projectDir}/bin/bulk_rnaseq_pipeline.R "\${ARGS[@]}"
    """
}
