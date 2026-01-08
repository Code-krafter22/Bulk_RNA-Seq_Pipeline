nextflow.enable.dsl=2

include { LIMMA_VOOM } from './modules/limma_voom.nf'

// ----------------------
// Inputs: single or batch
// ----------------------
workflow {

  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
        def id     = row.id ?: row.run ?: row.name ?: "run"
        def counts = file(row.counts)
        def meta   = file(row.meta)
        tuple(id, counts, meta)
      }
      | LIMMA_VOOM
  }
  else {
    if (!params.counts || !params.meta) {
      error "Provide either:\n  (1) --counts + --meta\n  OR\n  (2) --samplesheet (CSV with columns: id,counts,meta)\n"
    }

    Channel
      .of( tuple(params.id ?: "run1", file(params.counts), file(params.meta)) )
      | LIMMA_VOOM
  }
}
