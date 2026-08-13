#!/usr/bin/env nextflow

/*
========================================================================================
    10x - modular Nextflow pipeline for 10x Genomics single-cell data
========================================================================================
    Implemented modalities : gex
    Planned                : atac, arc (arc_gex + arc_atac)
    Github                 : https://github.com/OWNER/10x
========================================================================================
*/

nextflow.enable.dsl = 2

include { TENX } from './workflows/tenx'

workflow {

    log.info """\
    ==========================================================
     1 0 x   v${workflow.manifest.version}
    ==========================================================
     input              : ${params.input}
     reference          : ${params.reference}
     outdir             : ${params.outdir}
     cellranger version : ${params.cellranger_version}
     container          : ${params.cellranger_container ?: "ghcr.io/${params.ghcr_owner}/cellranger:${params.cellranger_version}"}
     profile            : ${workflow.profile}
    ==========================================================
    """.stripIndent()

    TENX ()
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\nDone. Results in ${params.outdir}\n"
        : "\nPipeline failed. See ${params.outdir}/pipeline_info for details.\n" )
}
