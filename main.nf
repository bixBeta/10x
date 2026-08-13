#!/usr/bin/env nextflow

/*
========================================================================================
    10x - modular Nextflow pipeline for 10x Genomics single-cell data
========================================================================================
    Implemented modalities : gex
    Planned                : atac, arc (arc_gex + arc_atac)
    Github                 : https://github.com/bixBeta/10x
========================================================================================
*/

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

    TENX()
}
