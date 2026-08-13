/*
========================================================================================
    TENX - main workflow: validate inputs, dispatch by modality
========================================================================================
*/

include { SAMPLESHEET } from '../subworkflows/local/samplesheet'
include { GEX         } from '../subworkflows/local/gex'

workflow TENX {

    main:

    if ( !params.input )
        error "No samplesheet given. Use --input samplesheet.csv (see assets/samplesheet_gex.csv)"

    SAMPLESHEET ( params.input )

    if ( params.validate_only ) {

        SAMPLESHEET.out.libraries.subscribe { meta, r1, r2, r3 ->
            log.info "[validate] ${meta.modality}  sample=${meta.sample}  library=${meta.id}  runs=${meta.n_runs}"
        }

    }
    else {

        def ch_versions = Channel.empty()

        // ---- GEX -------------------------------------------------------------
        // The reference is only required when the samplesheet actually contains
        // GEX libraries, so the check rides along the channel rather than
        // firing upfront.
        def ch_gex = SAMPLESHEET.out.gex.map { meta, r1, r2, r3 ->
            if ( !params.reference )
                error "Library '${meta.id}' is GEX but --reference (Cell Ranger transcriptome directory) was not provided"
            [ meta, r1, r2, r3 ]
        }

        def ch_reference = Channel.value(
            params.reference ? file(params.reference, checkIfExists: true) : []
        )

        GEX ( ch_gex, ch_reference )
        ch_versions = ch_versions.mix( GEX.out.versions )

        // ---- not yet implemented ---------------------------------------------
        SAMPLESHEET.out.atac
            .map { meta, r1, r2, r3 -> meta.id }
            .collect()
            .subscribe { ids ->
                log.warn "ATAC support lands in v0.2. Skipping ${ids.size()} librar${ids.size() == 1 ? 'y' : 'ies'}: ${ids.join(', ')}"
            }

        SAMPLESHEET.out.arc
            .map { meta, r1, r2, r3 -> meta.id }
            .collect()
            .subscribe { ids ->
                log.warn "ARC/multiome support lands in v0.3. Skipping ${ids.size()} librar${ids.size() == 1 ? 'y' : 'ies'}: ${ids.join(', ')}"
            }

        // ---- provenance -------------------------------------------------------
        ch_versions
            .map { it.text }
            .unique()
            .collectFile(
                name: 'software_versions.yml',
                storeDir: "${params.outdir}/pipeline_info",
                sort: true,
                cache: false
            )
    }
}
