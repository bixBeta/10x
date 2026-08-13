/*
========================================================================================
    GEX - 3'/5' gene expression quantification
========================================================================================
    One CELLRANGER_COUNT run per library. Libraries belonging to the same
    biological sample stay linked through meta.sample for later aggregation.
----------------------------------------------------------------------------------------
*/

include { CELLRANGER_COUNT } from '../../modules/local/cellranger/count/main'

workflow GEX {

    take:
    ch_libraries    // [ meta, [r1...], [r2...], [] ]
    ch_reference    // transcriptome dir

    main:

    def ch_versions = Channel.empty()

    CELLRANGER_COUNT ( ch_libraries, ch_reference.first() )
    ch_versions = ch_versions.mix( CELLRANGER_COUNT.out.versions.first() )

    emit:
    outs        = CELLRANGER_COUNT.out.outs
    filtered_h5 = CELLRANGER_COUNT.out.filtered_h5
    metrics     = CELLRANGER_COUNT.out.metrics
    versions    = ch_versions
}
