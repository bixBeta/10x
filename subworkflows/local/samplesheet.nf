/*
========================================================================================
    SAMPLESHEET - parse, validate and group the input samplesheet
========================================================================================

    Columns
    -------
      sample    (required)  biological sample / aggregation key
      library   (optional)  one physical 10x library = one GEM well.
                            Defaults to `sample` when absent.
      fastq_1   (required)
      fastq_2   (required)
      fastq_3   (required for atac-family, forbidden for gex-family)
      modality  (optional)  gex | atac | arc_gex | arc_atac. Defaults to gex.

    Grouping rules
    --------------
      Rows sharing (sample, library, modality) are the SAME physical library
      sequenced more than once (top-up runs, extra flowcells). They are merged
      into a single Cell Ranger run as consecutive pseudo-lanes.

      Rows differing in `library` are DIFFERENT GEM wells. Barcodes are not
      comparable across them, so each gets its own Cell Ranger run. They stay
      linked through `meta.sample` for downstream aggregation.

----------------------------------------------------------------------------------------
*/

def VALID_MODALITIES = ['gex', 'atac', 'arc_gex', 'arc_atac']
def ATAC_LIKE        = ['atac', 'arc_atac']

// Turn one CSV row into [ meta, [fastq_1, fastq_2, fastq_3 or null] ]
def parse_row(row, row_num) {

    def err = { msg -> error("[samplesheet] line ${row_num}: ${msg}") }

    def sample = row.sample?.trim()
    if ( !sample ) err("'sample' is required")
    if ( sample =~ /\s/ ) err("'sample' must not contain spaces: '${sample}'")

    def modality = (row.modality?.trim() ?: 'gex').toLowerCase()
    if ( !(modality in VALID_MODALITIES) )
        err("unknown modality '${modality}'. Valid: ${VALID_MODALITIES.join(', ')}")

    // library defaults to sample -> the simple one-library-per-sample case
    def library = row.library?.trim() ?: sample
    if ( library =~ /\s/ ) err("'library' must not contain spaces: '${library}'")

    if ( !row.fastq_1?.trim() ) err("'fastq_1' is required")
    if ( !row.fastq_2?.trim() ) err("'fastq_2' is required (10x reads are always paired)")

    def fq1 = file(row.fastq_1.trim(), checkIfExists: true)
    def fq2 = file(row.fastq_2.trim(), checkIfExists: true)
    def fq3 = null

    def has_fq3 = row.fastq_3?.trim() as Boolean

    if ( modality in ATAC_LIKE ) {
        if ( !has_fq3 )
            err("modality '${modality}' requires 'fastq_3' (ATAC reads are R1 + R2[barcode] + R3)")
        fq3 = file(row.fastq_3.trim(), checkIfExists: true)
    }
    else if ( has_fq3 ) {
        err("modality '${modality}' must not set 'fastq_3'")
    }

    [ fq1, fq2, fq3 ].findAll { it }.each { f ->
        if ( !(f.name ==~ /.*\.f(ast)?q\.gz$/) )
            log.warn "[samplesheet] line ${row_num}: '${f.name}' does not look like a gzipped FASTQ"
    }

    def meta = [
        id       : library,      // unique key for one Cell Ranger run
        sample   : sample,       // biological unit, used for aggregation
        library  : library,
        modality : modality,
        single_library : (library == sample)
    ]

    [ meta, [ fq1, fq2, fq3 ] ]
}


workflow SAMPLESHEET {

    take:
    samplesheet     // path to csv

    main:

    def ch_rows = Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        // keep a 1-based line number (+1 for the header) for error messages
        .map { row -> [ row, 0 ] }
        .toList()
        .flatMap { rows ->
            rows.withIndex().collect { pair, idx -> parse_row(pair[0], idx + 2) }
        }

    // ---- merge re-sequencing runs of the same physical library --------------
    def ch_libraries = ch_rows
        .map { meta, fqs -> [ [ meta.sample, meta.library, meta.modality ], meta, fqs ] }
        .groupTuple()
        .map { key, metas, fqsets ->

            // deterministic lane order regardless of channel arrival order
            def sorted = fqsets.sort { a, b -> a[0].toString() <=> b[0].toString() }

            def r1 = sorted.collect { it[0] }
            def r2 = sorted.collect { it[1] }
            def r3 = sorted.collect { it[2] }.findAll { it }

            def meta = metas[0] + [ n_runs: sorted.size() ]

            if ( meta.n_runs > 1 )
                log.info "[samplesheet] ${meta.id} (${meta.modality}): merging ${meta.n_runs} sequencing runs as pseudo-lanes"

            [ meta, r1, r2, r3 ]
        }

    // ---- report multi-library samples ---------------------------------------
    ch_libraries
        .map { meta, r1, r2, r3 -> [ meta.sample, meta.id ] }
        .groupTuple()
        .filter { sample, libs -> libs.size() > 1 }
        .subscribe { sample, libs ->
            log.info "[samplesheet] sample '${sample}' has ${libs.size()} libraries (${libs.sort().join(', ')}); each is quantified separately"
        }

    // ---- split by modality ---------------------------------------------------
    def ch_branched = ch_libraries.branch { meta, r1, r2, r3 ->
        gex : meta.modality == 'gex'
        atac: meta.modality == 'atac'
        arc : meta.modality in ['arc_gex', 'arc_atac']
    }

    emit:
    gex       = ch_branched.gex
    atac      = ch_branched.atac
    arc       = ch_branched.arc
    libraries = ch_libraries
}
