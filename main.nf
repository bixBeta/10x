nextflow.enable.dsl=2

// Project Params:
params.sheet            = "sample-sheet.csv"

// Module Params:
params.help             = false
params.listRefs         = false
params.fastqs           = null

// Default Params:
params.mode             = "gex"
params.id               = "TREx_ID"
params.ref              = null
params.chemistry        = "auto"
params.expectCells      = null
params.forceCells       = null
params.createBam        = false
params.introns          = null

// Command Line Channels     ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~

ch_pin          =    channel.value(params.id)
ch_sheet        =    channel.fromPath(params.sheet)


// ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~

if( params.help ) {

log.info """
1  0  x      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixBeta/10x -r main < args ... >

Args:
    * --help           : Prints this help documentation
    * --listRefs       : Get extended list of 10x references available for this pipeline
    * --id             : TREx Project ID
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -----------------------------------------------------------
        Sample Sheet Example: ( comma delimited file )
        |-------|-----------------|-----------------|
        | label | fastq1          | fastq2          |
        |-------|-----------------|-----------------|
        | SS1   | SS1_R1.fastq.gz | SS1_R2.fastq.gz |
        | SS2   | SS2_R1.fastq.gz | SS2_R2.fastq.gz |
        |-------|-----------------|-----------------|

        One row per fastq pair. Fastq names are arbitrary; they get renamed
        into the bcl2fastq convention Cell Ranger expects.

        Same library sequenced more than once ( top-ups, extra flow cells ):
        repeat the label. The rows are pooled into ONE cellranger count as
        consecutive lanes.

        |-------|-----------------------|-----------------------|
        | label | fastq1                | fastq2                |
        |-------|-----------------------|-----------------------|
        | SS1   | run1/SS1_R1.fastq.gz  | run1/SS1_R2.fastq.gz  |
        | SS1   | run2/SS1_R1.fastq.gz  | run2/SS1_R2.fastq.gz  |
        |-------|-----------------------|-----------------------|

        Several libraries from one sample ( separate GEM wells ): add the
        optional library column. Barcodes are not comparable across GEM wells,
        so each library is counted separately and stays tied to its label.

        |-------|-----------|-----------------|-----------------|
        | label | library   | fastq1          | fastq2          |
        |-------|-----------|-----------------|-----------------|
        | SS1   | SS1_wellA | A_R1.fastq.gz   | A_R2.fastq.gz   |
        | SS1   | SS1_wellB | B_R1.fastq.gz   | B_R2.fastq.gz   |
        |-------|-----------|-----------------|-----------------|
        -----------------------------------------------------------

    * --mode           : use 'gex'  for 3'/5' gene expression; default <gex>
                       : use 'atac' for scATAC                 ( v0.2, not yet implemented )
                       : use 'arc'  for multiome GEX + ATAC    ( v0.3, not yet implemented )
    * --ref            : 10x reference. Use --listRefs to see all available references.
                         Also supports a path value for a cellranger transcriptome dir.
    * --fastqs         : Use this param if fastq files are in the fastqs folder in the project directory;
                         If --fastqs is not specified, the fastqs must be supplied with absolute paths in the sample-sheet.csv
    * --chemistry      : Cell Ranger chemistry; default <auto>
    * --expectCells    : Expected number of recovered cells; default <null> ( Cell Ranger estimates )
    * --forceCells     : Force pipeline to use this number of cells
    * --createBam      : Emit the position-sorted BAM; default <false>
    * --introns        : Set true/false to override --include-introns; default <null> ( tool default )
    * --crversion      : Cell Ranger version to run, selects the container tag; default <9.0.1>
    * --container      : Full container override e.g. a local .sif path

"""

    exit 0
}


log.info """
1  0  x      W  O  R  K  F  L  O  W  -  @bixBeta
=========================================================================================================================
trexID       : ${params.id}
sheet        : ${params.sheet}
mode         : ${params.mode}
ref          : ${params.ref}
chemistry    : ${params.chemistry}
expectCells  : ${params.expectCells}
createBam    : ${params.createBam}
crversion    : ${params.crversion}
container    : ${params.container ?: "docker://ghcr.io/bixbeta/cellranger:" + params.crversion}
"""

// 10x reference MAP
refDir = [
GRCh38              :"/workdir/genomes/Homo_sapiens/hg38/10x/refdata-gex-GRCh38-2024-A",
GRCm39              :"/workdir/genomes/Mus_musculus/GRCm39/10x/refdata-gex-GRCm39-2024-A",
GRCh38_GRCm39       :"/workdir/genomes/Multi/10x/refdata-gex-GRCh38_and_GRCm39-2024-A"]

// ATAC / ARC references, wired up when those modes land
refDirAtac = [
GRCh38              :"/workdir/genomes/Homo_sapiens/hg38/10x/refdata-cellranger-atac-GRCh38-2020-A-2.0.0",
GRCm39              :"/workdir/genomes/Mus_musculus/GRCm39/10x/refdata-cellranger-atac-mm10-2020-A-2.0.0"]

refDirArc = [
GRCh38              :"/workdir/genomes/Homo_sapiens/hg38/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0",
GRCm39              :"/workdir/genomes/Mus_musculus/GRCm39/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0"]


if( params.listRefs ) {

    println("")
    log.info """
    Available 10x GEX References
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    refDir.each(printMap)

    log.info """
    Available 10x ATAC References
    =========================================================================================================================
    """
    .stripIndent()

    refDirAtac.each(printMap)

    log.info """
    Available 10x ARC References
    =========================================================================================================================
    """
    .stripIndent()

    refDirArc.each(printMap)

    exit 0
}

include {   CELLRANGER_COUNT         } from './modules/cellranger'
include {   DUMP_VERSIONS            } from './modules/versions'


// allows a user to pass a transcriptome path via --ref, same as --genome elsewhere
if( refDir.containsKey(params.ref) ){

    ref = refDir[params.ref]

} else {

    ref = params.ref
}

ref_ch = channel.value(ref)


/* ---------------------------------------------------------------------------------------------------------
Sample sheet -> one entry per physical 10x library

Rows sharing ( label, library ) are the same library sequenced more than once and
are pooled as consecutive lanes. Rows differing in library are separate GEM wells
and are counted separately.
------------------------------------------------------------------------------------------------------------ */

def readSheet() {

    def dir = params.fastqs ? "fastqs/" : ""

    ch_sheet
        | splitCsv( header:true )
        | map { row ->
              def label   = row.label?.trim()
              def library = row.library?.trim() ?: label

              if( !label )       error "sample-sheet: every row needs a label"
              if( !row.fastq1 )  error "sample-sheet: ${label} is missing fastq1"
              if( !row.fastq2 )  error "sample-sheet: ${label} is missing fastq2 ( 10x reads are always paired )"

              [ [label, library], [ file(dir + row.fastq1.trim()), file(dir + row.fastq2.trim()) ] ]
          }
        | groupTuple
        | map { key, pairs ->
              // sort so lane order does not depend on channel arrival order
              def sorted = pairs.sort { a, b -> a[0].toString() <=> b[0].toString() }
              [ key[1], key[0], sorted.collect{ it[0] }, sorted.collect{ it[1] } ]
          }
        | view { library, label, r1, r2 -> "LIBRARY >> ${library}  ( sample: ${label}, runs: ${r1.size()} )" }
}


/* ---------------------------------------------------------------------------------------------------------
GENE EXPRESSION Workflow
------------------------------------------------------------------------------------------------------------ */

workflow GEX {

    if( ref == null ){
        error "No reference provided. Use --ref < key or path >, see --listRefs"
    }

    meta_ch = readSheet()

    CELLRANGER_COUNT(meta_ch, ref_ch)

    ch_versions = CELLRANGER_COUNT.out.versions.first()

    DUMP_VERSIONS(ch_versions.collect())
}


workflow {

    if( params.mode == "gex" ){

        GEX()
    }

    else if( params.mode == "atac" ){

        error "atac mode is not implemented yet ( v0.2 )"
    }

    else if( params.mode == "arc" ){

        error "arc mode is not implemented yet ( v0.3 )"
    }

    else {

        error "Invalid mode provided: ${params.mode}"
    }

}
