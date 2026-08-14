nextflow.enable.dsl=2

// Project Params:
params.sheet            = "sample-sheet.csv"

// Module Params:
params.help             = false
params.listRefs         = false

// Default Params:
params.mode             = "gex"
params.id               = "TREx_ID"
params.ref              = null
params.chemistry        = "auto"
params.expectCells      = null
params.forceCells       = null
params.createBam        = false
params.introns          = null
params.r1length         = 28
params.r2length         = null

// --crversion, --container, --crpath, --maxforks, --localcores and --localmem
// are declared in nextflow.config, since process directives are resolved before
// this script is parsed.

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
        |-------|------------------------------------------------------------|
        | label | fastqs                                                     |
        |-------|------------------------------------------------------------|
        | JS4   | /local/Illumina/DRV/<run>/Unaligned/Project_10488923/Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2 |
        | JS5   | /local/Illumina/DRV/<run>/Unaligned/Project_10488923/Sample_SC2620_JS5_G3_Sofia_10488923_253GGLLT4_L2 |
        |-------|------------------------------------------------------------|

        fastqs is either the 10x delivery directory or ANY fastq inside it,
        whichever is easier to paste. A file resolves to its parent dir, since
        --fastqs always takes a directory. Nothing is copied, staged or renamed.
        Both of these are the same library:

            /local/.../Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2
            /local/.../Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2/SC2620_JS4_G3_Reign_10488923_253GGLLT4_S3_L002_R1_001.fastq.gz

        label is yours: short, user defined, and used for --id and for the
        output folder. It is never derived from the path.

        --sample is read off the files in the directory, since cellranger has
        to be given the fastq prefix:

            Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2/
              SC2620_JS4_G3_Reign_10488923_253GGLLT4_S3_L002_R1_001.fastq.gz
              -> --sample=SC2620_JS4_G3_Reign_10488923_253GGLLT4

        Same library sequenced more than once ( top-ups, extra flow cells ):
        repeat the label. The dirs are passed as ONE comma separated --fastqs,
        which is how cellranger pools them.

        |-------|--------------------------------------------|
        | label | fastqs                                     |
        |-------|--------------------------------------------|
        | JS4   | /local/Illumina/DRV/run1/.../Sample_SC2620_JS4_..._L2 |
        | JS4   | /local/Illumina/DRV/run2/.../Sample_SC2620_JS4_..._L3 |
        |-------|--------------------------------------------|

        Several libraries from one sample ( separate GEM wells ): add the
        optional library column. Barcodes are not comparable across GEM wells,
        so each library is counted separately and stays tied to its label.

        |-------|-----------|------------------------------------|
        | label | library   | fastqs                             |
        |-------|-----------|------------------------------------|
        | JS4   | JS4_wellA | /local/.../Sample_SC2620_JS4A_..._L2 |
        | JS4   | JS4_wellB | /local/.../Sample_SC2620_JS4B_..._L2 |
        |-------|-----------|------------------------------------|

        ATAC ( --mode atac, v0.2 ) and MULTIOME ( --mode arc, v0.3 ) use the
        same directory based sheet. R1 / R2 / R3 live inside the delivery dir,
        so there is no fastq3 column to fill in. arc adds a type column
        ( gex | atac ) so one label can carry both libraries.

        |-------|-----------|------|------------------------------|
        | label | library   | type | fastqs                       |
        |-------|-----------|------|------------------------------|
        | JS4   | JS4_gex   | gex  | /local/.../Sample_..._GEX_L2  |
        | JS4   | JS4_atac  | atac | /local/.../Sample_..._ATAC_L2 |
        |-------|-----------|------|------------------------------|
        -----------------------------------------------------------

    * --mode           : use 'gex'  for 3'/5' gene expression; default <gex>
                       : use 'atac' for scATAC                 ( v0.2, not yet implemented )
                       : use 'arc'  for multiome GEX + ATAC    ( v0.3, not yet implemented )

    * --ref            : 10x reference. Use --listRefs to see all available references.
                         Also supports a path value for a cellranger transcriptome dir.
    * --chemistry      : Cell Ranger chemistry; default <auto>
    * --expectCells    : Expected number of recovered cells; default <null> ( Cell Ranger estimates )
    * --forceCells     : Force pipeline to use this number of cells
    * --createBam      : Emit the position-sorted BAM; default <false>
    * --introns        : Set true/false to override --include-introns; default <null> ( tool default )
    * --r1length       : Trim R1 to this length; default <28>. Pass 0 to omit the flag.
    * --r2length       : Trim R2 to this length; default <null>

  Runtime / resources:
    * --localcores     : cellranger --localcores, also the cpus reserved; default <32>
    * --localmem       : cellranger --localmem in GB, also the memory reserved; default <180>
    * --maxforks       : how many processes run at once, pipeline wide; default <2>
    * --crversion      : Cell Ranger version, selects the container tag; default <9.0.1>
    * --crpath         : Run a native install instead of the container
                         e.g. --crpath /programs/cellranger-9.0.1/cellranger
    * --container      : Full container override e.g. a local .sif path

  The command this builds per library:

    <crpath|cellranger> count --id=<label> \\
      --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam> --r1-length=<r1length> \\
      --transcriptome=<ref> \\
      --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]>

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
r1length     : ${params.r1length}
localcores   : ${params.localcores}
localmem     : ${params.localmem}
maxforks     : ${params.maxforks}
crversion    : ${params.crversion}
cellranger   : ${params.crpath ?: (params.container ?: "docker://ghcr.io/bixbeta/cellranger:" + params.crversion)}
"""

// 10x reference MAP  ( /local/workdir/10x_analysis/REFS )
refDir = [
CanFam3_1           :"/local/workdir/10x_analysis/REFS/Canine/CanFam3_1",
GRCh38              :"/local/workdir/10x_analysis/REFS/Human/GRCh38",
GRCm39              :"/local/workdir/10x_analysis/REFS/Mouse/GRCm39"]

// ATAC / ARC references, wired up when those modes land
refDirAtac = [
GRCh38              :"/local/workdir/10x_analysis/REFS/ATAC/Human/GRCh38",
GRCm39              :"/local/workdir/10x_analysis/REFS/ATAC/Mouse/GRCm39"]

refDirArc = [
GRCh38              :"/local/workdir/10x_analysis/REFS/ARC/Human/GRCh38",
GRCm39              :"/local/workdir/10x_analysis/REFS/ARC/Mouse/GRCm39"]


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

// cellranger --sample must match the fastq prefix INSIDE the delivery dir, which
// is the dir name without the Sample_ prefix and the trailing lane suffix:
//
//   Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2/
//     SC2620_JS4_G3_Reign_10488923_253GGLLT4_S3_L002_R1_001.fastq.gz
//     -> SC2620_JS4_G3_Reign_10488923_253GGLLT4
//
// It is read off the files rather than parsed out of the dir name, so an
// unexpected naming variant fails loudly here instead of inside cellranger.
def sampleNames(fq) {

    def names = file(fq).list()
                    .findAll { it ==~ /.*_S\d+_L\d+_R1(_001)?\.f(ast)?q\.gz$/ }
                    .collect { it.replaceFirst(/_S\d+_L\d+_R1(_001)?\.f(ast)?q\.gz$/, '') }
                    .unique()
                    .sort()

    if( !names )
        error "sample-sheet: no *_S<n>_L<n>_R1_001.fastq.gz files found in ${fq}"

    names
}

// the same prefix, taken from a single file name
def sampleFromFile(name) {

    def m = name =~ /^(.+)_S\d+_L\d+_R\d(_001)?\.f(ast)?q\.gz$/
    m ? m[0][1] : null
}

// The sheet accepts either the delivery dir or any fastq inside it, whichever
// is easier to paste. A file resolves to its parent dir, since --fastqs always
// takes a directory.
def resolveFastqs(fq, label) {

    def d = file(fq)

    if( !d.exists() )
        error "sample-sheet: ${label} path does not exist: ${fq}"

    if( d.isDirectory() )
        return [ fq, sampleNames(fq) ]

    def pfx = sampleFromFile(d.name)
    if( !pfx )
        error "sample-sheet: ${label}: cannot read a 10x sample prefix from '${d.name}'. Expected <prefix>_S<n>_L<n>_R1_001.fastq.gz, or give the delivery directory instead."

    [ d.parent.toString(), [ pfx ] ]
}


def readSheet() {

    ch_sheet
        | splitCsv( header:true )
        | map { row ->
              def label   = row.label?.trim()
              def library = row.library?.trim() ?: label
              def fq      = row.fastqs?.trim()

              if( !label ) error "sample-sheet: every row needs a label"
              if( !fq )    error "sample-sheet: ${label} is missing fastqs ( a 10x delivery dir, or any fastq in it )"

              [ [label, library], resolveFastqs(fq, label) ]
          }
        | groupTuple
        | map { key, entries ->
              // several rows for one library = the same library sequenced more
              // than once. cellranger takes both as comma separated lists, and
              // the prefixes differ per flow cell so every one has to be listed.
              // dirs are deduped: two rows may name two files in one dir.
              def dirs    = entries.collect { it[0] }.unique().sort()
              def samples = entries.collectMany { it[1] }.unique().sort()

              [ key[1], key[0], dirs, samples ]
          }
        | view { library, label, dirs, samples ->
              "LIBRARY >> ${library}  ( label: ${label}, fastq dirs: ${dirs.size()}, sample: ${samples.join(',')} )" }
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
