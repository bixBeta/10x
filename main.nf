nextflow.enable.dsl=2

// Project Params:
params.sheet            = "sample-sheet.csv"

// Module Params:
params.help             = false
params.listRefs         = false
params.listPrograms     = false

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
    * --listPrograms   : List the 10x software installed under --programs
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

        MULTIOME ( --mode arc ): one label carries a Gene Expression library
        AND a Chromatin Accessibility library, told apart by library_type.
        They are not counted separately - cellranger-arc takes a libraries.csv
        naming both, which this pipeline writes for you.

        |-------|--------------------------------|-------------------------|
        | label | fastqs                         | library_type            |
        |-------|--------------------------------|-------------------------|
        | JS4   | /local/.../Project_10488522    | Gene Expression         |
        | JS4   | /local/.../Project_10488522    | Chromatin Accessibility |
        |-------|--------------------------------|-------------------------|

        becomes JS4_libraries.csv:

            fastqs,sample,library_type
            /local/.../Project_10488522,SC2619_JS4_BC_MG3_...,Gene Expression
            /local/.../Project_10488522,SC2619_JS4_MA_...,Chromatin Accessibility

        sample is optional. It is the fastq prefix, and when given it is used
        verbatim for --sample and for the sample column of libraries.csv. It is
        only REQUIRED when the fastqs path holds more than one prefix, e.g. a
        Project_* dir carrying both halves of a multiome pair - one row is one
        library, so that case is rejected rather than guessed at.

        library_type is optional and defaults to Gene Expression. It takes 10x's
        own vocabulary, so the value lands in libraries.csv verbatim:

            Gene Expression | Chromatin Accessibility | Antibody Capture
            CRISPR Guide Capture | Multiplexing Capture | VDJ-T | VDJ-B

        Shorthands are accepted: gex, rna, atac, adt, citeseq, hto, crispr,
        cmo, cellplex.
        -----------------------------------------------------------

    * --mode           : use 'gex'  for 3'/5' gene expression; default <gex>
                       : use 'arc'  for multiome GEX + ATAC
                       : use 'atac' for scATAC                 ( v0.2, not yet implemented )

    * --ref            : 10x reference. Use --listRefs to see all available references.
                         Also supports a path value for a cellranger transcriptome dir.
    * --chemistry      : Cell Ranger chemistry; default <auto>
    * --expectCells    : Expected number of recovered cells; default <null> ( Cell Ranger estimates )
    * --forceCells     : Force pipeline to use this number of cells
    * --createBam      : Emit the position-sorted BAM; default <false>
    * --introns        : Set true/false to override --include-introns; default <null> ( tool default )
    * --r1length       : Trim R1 to this length; default <28>. Pass 0 to omit the flag.
    * --r2length       : Trim R2 to this length; default <null>

  10x software ( local install, see --listPrograms ):
    * --crversion      : Cell Ranger version; default <9.0.1>
                         -> /programs/cellranger-<crversion>/cellranger
    * --arcversion     : cellranger-arc version; default <2.2.0>
                         -> /programs/cellranger-arc-<arcversion>/bin/cellranger-arc
    * --programs       : Where the installs live; default </programs>
    * --crpath         : Full path to the cellranger binary. Overrides crversion.
    * --arcpath        : Full path to the cellranger-arc binary. Overrides arcversion.

                         The resolved binary is checked before anything runs, so a
                         version that is not installed fails immediately and lists
                         what is.

  Containers ( optional, off by default ):
    * --container      : Run gex inside this image e.g. a .sif path or docker:// uri
    * --arccontainer   : Same for arc mode

  Runtime / resources:
    * --localcores     : cellranger --localcores, also the cpus reserved; default <32>
    * --localmem       : cellranger --localmem in GB, also the memory reserved; default <180>
    * --maxforks       : concurrent tasks PER PROCESS, so 2 means up to 2 cellranger
                         runs at once, each taking localcores / localmem; default <2>

  The command this builds, gex mode, one per library:

    <crbin> count --id=<label> \\
      --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam> --r1-length=<r1length> \\
      --transcriptome=<ref> \\
      --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]>

  arc mode, one per label:

    <arcbin> count --id=<label> \\
      --reference=<ref> --libraries=<label>_libraries.csv \\
      --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam>

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
cellranger   : ${params.mode == "arc" ? params.arcbin : params.crbin}
container    : ${params.mode == "arc" ? (params.arccontainer ?: "none, local install") : (params.container ?: "none, local install")}
"""

// The 10x binary is a local install. Fail here, with the resolved path and what
// is actually installed, rather than a task dying much later. Skipped for a
// container run ( the binary lives in the image ) and for -stub-run.
def checkProgram(bin, containerParam, what) {

    if( containerParam )       return
    if( workflow.stubRun )     return

    def f = file(bin)
    if( f.exists() ) return

    def installed = file("${params.programs}/cellranger*")
    def avail = ( installed instanceof List ? installed : ( installed ? [installed] : [] ) )
                    .collect { it.name }.sort()

    error """No ${what} found at: ${bin}
    Installed under ${params.programs}: ${ avail ? avail.join(', ') : 'nothing matching cellranger*' }
    Set --${what == 'cellranger-arc' ? 'arcversion' : 'crversion'} to an installed version, or give the full path with --${what == 'cellranger-arc' ? 'arcpath' : 'crpath'}."""
}


if( params.listPrograms ) {

    println("")
    log.info """
    10x software installed under ${params.programs}
    =========================================================================================================================
    """
    .stripIndent()

    def installed = file("${params.programs}/cellranger*")
    def dirs = ( installed instanceof List ? installed : ( installed ? [installed] : [] ) ).sort()

    if( !dirs ) {
        println "nothing matching ${params.programs}/cellranger*"
    }
    else {
        dirs.each { d ->
            def name = d.name
            def bin  = file("${d}/${name.replaceFirst(/-[0-9].*$/, '')}")
            def bbin = file("${d}/bin/${name.replaceFirst(/-[0-9].*$/, '')}")
            def hit  = bin.exists() ? bin : ( bbin.exists() ? bbin : null )
            println "${name} ----------- ${ hit ?: 'no executable found' }"
        }
    }

    exit 0
}

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
CanFam3             :"/local/workdir/10x_analysis/REFS/Canine/arc/CanFam3_Ensembl101annot",
GRCh38              :"/local/workdir/10x_analysis/REFS/Human/arc/GRCh38",
GRCm39              :"/local/workdir/10x_analysis/REFS/Mouse/arc/GRCm39"]


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
include {   CELLRANGER_ARC_COUNT     } from './modules/cellrangerarc'
include {   DUMP_VERSIONS            } from './modules/versions'


// allows a user to pass a reference path via --ref, same as --genome elsewhere.
// each mode has its own reference set, since gex / atac / arc references differ.
refMap = params.mode == "arc"  ? refDirArc
       : params.mode == "atac" ? refDirAtac
       :                         refDir

if( refMap.containsKey(params.ref) ){

    ref = refMap[params.ref]

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

    // A path may be the Sample_* dir itself, or the Project_* dir above it that
    // holds Sample_*/ subdirs - both are valid --fastqs values for cellranger,
    // so look one level down as well.
    def hits = []
    [ "${fq}/*_R1*.f*q.gz", "${fq}/*/*_R1*.f*q.gz" ].each { pattern ->
        def found = file(pattern)
        if( found instanceof List ) hits.addAll(found)
        else if( found )            hits.add(found)
    }

    def names = hits.collect { it.name }
                    .findAll { it ==~ /.*_S\d+_L\d+_R1(_001)?\.f(ast)?q\.gz$/ }
                    .collect { it.replaceFirst(/_S\d+_L\d+_R1(_001)?\.f(ast)?q\.gz$/, '') }
                    .unique()
                    .sort()

    if( !names )
        error "sample-sheet: no *_S<n>_L<n>_R1_001.fastq.gz files found in ${fq} or its immediate subdirectories"

    names
}

// library_type uses 10x's own vocabulary, so the value can be written straight
// into the libraries.csv that cellranger-arc / cellranger multi consume. Common
// shorthands are accepted and normalised to the canonical string.
def libraryType(raw, label) {

    def canonical = [
        "gene expression"       : "Gene Expression",
        "gex"                   : "Gene Expression",
        "rna"                   : "Gene Expression",
        "antibody capture"      : "Antibody Capture",
        "antibody"              : "Antibody Capture",
        "adt"                   : "Antibody Capture",
        "citeseq"               : "Antibody Capture",
        "cite-seq"              : "Antibody Capture",
        "hto"                   : "Antibody Capture",
        "crispr guide capture"  : "CRISPR Guide Capture",
        "crispr"                : "CRISPR Guide Capture",
        "multiplexing capture"  : "Multiplexing Capture",
        "cmo"                   : "Multiplexing Capture",
        "cellplex"              : "Multiplexing Capture",
        "chromatin accessibility": "Chromatin Accessibility",
        "atac"                  : "Chromatin Accessibility",
        "vdj-t"                 : "VDJ-T",
        "vdj-b"                 : "VDJ-B",
        "vdj"                   : "VDJ"
    ]

    if( !raw ) return "Gene Expression"

    def hit = canonical[ raw.toLowerCase() ]
    if( !hit )
        error "sample-sheet: ${label} has unknown library_type '${raw}'. Use one of: ${canonical.values().unique().sort().join(', ')}"

    hit
}


// the same prefix, taken from a single file name
def sampleFromFile(name) {

    def m = name =~ /^(.+)_S\d+_L\d+_R\d(_001)?\.f(ast)?q\.gz$/
    m ? m[0][1] : null
}

// The sheet accepts either the delivery dir or any fastq inside it, whichever
// is easier to paste. A file resolves to its parent dir, since --fastqs always
// takes a directory.
def resolveFastqs(fq, label, given) {

    def d = file(fq)

    if( !d.exists() )
        error "sample-sheet: ${label} path does not exist: ${fq}"

    // a file resolves to its parent, since --fastqs takes a directory
    def dir = d.isDirectory() ? fq : d.parent.toString()

    // An explicit sample column wins: it is what lands in --sample and in
    // libraries.csv, so the sheet can name it exactly like 10x does.
    if( given ) {
        def found = file("${dir}/${given}_S*_L*_R1*.f*q.gz")
        def any   = (found instanceof List) ? found : (found ? [found] : [])
        if( !any ) {
            def sub = file("${dir}/*/${given}_S*_L*_R1*.f*q.gz")
            any = (sub instanceof List) ? sub : (sub ? [sub] : [])
        }
        if( !any )
            error "sample-sheet: ${label}: no reads for sample '${given}' under ${dir}"
        return [ dir, [ given ] ]
    }

    if( !d.isDirectory() ) {
        def pfx = sampleFromFile(d.name)
        if( !pfx )
            error "sample-sheet: ${label}: cannot read a 10x sample prefix from '${d.name}'. Expected <prefix>_S<n>_L<n>_R1_001.fastq.gz, or give the delivery directory instead."
        return [ dir, [ pfx ] ]
    }

    // One row is one library, so one prefix. A Project_* dir holding several
    // samples is ambiguous and has to be pinned down with the sample column,
    // rather than quietly counting the wrong reads together.
    def names = sampleNames(dir)
    if( names.size() > 1 )
        error "sample-sheet: ${label}: ${dir} holds ${names.size()} sample prefixes (${names.join(', ')}). Add a sample column naming the one this row means."

    [ dir, names ]
}


def readRows() {

    ch_sheet
        | splitCsv( header:true )
        | map { row ->
              def label   = row.label?.trim()
              def library = row.library?.trim() ?: label
              def fq      = row.fastqs?.trim()

              if( !label ) error "sample-sheet: every row needs a label"
              if( !fq )    error "sample-sheet: ${label} is missing fastqs ( a 10x delivery dir, or any fastq in it )"

              def ltype = libraryType(row.library_type?.trim(), label)
              def res   = resolveFastqs(fq, label, row.sample?.trim())

              [ label, library, res[0], res[1], ltype ]
          }
}


/* ---------------------------------------------------------------------------------------------------------
GENE EXPRESSION Workflow
------------------------------------------------------------------------------------------------------------ */

workflow GEX {

    checkProgram(params.crbin, params.container, "cellranger")

    if( ref == null ){
        error "No reference provided. Use --ref < key or path >, see --listRefs"
    }

    meta_ch = readRows()
        | map { label, library, dir, samples, ltype ->
              if( ltype != "Gene Expression" )
                  error "sample-sheet: ${label} is '${ltype}'. --mode gex counts Gene Expression only; use --mode arc for multiome."
              [ [label, library], [ dir, samples ] ]
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

    CELLRANGER_COUNT(meta_ch, ref_ch)

    ch_versions = CELLRANGER_COUNT.out.versions.first()

    DUMP_VERSIONS(ch_versions.collect())
}


/* ---------------------------------------------------------------------------------------------------------
MULTIOME Workflow

One label carries a Gene Expression library and a Chromatin Accessibility
library. They are NOT counted separately: cellranger-arc takes a libraries.csv
listing both, which is generated per label from the sheet.
------------------------------------------------------------------------------------------------------------ */

workflow ARC {

    checkProgram(params.arcbin, params.arccontainer, "cellranger-arc")

    if( ref == null ){
        error "No reference provided. Use --ref < key or path >, see --listRefs"
    }

    arc_ch = readRows()
        | map { label, library, dir, samples, ltype -> [ label, [ dir, samples, ltype ] ] }
        | groupTuple
        | map { label, entries ->

              def types = entries.collect { it[2] }.unique()

              def unsupported = types.findAll { !( it in ["Gene Expression", "Chromatin Accessibility"] ) }
              if( unsupported )
                  error "arc: ${label} has library_type(s) cellranger-arc cannot take: ${unsupported.join(', ')}"

              if( !types.contains("Gene Expression") )
                  error "arc: ${label} has no Gene Expression library. Multiome needs both, see --help"
              if( !types.contains("Chromatin Accessibility") )
                  error "arc: ${label} has no Chromatin Accessibility library. Multiome needs both, see --help"

              // one libraries.csv row per fastq dir + sample prefix
              def rows = entries.collectMany { e -> e[1].collect { s -> [ e[0], s, e[2] ] } }
                                .unique()
                                .sort { a1, b1 -> a1[2] <=> b1[2] ?: a1[1] <=> b1[1] }

              [ label, rows ]
          }
        | view { label, rows ->
              "MULTIOME >> ${label}  ( libraries: " + rows.collect{ "${it[1]} [${it[2]}]" }.join(' + ') + " )" }

    CELLRANGER_ARC_COUNT(arc_ch, ref_ch)

    DUMP_VERSIONS( CELLRANGER_ARC_COUNT.out.versions.first().collect() )
}


workflow {

    if( params.mode == "gex" ){

        GEX()
    }

    else if( params.mode == "atac" ){

        error "atac mode is not implemented yet ( v0.2 )"
    }

    else if( params.mode == "arc" ){

        ARC()
    }

    else {

        error "Invalid mode provided: ${params.mode}"
    }

}
