nextflow.enable.dsl=2

// Project Params:
params.sheet            = "sample-sheet.csv"

// Module Params:
params.help             = false
params.listRefs         = false
params.listPrograms     = false

// Default Params:
params.mode             = "gex"
params.id               = "BRC_ID"
params.ref              = null
params.chemistry        = "auto"
params.expectCells      = null
params.forceCells       = null
params.createBam        = false
params.introns          = null
params.r1length         = 28
params.r2length         = null

// atac only
params.peaks            = null
params.dimReduce        = null

// --crversion, --container, --crpath, --maxforks, --localcores and --localmem
// are declared in nextflow.config, since process directives are resolved before
// this script is parsed.

// Command Line Channels     ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~

ch_pin          =    channel.value(params.id)
ch_sheet        =    channel.fromPath(params.sheet)
ch_mqc_config   =    channel.value( file( params.multiqcconfig ?: "${projectDir}/assets/multiqc_config.yaml",
                                          checkIfExists: true ) )


// ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~

if( params.help ) {

log.info """
1  0  x      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run bixBeta/10x -r main -params-file params.yaml
    nextflow run bixBeta/10x -r main < args ... >

    params.yaml ships with the repo and carries every default, commented.
    Anything on the command line overrides it.

Args:
    * --help           : Prints this help documentation
    * --listRefs       : Get extended list of 10x references available for this pipeline
    * --listPrograms   : List the 10x software installed under --programs
    * --id             : BRC Project ID
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -----------------------------------------------------------
        Sample Sheet Example: ( comma delimited file )
        |-------|------------------------------------------------------------|
        | label | fastqs                                                     |
        |-------|------------------------------------------------------------|
        | JS4   | /local/Illumina/DRV/<run>/Unaligned/Project_104***/Sample_SC2620_JS4_G3_***_104***_253***_L2 |
        | JS5   | /local/Illumina/DRV/<run>/Unaligned/Project_104***/Sample_SC2620_JS5_G3_***_104***_253***_L2 |
        |-------|------------------------------------------------------------|

        fastqs is either the 10x delivery directory or ANY fastq inside it,
        whichever is easier to paste. A file resolves to its parent dir, since
        --fastqs always takes a directory. Nothing is copied, staged or renamed.
        Both of these are the same library:

            /local/.../Sample_SC2620_JS4_G3_***_104***_253***_L2
            /local/.../Sample_SC2620_JS4_G3_***_104***_253***_L2/SC2620_JS4_G3_***_104***_253***_S3_L002_R1_001.fastq.gz

        label is yours: short, user defined, and used for --id and for the
        output folder. It is never derived from the path.

        --sample is read off the files in the directory, since cellranger has
        to be given the fastq prefix:

            Sample_SC2620_JS4_G3_***_104***_253***_L2/
              SC2620_JS4_G3_***_104***_253***_S3_L002_R1_001.fastq.gz
              -> --sample=SC2620_JS4_G3_***_104***_253***

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
        | JS4   | /local/.../Project_104***    | Gene Expression         |
        | JS4   | /local/.../Project_104***    | Chromatin Accessibility |
        |-------|--------------------------------|-------------------------|

        becomes JS4_libraries.csv:

            fastqs,sample,library_type
            /local/.../Project_104***,SC2619_JS4_BC_MG3_...,Gene Expression
            /local/.../Project_104***,SC2619_JS4_MA_...,Chromatin Accessibility

        scATAC ( --mode atac ): one row per library, the same shape as gex -
        cellranger-atac takes --fastqs / --sample, not a libraries.csv. Only
        Chromatin Accessibility is accepted, and library_type can be left out
        since that is the only thing this mode counts.

        |-------|----------------------------------------|-------------------------|
        | label | fastqs                                 | library_type            |
        |-------|----------------------------------------|-------------------------|
        | CK1   | /local/.../Sample_SC2617_CK1_MA_..._L1  | Chromatin Accessibility |
        |-------|----------------------------------------|-------------------------|

        To count the ATAC half of a MULTIOME run on its own, use --mode atac
        with --chemistry ARC-v1 and an arc reference. --mode arc cannot do it:
        cellranger-arc needs both libraries and has no single assay mode.

        sample is optional. It is the fastq prefix - the part BEFORE the
        _S<n>_L<n>_R<n>_001 tail, not the file name - and when given it is used
        verbatim for --sample and for the sample column of libraries.csv. It is
        only REQUIRED when the fastqs path holds more than one prefix, e.g. a
        Project_* dir carrying both halves of a multiome pair - one row is one
        library, so that case is rejected rather than guessed at.

        library_type is optional and defaults to whatever the mode counts:
        Gene Expression, or Chromatin Accessibility under --mode atac. It takes
        10x's own vocabulary, so the value lands in libraries.csv verbatim:

            Gene Expression | Chromatin Accessibility | Antibody Capture
            CRISPR Guide Capture | Multiplexing Capture | VDJ-T | VDJ-B

        Shorthands are accepted: gex, rna, atac, adt, citeseq, hto, crispr,
        cmo, cellplex.
        -----------------------------------------------------------

    * --mode           : use 'gex'  for 3'/5' gene expression; default <gex>
                       : use 'arc'  for multiome GEX + ATAC
                       : use 'atac' for scATAC ( also the ATAC half of a multiome
                         run on its own, with --chemistry ARC-v1 )

    * --ref            : 10x reference. Use --listRefs to see all available references.
                         Also supports a path value for a cellranger transcriptome dir.
    * --chemistry      : Cell Ranger chemistry; default <auto>
                         atac: left to the tool unless set. 'ARC-v1' counts the
                         ATAC library of a multiome experiment on its own.
    * --expectCells    : Expected number of recovered cells; default <null> ( Cell Ranger estimates )
    * --forceCells     : Force pipeline to use this number of cells
    * --createBam      : Emit the position-sorted BAM; default <false>
    * --introns        : Set true/false to override --include-introns; default <null> ( tool default )
    * --r1length       : Trim R1 to this length; default <28>. Pass 0 to omit the flag.
    * --r2length       : Trim R2 to this length; default <null>
    * --peaks          : atac only. BED of peaks to use instead of calling them;
                         default <null> ( cellranger-atac calls its own )
    * --dimReduce      : atac only. cellranger-atac --dim-reduce: lsa | plsa | pca
                         default <null> ( tool default, lsa )
    * --multiqc        : Aggregate the web summaries into one MultiQC report; default <true>
                         MultiQC reads cellranger / cellranger-arc web_summary.html,
                         so the report covers whatever was counted this run.
    * --multiqcversion : MultiQC version, selects the image; default <1.35>
    * --multiqcsif     : Full path to a specific multiqc image
    * --multiqcpath    : Run a multiqc binary instead of an image ( default: multiqc on PATH )
    * --multiqcconfig  : Custom multiqc config yaml
    * --aggr           : Run cellranger aggr for labels that carry several libraries
                         ( separate GEM wells ); default <true>. A label with one
                         library is left alone: re-sequencing runs were already
                         pooled into its single count.
    * --normalize      : cellranger aggr --normalize; default <mapped>
                       : 'mapped' subsamples deeper libraries to match the shallowest
                       : 'none'   keeps every read

  10x software:
    * --crversion      : Cell Ranger version; default <9.0.1>
    * --arcversion     : cellranger-arc version; default <2.2.0>
    * --atacversion    : cellranger-atac version; default <2.2.0>
    * --engine         : 'singularity' runs a locally built .sif; default <singularity>
                       : 'local'       runs an install under --programs

      singularity ( default ):
        * --sifdir     : Where the images live; default </local/workdir/singularity>
                         -> <sifdir>/cellranger-<crversion>.sif
                         -> <sifdir>/cellranger-arc-<arcversion>.sif
                         -> <sifdir>/cellranger-atac-<atacversion>.sif
        * --crsif      : Full path to the gex image. Overrides sifdir + crversion.
        * --arcsif     : Full path to the arc image.
        * --atacsif    : Full path to the atac image.

                         Build one with:
                         containers/build-sif.sh cellranger      10.1.0 <tarball> <sifdir>
                         containers/build-sif.sh cellranger-arc  2.2.0  <tarball> <sifdir>
                         containers/build-sif.sh cellranger-atac 2.2.0  <tarball> <sifdir>

      local:
        * --programs   : Where the installs live; default </programs>
                         -> /programs/cellranger-<crversion>/cellranger
                         -> /programs/cellranger-arc-<arcversion>/bin/cellranger-arc
                         -> /programs/cellranger-atac-<atacversion>/cellranger-atac
        * --crpath     : Full path to the cellranger binary. Implies --engine local.
        * --arcpath    : Full path to the cellranger-arc binary. Implies --engine local.
        * --atacpath   : Full path to the cellranger-atac binary. Implies --engine local.
        * --listPrograms : List what is installed under --programs

      Either way the resolved image or binary is checked before anything runs, so
      a version that is not there fails immediately and lists what is.

  Runtime / resources:
    * --localcores     : cellranger --localcores, also the cpus reserved; default <32>
    * --localmem       : cellranger --localmem in GB, also the memory reserved; default <128>
    * --maxforks       : concurrent tasks PER PROCESS, so 2 would run 2 cellranger
                         jobs at once, each taking localcores / localmem; default <1>

  The command this builds, gex mode, one per library:

    <crbin> count --id=<label> \\
      --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam> --r1-length=<r1length> \\
      --transcriptome=<ref> \\
      --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]>

  arc mode, one per label:

    <arcbin> count --id=<label> \\
      --reference=<ref> --libraries=<label>_libraries.csv \\
      --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam>

  atac mode, one per library:

    <atacbin> count --id=<label> \\
      --localcores=<localcores> --localmem=<localmem> \\
      --reference=<ref> \\
      --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]> \\
      [ --chemistry=<chemistry> --force-cells=<forceCells> --peaks=<peaks> --dim-reduce=<dimReduce> ]

"""

    exit 0
}


// Each mode drives a different tool, so the banner reports the one that will
// actually run rather than always naming cellranger.
modeTool    = params.mode == "arc" ? "cellranger-arc" : params.mode == "atac" ? "cellranger-atac" : "cellranger"
modeVersion = params.mode == "arc" ? params.arcversion : params.mode == "atac" ? params.atacversion : params.crversion
modeEngine  = params.mode == "arc" ? params.arcengine  : params.mode == "atac" ? params.atacengine  : params.crengine
modeImage   = params.mode == "arc" ? ( params.arcimage  ?: params.arcbin  )
            : params.mode == "atac" ? ( params.atacimage ?: params.atacbin )
            :                         ( params.crimage   ?: params.crbin  )

log.info """
1  0  x      W  O  R  K  F  L  O  W  -  @bixBeta
=========================================================================================================================
brcID        : ${params.id}
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
version      : ${modeTool} ${modeVersion}
engine       : ${modeEngine}
cellranger   : ${modeImage}
"""

// Fail here, with the resolved path, rather than a task dying much later.
// singularity: the .sif must exist. local: the binary must exist.
// Skipped for -stub-run, and for a remote image uri that is not a local file.
def checkEngine(engine, bin, image, what) {

    if( workflow.stubRun ) return

    if( engine == "singularity" ) {

        def local = image.replaceFirst(/^file:\/\//, '')

        // docker:// and friends are resolved by singularity, not by us
        if( local ==~ /^[a-z]+:\/\/.*/ ) return

        if( file(local).exists() ) return

        def built = file("${params.sifdir}/*.sif")
        def avail = ( built instanceof List ? built : ( built ? [built] : [] ) )
                        .collect { it.name }.sort()

        error """No ${what} image at: ${image}
    Present in ${params.sifdir}: ${ avail ? avail.join(', ') : 'no .sif files' }
    Build one with containers/build-sif.sh, point at another dir with --sifdir,
    or run the local install instead with --engine local."""
    }

    // a bare name is resolved on PATH, not a path we can check
    if( !bin.contains('/') ) return

    if( file(bin).exists() ) return

    def installed = file("${params.programs}/cellranger*")
    def avail = ( installed instanceof List ? installed : ( installed ? [installed] : [] ) )
                    .collect { it.name }.sort()

    error """No ${what} found at: ${bin}
    Installed under ${params.programs}: ${ avail ? avail.join(', ') : 'nothing matching cellranger*' }
    Set a version that is installed, or give the full path."""
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

// 10x reference MAP
refsBase = "/local/workdir/10x_analysis/REFS"

refDir = [
Apoculata_stonyCoral    :"${refsBase}/Apoculata_stonyCoral/240514_fromSarahArnold/apoculata",
Canine                  :"${refsBase}/Canine/CanFam3_1",
Combo_Human_Mouse       :"${refsBase}/Combo_Human_Mouse/refdata-gex-GRCh38_and_GRCm39-2024-A",
Combo_Human_Mouse_2020  :"${refsBase}/Combo_Human_Mouse/refdata-gex-GRCh38-and-mm10-2020-A",
Feline                  :"${refsBase}/Feline/Fca126",
Horse                   :"${refsBase}/Horse/EquCab3/ENSEMBL_annot116/spaceranger/EquCab3_ENS116",
Horse_ENS112            :"${refsBase}/Horse/EquCab3/ENSEMBL_annot112/spaceranger/EquCab3_ENS112",
Human                   :"${refsBase}/Human/refdata-gex-GRCh38-2024-A",
MAIZE                   :"${refsBase}/MAIZE/MAIZE_CellRanger",
Mouse                   :"${refsBase}/Mouse/refdata-gex-GRCm39-2024-A",
Nematostella            :"${refsBase}/Nematostella/Nematostella2_2/Nematostella2_2_standard"]

// ARC references carry a regions/ dir alongside fasta / genes / star
refDirArc = [
Canine                  :"${refsBase}/Canine/arc/CanFam3_Ensembl101annot",
Human                   :"${refsBase}/Human/refdata-cellranger-arc-GRCh38-2024-A",
Mouse                   :"${refsBase}/Mouse/refdata-cellranger-arc-GRCm39-2024-A",
Nematostella            :"${refsBase}/Nematostella/Nematostella2_2/Nematostella2_2",
Nematostella_jaNemVect1 :"${refsBase}/Nematostella/jaNemVect1.1/NCBI/jaNemVect1_1"]

// cellranger-atac --reference takes either an ATAC or an ARC reference, so atac
// mode reuses the arc set rather than keeping a duplicate list of the same dirs.
refDirAtac = refDirArc

// Not wired to a mode: VDJ needs cellranger multi / vdj, which this pipeline
// does not run yet.
//   ${refsBase}/Mouse/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0


// A cellranger reference is a directory holding reference.json. A key that
// points at a parent - a species dir holding several builds - would only fail
// much later inside cellranger, so say so here and name the builds found.
def checkReference(ref, what) {

    if( workflow.stubRun ) return

    def d = file(ref)

    if( !d.exists() )
        error "No ${what} reference at: ${ref}\n    See --listRefs, or give a full path with --ref."

    if( file("${ref}/reference.json").exists() ) return

    def kids = file("${ref}/*/reference.json")
    def builds = ( kids instanceof List ? kids : ( kids ? [kids] : [] ) )
                    .collect { it.parent.name }.sort()

    if( builds )
        error """${ref} is not a reference itself, it holds ${builds.size()}:
    ${builds.join(', ')}
    Point --ref at one of them, e.g. --ref ${ref}/${builds[0]}"""

    log.warn "[ref] no reference.json under ${ref} - passing it to ${what} as given"
}


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
include {   CELLRANGER_ATAC_COUNT    } from './modules/cellrangeratac'
include {   WEBSUMMARY_MQC           } from './modules/websummary'
include {   CELLRANGER_AGGR          } from './modules/cellrangeraggr'
include {   COMBINE_METRICS          } from './modules/metrics'
include {   MULTIQC                  } from './modules/multiqc'
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

// R1s for one prefix, or for all of them with prefix "*". A path may be the
// Sample_* dir itself, or the Project_* dir above it that holds Sample_*/
// subdirs - both are valid --fastqs values for cellranger, so every lookup
// looks one level down as well.
def readsUnder(dir, prefix) {

    def hits = []
    [ "${dir}/${prefix}_S*_L*_R1*.f*q.gz", "${dir}/*/${prefix}_S*_L*_R1*.f*q.gz" ].each { pattern ->
        def found = file(pattern)
        if( found instanceof List ) hits.addAll(found)
        else if( found )            hits.add(found)
    }

    hits
}

// cellranger --sample must match the fastq prefix INSIDE the delivery dir, which
// is the dir name without the Sample_ prefix and the trailing lane suffix:
//
//   Sample_SC2620_JS4_G3_***_104***_253***_L2/
//     SC2620_JS4_G3_***_104***_253***_S3_L002_R1_001.fastq.gz
//     -> SC2620_JS4_G3_***_104***_253***
//
// It is read off the files rather than parsed out of the dir name, so an
// unexpected naming variant fails loudly here instead of inside cellranger.
def sampleNames(fq) {

    def names = readsUnder(fq, "*").collect { it.name }
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

    // No column at all means the only thing the mode counts, so an atac sheet
    // does not have to repeat Chromatin Accessibility on every row.
    if( !raw ) return params.mode == "atac" ? "Chromatin Accessibility" : "Gene Expression"

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

// What a sample column probably meant, when what it says cannot be found: the
// illumina tail trimmed back off, whether a whole file name was pasted or just
// the _S1_L001 that ends the directory listing.
def samplePrefixHint(given) {

    def stripped = given.replaceFirst(/\.f(ast)?q(\.gz)?$/, '')
                        .replaceFirst(/_R\d(_001)?$/, '')
                        .replaceFirst(/_L\d+$/, '')
                        .replaceFirst(/_S\d+$/, '')

    stripped != given && stripped ? stripped : null
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
        if( readsUnder(dir, given) )
            return [ dir, [ given ] ]

        // The commonest way to get this column wrong is to paste the file name
        // stem, tail and all. Say so, rather than only reporting nothing found.
        def hint = samplePrefixHint(given)
        if( hint && readsUnder(dir, hint) ) {
            def drop = sampleNames(dir).size() == 1
                     ? " Or drop the sample column: ${dir} holds this one prefix only, so the pipeline reads it off the files."
                     : ""
            error "sample-sheet: ${label}: sample '${given}' is a fastq file name, not a 10x sample prefix. Use '${hint}', the part before _S<n>_L<n>_R<n>_001.${drop}"
        }

        error "sample-sheet: ${label}: no reads for sample '${given}' under ${dir}. Expected ${given}_S<n>_L<n>_R1_001.fastq.gz."
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

    checkEngine(params.crengine, params.crbin, params.crimage, "cellranger")
    checkReference(ref, "cellranger")

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

    // one table for the whole run, header once
    COMBINE_METRICS(
        CELLRANGER_COUNT.out.run_metrics
            | collect
            | filter { it }
            | map { csvs -> [ "ALL_metrics_summary.csv", "", csvs ] }
    )

    // A label with several libraries means several GEM wells, counted
    // separately. Combining those is what cellranger aggr is for. A label with
    // one library needs nothing: re-sequencing runs were already pooled into
    // that single count.
    if( params.aggr ) {

        aggr_ch = CELLRANGER_COUNT.out.molecule_info
            | groupTuple
            | filter { label, libraries, h5s -> libraries.size() > 1 }
            | map { label, libraries, h5s ->
                  // pair them up and sort by library, so the csv rows and the
                  // staged files stay in step and the order is deterministic
                  def pairs = [ libraries, h5s ].transpose().sort { it[0] }
                  [ label, pairs.collect { it[0] }, pairs.collect { it[1] } ]
              }
            | view { label, libraries, h5s ->
                  "AGGR >> ${label}  ( ${libraries.size()} libraries: ${libraries.join(', ')} )" }

        CELLRANGER_AGGR(aggr_ch)

        ch_versions = ch_versions.mix( CELLRANGER_AGGR.out.versions.first() )
    }

    // DUMP_VERSIONS first: it writes the Software Versions section MultiQC
    // renders, so the table is an input to the report rather than an output of
    // it. MultiQC's own version is therefore not in the table - it cannot be,
    // since the table has to exist before multiqc runs.
    DUMP_VERSIONS(ch_versions.collect())

    // MultiQC reads the web summaries, the aggregated one included
    if( params.multiqc ) {

        checkEngine(params.multiqcengine, params.multiqcbin, params.multiqcimage, "multiqc")

        mqc_ch = params.aggr
            ? CELLRANGER_COUNT.out.run_web_summary.mix( CELLRANGER_AGGR.out.run_web_summary )
            : CELLRANGER_COUNT.out.run_web_summary

        MULTIQC( mqc_ch.collect(), ch_mqc_config, DUMP_VERSIONS.out.mqc_yml )
    }
}


/* ---------------------------------------------------------------------------------------------------------
MULTIOME Workflow

One label carries a Gene Expression library and a Chromatin Accessibility
library. They are NOT counted separately: cellranger-arc takes a libraries.csv
listing both, which is generated per label from the sheet.
------------------------------------------------------------------------------------------------------------ */

workflow ARC {

    checkEngine(params.arcengine, params.arcbin, params.arcimage, "cellranger-arc")
    checkReference(ref, "cellranger-arc")

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

    // arc reports summary.csv, whose columns differ from the gex table, so it
    // gets its own combined file rather than being stacked onto that one
    COMBINE_METRICS(
        CELLRANGER_ARC_COUNT.out.run_metrics
            | collect
            | filter { it }
            | map { csvs -> [ "ALL_summary.csv", "cellranger_arc_metrics", csvs ] }
    )

    DUMP_VERSIONS( CELLRANGER_ARC_COUNT.out.versions.first().collect() )

    if( params.multiqc ) {

        checkEngine(params.multiqcengine, params.multiqcbin, params.multiqcimage, "multiqc")

        MULTIQC(
            CELLRANGER_ARC_COUNT.out.run_web_summary.mix( COMBINE_METRICS.out.mqc ).collect(),
            ch_mqc_config,
            DUMP_VERSIONS.out.mqc_yml
        )
    }
}


/* ---------------------------------------------------------------------------------------------------------
scATAC Workflow

cellranger-atac counts one library at a time from --fastqs / --sample, like the
gex path and unlike arc. Nothing here is specific to a standalone ATAC kit: the
ATAC library of a multiome run is counted the same way, with --chemistry ARC-v1.
------------------------------------------------------------------------------------------------------------ */

workflow ATAC {

    checkEngine(params.atacengine, params.atacbin, params.atacimage, "cellranger-atac")
    checkReference(ref, "cellranger-atac")

    if( ref == null ){
        error "No reference provided. Use --ref < key or path >, see --listRefs"
    }

    meta_ch = readRows()
        | map { label, library, dir, samples, ltype ->
              if( ltype != "Chromatin Accessibility" )
                  error "sample-sheet: ${label} is '${ltype}'. --mode atac counts Chromatin Accessibility only; use --mode arc for multiome."
              [ [label, library], [ dir, samples ] ]
          }
        | groupTuple
        | map { key, entries ->
              // same pooling rule as gex: repeated rows for one library are
              // re-sequencing runs, passed as comma separated lists
              def dirs    = entries.collect { it[0] }.unique().sort()
              def samples = entries.collectMany { it[1] }.unique().sort()

              [ key[1], key[0], dirs, samples ]
          }
        | view { library, label, dirs, samples ->
              "ATAC >> ${library}  ( label: ${label}, fastq dirs: ${dirs.size()}, sample: ${samples.join(',')} )" }

    CELLRANGER_ATAC_COUNT(meta_ch, ref_ch)

    // cellranger-atac reports summary.csv, with its own columns again, so it
    // gets its own combined table. cellranger-atac aggr is not wired up yet, so
    // several libraries under one label stay as separate counts.
    COMBINE_METRICS(
        CELLRANGER_ATAC_COUNT.out.run_metrics
            | collect
            | filter { it }
            | map { csvs -> [ "ALL_summary.csv", "cellranger_atac_metrics", csvs ] }
    )

    // Everything the web summary shows - the grouped metric tables, the knee,
    // TSS and insert size curves - parsed out of its embedded json, since no
    // multiqc module will do it. summary.csv only carries the flat metrics.
    WEBSUMMARY_MQC(
        CELLRANGER_ATAC_COUNT.out.run_web_summary
            | collect
            | map { htmls -> [ "cellranger_atac_metrics", htmls ] },
        channel.value(params.mqcpoints)
    )

    DUMP_VERSIONS( CELLRANGER_ATAC_COUNT.out.versions.first().collect() )

    if( params.multiqc ) {

        checkEngine(params.multiqcengine, params.multiqcbin, params.multiqcimage, "multiqc")

        // The html itself is published for humans but NOT handed to multiqc:
        // nothing would parse it, and the cellranger_arc module trying and
        // failing would take the whole report down ( MultiQC issue #3609 ).
        // What multiqc gets is what this pipeline rendered from it.
        MULTIQC(
            COMBINE_METRICS.out.mqc.mix( WEBSUMMARY_MQC.out.mqc ).collect(),
            ch_mqc_config,
            DUMP_VERSIONS.out.mqc_yml
        )
    }
}


workflow {

    if( params.mode == "gex" ){

        GEX()
    }

    else if( params.mode == "atac" ){

        ATAC()
    }

    else if( params.mode == "arc" ){

        ARC()
    }

    else {

        error "Invalid mode provided: ${params.mode}"
    }

}
