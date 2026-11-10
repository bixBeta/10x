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

// spatial only: spaceranger counts one Visium capture area, so it needs the
// microscope image and the slide identity alongside the reads.
// The read lengths are their own params rather than --r1length / --r2length:
// a Visium library is not a gex one, and gex's default of 28 would silently
// trim a 43 nt CytAssist R1.
params.imagedir         = "visium_images"
params.slide            = null
params.probeset         = null
params.nosecondary      = false
params.spatialR1length  = null
params.spatialR2length  = null

// multi only: cellranger multi needs a reference per assay present in the sheet
params.vdjref           = null
params.featureref       = null

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

        VISIUM ( --mode spatial ): one row per capture area. spaceranger takes
        --fastqs / --sample like the gex path, plus two things only spatial
        has: WHICH slide and area the section sat on, and an IMAGE of the
        tissue to align the spots to.

        |-------|-----------------------------|------------|------|------------------------|
        | label | fastqs                      | slide      | area | cytaimage              |
        |-------|-----------------------------|------------|------|------------------------|
        | YF1   | /local/.../Sample_..._A1_..._L1 | H1-DYCH77W | A1   | CAVG..._A1_yf1.tif |
        | YF2   | /local/.../Sample_..._B1_..._L1 | H1-DYCH77W | B1   | CAVG..._B1_yf2.tif |
        |-------|-----------------------------|------------|------|------------------------|

        slide is the serial printed on the slide and area is the capture area
        the section sits in ( A1 | B1 | C1 | D1 ). One slide normally carries
        a whole run, so --slide sets it for every row and only area has to be
        in the sheet. If the serial is lost, put one of spaceranger's keywords
        in slide instead - visium-1, visium-2, visium-2-large - and leave area
        out; that becomes --unknown-slide.

        An image is REQUIRED, in whichever column matches what was taken:

            cytaimage        the tif the CytAssist instrument wrote
                             ( CytAssist assays, Visium HD included )
            image            a brightfield H&E image
            darkimage        a fluorescence image
            colorizedimage   a composite colour fluorescence image

        A bare file name is looked up under --imagedir, which defaults to
        visium_images next to where you launch - so drop the tifs the lab sent
        in there and the sheet only has to name them. A full path works too.

        Three more optional columns, all files, all resolved the same way:

            probeset     --probe-set, for a probe based ( FFPE / CytAssist )
                         assay that needs a panel other than the reference's
            alignment    --loupe-alignment, the json Loupe Browser exports
                         when the fiducial frame had to be aligned by hand.
                         Needed whenever spaceranger's own alignment missed,
                         which is common for fluorescence and folded sections
            slidefile    --slidefile, the slide layout .gpr. spaceranger
                         normally fetches the layout from 10x, so this is what
                         an offline compute node needs

        A capture area sequenced over several lanes or runs: repeat the label,
        exactly as in gex mode. The slide, area and image columns must match on
        every row of that label - they describe one physical section.

        MULTI ( --mode multi ): one label carries every library of one sample,
        told apart by library_type. cellranger multi takes them together
        through a config csv, which this pipeline writes for you.

        |-------|----------------------------|-------|------------------|
        | label | fastqs                     | lanes | library_type     |
        |-------|----------------------------|-------|------------------|
        | EW1   | /local/.../Project_104***  | 2     | Gene Expression  |
        | EW1   | /local/.../Project_104***  | 2     | Antibody Capture |
        | EW1   | /local/.../Project_104***  | 5     | VDJ              |
        |-------|----------------------------|-------|------------------|

        becomes EW1_multi_config.csv:

            [gene-expression]
            reference,<--ref>
            create-bam,<--createBam>

            [vdj]
            reference,<--vdjref>

            [feature]
            reference,<--featureref>

            [libraries]
            fastq_id,fastqs,lanes,feature_types
            SC2622_EW1_G5_***,/local/.../Project_104***,2,Gene Expression
            SC2622_EW1_CSP_***,/local/.../Project_104***,2,Antibody Capture
            SC2622_EW1_VDJ_***,/local/.../Project_104***,5,VDJ

        lanes is OPTIONAL and multi only. It restricts the library to those
        flow cell lanes; leave it out and 'any' is used. Several lanes are
        separated with a PIPE, since the comma is the csv delimiter: 1|2|3.
        A range works too: 1-4.

        The [vdj] and [feature] sections are only written when a library needs
        them, since an empty section is an error rather than a no-op.

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

    * --mode           : use 'gex'   for 3'/5' gene expression; default <gex>
                       : use 'arc'   for multiome GEX + ATAC
                       : use 'atac'  for scATAC ( also the ATAC half of a multiome
                         run on its own, with --chemistry ARC-v1 )
                       : use 'multi' for cellranger multi: GEX plus any of
                         Antibody Capture ( CITE-seq ), CRISPR and VDJ
                       : use 'spatial' for Visium / Visium HD through spaceranger

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
    * --imagedir       : spatial only. Where a bare image file name in the sheet
                         is looked up; default <visium_images>, relative to the
                         launch dir. A full path in the sheet skips this.
    * --slide          : spatial only. Slide serial for every row that does not
                         give its own, e.g. H1-DYCH77W. One slide normally
                         carries a whole run.
    * --probeset       : spatial only. --probe-set panel csv for every row that
                         does not give its own; default <null>
    * --nosecondary    : spatial only. Skip spaceranger's secondary analysis
                         ( clustering, t-SNE, UMAP ); default <false>
    * --spatialR1length: spatial only. --r1-length; default <null> ( spaceranger
                         reads the length itself ). NOT --r1length: a CytAssist
                         library carries a 43 nt R1, so the gex default of 28
                         would trim every read.
    * --spatialR2length: spatial only. --r2-length; default <null>
    * --vdjref         : multi only. VDJ reference for the [vdj] section; a key
                         from --listRefs or a full path. Required when the sheet
                         carries a VDJ library.
    * --featureref     : multi only. The feature barcode panel csv for the
                         [feature] section. Required when the sheet carries an
                         Antibody Capture or CRISPR library.
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
    * --srversion      : spaceranger version; default <4.0.1>
    * --engine         : 'singularity' runs a locally built .sif; default <singularity>
                       : 'local'       runs an install under --programs

      singularity ( default ):
        * --sifdir     : Where the images live; default </local/workdir/singularity>
                         -> <sifdir>/cellranger-<crversion>.sif
                         -> <sifdir>/cellranger-arc-<arcversion>.sif
                         -> <sifdir>/cellranger-atac-<atacversion>.sif
                         -> <sifdir>/spaceranger-<srversion>.sif
        * --crsif      : Full path to the gex image. Overrides sifdir + crversion.
        * --arcsif     : Full path to the arc image.
        * --atacsif    : Full path to the atac image.
        * --srsif      : Full path to the spaceranger image.

                         Build one with:
                         containers/build-sif.sh cellranger      10.1.0 <tarball> <sifdir>
                         containers/build-sif.sh cellranger-arc  2.2.0  <tarball> <sifdir>
                         containers/build-sif.sh cellranger-atac 2.2.0  <tarball> <sifdir>
                         containers/build-sif.sh spaceranger     4.0.1  <tarball> <sifdir>

      local:
        * --programs   : Where the installs live; default </programs>
                         -> /programs/cellranger-<crversion>/cellranger
                         -> /programs/cellranger-arc-<arcversion>/bin/cellranger-arc
                         -> /programs/cellranger-atac-<atacversion>/cellranger-atac
                         -> /programs/spaceranger-<srversion>/spaceranger
        * --crpath     : Full path to the cellranger binary. Implies --engine local.
        * --arcpath    : Full path to the cellranger-arc binary. Implies --engine local.
        * --atacpath   : Full path to the cellranger-atac binary. Implies --engine local.
        * --srpath     : Full path to the spaceranger binary. Implies --engine local.
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

  spatial mode, one per capture area:

    <srbin> count --id=<label> \\
      --localcores=<localcores> --localmem=<localmem> \\
      --transcriptome=<ref> \\
      --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]> \\
      --create-bam=<createBam> \\
      --slide=<slide> --area=<area>   (  or --unknown-slide=<visium-N>  ) \\
      --cytaimage=<tif> [ --image / --darkimage / --colorizedimage ] \\
      [ --probe-set=<csv> --loupe-alignment=<json> --slidefile=<gpr> \\
        --nosecondary --r1-length=<spatialR1length> --r2-length=<spatialR2length> ]

"""

    exit 0
}


// Each mode drives a different tool, so the banner reports the one that will
// actually run rather than always naming cellranger.
modeTool    = params.mode == "arc"     ? "cellranger-arc"
            : params.mode == "atac"    ? "cellranger-atac"
            : params.mode == "spatial" ? "spaceranger"
            :                            "cellranger"
modeVersion = params.mode == "arc"     ? params.arcversion
            : params.mode == "atac"    ? params.atacversion
            : params.mode == "spatial" ? params.srversion
            :                            params.crversion
modeEngine  = params.mode == "arc"     ? params.arcengine
            : params.mode == "atac"    ? params.atacengine
            : params.mode == "spatial" ? params.srengine
            :                            params.crengine
modeImage   = params.mode == "arc"     ? ( params.arcimage  ?: params.arcbin  )
            : params.mode == "atac"    ? ( params.atacimage ?: params.atacbin )
            : params.mode == "spatial" ? ( params.srimage    ?: params.srbin   )
            :                            ( params.crimage   ?: params.crbin  )

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

if( params.mode == "spatial" ) {
    log.info """imagedir     : ${params.imagedir}
slide        : ${params.slide ?: '< from the sheet >'}
probeset     : ${params.probeset ?: '< none >'}
nosecondary  : ${params.nosecondary}
r1length     : ${params.spatialR1length ?: '< left to spaceranger >'}
r2length     : ${params.spatialR2length ?: '< left to spaceranger >'}
"""
}

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

    def installed = file("${params.programs}/{cellranger,spaceranger}*")
    def avail = ( installed instanceof List ? installed : ( installed ? [installed] : [] ) )
                    .collect { it.name }.sort()

    error """No ${what} found at: ${bin}
    Installed under ${params.programs}: ${ avail ? avail.join(', ') : 'nothing matching cellranger* or spaceranger*' }
    Set a version that is installed, or give the full path."""
}


if( params.listPrograms ) {

    println("")
    log.info """
    10x software installed under ${params.programs}
    =========================================================================================================================
    """
    .stripIndent()

    def installed = file("${params.programs}/{cellranger,spaceranger}*")
    def dirs = ( installed instanceof List ? installed : ( installed ? [installed] : [] ) ).sort()

    if( !dirs ) {
        println "nothing matching ${params.programs}/cellranger* or ${params.programs}/spaceranger*"
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

// VDJ references, for --mode multi. These are a different build from the gex
// ones, so --vdjref is its own key rather than being derived from --ref.
refDirVdj = [
Human                   :"${refsBase}/Human/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0",
Mouse                   :"${refsBase}/Mouse/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0"]


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

    log.info """
    Available 10x VDJ References   ( --vdjref, for --mode multi )
    =========================================================================================================================
    """
    .stripIndent()

    refDirVdj.each(printMap)

    exit 0
}

include {   CELLRANGER_COUNT         } from './modules/cellranger'
include {   CELLRANGER_ARC_COUNT     } from './modules/cellrangerarc'
include {   CELLRANGER_ATAC_COUNT    } from './modules/cellrangeratac'
include {   SPACERANGER_COUNT        } from './modules/spaceranger'
include {   WEBSUMMARY_MQC           } from './modules/websummary'
include {   CELLRANGER_MULTI         } from './modules/cellrangermulti'
include {   MULTI_METRICS_MQC        } from './modules/multimetrics'
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

// --vdjref takes a key from refDirVdj or a full path, like --ref
vdjRef = refDirVdj.containsKey(params.vdjref) ? refDirVdj[params.vdjref] : params.vdjref


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


// A spatial input that is a file: the microscope image, the probe set csv, a
// Loupe manual alignment json, a slide layout gpr. A bare name is looked up
// under --imagedir, so the sheet can name the tif the lab delivered rather than
// repeating a long path on every row; an absolute path is taken as given.
// Absent returns [], which is what an optional path input in a module takes.
def spatialFile(given, label, what) {

    if( !given ) return []

    // file() resolves a relative path against the launch dir and returns an
    // absolute Path either way, so a sheet value and --imagedir are handled by
    // the same call rather than by string surgery.
    def named = file(given)
    def base  = file(params.imagedir)
    def under = file("${base}/${given}")

    if( under.exists() )   return under
    if( named.exists() )   return named
    // -stub-run has no real inputs, so the name is passed through unchecked
    if( workflow.stubRun ) return under

    error """sample-sheet: ${label}: no ${what} named '${given}'
    looked for  ${under}
    and for     ${named}
    Put the images under ${base} and name them in the sheet, or give a full path."""
}


// spaceranger counts one capture area at a time from --fastqs / --sample, like
// the gex path. What is different is that every area also carries its slide
// identity and a picture of the tissue, so those columns ride along with the
// reads and the whole sheet is validated before anything is submitted - an
// image typo should not surface an hour into a run.
def readRowsSpatial() {

    ch_sheet
        | splitCsv( header:true )
        | map { row ->
              def label   = row.label?.trim()
              def library = row.library?.trim() ?: label
              def fq      = row.fastqs?.trim()

              if( !label ) error "sample-sheet: every row needs a label"
              if( !fq )    error "sample-sheet: ${label} is missing fastqs ( a 10x delivery dir, or any fastq in it )"

              def ltype = libraryType(row.library_type?.trim(), label)
              if( ltype != "Gene Expression" )
                  error "sample-sheet: ${label} is '${ltype}'. --mode spatial counts Gene Expression only."

              def res = resolveFastqs(fq, label, row.sample?.trim())

              // one slide usually carries every area of a run, so --slide is
              // taken as a default and only the area has to be per row
              def slide = row.slide?.trim() ?: params.slide
              def area  = row.area?.trim()

              if( !slide )
                  error """sample-sheet: ${label} has no slide. It is the serial printed on the Visium slide, e.g. H1-DYCH77W.
    Put it in a slide column, or pass --slide when the whole run is one slide.
    If the serial is lost, use one of spaceranger's keywords instead -
    visium-1, visium-2 or visium-2-large - and leave area out."""

              if( slide.startsWith("visium-") ) {
                  if( area )
                      error "sample-sheet: ${label}: slide '${slide}' is an unknown-slide keyword, which has no capture area. Drop the area, or give the real slide serial."
              }
              else if( !area ) {
                  error "sample-sheet: ${label} has no area. It is the capture area the section sits in: A1, B1, C1 or D1."
              }

              def images = [ row.cytaimage?.trim(), row.image?.trim(),
                             row.darkimage?.trim(), row.colorizedimage?.trim() ]

              if( !images.any { it } )
                  error """sample-sheet: ${label} has no image. spaceranger aligns the reads to a picture of the tissue, so one is required:
    cytaimage        the tif the CytAssist instrument wrote ( CytAssist assays, Visium HD included )
    image            a brightfield H&E image
    darkimage        a fluorescence image
    colorizedimage   a composite colour fluorescence image
    A bare file name is looked up under ${params.imagedir}; a full path works too."""

              [ label, library, res[0], res[1], slide, area ?: "",
                spatialFile(images[0], label, "cytaimage"),
                spatialFile(images[1], label, "image"),
                spatialFile(images[2], label, "darkimage"),
                spatialFile(images[3], label, "colorizedimage"),
                spatialFile(row.probeset?.trim() ?: params.probeset, label, "probe set csv"),
                spatialFile(row.alignment?.trim(), label, "loupe alignment json"),
                spatialFile(row.slidefile?.trim(), label, "slide file") ]
          }
}


// cellranger multi names every library individually in its config, so unlike
// the gex path nothing is pooled into comma separated lists: one sheet row with
// several prefixes under it becomes several [libraries] rows.
def readRowsMulti() {

    ch_sheet
        | splitCsv( header:true )
        | flatMap { row ->
              def label = row.label?.trim()
              def fq    = row.fastqs?.trim()

              if( !label ) error "sample-sheet: every row needs a label"
              if( !fq )    error "sample-sheet: ${label} is missing fastqs ( a 10x delivery dir, or any fastq in it )"

              def ltype = libraryType(row.library_type?.trim(), label)
              def res   = resolveFastqs(fq, label, row.sample?.trim())
              def lanes = row.lanes?.trim()

              res[1].collect { s -> [ label, [ res[0], s, ltype, lanes ] ] }
          }
        | groupTuple
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

    // The web summary's grouped metric tables AND its plotly figures, lifted
    // out of the embedded json since no multiqc module will do it. The tables
    // are richer than summary.csv - its nucleosome fractions are not in the
    // csv at all - and the figures are handed to plotly as cellranger built
    // them, so they look like the web summary because they are its plots.
    WEBSUMMARY_MQC(
        CELLRANGER_ATAC_COUNT.out.run_web_summary
            | collect
            | map { htmls -> [ "cellranger_atac_metrics", htmls ] }
    )

    DUMP_VERSIONS( CELLRANGER_ATAC_COUNT.out.versions.first().collect() )

    if( params.multiqc ) {

        checkEngine(params.multiqcengine, params.multiqcbin, params.multiqcimage, "multiqc")

        // The html itself is published for humans but NOT handed to multiqc:
        // nothing would parse it, and the cellranger_arc module trying and
        // failing would take the whole report down ( MultiQC issue #3609 ).
        // What multiqc gets is what this pipeline rendered from it.
        MULTIQC(
            COMBINE_METRICS.out.mqc
                .mix( WEBSUMMARY_MQC.out.mqc )
                .mix( WEBSUMMARY_MQC.out.mqc_html )
                .collect(),
            ch_mqc_config,
            DUMP_VERSIONS.out.mqc_yml
        )
    }
}


/* ---------------------------------------------------------------------------------------------------------
VISIUM SPATIAL Workflow

spaceranger counts one capture area at a time from --fastqs / --sample, like the
gex path. Everything specific to spatial is an extra INPUT rather than an extra
step: the slide serial and capture area, and a microscope image to align the
reads to. Those are per area, so they live in the sheet.
------------------------------------------------------------------------------------------------------------ */

workflow SPATIAL {

    checkEngine(params.srengine, params.srbin, params.srimage, "spaceranger")

    if( ref == null ){
        error "No reference provided. Use --ref < key or path >, see --listRefs"
    }

    checkReference(ref, "spaceranger")

    // --r1length is the gex flag, and its default of 28 would silently trim a
    // 43 nt CytAssist R1 down to a gex barcode. Rather than quietly ignore a
    // value someone meant, say where it should go.
    if( params.r1length != 28 || workflow.commandLine =~ /--r1length/ )
        error """spatial: --r1length / --r2length are the gene expression flags and are not used here.
    A Visium read is not a gex read - a CytAssist library carries a 43 nt R1 -
    so this mode has --spatialR1length / --spatialR2length instead, and passes
    neither unless asked, leaving spaceranger to read the length itself."""

    meta_ch = readRowsSpatial()
        | map { r -> [ [ r[0], r[1] ], r[2..-1] ] }
        | groupTuple
        | map { key, entries ->

              // same pooling rule as gex: repeated rows for one library are
              // re-sequencing runs, passed as comma separated lists
              def dirs    = entries.collect { it[0] }.unique().sort()
              def samples = entries.collectMany { it[1] }.unique().sort()

              // a capture area is one physical thing, so everything after the
              // reads has to agree across its rows. Two rows naming different
              // images would otherwise have one of them picked at random.
              def rest = entries.collect { it[2..-1] }.unique()
              if( rest.size() > 1 )
                  error "sample-sheet: ${key[1]} appears more than once with different slide / area / image values. Repeat a label only for re-sequencing runs of the same capture area."

              [ key[1], key[0], dirs, samples ] + rest[0]
          }
        | view { r ->
              "VISIUM >> ${r[0]}  ( label: ${r[1]}, fastq dirs: ${r[2].size()}, sample: ${r[3].join(',')}, slide: ${r[4]}${ r[5] ? ' ' + r[5] : '' } )" }

    SPACERANGER_COUNT(meta_ch, ref_ch)

    // spaceranger reports a wide metrics_summary.csv like gex does, so the rows
    // stack straight into one table. It gets its own section id because the
    // columns are spots and slides rather than cells.
    COMBINE_METRICS(
        SPACERANGER_COUNT.out.run_metrics
            | collect
            | filter { it }
            | map { csvs -> [ "ALL_metrics_summary.csv", "spaceranger_metrics", csvs ] }
    )

    DUMP_VERSIONS( SPACERANGER_COUNT.out.versions.first().collect() )

    if( params.multiqc ) {

        checkEngine(params.multiqcengine, params.multiqcbin, params.multiqcimage, "multiqc")

        // MultiQC ships no spaceranger module, and its cellranger module only
        // matches a count web summary, so the html is published for a human to
        // open but NOT handed over: a module that tried and failed would take
        // the whole report down ( MultiQC issue #3609 ). What multiqc gets is
        // what this pipeline rendered from the metrics.
        MULTIQC(
            COMBINE_METRICS.out.mqc.collect(),
            ch_mqc_config,
            DUMP_VERSIONS.out.mqc_yml
        )
    }
}


/* ---------------------------------------------------------------------------------------------------------
MULTI Workflow

One label carries every library of one sample - Gene Expression plus any of
Antibody Capture, CRISPR and VDJ - and cellranger multi takes them together
through a config csv, which is generated per label from the sheet.
------------------------------------------------------------------------------------------------------------ */

workflow MULTI {

    checkEngine(params.crengine, params.crbin, params.crimage, "cellranger")
    checkReference(ref, "cellranger")

    if( ref == null ){
        error "No reference provided. Use --ref < key or path >, see --listRefs"
    }

    multi_ch = readRowsMulti()
        | map { label, rows ->

              def supported = ["Gene Expression", "Antibody Capture",
                               "CRISPR Guide Capture", "Multiplexing Capture",
                               "VDJ", "VDJ-T", "VDJ-B"]

              def types = rows.collect { it[2] }.unique()

              def unsupported = types.findAll { !( it in supported ) }
              if( unsupported )
                  error "multi: ${label} has library_type(s) cellranger multi cannot take: ${unsupported.join(', ')}"

              if( !types.contains("Gene Expression") )
                  error "multi: ${label} has no Gene Expression library. cellranger multi needs one, see --help"

              // a section with no reference is an error inside cellranger, so
              // the reference is required as soon as a library needs it
              if( types.any { it.startsWith("VDJ") } && !vdjRef )
                  error "multi: ${label} has a VDJ library but no --vdjref was given. See --listRefs."

              if( types.any { it in ["Antibody Capture", "CRISPR Guide Capture"] } && !params.featureref )
                  error "multi: ${label} has a ${types.find{ it in ['Antibody Capture','CRISPR Guide Capture'] }} library but no --featureref was given.\n    It is the panel csv, e.g. --featureref /local/.../TotalSeq_C_..._Antibody_reference.csv"

              // Deterministic order, and the one 10x writes by hand: gene
              // expression first, then the feature barcode assays, VDJ last.
              // Sorting the library_type strings alphabetically instead would
              // put Antibody Capture above Gene Expression, which is valid but
              // does not match a hand written config.
              // A list and indexOf, not a map of ranks: the first rank would
              // be 0, and groovy's elvis treats 0 as false, so `rank[t] ?: 99`
              // silently sorted Gene Expression LAST. The elvis on the second
              // comparison below is correct - there 0 means "equal, fall
              // through to the tiebreak", which is the idiom.
              def assays = [ "Gene Expression", "Antibody Capture",
                             "CRISPR Guide Capture", "Multiplexing Capture",
                             "VDJ", "VDJ-T", "VDJ-B" ]

              def ordered = rows.unique().sort { a1, b1 ->
                  def ra = assays.indexOf(a1[2])
                  def rb = assays.indexOf(b1[2])
                  ( ra < 0 ? 99 : ra ) <=> ( rb < 0 ? 99 : rb ) ?: a1[1] <=> b1[1] }

              [ label, ordered ]
          }
        | view { label, rows ->
              "MULTI >> ${label}  ( libraries: " + rows.collect{ "${it[1]} [${it[2]}]" }.join(' + ') + " )" }

    if( vdjRef )           checkReference(vdjRef, "cellranger vdj")
    if( params.featureref && !workflow.stubRun && !file(params.featureref).exists() )
        error "No feature reference csv at: ${params.featureref}"

    CELLRANGER_MULTI(
        multi_ch,
        ref_ch,
        channel.value( vdjRef ?: "" ),
        channel.value( params.featureref ?: "" )
    )

    // One table for the whole run, header once, each row carrying the sample it
    // came from. The section is empty, so nothing is rendered for multiqc from
    // it: a multi metrics_summary.csv is long form ( one metric per row ) and
    // needs its own layout rather than the wide one the other modes stack.
    COMBINE_METRICS(
        CELLRANGER_MULTI.out.run_metrics
            | collect
            | filter { it }
            | map { csvs -> [ "ALL_metrics_summary.csv", "", csvs ] }
    )

    // The same csvs pivoted into one multiqc table per ( Category, Library
    // Type ): rows are samples, columns are metric names.
    MULTI_METRICS_MQC(
        CELLRANGER_MULTI.out.run_metrics
            | collect
            | filter { it }
            | map { csvs -> [ "cellranger_multi_metrics", csvs ] }
    )

    DUMP_VERSIONS( CELLRANGER_MULTI.out.versions.first().collect() )

    if( params.multiqc ) {

        checkEngine(params.multiqcengine, params.multiqcbin, params.multiqcimage, "multiqc")

        // The web summaries go in as well. Nothing parses a multi summary
        // today - multiqc's cellranger module matches on subcommand "count" -
        // but they cost nothing and start working if that changes.
        MULTIQC(
            CELLRANGER_MULTI.out.run_web_summary.mix( MULTI_METRICS_MQC.out.mqc ).collect(),
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

    else if( params.mode == "spatial" ){

        SPATIAL()
    }

    else if( params.mode == "multi" ){

        MULTI()
    }

    else {

        error "Invalid mode provided: ${params.mode}"
    }

}
