chemistry       = params.chemistry
expectCells     = params.expectCells
forceCells      = params.forceCells
createBam       = params.createBam
introns         = params.introns
r1length        = params.r1length
r2length        = params.r2length

// resolved in nextflow.config from --crpath, or --crversion + --programs
crBin           = params.crbin

// Cell Ranger 8 dropped --no-bam in favour of a mandatory --create-bam
crMajor         = params.crversion.tokenize('.')[0] as Integer


process CELLRANGER_COUNT {

    tag "$library"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER/${library}" , mode: "symlink", overwrite: true , pattern: "outs/**"

    input:
        tuple val(library), val(label), val(fastqs), val(samples)
        val ref

    output:
        tuple val(library), path("outs")                        , emit: outs
        path "outs/web_summary.html"                            , emit: web_summary
        path "outs/metrics_summary.csv"                         , emit: metrics
        path "outs/*filtered_feature_bc_matrix.h5"              , emit: filtered_h5   , optional: true
        path "outs/*.bam*"                                      , emit: bam           , optional: true
        path "versions.yml"                                     , emit: versions

    script:

    // fastqs are 10x delivery dirs, passed through untouched - no staging, no
    // renaming. Several runs of the same library become the comma separated
    // lists cellranger expects, for both the dirs and the sample prefixes.
    def fastqArg    = fastqs instanceof List ? fastqs.join(',') : "${fastqs}"
    def sampleArg   = samples instanceof List ? samples.join(',') : "${samples}"

    def bamArg      = crMajor >= 8 ? "--create-bam=${createBam}" : ( createBam ? "" : "--no-bam" )
    def r1Arg       = r1length    ? "--r1-length=${r1length}"        : ""
    def r2Arg       = r2length    ? "--r2-length=${r2length}"        : ""
    def chemArg     = chemistry   ? "--chemistry=${chemistry}"       : ""
    def cellsArg    = expectCells ? "--expect-cells=${expectCells}"  : ""
    def forceArg    = forceCells  ? "--force-cells=${forceCells}"    : ""
    def intronArg   = introns != null ? "--include-introns=${introns}" : ""

    """
    ${crBin} count \
        --id=${library} \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        ${bamArg} ${r1Arg} ${r2Arg} \
        --transcriptome=${ref} \
        --fastqs=${fastqArg} \
        --sample=${sampleArg} \
        ${chemArg} ${cellsArg} ${forceArg} ${intronArg}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${library}/outs outs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$( ${crBin} --version 2>&1 | head -1 | sed 's/.*cellranger[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p outs
    touch outs/web_summary.html
    touch outs/metrics_summary.csv
    touch outs/filtered_feature_bc_matrix.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: ${params.crversion}
    END_VERSIONS
    """
}
