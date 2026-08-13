chemistry       = params.chemistry
expectCells     = params.expectCells
forceCells      = params.forceCells
createBam       = params.createBam
introns          = params.introns
r1length        = params.r1length
r2length        = params.r2length

// binary: a native install via --crpath, otherwise whatever is on PATH in the container
crBin           = params.crpath ?: "cellranger"

// Cell Ranger 8 dropped --no-bam in favour of a mandatory --create-bam
crMajor         = params.crversion.tokenize('.')[0] as Integer


process CELLRANGER_COUNT {

    tag "$library"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER/${library}" , mode: "symlink", overwrite: true , pattern: "outs/**"

    input:
        tuple val(library), val(label), path(r1, stageAs: "R1_??.fastq.gz"), path(r2, stageAs: "R2_??.fastq.gz")
        val ref

    output:
        tuple val(library), path("outs")                        , emit: outs
        path "outs/web_summary.html"                            , emit: web_summary
        path "outs/metrics_summary.csv"                         , emit: metrics
        path "outs/*filtered_feature_bc_matrix.h5"              , emit: filtered_h5   , optional: true
        path "outs/*.bam*"                                      , emit: bam           , optional: true
        path "versions.yml"                                     , emit: versions

    script:

    def bamArg      = crMajor >= 8 ? "--create-bam=${createBam}" : ( createBam ? "" : "--no-bam" )
    def r1Arg       = r1length    ? "--r1-length=${r1length}"        : ""
    def r2Arg       = r2length    ? "--r2-length=${r2length}"        : ""
    def chemArg     = chemistry   ? "--chemistry=${chemistry}"       : ""
    def cellsArg    = expectCells ? "--expect-cells=${expectCells}"  : ""
    def forceArg    = forceCells  ? "--force-cells=${forceCells}"    : ""
    def intronArg   = introns != null ? "--include-introns=${introns}" : ""

    """
    # Cell Ranger only accepts bcl2fastq style names, so each input row is
    # symlinked in as one lane: <library>_S1_L00N_R1_001.fastq.gz
    mkdir -p fqs
    lane=0
    for f in R1_*.gz ; do
        lane=\$(( lane + 1 ))
        laneid=\$( printf "L%03d" \$lane )
        suffix=\${f#R1_}

        ln -s ../"\$f"            "fqs/${library}_S1_\${laneid}_R1_001.fastq.gz"
        ln -s ../"R2_\${suffix}"  "fqs/${library}_S1_\${laneid}_R2_001.fastq.gz"
    done

    echo "staged \$lane sequencing run(s) for ${library}"

    ${crBin} count \
        --id=${library} \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        ${bamArg} ${r1Arg} ${r2Arg} \
        --transcriptome=${ref} \
        --fastqs=fqs \
        --sample=${library} \
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
