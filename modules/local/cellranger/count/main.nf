/*
========================================================================================
    CELLRANGER_COUNT
========================================================================================
    One run per physical 10x library (one GEM well).

    FASTQ staging
    -------------
    Cell Ranger only accepts bcl2fastq-style filenames, so the samplesheet's
    arbitrarily-named files are staged as R1_01/R2_01/R1_02/... by Nextflow and
    then symlinked into

        cr_fastqs/<library>_S1_L001_R1_001.fastq.gz
        cr_fastqs/<library>_S1_L001_R2_001.fastq.gz
        cr_fastqs/<library>_S1_L002_R1_001.fastq.gz   <- 2nd sequencing run
        ...

    Each input row becomes one pseudo-lane, which is how re-sequencing runs of
    the same library get pooled into a single count.

    Version pinning
    ---------------
    The image tag follows --cellranger_version, so a run can be pinned to any
    built version without touching the code:

        nextflow run . --cellranger_version 7.2.0
        nextflow run . --cellranger_container /shared/sif/cellranger-9.0.1.sif
----------------------------------------------------------------------------------------
*/

process CELLRANGER_COUNT {

    tag   "${meta.id}"
    label 'process_cellranger'

    container { params.cellranger_container ?: "ghcr.io/${params.ghcr_owner}/cellranger:${params.cellranger_version}" }

    input:
    tuple val(meta), path(r1, stageAs: 'R1_??.fastq.gz'), path(r2, stageAs: 'R2_??.fastq.gz'), path(r3, stageAs: 'R3_??.fastq.gz')
    path  reference

    output:
    tuple val(meta), path("${meta.id}/outs")                                , emit: outs
    tuple val(meta), path("${meta.id}/outs/*filtered_feature_bc_matrix.h5") , emit: filtered_h5 , optional: true
    tuple val(meta), path("${meta.id}/outs/*raw_feature_bc_matrix.h5")      , emit: raw_h5      , optional: true
    tuple val(meta), path("${meta.id}/outs/metrics_summary.csv")            , emit: metrics     , optional: true
    tuple val(meta), path("${meta.id}/outs/web_summary.html")               , emit: web_summary , optional: true
    tuple val(meta), path("${meta.id}/outs/*.bam*")                         , emit: bam         , optional: true
    path  "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: meta.id
    def mem_gb   = Math.max( (task.memory.toGiga() as int) - 4, 4 )

    // Cell Ranger 8 dropped --no-bam in favour of a mandatory --create-bam
    def cr_major = (params.cellranger_version =~ /^(\d+)/) ? ((params.cellranger_version =~ /^(\d+)/)[0][1] as int) : 9
    def bam_arg  = cr_major >= 8
                    ? "--create-bam=${params.create_bam}"
                    : ( params.create_bam ? '' : '--no-bam' )

    def chem_arg    = params.chemistry       ? "--chemistry=${params.chemistry}"             : ''
    def cells_arg   = params.expect_cells    ? "--expect-cells=${params.expect_cells}"       : ''
    def force_arg   = params.force_cells     ? "--force-cells=${params.force_cells}"         : ''
    def intron_arg  = params.include_introns != null ? "--include-introns=${params.include_introns}" : ''

    """
    set -euo pipefail
    shopt -s nullglob

    # ---- stage FASTQs under Cell Ranger's required naming convention --------
    mkdir -p cr_fastqs
    lane=0
    for f in R1_*.gz ; do
        lane=\$(( lane + 1 ))
        lane_id=\$( printf "L%03d" \$lane )
        suffix=\${f#R1_}

        ln -s ../"\$f" "cr_fastqs/${prefix}_S1_\${lane_id}_R1_001.fastq.gz"

        if [ -e "R2_\${suffix}" ] ; then
            ln -s ../"R2_\${suffix}" "cr_fastqs/${prefix}_S1_\${lane_id}_R2_001.fastq.gz"
        else
            echo "ERROR: no R2 mate found for \$f" >&2 ; exit 1
        fi

        if [ -e "R3_\${suffix}" ] ; then
            ln -s ../"R3_\${suffix}" "cr_fastqs/${prefix}_S1_\${lane_id}_R3_001.fastq.gz"
        fi
    done

    if [ "\$lane" -eq 0 ] ; then
        echo "ERROR: no FASTQ files were staged for ${prefix}" >&2 ; exit 1
    fi
    echo "Staged \$lane sequencing run(s) as pseudo-lanes for ${prefix}"

    # ---- quantify ------------------------------------------------------------
    cellranger count \\
        --id="${prefix}" \\
        --transcriptome="${reference}" \\
        --fastqs=cr_fastqs \\
        --sample="${prefix}" \\
        --localcores=${task.cpus} \\
        --localmem=${mem_gb} \\
        ${bam_arg} \\
        ${chem_arg} \\
        ${cells_arg} \\
        ${force_arg} \\
        ${intron_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$( cellranger --version 2>&1 | head -n1 | sed 's/.*cellranger[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p ${prefix}/outs/filtered_feature_bc_matrix
    touch ${prefix}/outs/web_summary.html
    touch ${prefix}/outs/metrics_summary.csv
    touch ${prefix}/outs/filtered_feature_bc_matrix.h5
    touch ${prefix}/outs/raw_feature_bc_matrix.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: ${params.cellranger_version}
    END_VERSIONS
    """
}
