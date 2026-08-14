createBam       = params.createBam

// binary: a native install via --arcpath, otherwise whatever is on PATH in the container
arcBin          = params.arcpath ?: "cellranger-arc"


process CELLRANGER_ARC_COUNT {

    tag "$label"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER_ARC/${label}" , mode: "symlink", overwrite: true , pattern: "outs/**"
    publishDir "CELLRANGER_ARC/${label}" , mode: "copy"   , overwrite: true , pattern: "*_libraries.csv"

    input:
        tuple val(label), val(rows)
        val ref

    output:
        tuple val(label), path("outs")                          , emit: outs
        path "${label}_libraries.csv"                           , emit: libraries
        path "outs/web_summary.html"                            , emit: web_summary
        path "outs/summary.csv"                                 , emit: metrics       , optional: true
        path "outs/*filtered_feature_bc_matrix.h5"              , emit: filtered_h5   , optional: true
        path "versions.yml"                                     , emit: versions

    script:

    // cellranger-arc takes a libraries.csv rather than --fastqs / --sample:
    //
    //   fastqs,sample,library_type
    //   /local/Illumina/.../Project_10488522,SC2619_JS4_BC_...,Gene Expression
    //   /local/Illumina/.../Project_10488522,SC2619_JS4_MA_...,Chromatin Accessibility
    //
    // It is written here rather than by hand, straight from the sample sheet.
    def csvRows = rows.collect { "echo '${it[0]},${it[1]},${it[2]}' >> ${label}_libraries.csv" }.join('\n    ')

    """
    echo 'fastqs,sample,library_type' > ${label}_libraries.csv
    ${csvRows}

    echo "--- ${label}_libraries.csv ---"
    cat ${label}_libraries.csv

    ${arcBin} count \
        --id=${label} \
        --reference=${ref} \
        --libraries=${label}_libraries.csv \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()} \
        --create-bam=${createBam}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${label}/outs outs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: \$( ${arcBin} --version 2>&1 | head -1 | sed 's/.*cellranger-arc[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    def csvRows = rows.collect { "echo '${it[0]},${it[1]},${it[2]}' >> ${label}_libraries.csv" }.join('\n    ')
    """
    echo 'fastqs,sample,library_type' > ${label}_libraries.csv
    ${csvRows}

    mkdir -p outs
    touch outs/web_summary.html
    touch outs/summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: ${params.arcversion}
    END_VERSIONS
    """
}
