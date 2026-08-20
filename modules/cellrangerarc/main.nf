createBam       = params.createBam
checkversion    = params.checkversion
wantversion     = params.arcversion

// resolved in nextflow.config from --arcpath, or --arcversion + --programs
arcBin          = params.arcbin


process CELLRANGER_ARC_COUNT {

    tag "$label"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER_ARC/${label}" , mode: "symlink", overwrite: true , pattern: "outs"
    publishDir "CELLRANGER_ARC/${label}" , mode: "copy"   , overwrite: true , pattern: "*_libraries.csv"

    // flat collections across labels, for a quick look over a whole run
    publishDir "filtered_counts/${label}" , mode: "symlink", overwrite: true , pattern: "filtered_feature_bc_matrix.h5"
    publishDir "web_summary_htmls"       , mode: "symlink", overwrite: true , pattern: "*_web_summary.html"
    publishDir "summary_metrics"         , mode: "symlink", overwrite: true , pattern: "*_summary.csv"

    input:
        tuple val(label), val(rows)
        val ref

    output:
        tuple val(label), path("outs")                          , emit: outs
        path "${label}_libraries.csv"                           , emit: libraries
        path "outs/web_summary.html"                            , emit: web_summary
        path "outs/summary.csv"                                 , emit: metrics       , optional: true
        path "outs/*filtered_feature_bc_matrix.h5"              , emit: filtered_h5   , optional: true
        path "${label}_web_summary.html"                        , emit: run_web_summary , optional: true
        path "${label}_summary.csv"                             , emit: run_metrics     , optional: true
        path "filtered_feature_bc_matrix.h5"                    , emit: filtered_link   , optional: true
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

    # The image or install is only labelled with a version - ask the tool itself
    # before counting anything, so a mislabelled sif cannot silently produce
    # results from the wrong cellranger-arc.
    if [ "${checkversion}" = "true" ] ; then
        have=\$( ${arcBin} --version 2>&1 | head -1 | sed 's/.*cellranger-arc[- ]//; s/[^0-9.].*\$//' )
        if [ -z "\$have" ] ; then
            echo "WARN: could not read a version from cellranger-arc --version, skipping the check" >&2
        elif [ "\$have" != "${wantversion}" ] ; then
            echo "ERROR: asked for cellranger-arc ${wantversion} but ${arcBin} reports \$have" >&2
            echo "       Set --arcversion to match, point at another image with --arcsif, or pass --checkversion false." >&2
            exit 1
        fi
    fi


    ${arcBin} count \\
        --id=${label} \\
        --reference=${ref} \\
        --libraries=${label}_libraries.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --create-bam=${createBam}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${label}/outs outs

    # label prefixed links for the run level collections
    # if/then, not &&: nextflow runs task scripts under bash -ue, so a false
    # test as the last command of a line would abort the task
    if [ -e outs/web_summary.html ] ; then
        ln -s outs/web_summary.html ${label}_web_summary.html
    fi
    if [ -e outs/summary.csv ] ; then
        ln -s outs/summary.csv ${label}_summary.csv
    fi
    if [ -e outs/filtered_feature_bc_matrix.h5 ] ; then
        ln -s outs/filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5
    fi

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

    touch outs/filtered_feature_bc_matrix.h5

    ln -s outs/web_summary.html ${label}_web_summary.html
    ln -s outs/summary.csv      ${label}_summary.csv
    ln -s outs/filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-arc: ${params.arcversion}
    END_VERSIONS
    """
}
