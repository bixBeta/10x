projectId       = params.id
checkversion    = params.checkversion
wantversion     = params.multiqcversion

// resolved in nextflow.config from --multiqcpath, or --multiqcversion + --sifdir
mqcBin          = params.multiqcbin


process MULTIQC {

    tag "$projectId"

    label 'process_low'

    // the report is small and is the thing people keep, so copy rather than link
    publishDir "multiqc" , mode: "copy", overwrite: true

    input:
        path web_summaries
        path mqc_config
        path mqc_versions

    output:
        path "${projectId}_multiqc_report.html"      , emit: report
        path "${projectId}_multiqc_report_data"      , emit: data   , optional: true
        path "versions.yml"                          , emit: versions

    script:

    // MultiQC finds cellranger and cellranger-arc runs by matching content
    // inside web_summary.html, so the labelled summaries are the whole input;
    // metrics_summary.csv is not one of its search patterns.
    //
    // mqc_versions is software_versions_mqc.yml. Anything named *_mqc.yml is
    // picked up as custom content, which is how the Software Versions table
    // gets in: MultiQC's own version detection is off in the config.

    """
    if [ "${checkversion}" = "true" ] ; then
        have=\$( ${mqcBin} --version 2>&1 | head -1 | sed 's/[^0-9]*//; s/[^0-9.].*\$//' )
        if [ -z "\$have" ] ; then
            echo "WARN: could not read a version from multiqc --version, skipping the check" >&2
        elif [ "\$have" != "${wantversion}" ] ; then
            echo "ERROR: asked for multiqc ${wantversion} but ${mqcBin} reports \$have" >&2
            echo "       Set --multiqcversion to match, point at another image with --multiqcsif, or pass --checkversion false." >&2
            exit 1
        fi
    fi

    ${mqcBin} . \\
        --title "${projectId}" \\
        --filename ${projectId}_multiqc_report.html \\
        --config ${mqc_config} \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( ${mqcBin} --version 2>&1 | head -1 | sed 's/[^0-9]*//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${projectId}_multiqc_report.html
    mkdir -p ${projectId}_multiqc_report_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: ${params.multiqcversion}
    END_VERSIONS
    """
}
