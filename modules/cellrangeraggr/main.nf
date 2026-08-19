normalize       = params.normalize
checkversion    = params.checkversion
wantversion     = params.crversion

// resolved in nextflow.config from --crpath, or --crversion + --programs
crBin           = params.crbin


process CELLRANGER_AGGR {

    tag "$label"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER_AGGR/${label}" , mode: "symlink", overwrite: true , pattern: "outs/**"
    publishDir "CELLRANGER_AGGR/${label}" , mode: "copy"   , overwrite: true , pattern: "*_aggr.csv"
    publishDir "web_summary_htmls"        , mode: "symlink", overwrite: true , pattern: "*_aggr_web_summary.html"

    input:
        tuple val(label), val(libraries), path(molecule_h5s)

    output:
        tuple val(label), path("outs")                          , emit: outs
        path "${label}_aggr.csv"                                , emit: csv
        path "${label}_aggr_web_summary.html"                   , emit: run_web_summary , optional: true
        path "versions.yml"                                     , emit: versions

    script:

    // cellranger aggr takes a csv of one row per library. The molecule files
    // arrive already named <library>_molecule_info.h5, so the rows can name
    // them directly and stay in step with the library list.
    def csvRows = libraries.collect { "echo '${it},${it}_molecule_info.h5' >> ${label}_aggr.csv" }.join('\n    ')
    def normArg = normalize ? "--normalize=${normalize}" : ""

    """
    echo 'sample_id,molecule_h5' > ${label}_aggr.csv
    ${csvRows}

    echo "--- ${label}_aggr.csv ---"
    cat ${label}_aggr.csv

    if [ "${checkversion}" = "true" ] ; then
        have=\$( ${crBin} --version 2>&1 | head -1 | sed 's/.*cellranger[- ]//; s/[^0-9.].*\$//' )
        if [ -n "\$have" ] && [ "\$have" != "${wantversion}" ] ; then
            echo "ERROR: asked for cellranger ${wantversion} but ${crBin} reports \$have" >&2
            exit 1
        fi
    fi

    ${crBin} aggr \\
        --id=${label} \\
        --csv=${label}_aggr.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        ${normArg}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${label}/outs outs

    if [ -e outs/web_summary.html ] ; then
        ln -s outs/web_summary.html ${label}_aggr_web_summary.html
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$( ${crBin} --version 2>&1 | head -1 | sed 's/.*cellranger[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    def csvRows = libraries.collect { "echo '${it},${it}_molecule_info.h5' >> ${label}_aggr.csv" }.join('\n    ')
    """
    echo 'sample_id,molecule_h5' > ${label}_aggr.csv
    ${csvRows}

    mkdir -p outs
    touch outs/web_summary.html
    ln -s outs/web_summary.html ${label}_aggr_web_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: ${params.crversion}
    END_VERSIONS
    """
}
