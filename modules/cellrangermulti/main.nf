createBam       = params.createBam
chemistry       = params.chemistry
expectCells     = params.expectCells
forceCells      = params.forceCells
introns         = params.introns
r1length        = params.r1length
r2length        = params.r2length
checkversion    = params.checkversion
wantversion     = params.crversion

// resolved in nextflow.config from --crpath, or --crversion + --programs
crBin           = params.crbin


process CELLRANGER_MULTI {

    tag "$label"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER_MULTI/${label}" , mode: "symlink", overwrite: true , pattern: "outs"
    publishDir "CELLRANGER_MULTI/${label}" , mode: "copy"   , overwrite: true , pattern: "*_multi_config.csv"

    // flat collections across labels, for a quick look over a whole run
    publishDir "filtered_counts/${label}"  , mode: "symlink", overwrite: true , pattern: "*_filtered_feature_bc_matrix.h5"
    publishDir "web_summary_htmls"         , mode: "symlink", overwrite: true , pattern: "*_web_summary.html"
    publishDir "summary_metrics"           , mode: "symlink", overwrite: true , pattern: "*_metrics_summary.csv"

    input:
        tuple val(label), val(rows)
        val gexref
        val vdjref
        val featureref

    output:
        tuple val(label), path("outs")                  , emit: outs
        path "${label}_multi_config.csv"                , emit: config
        path "*_web_summary.html"                       , emit: run_web_summary , optional: true
        path "*_metrics_summary.csv"                    , emit: run_metrics     , optional: true
        path "*_filtered_feature_bc_matrix.h5"          , emit: filtered_link   , optional: true
        path "versions.yml"                             , emit: versions

    script:

    // cellranger multi takes one config csv per run, which this builds from the
    // sample sheet. Sections are only emitted when a library needs them: [vdj]
    // for a VDJ library, [feature] for Antibody Capture or CRISPR, since an
    // empty section is an error rather than a no-op.
    def types    = rows.collect { it[2] }.unique()
    def hasVdj   = types.any { it.startsWith("VDJ") }
    def hasFeat  = types.any { it in ["Antibody Capture", "CRISPR Guide Capture"] }

    def lines = []
    lines << "[gene-expression]"
    lines << "reference,${gexref}"
    lines << "create-bam,${createBam}"
    if( chemistry && chemistry != "auto" ) lines << "chemistry,${chemistry}"
    if( expectCells )                      lines << "expect-cells,${expectCells}"
    if( forceCells )                       lines << "force-cells,${forceCells}"
    if( introns != null )                  lines << "include-introns,${introns}"
    if( r1length )                         lines << "r1-length,${r1length}"
    if( r2length )                         lines << "r2-length,${r2length}"

    if( hasVdj ) {
        lines << ""
        lines << "[vdj]"
        lines << "reference,${vdjref}"
    }

    if( hasFeat ) {
        lines << ""
        lines << "[feature]"
        lines << "reference,${featureref}"
    }

    lines << ""
    lines << "[libraries]"
    lines << "fastq_id,fastqs,lanes,feature_types"
    rows.each { fq, sample, ltype, lanes ->
        lines << "${sample},${fq},${lanes ?: 'any'},${ltype}"
    }

    def writeConfig = lines.collect { "echo '${it}' >> ${label}_multi_config.csv" }.join('\n    ')

    """
    : > ${label}_multi_config.csv
    ${writeConfig}

    echo "--- ${label}_multi_config.csv ---"
    cat ${label}_multi_config.csv

    # The image or install is only labelled with a version - ask the tool itself
    # before counting anything, so a mislabelled sif cannot silently produce
    # results from the wrong cellranger.
    if [ "${checkversion}" = "true" ] ; then
        have=\$( ${crBin} --version 2>&1 | head -1 | sed 's/.*cellranger[- ]//; s/[^0-9.].*\$//' )
        if [ -z "\$have" ] ; then
            echo "WARN: could not read a version from cellranger --version, skipping the check" >&2
        elif [ "\$have" != "${wantversion}" ] ; then
            echo "ERROR: asked for cellranger ${wantversion} but ${crBin} reports \$have" >&2
            echo "       Set --crversion to match, point at another image with --crsif, or pass --checkversion false." >&2
            exit 1
        fi
    fi

    ${crBin} multi \\
        --id=${label} \\
        --csv=${label}_multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${label}/outs outs

    # multi writes one directory per sample under per_sample_outs: one for a
    # plain run, several when it is multiplexed. Each is linked with the label
    # in front so the flat collections cannot collide.
    # if/then, not &&: nextflow runs task scripts under bash -ue, so a false
    # test as the last command of a line would abort the task
    for d in outs/per_sample_outs/*/ ; do
        if [ -d "\$d" ] ; then
            s=\$( basename "\$d" )
            if [ -e "\${d}web_summary.html" ] ; then
                ln -s "\${d}web_summary.html" ${label}_\${s}_web_summary.html
            fi
            if [ -e "\${d}metrics_summary.csv" ] ; then
                ln -s "\${d}metrics_summary.csv" ${label}_\${s}_metrics_summary.csv
            fi
            if [ -e "\${d}count/sample_filtered_feature_bc_matrix.h5" ] ; then
                ln -s "\${d}count/sample_filtered_feature_bc_matrix.h5" ${label}_\${s}_filtered_feature_bc_matrix.h5
            fi
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$( ${crBin} --version 2>&1 | head -1 | sed 's/.*cellranger[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    def types    = rows.collect { it[2] }.unique()
    def hasVdj   = types.any { it.startsWith("VDJ") }
    def hasFeat  = types.any { it in ["Antibody Capture", "CRISPR Guide Capture"] }
    def lines = [ "[gene-expression]", "reference,${gexref}", "create-bam,${createBam}" ]
    if( hasVdj )  lines.addAll([ "", "[vdj]", "reference,${vdjref}" ])
    if( hasFeat ) lines.addAll([ "", "[feature]", "reference,${featureref}" ])
    lines.addAll([ "", "[libraries]", "fastq_id,fastqs,lanes,feature_types" ])
    rows.each { fq, sample, ltype, lanes -> lines << "${sample},${fq},${lanes ?: 'any'},${ltype}" }
    def writeConfig = lines.collect { "echo '${it}' >> ${label}_multi_config.csv" }.join('\n    ')
    """
    : > ${label}_multi_config.csv
    ${writeConfig}

    mkdir -p outs/per_sample_outs/${label}/count
    touch outs/per_sample_outs/${label}/web_summary.html
    touch outs/per_sample_outs/${label}/metrics_summary.csv
    touch outs/per_sample_outs/${label}/count/sample_filtered_feature_bc_matrix.h5

    ln -s outs/per_sample_outs/${label}/web_summary.html    ${label}_${label}_web_summary.html
    ln -s outs/per_sample_outs/${label}/metrics_summary.csv ${label}_${label}_metrics_summary.csv
    ln -s outs/per_sample_outs/${label}/count/sample_filtered_feature_bc_matrix.h5 ${label}_${label}_filtered_feature_bc_matrix.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: ${params.crversion}
    END_VERSIONS
    """
}
