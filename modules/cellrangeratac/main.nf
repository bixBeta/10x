chemistry       = params.chemistry
forceCells      = params.forceCells
peaks           = params.peaks
dimReduce       = params.dimReduce
checkversion    = params.checkversion
wantversion     = params.atacversion

// resolved in nextflow.config from --atacpath, or --atacversion + --programs
atacBin         = params.atacbin


process CELLRANGER_ATAC_COUNT {

    tag "$library"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "CELLRANGER_ATAC/${library}" , mode: "symlink", overwrite: true , pattern: "outs"

    // flat collections across libraries, for a quick look over a whole run.
    // Signac / ArchR want the fragments file and singlecell.csv as much as the
    // peak matrix, so all three are linked into one folder per library rather
    // than just the matrix the gex path publishes.
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "filtered_peak_bc_matrix.h5"
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "fragments.tsv.gz*"
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "singlecell.csv"

    publishDir "web_summary_htmls"      , mode: "symlink", overwrite: true , pattern: "*_web_summary.html"
    publishDir "summary_metrics"        , mode: "symlink", overwrite: true , pattern: "*_summary.csv"

    input:
        tuple val(library), val(label), val(fastqs), val(samples)
        val ref

    output:
        tuple val(library), path("outs")                        , emit: outs
        path "outs/web_summary.html"                            , emit: web_summary
        path "outs/summary.csv"                                 , emit: metrics
        path "outs/filtered_peak_bc_matrix.h5"                  , emit: filtered_h5   , optional: true
        path "outs/fragments.tsv.gz"                            , emit: fragments     , optional: true
        path "${library}_web_summary.html"                      , emit: run_web_summary , optional: true
        path "${library}_summary.csv"                           , emit: run_metrics     , optional: true
        path "filtered_peak_bc_matrix.h5"                       , emit: filtered_link   , optional: true
        path "fragments.tsv.gz*"                                , emit: fragment_links  , optional: true
        path "singlecell.csv"                                   , emit: singlecell_link , optional: true
        path "versions.yml"                                     , emit: versions

    script:

    // same directory-in, nothing-renamed contract as the gex path: several runs
    // of one library become comma separated lists of dirs and prefixes
    def fastqArg    = fastqs instanceof List ? fastqs.join(',') : "${fastqs}"
    def sampleArg   = samples instanceof List ? samples.join(',') : "${samples}"

    // cellranger-atac detects the assay configuration itself and takes no
    // 'auto' literal, unlike cellranger - so the flag is only passed when it
    // was asked for. The value that matters is ARC-v1: the ATAC library of a
    // Multiome experiment, counted on its own.
    def chemArg     = ( chemistry && chemistry != "auto" ) ? "--chemistry=${chemistry}" : ""
    def forceArg    = forceCells  ? "--force-cells=${forceCells}"  : ""
    def peaksArg    = peaks       ? "--peaks=${peaks}"             : ""
    def dimArg      = dimReduce   ? "--dim-reduce=${dimReduce}"    : ""

    """
    # The image or install is only labelled with a version - ask the tool itself
    # before counting anything, so a mislabelled sif cannot silently produce
    # results from the wrong cellranger-atac.
    if [ "${checkversion}" = "true" ] ; then
        have=\$( ${atacBin} --version 2>&1 | head -1 | sed 's/.*cellranger-atac[- ]//; s/[^0-9.].*\$//' )
        if [ -z "\$have" ] ; then
            echo "WARN: could not read a version from cellranger-atac --version, skipping the check" >&2
        elif [ "\$have" != "${wantversion}" ] ; then
            echo "ERROR: asked for cellranger-atac ${wantversion} but ${atacBin} reports \$have" >&2
            echo "       Set --atacversion to match, point at another image with --atacsif, or pass --checkversion false." >&2
            exit 1
        fi
    fi

    ${atacBin} count \\
        --id=${library} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --reference=${ref} \\
        --fastqs=${fastqArg} \\
        --sample=${sampleArg} \\
        ${chemArg} ${forceArg} ${peaksArg} ${dimArg}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${library}/outs outs

    # label prefixed links for the run level collections
    # if/then, not &&: nextflow runs task scripts under bash -ue, so a false
    # test as the last command of a line would abort the task
    if [ -e outs/web_summary.html ] ; then
        ln -s outs/web_summary.html ${library}_web_summary.html
    fi
    if [ -e outs/summary.csv ] ; then
        ln -s outs/summary.csv ${library}_summary.csv
    fi

    # a Signac ready folder: peak matrix, fragments + index, per barcode table
    if [ -e outs/filtered_peak_bc_matrix.h5 ] ; then
        ln -s outs/filtered_peak_bc_matrix.h5 filtered_peak_bc_matrix.h5
    fi
    if [ -e outs/fragments.tsv.gz ] ; then
        ln -s outs/fragments.tsv.gz fragments.tsv.gz
    fi
    if [ -e outs/fragments.tsv.gz.tbi ] ; then
        ln -s outs/fragments.tsv.gz.tbi fragments.tsv.gz.tbi
    fi
    if [ -e outs/singlecell.csv ] ; then
        ln -s outs/singlecell.csv singlecell.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-atac: \$( ${atacBin} --version 2>&1 | head -1 | sed 's/.*cellranger-atac[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p outs
    touch outs/web_summary.html
    touch outs/summary.csv
    touch outs/filtered_peak_bc_matrix.h5
    touch outs/fragments.tsv.gz outs/fragments.tsv.gz.tbi
    touch outs/singlecell.csv

    ln -s outs/web_summary.html ${library}_web_summary.html
    ln -s outs/summary.csv      ${library}_summary.csv
    ln -s outs/filtered_peak_bc_matrix.h5 filtered_peak_bc_matrix.h5
    ln -s outs/fragments.tsv.gz           fragments.tsv.gz
    ln -s outs/fragments.tsv.gz.tbi       fragments.tsv.gz.tbi
    ln -s outs/singlecell.csv             singlecell.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger-atac: ${params.atacversion}
    END_VERSIONS
    """
}
