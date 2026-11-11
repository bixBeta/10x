createBam       = params.createBam
nosecondary     = params.nosecondary
r1length        = params.spatialR1length
r2length        = params.spatialR2length
checkversion    = params.checkversion
wantversion     = params.srversion

// resolved in nextflow.config from --srpath, or --srversion + --programs
srBin           = params.srbin


process SPACERANGER_COUNT {

    tag "$library"

    // cpus / memory come from --localcores and --localmem, see nextflow.config

    publishDir "SPACERANGER/${library}" , mode: "symlink", overwrite: true , pattern: "outs"

    // flat collections across capture areas, for a quick look over a whole run.
    // Seurat's Load10X_Spatial and squidpy both want the spatial/ folder -
    // tissue_positions.csv and scalefactors_json.json - as much as the matrix,
    // so it is linked alongside rather than left inside outs/ only.
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "filtered_feature_bc_matrix.h5"
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "spatial"
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "binned_outputs"
    publishDir "filtered_counts/${library}" , mode: "symlink", overwrite: true , pattern: "cloupe*.cloupe"

    publishDir "web_summary_htmls"      , mode: "symlink", overwrite: true , pattern: "*_web_summary.html"
    publishDir "summary_metrics"        , mode: "symlink", overwrite: true , pattern: "*_metrics_summary.csv"

    input:
        tuple val(library), val(label), val(fastqs), val(samples),
              val(slide), val(area),
              path(cytaimage), path(image), path(darkimage), path(colorizedimage),
              path(probeset), path(alignment), path(slidefile)
        val ref

    output:
        tuple val(library), path("outs")                        , emit: outs
        path "outs/web_summary.html"                            , emit: web_summary
        path "outs/metrics_summary.csv"                         , emit: metrics
        path "${library}_web_summary.html"                      , emit: run_web_summary , optional: true
        path "${library}_metrics_summary.csv"                   , emit: run_metrics     , optional: true
        path "filtered_feature_bc_matrix.h5"                    , emit: filtered_link   , optional: true
        path "spatial"                                          , emit: spatial_link    , optional: true
        path "binned_outputs"                                   , emit: binned_link     , optional: true
        path "cloupe*.cloupe"                                   , emit: cloupe_link     , optional: true
        path "versions.yml"                                     , emit: versions

    script:

    // --create-bam is spaceranger 3.0 and later; 2.x took --no-bam instead, so
    // an older --srversion would need that line changed as well as the number.
    //
    // same directory-in, nothing-renamed contract as the gex path: a capture
    // area sequenced more than once becomes comma separated dirs and prefixes
    def fastqArg    = fastqs instanceof List ? fastqs.join(',') : "${fastqs}"
    def sampleArg   = samples instanceof List ? samples.join(',') : "${samples}"

    // Either a real slide serial plus its capture area, or one of spaceranger's
    // unknown-slide keywords when the serial was lost. The two are mutually
    // exclusive flags, and main.nf has already rejected the combinations that
    // are not, so this only has to pick the right one.
    def slideArg    = slide?.startsWith("visium-")
                    ? "--unknown-slide=${slide}"
                    : "--slide=${slide} --area=${area}"

    // The images are staged, so these are bare names in the task dir. An absent
    // optional input arrives as an empty list, which is falsy in groovy.
    def imageArgs   = [
        cytaimage       ? "--cytaimage=${cytaimage}"             : "",
        image           ? "--image=${image}"                     : "",
        darkimage       ? "--darkimage=${darkimage}"             : "",
        colorizedimage  ? "--colorizedimage=${colorizedimage}"   : "",
    ].findAll{ it }.join(' ')

    def probeArg    = probeset  ? "--probe-set=${probeset}"        : ""
    def alignArg    = alignment ? "--loupe-alignment=${alignment}" : ""
    def slidefileArg = slidefile ? "--slidefile=${slidefile}"      : ""
    // toString, not plain truthiness: --nosecondary false arrives as the STRING
    // "false", which groovy counts as true
    def nosecArg    = nosecondary.toString() == "true" ? "--nosecondary" : ""

    // Read lengths are their own params here, NOT --r1length: a Visium library
    // is not a gex one. gex's default of 28 would silently trim a 43 nt
    // CytAssist R1, so the flag is left out altogether unless it was asked for
    // and spaceranger reads the length off the reads itself.
    def r1Arg       = r1length ? "--r1-length=${r1length}" : ""
    def r2Arg       = r2length ? "--r2-length=${r2length}" : ""

    """
    # The image or install is only labelled with a version - ask the tool itself
    # before counting anything, so a mislabelled sif cannot silently produce
    # results from the wrong spaceranger.
    if [ "${checkversion}" = "true" ] ; then
        have=\$( ${srBin} --version 2>&1 | head -1 | sed 's/.*spaceranger[- ]//; s/[^0-9.].*\$//' )
        if [ -z "\$have" ] ; then
            echo "WARN: could not read a version from spaceranger --version, skipping the check" >&2
        elif [ "\$have" != "${wantversion}" ] ; then
            echo "ERROR: asked for spaceranger ${wantversion} but ${srBin} reports \$have" >&2
            echo "       Set --srversion to match, point at another image with --srsif, or pass --checkversion false." >&2
            exit 1
        fi
    fi

    ${srBin} count \\
        --id=${library} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --transcriptome=${ref} \\
        --fastqs=${fastqArg} \\
        --sample=${sampleArg} \\
        --create-bam=${createBam} \\
        ${slideArg} \\
        ${imageArgs} \\
        ${probeArg} ${alignArg} ${slidefileArg} ${nosecArg} ${r1Arg} ${r2Arg}

    # lift outs/ up so publishDir and downstream modules see a stable path
    mv ${library}/outs outs

    # label prefixed links for the run level collections
    # if/then, not &&: nextflow runs task scripts under bash -ue, so a false
    # test as the last command of a line would abort the task
    if [ -e outs/web_summary.html ] ; then
        ln -s outs/web_summary.html ${library}_web_summary.html
    fi
    if [ -e outs/metrics_summary.csv ] ; then
        ln -s outs/metrics_summary.csv ${library}_metrics_summary.csv
    fi

    # a Seurat ready folder: the matrix, the spatial/ metadata beside it, and
    # the loupe file a collaborator will actually open
    if [ -e outs/filtered_feature_bc_matrix.h5 ] ; then
        ln -s outs/filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5
    fi
    if [ -d outs/spatial ] ; then
        ln -s outs/spatial spatial
    fi
    # Visium HD suffixes it with the bin size - cloupe_008um.cloupe - and can
    # write one per bin, so this is a glob rather than the plain name a
    # standard Visium run produces.
    for f in outs/cloupe*.cloupe ; do
        if [ -e "\$f" ] ; then
            ln -s "\$f" "\$( basename "\$f" )"
        fi
    done

    # Visium HD writes no top level matrix: every bin size gets its own
    # directory under binned_outputs/, each with its own spatial/ folder.
    if [ -d outs/binned_outputs ] ; then
        ln -s outs/binned_outputs binned_outputs
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: \$( ${srBin} --version 2>&1 | head -1 | sed 's/.*spaceranger[- ]//; s/[^0-9.].*\$//' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p outs/spatial
    touch outs/web_summary.html
    touch outs/metrics_summary.csv
    touch outs/filtered_feature_bc_matrix.h5
    touch outs/cloupe_008um.cloupe
    touch outs/spatial/tissue_positions.csv outs/spatial/scalefactors_json.json

    ln -s outs/web_summary.html    ${library}_web_summary.html
    ln -s outs/metrics_summary.csv ${library}_metrics_summary.csv
    ln -s outs/filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5
    ln -s outs/spatial                       spatial
    ln -s outs/cloupe_008um.cloupe           cloupe_008um.cloupe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: ${params.srversion}
    END_VERSIONS
    """
}
