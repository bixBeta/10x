process COMBINE_METRICS {

    label 'process_low'

    // a generated summary, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(outname), path(csvs)

    output:
        path outname , emit: combined

    script:
    // header from the first file, rows from all of them
    def files = csvs.collect { it.name }.sort().join(' ')
    """
    awk '(NR == 1) || (FNR > 1)' ${files} > ${outname}
    """

    stub:
    def files = csvs.collect { it.name }.sort().join(' ')
    """
    awk '(NR == 1) || (FNR > 1)' ${files} > ${outname}
    """
}
