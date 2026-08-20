process COMBINE_METRICS {

    label 'process_low'

    // a generated summary, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(outname), path(csvs)

    output:
        path outname , emit: combined

    script:
    // Header from the first file, rows from all of them, each prefixed with the
    // library it came from. cellranger's own metrics_summary.csv carries no
    // sample column, so without this the stacked rows cannot be told apart.
    // The label is the file name minus the suffix, which is how the per library
    // links were named.
    def files = csvs.collect { it.name }.sort().join(' ')
    """
    awk -F, '
        FNR == 1 {
            lab = FILENAME
            sub(/.*\\//, "", lab)
            sub(/_metrics_summary\\.csv\$/, "", lab)
            sub(/_summary\\.csv\$/, "", lab)
            if ( NR == 1 ) print "label," \$0
            next
        }
        { print lab "," \$0 }
    ' ${files} > ${outname}
    """

    stub:
    def files = csvs.collect { it.name }.sort().join(' ')
    """
    awk -F, '
        FNR == 1 {
            lab = FILENAME
            sub(/.*\\//, "", lab)
            sub(/_metrics_summary\\.csv\$/, "", lab)
            sub(/_summary\\.csv\$/, "", lab)
            if ( NR == 1 ) print "label," \$0
            next
        }
        { print lab "," \$0 }
    ' ${files} > ${outname}
    """
}
