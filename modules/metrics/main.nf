process COMBINE_METRICS {

    label 'process_low'

    // a generated summary, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(outname), val(section), path(csvs)

    output:
        path outname          , emit: combined
        path "*_mqc.yml"      , emit: mqc      , optional: true

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

    # MultiQC's cellranger_arc module cannot parse cellranger-arc 2.2 web
    # summaries ( MultiQC issue #3609, open ), so for arc the pipeline renders
    # the metrics itself as custom content. csv, not awk: these tables quote
    # fields that contain commas.
    if [ -n "${section}" ] ; then
        if command -v python3 > /dev/null 2>&1 ; then
            python3 - "${outname}" "${section}" <<'PYEOF'
import csv, html, sys

table, section = sys.argv[1], sys.argv[2]
with open(table, newline="") as fh:
    rows = list(csv.reader(fh))

if rows:
    head, body = rows[0], rows[1:]
    cells = "".join(f"<th>{html.escape(c)}</th>" for c in head)
    parts = ['<table class="table table-condensed">', f"<thead><tr>{cells}</tr></thead>", "<tbody>"]
    for r in body:
        parts.append("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in r) + "</tr>")
    parts.append("</tbody></table>")

    # ONE line, indented: every line of a "data: |" block scalar has to be
    # indented past the key, and an unindented line silently ends the block and
    # makes the file invalid YAML.
    with open(f"{section}_mqc.yml", "w") as fh:
        fh.write(f"id: {section}\\n")
        fh.write("section_name: Cell Ranger ARC metrics\\n")
        fh.write("description: Read from each summary.csv, since MultiQC cannot parse cellranger-arc 2.2 web summaries (MultiQC issue 3609).\\n")
        fh.write("plot_type: html\\n")
        fh.write("data: |\\n")
        fh.write("  " + "".join(parts) + "\\n")
PYEOF
        else
            echo "WARN: no python3, skipping the ${section} MultiQC section" >&2
        fi
    fi
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
