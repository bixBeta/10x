process COMBINE_METRICS {

    label 'process_low'

    // a generated summary, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(outname), val(section), path(csvs)

    output:
        path outname          , emit: combined
        path "*_mqc.json"     , emit: mqc      , optional: true

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
import csv, json, sys

table, section = sys.argv[1], sys.argv[2]
with open(table, newline="") as fh:
    rows = list(csv.DictReader(fh))

# Cell Ranger ARC prefixes its columns by assay, so the prefix is the grouping.
# A handful of headline metrics stay visible; the rest are one click away under
# "Configure columns" rather than 30 columns wide by default.
VISIBLE = {
    "Estimated number of cells",
    "ATAC Median high-quality fragments per cell",
    "ATAC TSS enrichment score",
    "ATAC Fraction of high-quality fragments in cells",
    "ATAC Fraction of transposition events in peaks in cells",
    "GEX Median genes per cell",
    "GEX Median UMI counts per cell",
    "GEX Mean raw reads per cell",
    "GEX Fraction of transcriptomic reads in cells",
}
SKIP = {"label", "Sample ID"}

def number(value):
    try:
        return int(value.replace(",", ""))
    except (ValueError, AttributeError):
        pass
    try:
        return float(value.replace(",", ""))
    except (ValueError, AttributeError):
        return value

data, headers = {}, {}
for row in rows:
    sample = row.get("label") or row.get("Sample ID") or "sample"
    data[sample] = {}
    for col, value in row.items():
        if col in SKIP or col is None:
            continue
        if col.startswith("ATAC "):
            namespace, title = "ATAC", col[5:]
        elif col.startswith("GEX "):
            namespace, title = "GEX", col[4:]
        else:
            namespace, title = "ARC", col
        data[sample][col] = number(value)
        headers[col] = {
            "title": title,
            "namespace": namespace,
            "description": col,
            "hidden": col not in VISIBLE,
        }

doc = {
    "id": section,
    "section_name": "Cell Ranger ARC metrics",
    "description": (
        "From each summary.csv. Columns are grouped by assay; use Configure "
        "columns for the rest. MultiQC cannot parse cellranger-arc 2.2 web "
        "summaries itself (MultiQC issue 3609)."
    ),
    "plot_type": "table",
    "pconfig": {"id": section + "_table", "namespace": "Cell Ranger ARC", "col1_header": "Library"},
    "headers": headers,
    "data": data,
}

# JSON, which MultiQC reads as *_mqc.json custom content. It is also valid
# YAML, and it removes the indentation traps of hand written block scalars.
with open(section + "_mqc.json", "w") as fh:
    json.dump(doc, fh, indent=1)
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
