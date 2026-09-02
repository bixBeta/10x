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

# Cell Ranger ARC prefixes its columns by assay. MultiQC's header "namespace"
# only shows in tooltips, so the assays are split into their own sections
# instead - each one narrow enough to read.
GROUPS = [
    ("library", "Cell Ranger ARC: library", None,    None),
    ("atac",    "Cell Ranger ARC: ATAC",    "ATAC ", {
        "Median high-quality fragments per cell",
        "TSS enrichment score",
        "Fraction of high-quality fragments in cells",
        "Fraction of transposition events in peaks in cells",
        "Number of peaks",
        "Sequenced read pairs",
    }),
    ("gex",     "Cell Ranger ARC: GEX",     "GEX ",  {
        "Median genes per cell",
        "Median UMI counts per cell",
        "Mean raw reads per cell",
        "Fraction of transcriptomic reads in cells",
        "Total genes detected",
        "Sequenced read pairs",
    }),
]
SKIP = {"label"}

def number(value):
    for cast in (int, float):
        try:
            return cast(value.replace(",", ""))
        except (ValueError, AttributeError):
            continue
    return value

for suffix, title, prefix, visible in GROUPS:
    data, headers = {}, {}
    for row in rows:
        sample = row.get("label") or row.get("Sample ID") or "sample"
        cells = {}
        for col, value in row.items():
            if col in SKIP or col is None:
                continue
            is_assay = col.startswith("ATAC ") or col.startswith("GEX ")
            if prefix is None:
                if is_assay:
                    continue
                name = col
            else:
                if not col.startswith(prefix):
                    continue
                name = col[len(prefix):]
            cells[name] = number(value)
            headers[name] = {
                "title": name,
                "description": col,
                "hidden": bool(visible) and name not in visible,
            }
        if cells:
            data[sample] = cells

    if not data:
        continue

    doc = {
        "id": section + "_" + suffix,
        "section_name": title,
        "description": (
            "From each summary.csv. MultiQC cannot parse cellranger-arc 2.2 web "
            "summaries itself (MultiQC issue 3609)."
        ),
        "plot_type": "table",
        "pconfig": {"id": section + "_" + suffix + "_table", "col1_header": "Library"},
        "headers": headers,
        "data": data,
    }
    with open(section + "_" + suffix + "_mqc.json", "w") as fh:
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
