process MULTI_METRICS_MQC {

    label 'process_low'

    // small generated json, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(section), path(csvs)

    output:
        path "*_mqc.json" , emit: mqc , optional: true

    script:
    // cellranger multi writes a LONG metrics_summary.csv - one metric per row,
    // with Category, Library Type, Grouped By, Group Name, Metric Name and
    // Metric Value - rather than the wide one count and atac produce. So it is
    // pivoted here into one multiqc table per ( Category, Library Type ):
    // rows are samples, columns are metric names.
    //
    // The csv is used rather than the web summary html on purpose. A multi
    // summary keeps everything under data.library_websummary / sample_websummary
    // in per assay tabs, and its barcode rank plot alone carries 30 to 50
    // traces per sample; the csv carries the same numbers in a form that
    // pivots without any of that.
    //
    // No backslashes or ${} below: this is a groovy triple quoted string, and
    // a regex would have to be double escaped to survive it.
    """
    python3 - ${section} ${csvs} <<'PYEOF'
import csv, json, os, sys

section = sys.argv[1]
paths   = sys.argv[2:]

# Category is Cells or Library; Library Type is Gene Expression, Antibody
# Capture, VDJ T and so on. One table per pair, ordered so the assay a reader
# cares about most comes first.
CATEGORY_ORDER = ["Cells", "Library"]
ASSAY_ORDER    = ["Gene Expression", "Antibody Capture", "CRISPR Guide Capture",
                  "Multiplexing Capture", "VDJ T", "VDJ B", "VDJ"]

# The Library tables run to twenty odd columns. These are shown; the rest are
# emitted but hidden, one click away under Configure columns. The Cells tables
# are small enough to show whole.
VISIBLE = {
    "Gene Expression": {
        "Estimated number of cells", "Mean reads per cell",
        "Number of reads in the library", "Sequencing saturation",
        "Confidently mapped to transcriptome", "Confidently mapped reads in cells",
        "Valid barcodes", "Q30 barcodes",
    },
    "Antibody Capture": {
        "Estimated number of cells", "Mean reads per cell",
        "Number of reads in the library", "Sequencing saturation",
        "Fraction antibody reads usable", "Antibody reads in cells",
        "Valid barcodes", "Q30 barcodes",
    },
    "VDJ T": {
        "Estimated number of cells", "Mean reads per cell",
        "Mean used reads per cell", "Fraction reads in cells",
        "Reads mapped to any V(D)J gene", "Valid barcodes", "Q30 barcodes",
    },
}
VISIBLE["VDJ B"] = VISIBLE["VDJ T"]


def tidy(name):
    # the links are <label>_<sample>_metrics_summary.csv, and for a run that is
    # not multiplexed the sample IS the label, giving EW4_SC2622_EW4_SC2622
    half = len(name) // 2
    if len(name) % 2 == 1 and name[half] == "_" and name[:half] == name[half + 1:]:
        return name[:half]
    return name


def number(value):
    if not isinstance(value, str):
        return value
    bare = value.replace(",", "").strip()
    pct = bare.endswith("%")
    if pct:
        bare = bare[:-1]
    for cast in (int, float):
        try:
            return cast(bare)
        except ValueError:
            continue
    return value


def rank(value, order):
    return order.index(value) if value in order else len(order)


# sample -> (category, library type) -> metric -> value
data = {}
pct_cols = {}

for path in paths:
    label = os.path.basename(path)
    for tail in ("_metrics_summary.csv", ".csv"):
        if label.endswith(tail):
            label = label[: -len(tail)]
            break
    label = tidy(label)

    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))

    if not rows or "Metric Name" not in rows[0]:
        sys.stderr.write("WARN: " + path + " is not a multi metrics_summary.csv, skipping\\n")
        continue

    # a metric can repeat within a section when the rows are grouped, e.g. per
    # fastq id. Qualify those with the group so columns cannot collide.
    seen = {}
    for row in rows:
        key = (row.get("Category") or "", row.get("Library Type") or "")
        seen.setdefault(key, {}).setdefault(row.get("Metric Name") or "", []).append(row)

    for key, metrics in seen.items():
        cells = data.setdefault(label, {}).setdefault(key, {})

        def put(name, hit):
            raw = hit.get("Metric Value") or ""
            cells[name] = number(raw)
            if raw.strip().endswith("%"):
                pct_cols.setdefault(key, set()).add(name)

        for metric, hits in metrics.items():
            if len(hits) == 1:
                put(metric, hits[0])
                continue

            # The same metric appears once per grouping, e.g. Number of reads
            # grouped by Fastq ID and again by Physical library ID. Qualify
            # with GROUPED BY, not Group Name: the group name is the fastq id
            # itself, which differs per sample and would give every sample its
            # own column instead of one shared one.
            by = {}
            for hit in hits:
                by.setdefault((hit.get("Grouped By") or "").strip(), []).append(hit)

            for grouped_by, bucket in by.items():
                if len(bucket) == 1:
                    put(metric + (" (" + grouped_by + ")" if grouped_by else ""), bucket[0])
                    continue
                # several rows within one grouping - a library sequenced more
                # than once - so the group name is what tells them apart
                for hit in bucket:
                    group = (hit.get("Group Name") or "").strip()
                    put(metric + " (" + (group or grouped_by) + ")", hit)

if not data:
    sys.stderr.write("no parsable metrics_summary.csv, nothing to render\\n")
    sys.exit(0)

sections = sorted(
    {k for per_label in data.values() for k in per_label},
    key=lambda k: (rank(k[0], CATEGORY_ORDER), rank(k[1], ASSAY_ORDER)))

for category, assay in sections:
    table, headers = {}, {}
    for label in sorted(data):
        cells = data[label].get((category, assay))
        if cells:
            table[label] = cells
            for name in cells:
                show = VISIBLE.get(assay)
                headers[name] = {
                    "title": name,
                    "description": name,
                    # the Cells tables are short, so nothing is hidden there
                    "hidden": category == "Library" and bool(show) and name not in show,
                }
    if not table:
        continue

    for name in pct_cols.get((category, assay), set()):
        if name in headers:
            headers[name]["suffix"] = "%"

    slug = (category + "_" + assay).lower().replace(" ", "_")
    sid  = section + "_" + slug
    title = "Cell Ranger multi: " + category + " - " + assay

    with open(sid + "_mqc.json", "w") as fh:
        json.dump({
            "id": sid,
            "section_name": title,
            "description": "From each per sample metrics_summary.csv.",
            "plot_type": "table",
            "pconfig": {"id": sid + "_table", "col1_header": "Sample"},
            "headers": headers,
            "data": table,
        }, fh, indent=1)
    print("wrote " + sid + "_mqc.json")
PYEOF
    """

    stub:
    """
    echo '{"id":"stub","section_name":"stub","plot_type":"table","data":{"a":{"b":1}}}' > ${section}_stub_mqc.json
    """
}
