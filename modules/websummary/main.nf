process WEBSUMMARY_MQC {

    label 'process_low'

    // small generated json, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(section), path(summaries)

    output:
        path "*_mqc.json" , emit: mqc , optional: true

    script:
    // A Cell Ranger web summary carries the whole report as one JSON blob on a
    // "const data = {...}" line - tables, plot traces, axis labels and all.
    // MultiQC's own cellranger module reads it exactly this way, but there is
    // no cellranger-atac module to do it, so the pipeline does.
    //
    // No backslashes or ${} below: this is a groovy triple quoted string, and
    // a regex would have to be double escaped to survive it. String methods
    // are used instead, which keeps the python readable here.
    """
    python3 - ${section} ${summaries} <<'PYEOF'
import json, os, sys

section = sys.argv[1]
paths   = sys.argv[2:]

# ---- what to lift out of the web summary -----------------------------------
# Keys confirmed against a real cellranger-atac 2.2.0 web_summary.html. The
# tables are the ones the report shows, kept grouped the same way, which is
# richer than summary.csv: the mapping table's nucleosome fractions are not in
# the csv at all.
TABLES = [
    ("sequencing", "Sequencing",          "sequencing_info_table"),
    ("cells",      "Cells",               "cell_metrics_table"),
    ("mapping",    "Mapping",             "mapping_table"),
    ("targeting",  "Targeting",           "targeting_table"),
    ("complexity", "Complexity",          "bulk_complexity_table"),
]

# ( suffix, title, key, one series per trace, log axes )
# The tSNE plots are deliberately left out: a 2D cluster scatter is not
# comparable between libraries, which is the whole point of aggregating.
PLOTS = [
    ("knee",       "Barcode Rank",          "barcode_knee_plot",        True,  True),
    ("fragments",  "Fragment Distribution", "fragments_per_cell_plot",  True,  True),
    ("tss",        "TSS Enrichment",        "tss_enrichment_plot",      False, False),
    ("insert",     "Insert Size",           "non_log_insert_size_plot", False, False),
    ("efficiency", "Assay Efficiency",      "assay_efficiency_plot",    False, False),
]

# a plot trace can carry 100k points; four libraries of that makes the report
# unusable, and MultiQC caps its own cellranger series similarly
MAX_POINTS = 1000


def load(path):
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("const data") and "=" in line:
                blob = line.split("=", 1)[1].strip().rstrip(";")
                try:
                    return json.loads(blob)
                except ValueError:
                    continue
    return None


def clean(text):
    # titles arrive as ' Cells<br><sup>HS2 - </sup>': drop the subtitle, strip
    # any tags, collapse whitespace
    if not isinstance(text, str):
        text = "" if text is None else str(text)
    text = text.split("<br>")[0]
    out, depth = [], 0
    for ch in text:
        if ch == "<":
            depth += 1
        elif ch == ">":
            depth = depth - 1 if depth > 0 else 0
        elif depth == 0:
            out.append(ch)
    return " ".join("".join(out).split())


def axis(layout, name):
    return clean((layout.get(name) or {}).get("title"))


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


def thin(xs, ys):
    pairs = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pairs) <= MAX_POINTS:
        return pairs
    step = len(pairs) / float(MAX_POINTS)
    keep = [pairs[int(i * step)] for i in range(MAX_POINTS)]
    keep[-1] = pairs[-1]
    return keep


def write(doc, name):
    with open(name, "w") as fh:
        json.dump(doc, fh, indent=1)
    print("wrote " + name)


# ---- read every summary once ------------------------------------------------
reports = []
for path in paths:
    doc = load(path)
    if doc is None:
        sys.stderr.write("WARN: no embedded json in " + path + ", skipping\\n")
        continue
    root = doc.get("summary", doc)
    if "cell_metrics_table" not in root:
        sys.stderr.write("WARN: " + path + " is not a cellranger-atac summary, skipping\\n")
        continue
    # the label is how every other section names this library
    label = os.path.basename(path)
    for tail in ("_web_summary.html", ".html"):
        if label.endswith(tail):
            label = label[: -len(tail)]
            break
    reports.append((label or (root.get("sample") or {}).get("id") or "sample", root))

if not reports:
    sys.stderr.write("no parsable web summaries, nothing to render\\n")
    sys.exit(0)

reports.sort()

# ---- tables ----------------------------------------------------------------
for suffix, title, key in TABLES:
    data, headers, pct_cols = {}, {}, set()
    for label, root in reports:
        node = root.get(key) or {}
        rows = node.get("rows") or []
        cells = {}
        for row in rows:
            if not isinstance(row, list) or len(row) < 2:
                continue
            name, raw = clean(row[0]), row[1]
            if not name:
                continue
            cells[name] = number(raw)
            if isinstance(raw, str) and raw.strip().endswith("%"):
                pct_cols.add(name)
            headers[name] = {"title": name, "description": name}
        if cells:
            data[label] = cells
    if not data:
        continue
    for name in pct_cols:
        headers[name]["suffix"] = "%"
    sid = section + "_ws_" + suffix
    write({
        "id": sid,
        "section_name": "Cell Ranger ATAC: " + title,
        "description": "From each web_summary.html. MultiQC ships no cellranger-atac module, so the pipeline parses it.",
        "plot_type": "table",
        "pconfig": {"id": sid + "_table", "col1_header": "Library"},
        "headers": headers,
        "data": data,
    }, sid + "_mqc.json")

# ---- plots -----------------------------------------------------------------
for suffix, title, key, per_trace, logaxes in PLOTS:
    series, xlab, ylab = {}, "", ""
    for label, root in reports:
        node = root.get(key) or {}
        layout = node.get("layout") or {}
        xlab = xlab or axis(layout, "xaxis")
        ylab = ylab or axis(layout, "yaxis")
        for trace in node.get("data") or []:
            if not isinstance(trace, dict):
                continue
            pairs = thin(trace.get("x") or [], trace.get("y") or [])
            if not pairs:
                continue
            tname = clean(trace.get("name"))
            name = label + " " + tname if (per_trace and tname) else label
            # json keys are strings; multiqc coerces them with float() for a
            # line plot, so the x axis stays numeric
            series[name] = dict((str(x), y) for x, y in pairs)
    if not series:
        continue
    sid = section + "_ws_" + suffix
    pconfig = {"id": sid + "_plot", "title": "Cell Ranger ATAC: " + title,
               "xlab": xlab, "ylab": ylab}
    if logaxes:
        pconfig["xlog"] = True
        pconfig["ylog"] = True
    write({
        "id": sid,
        "section_name": "Cell Ranger ATAC: " + title,
        "description": "From each web_summary.html, thinned to " + str(MAX_POINTS) + " points per series.",
        "plot_type": "linegraph",
        "pconfig": pconfig,
        "data": series,
    }, sid + "_mqc.json")
PYEOF
    """

    stub:
    """
    echo '{"id":"stub","section_name":"stub","plot_type":"table","data":{"a":{"b":1}}}' > ${section}_ws_stub_mqc.json
    """
}
