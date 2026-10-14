process WEBSUMMARY_MQC {

    label 'process_low'

    // small generated json, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(section), path(summaries)

    output:
        path "*_mqc.json" , emit: mqc      , optional: true
        path "*_mqc.html" , emit: mqc_html , optional: true

    script:
    // A Cell Ranger web summary carries the whole report as one JSON blob on a
    // "const data = {...}" line. MultiQC's own cellranger module reads it this
    // way, but there is no cellranger-atac module, so the pipeline does.
    //
    // The metric tables become multiqc tables. The PLOTS are not redrawn:
    // each one is already a plotly figure ( data + layout ), and multiqc
    // bundles plotly, so the figure is handed to Plotly.newPlot as it stands.
    // Rebuilding them as multiqc line graphs was tried and looked wrong -
    // cellranger's own styling is lost and the per base noise it hides in a
    // small panel reads as a hairy band at report width.
    //
    // No backslashes or ${} below: this is a groovy triple quoted string, and
    // a regex would have to be double escaped to survive it. String methods
    // are used instead, which keeps the python readable here.
    """
    python3 - ${section} ${summaries} <<'PYEOF'
import json, os, sys

section = sys.argv[1]
paths   = sys.argv[2:]

# Keys confirmed against a real cellranger-atac 2.2.0 web_summary.html. The
# tables are kept grouped the way the report groups them, which is richer than
# summary.csv: the mapping table's nucleosome fractions are not in the csv.
TABLES = [
    ("sequencing", "Sequencing",          "sequencing_info_table"),
    ("cells",      "Cells",               "cell_metrics_table"),
    ("mapping",    "Mapping",             "mapping_table"),
    ("targeting",  "Targeting",           "targeting_table"),
    ("complexity", "Complexity",          "bulk_complexity_table"),
]

# Each of these is a plotly figure in the web summary and is passed through
# untouched apart from naming, colour and the per library subtitle. The tSNE
# scatters are left out: a per library 2D embedding is not comparable between
# libraries, which is the whole point of putting them on one axis.
PLOTS = [
    ("knee",       "Barcode Rank",          "barcode_knee_plot"),
    ("fragments",  "Fragment Distribution", "fragments_per_cell_plot"),
    ("tss",        "TSS Enrichment",        "tss_enrichment_plot"),
    ("insert",     "Insert Size",           "non_log_insert_size_plot"),
    ("efficiency", "Assay Efficiency",      "assay_efficiency_plot"),
]

# one colour per library, so a library keeps its colour across every plot
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
DASHES = ["solid", "dot", "dash", "dashdot"]


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


def number(value):
    if not isinstance(value, str):
        return value
    bare = value.replace(",", "").strip()
    if bare.endswith("%"):
        bare = bare[:-1]
    for cast in (int, float):
        try:
            return cast(bare)
        except ValueError:
            continue
    return value


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

for suffix, title, key in TABLES:
    data, headers, pct_cols = {}, {}, set()
    for label, root in reports:
        rows = (root.get(key) or {}).get("rows") or []
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
    with open(sid + "_mqc.json", "w") as fh:
        json.dump({
            "id": sid,
            "section_name": "Cell Ranger ATAC: " + title,
            "description": "From each web_summary.html. MultiQC ships no cellranger-atac module, so the pipeline parses it.",
            "plot_type": "table",
            "pconfig": {"id": sid + "_table", "col1_header": "Library"},
            "headers": headers,
            "data": data,
        }, fh, indent=1)
    print("wrote " + sid + "_mqc.json")


def plain_title(layout):
    # ' Cells<br><sup>HS2 - </sup>' names one library, which is wrong once
    # several are overlaid. Keep the layout otherwise exactly as cellranger
    # built it.
    title = layout.get("title")
    if isinstance(title, dict):
        title["text"] = clean(title.get("text"))
    else:
        layout["title"] = clean(title)


for suffix, title, key in PLOTS:
    traces, layout = [], None
    for i, (label, root) in enumerate(reports):
        fig = root.get(key) or {}
        if not fig.get("data"):
            continue
        if layout is None:
            layout = dict(fig.get("layout") or {})
        for j, trace in enumerate(fig["data"]):
            if not isinstance(trace, dict):
                continue
            trace = dict(trace)
            tname = clean(trace.get("name"))
            trace["name"] = label + (" " + tname if tname else "")
            # the library picks the colour, the trace within it the dash, so
            # Cells and Non-cells stay apart without inventing a second hue
            line = dict(trace.get("line") or {})
            line["color"] = PALETTE[i % len(PALETTE)]
            line["dash"] = DASHES[j % len(DASHES)]
            trace["line"] = line
            trace["legendgroup"] = label
            traces.append(trace)

    if not traces or layout is None:
        continue

    plain_title(layout)
    layout["showlegend"] = True
    layout["height"] = 450
    layout["margin"] = {"l": 70, "r": 30, "t": 50, "b": 60}

    sid = section + "_ws_" + suffix
    div = "div_" + sid
    with open(sid + "_mqc.html", "w") as fh:
        fh.write("<!--\n")
        fh.write("id: '" + sid + "'\n")
        fh.write("section_name: 'Cell Ranger ATAC: " + title + "'\n")
        fh.write("description: 'Cell Ranger own figure from each web_summary.html, every library on one axis.'\n")
        fh.write("-->\n")
        fh.write('<div id="' + div + '" style="width:100%;height:450px"></div>\n')
        fh.write("<script>\n(function () {\n")
        fh.write("  var figure = " + json.dumps({"data": traces, "layout": layout}) + ";\n")
        fh.write("  function draw() {\n")
        fh.write('    var el = document.getElementById("' + div + '");\n')
        fh.write("    if (!el || !window.Plotly) { return setTimeout(draw, 100); }\n")
        fh.write("    window.Plotly.newPlot(el, figure.data, figure.layout, "
                 "{responsive: true, displaylogo: false});\n")
        fh.write('    el.setAttribute("data-rendered", "yes");\n')
        fh.write("  }\n")
        fh.write('  if (document.readyState === "loading") {\n')
        fh.write("    document.addEventListener(\"DOMContentLoaded\", draw);\n")
        fh.write("  } else { draw(); }\n")
        fh.write("})();\n</script>\n")
    print("wrote " + sid + "_mqc.html")
PYEOF
    """

    stub:
    """
    echo '{"id":"stub","section_name":"stub","plot_type":"table","data":{"a":{"b":1}}}' > ${section}_ws_stub_mqc.json
    """
}
