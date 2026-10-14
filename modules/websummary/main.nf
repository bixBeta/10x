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
    // "const data = {...}" line. MultiQC's own cellranger module reads it this
    // way, but there is no cellranger-atac module, so the pipeline does.
    //
    // Only the metric TABLES are lifted out. Redrawing the curves as multiqc
    // line plots was tried and dropped: cellranger's own plotly figures carry
    // styling a redraw loses, and the per base noise the web summary hides in
    // a small panel reads as a hairy band at report width.
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
PYEOF
    """

    stub:
    """
    echo '{"id":"stub","section_name":"stub","plot_type":"table","data":{"a":{"b":1}}}' > ${section}_ws_stub_mqc.json
    """
}
