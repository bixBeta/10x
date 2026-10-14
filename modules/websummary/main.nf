process WEBSUMMARY_MQC {

    label 'process_low'

    // small generated json, so copied rather than linked: it should outlive work/
    publishDir "summary_metrics" , mode: "copy", overwrite: true

    input:
        tuple val(section), path(summaries)
        val maxpoints

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
    python3 - ${section} ${maxpoints} ${summaries} <<'PYEOF'
import json, math, os, sys

section = sys.argv[1]
# Target points per series. The curves are noisy at full resolution - a TSS
# profile swings a whole enrichment unit between adjacent bases - so y is
# AVERAGED within x bins rather than sampled. Sampling keeps every spike and
# draws a comb; averaging gives the smooth shape the web summary shows, which
# only looks clean there because it is drawn into a small panel. 0 keeps every
# raw point.
NBINS = int(sys.argv[2])
paths   = sys.argv[3:]

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
# log axes also decide how the binning is spaced: a knee runs over six decades,
# so equal width linear bins would put almost every point in the first one.
# The tSNE plots are deliberately left out: a 2D cluster scatter is not
# comparable between libraries, which is the whole point of aggregating.
PLOTS = [
    ("knee",       "Barcode Rank",          "barcode_knee_plot",        True,  True),
    ("fragments",  "Fragment Distribution", "fragments_per_cell_plot",  True,  True),
    ("tss",        "TSS Enrichment",        "tss_enrichment_plot",      False, False),
    ("insert",     "Insert Size",           "non_log_insert_size_plot", False, False),
    ("efficiency", "Assay Efficiency",      "assay_efficiency_plot",    False, False),
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


def series_of(xs, ys, logx):
    pairs = []
    for x, y in zip(xs, ys):
        if x is None or y is None:
            continue
        try:
            fx, fy = float(x), float(y)
        except (TypeError, ValueError):
            continue
        if fx != fx or fy != fy:        # NaN
            continue
        if logx and fx <= 0:            # cannot sit on a log axis
            continue
        pairs.append((fx, fy))

    # sorted by x: multiqc keys a line series by x, so a trace emitted in some
    # other order would otherwise be drawn as a zigzag
    pairs.sort()

    if not NBINS or len(pairs) <= NBINS:
        return pairs

    # bin on log10(x) for a log axis, so the bins are even on screen
    def pos(v):
        return math.log10(v) if logx else v

    lo, hi = pos(pairs[0][0]), pos(pairs[-1][0])
    if hi <= lo:
        return pairs
    width = (hi - lo) / NBINS

    acc = {}
    for fx, fy in pairs:
        idx = int((pos(fx) - lo) / width)
        if idx >= NBINS:
            idx = NBINS - 1
        sx, sy, n = acc.get(idx, (0.0, 0.0, 0))
        acc[idx] = (sx + fx, sy + fy, n + 1)

    # the mean x of each bin, not its midpoint, so the curve stays where the
    # data actually is
    return [(acc[i][0] / acc[i][2], acc[i][1] / acc[i][2]) for i in sorted(acc)]


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
            pairs = series_of(trace.get("x") or [], trace.get("y") or [], logaxes)
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
        "description": ("From each web_summary.html." if not NBINS else
                        "From each web_summary.html, averaged into " + str(NBINS) + " bins per series."),
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
