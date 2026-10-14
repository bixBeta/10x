#!/usr/bin/env python3
"""Assert the shape of the web summary sections modules/websummary renders.

check_mqc_yaml.py proves the files are valid MultiQC custom content. This
proves they say the right thing: one row per library, the metrics the csv does
not carry, and no html left in the titles the web summary wraps them in.

    python3 test/check_ws_sections.py summary_metrics CK1 CK2
"""
import glob
import json
import os
import sys

# Each of these is a plotly figure passed through from the web summary, one
# section per plot, with every library on one axis.
PLOTS = ["knee", "fragments", "tss", "insert", "efficiency"]

WANT = {
    "sequencing": ["Sequenced read pairs", "Valid barcodes"],
    "cells":      ["Estimated number of cells", "Median high-quality fragments per cell"],
    # neither of these is in summary.csv, which is why the tables are worth having
    "mapping":    ["Fragments in nucleosome-free regions",
                   "Fragments flanking a single nucleosome"],
    "targeting":  ["Number of peaks", "TSS enrichment score"],
    "complexity": ["Percent duplicates"],
}


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    outdir, labels = sys.argv[1], sorted(sys.argv[2:])
    problems, seen = [], 0

    for part, needed in WANT.items():
        path = os.path.join(outdir, f"cellranger_atac_metrics_ws_{part}_mqc.json")
        if not os.path.exists(path):
            problems.append(f"{part}: not rendered")
            continue
        with open(path) as fh:
            doc = json.load(fh)

        if sorted(doc["data"]) != labels:
            problems.append(f"{part}: rows are {sorted(doc['data'])}, expected {labels}")

        cols = set(next(iter(doc["data"].values())))
        for col in needed:
            if col not in cols:
                problems.append(f"{part}: missing column {col!r}")
        seen += len(cols)

        if "<" in json.dumps(doc["headers"]) or "<" in doc["section_name"]:
            problems.append(f"{part}: html left in the titles")

    # the plots are cellranger's own plotly figures, embedded as html custom
    # content - multiqc bundles plotly, so the figure renders as it stands
    for part in PLOTS:
        path = os.path.join(outdir, f"cellranger_atac_metrics_ws_{part}_mqc.html")
        if not os.path.exists(path):
            problems.append(f"{part}: figure not rendered")
            continue
        with open(path) as fh:
            body = fh.read()

        for needed in ("Plotly.newPlot", "data-rendered", f"id: 'cellranger_atac_metrics_ws_{part}'"):
            if needed not in body:
                problems.append(f"{part}: html is missing {needed!r}")

        fig = None
        for line in body.splitlines():
            line = line.strip()
            if line.startswith("var figure ="):
                try:
                    fig = json.loads(line[len("var figure ="):].strip().rstrip(";"))
                except ValueError as exc:
                    problems.append(f"{part}: figure json does not parse: {exc}")
                break
        if fig is None:
            problems.append(f"{part}: no figure json found")
            continue

        names = [t.get("name", "") for t in fig["data"]]
        missing = [lab for lab in labels if not any(n.startswith(lab) for n in names)]
        if missing:
            problems.append(f"{part}: no trace for {missing}, names are {names}")

        # every trace needs points, or the figure draws empty
        empty = [n for t, n in zip(fig["data"], names) if not t.get("x")]
        if empty:
            problems.append(f"{part}: traces with no points: {empty}")

        # a library keeps one colour across every plot, dashes separate the
        # traces within it
        colours = {}
        for t, n in zip(fig["data"], names):
            lab = n.split(" ")[0]
            colours.setdefault(lab, set()).add((t.get("line") or {}).get("color"))
        multi = {k: v for k, v in colours.items() if len(v) > 1}
        if multi:
            problems.append(f"{part}: a library has more than one colour: {multi}")

        # the per library subtitle has to go once several are overlaid
        title = fig["layout"].get("title")
        text = title.get("text") if isinstance(title, dict) else title
        if text and ("<" in text or any(lab in text for lab in labels)):
            problems.append(f"{part}: title still names one library: {text!r}")

    if problems:
        print("web summary sections are wrong:")
        for p in problems:
            print("  " + p)
        return 1

    print(f"web summary sections OK ({len(WANT)} tables, {seen} columns, "
          f"{len(PLOTS)} embedded figures, {len(labels)} libraries)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
