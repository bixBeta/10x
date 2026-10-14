#!/usr/bin/env python3
"""Assert the shape of the web summary tables modules/websummary renders.

check_mqc_yaml.py proves the files are valid MultiQC custom content. This
proves they say the right thing: one row per library, the metrics the csv does
not carry, and no html left in the titles the web summary wraps them in.

    python3 test/check_ws_sections.py summary_metrics CK1 CK2
"""
import glob
import json
import os
import sys

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

    # the curves are deliberately not redrawn any more
    stray = glob.glob(os.path.join(outdir, "*_ws_knee_*")) + \
            glob.glob(os.path.join(outdir, "*_ws_tss_*"))
    if stray:
        problems.append(f"line plots should not be rendered: {stray}")

    if problems:
        print("web summary sections are wrong:")
        for p in problems:
            print("  " + p)
        return 1

    print(f"web summary tables OK ({len(WANT)} sections, {seen} columns, "
          f"{len(labels)} libraries)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
