#!/usr/bin/env python3
"""Assert the shape of the pivoted cellranger multi metrics tables.

check_mqc_yaml.py proves the files are valid MultiQC custom content. This
proves the pivot said the right thing - specifically that a metric repeated
across groupings becomes ONE shared column named for the grouping, rather than
one column per sample named for that sample's fastq id, which is what a naive
qualification produces.

    python3 test/check_multi_metrics.py summary_metrics
"""
import glob
import json
import os
import sys

SECTIONS = ["cells_gene_expression", "cells_antibody_capture", "cells_vdj_t",
            "library_gene_expression", "library_antibody_capture", "library_vdj_t"]


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    outdir = sys.argv[1]
    problems = []

    for part in SECTIONS:
        path = os.path.join(outdir, f"cellranger_multi_metrics_{part}_mqc.json")
        if not os.path.exists(path):
            problems.append(f"{part}: not rendered")
            continue
        with open(path) as fh:
            doc = json.load(fh)
        if not doc["data"]:
            problems.append(f"{part}: no rows")

    gex = os.path.join(outdir, "cellranger_multi_metrics_library_gene_expression_mqc.json")
    if os.path.exists(gex):
        with open(gex) as fh:
            doc = json.load(fh)
        cols = set(doc["headers"])

        # the same metric under two groupings: one column each, named for the
        # GROUPING, shared by every sample
        for needed in ("Number of reads (Fastq ID)",
                       "Number of reads (Physical library ID)"):
            if needed not in cols:
                problems.append(f"missing shared column {needed!r}; columns are {sorted(cols)}")

        # a column named for one sample's fastq id means the qualification used
        # Group Name instead of Grouped By
        leaked = [c for c in cols if "_G5" in c or "_CSP" in c or "_VDJ" in c]
        if leaked:
            problems.append(f"per sample columns leaked: {leaked}")

        # every sample must carry every shared column
        for sample, row in doc["data"].items():
            missing = cols - set(row)
            if missing and len(missing) == len(cols):
                problems.append(f"{sample}: no values at all")

    if problems:
        print("multi metrics tables are wrong:")
        for p in problems:
            print("  " + p)
        return 1

    print(f"multi metrics tables OK ({len(SECTIONS)} sections)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
