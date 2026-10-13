#!/usr/bin/env python3
"""Assert the shape of the web summary sections modules/websummary renders.

check_mqc_yaml.py proves the files are valid MultiQC custom content. This
proves they say the right thing: one series per library per trace, thinned,
log axes on the knee, titles stripped of the html the web summary wraps them
in, and the 2D cluster scatters left out.

    python3 test/check_ws_sections.py summary_metrics CK1 CK2
"""
import glob
import json
import os
import sys

MAX_POINTS = 1000


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    outdir, labels = sys.argv[1], sorted(sys.argv[2:])

    knee = os.path.join(outdir, "cellranger_atac_metrics_ws_knee_mqc.json")
    with open(knee) as fh:
        doc = json.load(fh)

    problems = []

    want = sorted(f"{lab} {trace}" for lab in labels for trace in ("Cells", "Non-cells"))
    got = sorted(doc["data"])
    if got != want:
        problems.append(f"knee series are {got}, expected {want}")

    pconfig = doc["pconfig"]
    if not (pconfig.get("xlog") and pconfig.get("ylog")):
        problems.append(f"knee should be log/log, pconfig is {pconfig}")

    if "<" in json.dumps(pconfig):
        problems.append(f"html left in the knee pconfig: {pconfig}")

    fat = {k: len(v) for k, v in doc["data"].items() if len(v) > MAX_POINTS}
    if fat:
        problems.append(f"series not thinned to {MAX_POINTS}: {fat}")

    # a single trace plot is named by the library alone, not "<label> Cells"
    tss = os.path.join(outdir, "cellranger_atac_metrics_ws_tss_mqc.json")
    with open(tss) as fh:
        tdoc = json.load(fh)
    if sorted(tdoc["data"]) != labels:
        problems.append(f"tss series are {sorted(tdoc['data'])}, expected {labels}")

    # the tSNE scatters are per library and not comparable, so they stay out
    stray = glob.glob(os.path.join(outdir, "*cluster*"))
    if stray:
        problems.append(f"cluster plots should not be imported: {stray}")

    # the richer metrics the csv does not carry
    mapping = os.path.join(outdir, "cellranger_atac_metrics_ws_mapping_mqc.json")
    with open(mapping) as fh:
        mdoc = json.load(fh)
    cols = set(next(iter(mdoc["data"].values())))
    for needed in ("Fragments in nucleosome-free regions",
                   "Fragments flanking a single nucleosome"):
        if needed not in cols:
            problems.append(f"mapping table is missing {needed!r}")

    if problems:
        print("web summary sections are wrong:")
        for p in problems:
            print("  " + p)
        return 1

    print(f"web summary sections OK ({len(doc['data'])} knee series, "
          f"{len(cols)} mapping columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
