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

# nextflow.config's params.mqcpoints. Series are averaged into this many bins.
NBINS = 200

# Mean second difference of the TSS curve. The fixture carries 1 bp noise like
# the real thing: unbinned it measures about 0.19, binned about 0.06. Anything
# above this means the averaging stopped happening and the plot is a comb again.
MAX_ROUGHNESS = 0.12


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

    fat = {k: len(v) for k, v in doc["data"].items() if len(v) > NBINS}
    if fat:
        problems.append(f"series not binned to {NBINS}: {fat}")

    # multiqc keys a line series by x, so an unsorted trace draws as a zigzag
    for name, series in doc["data"].items():
        xs = [float(x) for x in series]
        if xs != sorted(xs):
            problems.append(f"{name} is not sorted by x")

    # a single trace plot is named by the library alone, not "<label> Cells"
    tss = os.path.join(outdir, "cellranger_atac_metrics_ws_tss_mqc.json")
    with open(tss) as fh:
        tdoc = json.load(fh)
    if sorted(tdoc["data"]) != labels:
        problems.append(f"tss series are {sorted(tdoc['data'])}, expected {labels}")

    # and it has to be smooth: averaging within bins is the whole point, since
    # a per base TSS profile swings a full enrichment unit between neighbours
    for name, series in tdoc["data"].items():
        ys = [series[x] for x in sorted(series, key=float)]
        if len(ys) < 3:
            problems.append(f"tss {name} has only {len(ys)} points")
            continue
        rough = [abs(ys[i - 1] - 2 * ys[i] + ys[i + 1]) for i in range(1, len(ys) - 1)]
        mean_rough = sum(rough) / len(rough)
        if mean_rough > MAX_ROUGHNESS:
            problems.append(f"tss {name} roughness {mean_rough:.3f} > {MAX_ROUGHNESS}")

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

    widest = max(len(s) for s in doc["data"].values())
    print(f"web summary sections OK ({len(doc['data'])} knee series, "
          f"{widest} points in the widest, tss roughness {mean_rough:.3f}, "
          f"{len(cols)} mapping columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
