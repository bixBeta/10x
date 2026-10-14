#!/usr/bin/env python3
"""Write a fake cellranger-atac web_summary.html for CI.

The key names, the ' Cells<br><sup>HS2 - </sup>' style titles and the leading
spaces in trace names are copied from a real cellranger-atac 2.2.0 summary, so
the parser in modules/websummary is tested against the shape it will actually
meet. Only the numbers are invented.

    python3 test/make_web_summary.py HS2 out/web_summary.html
"""
import json
import math
import os
import random
import sys


def table(rows):
    return {"rows": [[name, value] for name, value in rows]}


def trace(name, xs, ys, mode="lines", kind="scatter"):
    return {"name": name, "type": kind, "mode": mode, "x": xs, "y": ys}


def plot(title, label, xlab, ylab, traces):
    return {
        "layout": {
            "title": title + "<br><sup>" + label + " - </sup>",
            "xaxis": {"title": xlab},
            "yaxis": {"title": ylab},
        },
        "data": traces,
    }


def wobble(seed):
    """Deterministic 1 bp noise, like a real TSS or insert size profile.

    Without it the fixture is a clean analytic curve and the binning in
    modules/websummary is never actually exercised.
    """
    rng = random.Random(seed)
    return lambda: 1.0 + rng.uniform(-0.09, 0.09)


def build(label):
    # a knee-ish curve, and deliberately more points than the parser keeps so
    # the thinning is exercised
    tss_noise = wobble(label + "tss")
    ins_noise = wobble(label + "ins")

    n = 2500
    rank = list(range(1, n + 1))
    frags = [max(1, int(20000 * math.exp(-i / 400.0))) for i in rank]
    cut = 1800

    return {
        "sample": {"id": label, "description": "", "pipeline": "Cell Ranger ATAC"},
        "alarms": {"alarms": []},
        "is_barnyard": 0,
        "has_secondary": 1,
        "pipeline_info_table": table([
            ("Sample ID", label),
            ("Sample description", ""),
            ("Pipeline version", "cellranger-atac-2.2.0"),
            ("Reference path", "/local/ref"),
            ("Chemistry", "ARC-v1"),
            ("Organism", "Mouse"),
        ]),
        "sequencing_info_table": table([
            ("Sequenced read pairs", "423,110,222"),
            ("Valid barcodes", "96.4%"),
            ("Q30 bases in barcode", "94.1%"),
            ("Q30 bases in read 1", "93.2%"),
            ("Q30 bases in read 2", "92.8%"),
            ("Q30 bases in sample index i1", "91.0%"),
        ]),
        "cell_metrics_table": table([
            ("Estimated number of cells", "5,062"),
            ("Mean raw read pairs per cell", "83,585"),
            ("Fraction of high-quality fragments in cells", "63.8%"),
            ("Fraction of transposition events in peaks in cells", "51.1%"),
            ("Median high-quality fragments per cell", "14,416"),
        ]),
        "mapping_table": table([
            ("Confidently mapped read pairs", "92.3%"),
            ("Unmapped read pairs", "2.1%"),
            ("Non-nuclear read pairs", "1.4%"),
            ("Fragments in nucleosome-free regions", "48.7%"),
            ("Fragments flanking a single nucleosome", "31.2%"),
        ]),
        "targeting_table": table([
            ("Number of peaks", "74,832"),
            ("Fraction of genome in peaks", "4.0%"),
            ("TSS enrichment score", "7.1334"),
            ("Fraction of high-quality fragments overlapping TSS", "28.9%"),
            ("Fraction of high-quality fragments overlapping peaks", "45.2%"),
        ]),
        "bulk_complexity_table": table([("Percent duplicates", "39.1%")]),
        "barcode_knee_plot": plot(
            " Cells", label, "Barcodes", " Fragments Overlapping Peaks",
            [trace("Non-cells", rank[cut:], frags[cut:]),
             trace(" Cells", rank[:cut], frags[:cut])],
        ),
        "fragments_per_cell_plot": plot(
            " Fragment Distribution", label, " Fragments Per Barcode", "Barcodes",
            [trace("Non-cells", rank[cut:], frags[cut:], mode=None),
             trace(" Cells", rank[:cut], frags[:cut], mode=None)],
        ),
        "tss_enrichment_plot": plot(
            "Enrichment around TSS", label,
            "Relative Position (bp from TSS)", "Relative Enrichment",
            [trace("Cells", list(range(-1000, 1001)),
                   [(1 + 6 * math.exp(-(x / 150.0) ** 2)) * tss_noise()
                    for x in range(-1000, 1001)])],
        ),
        "non_log_insert_size_plot": plot(
            "Insert Size Distribution", label,
            "Insert Size", "Fragment Count (linear scale)",
            [trace("Cells", list(range(775)),
                   [int(1e6 * math.exp(-i / 120.0) * ins_noise()) for i in range(775)])],
        ),
        "assay_efficiency_plot": plot(
            "Median Unique Fragments per Cell", label,
            "Mean Read Pairs Per Cell", "Median Unique Fragments Per Cell",
            [trace("Observed", list(range(0, 101)), [i * 150 for i in range(0, 101)])],
        ),
        # 2D cluster scatters, which the parser must ignore
        "graph_cluster_plot": plot(
            "Cell Clustering (By Cluster)", label, "tSNE axis 1", "tSNE axis 2",
            [trace("Cluster 1 (1517)", [1, 2, 3], [1, 2, 3], mode="markers", kind="scattergl")],
        ),
        "depth_cluster_plot": plot(
            "Cell Clustering (By Depth)", label, "tSNE axis 1", "tSNE axis2",
            [trace("Cells", [1, 2, 3], [1, 2, 3], mode="markers", kind="scattergl")],
        ),
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    label, out = sys.argv[1], sys.argv[2]
    d = os.path.dirname(out)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("<html><head><script>\n")
        fh.write("const data = " + json.dumps(build(label)) + ";\n")
        fh.write("</script></head><body>fake</body></html>\n")
    print("wrote " + out)


if __name__ == "__main__":
    main()
