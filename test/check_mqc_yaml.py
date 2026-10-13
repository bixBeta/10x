#!/usr/bin/env python3
"""Parse the MultiQC custom content the pipeline writes.

Covers both shapes: the html Software Versions section and the grouped ARC
metrics table. JSON is valid YAML, so one loader reads either file.

Grepping these files is not enough: an unindented line silently ends a
"data: |" block scalar, so the file still contains the strings you would grep
for while being invalid YAML. MultiQC then drops the section, and with it the
rest of the custom content.
"""
import sys

import yaml

REQUIRED = ("id", "section_name", "plot_type", "data")

def main(paths):
    bad = False
    for path in paths:
        try:
            with open(path) as fh:
                doc = yaml.safe_load(fh)
        except yaml.YAMLError as exc:
            print(f"{path}: invalid YAML: {exc}")
            bad = True
            continue

        if not isinstance(doc, dict):
            print(f"{path}: not a mapping")
            bad = True
            continue

        missing = [k for k in REQUIRED if k not in doc]
        if missing:
            print(f"{path}: missing {', '.join(missing)}")
            bad = True
            continue

        kind = doc["plot_type"]
        if kind == "html":
            if "<table" not in doc["data"]:
                print(f"{path}: html section with no table")
                bad = True
                continue
            print(f"{path}: valid html, {doc['data'].count('<tr>')} rows")

        elif kind == "table":
            if not isinstance(doc["data"], dict) or not doc["data"]:
                print(f"{path}: table section with no rows")
                bad = True
                continue
            headers = doc.get("headers", {})
            missing = {c for row in doc["data"].values() for c in row} - set(headers)
            if missing:
                print(f"{path}: columns without a header: {sorted(missing)[:3]}")
                bad = True
                continue
            hidden = [h for h in headers.values() if h.get("hidden")]
            shown = [h for h in headers.values() if not h.get("hidden")]
            # a wide table with nothing hidden is the unreadable one we are
            # avoiding; a table with nothing shown is empty. A narrow table
            # needs no hiding at all.
            if not shown or (len(headers) > 12 and not hidden):
                print(f"{path}: {len(shown)} visible / {len(hidden)} hidden of {len(headers)}")
                bad = True
                continue
            print(f"{path}: valid table '{doc['section_name']}', "
                  f"{len(doc['data'])} rows, {len(shown)} visible / {len(hidden)} hidden")

        elif kind == "linegraph":
            if not isinstance(doc["data"], dict) or not doc["data"]:
                print(f"{path}: line plot with no series")
                bad = True
                continue
            # MultiQC coerces x keys with float() for a line plot. A key that
            # does not convert leaves it a categorical axis instead, silently.
            bad_x = []
            for sname, series in doc["data"].items():
                if not isinstance(series, dict) or not series:
                    bad_x.append(f"{sname}: empty")
                    continue
                for x in series:
                    try:
                        float(x)
                    except (TypeError, ValueError):
                        bad_x.append(f"{sname}: x={x!r}")
                        break
            if bad_x:
                print(f"{path}: non numeric x values: {bad_x[:3]}")
                bad = True
                continue
            widest = max(len(s) for s in doc["data"].values())
            print(f"{path}: valid line plot '{doc['section_name']}', "
                  f"{len(doc['data'])} series, up to {widest} points")

        else:
            print(f"{path}: unexpected plot_type {kind}")
            bad = True

    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
