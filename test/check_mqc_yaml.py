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
            spaces = sorted({h.get("namespace", "") for h in headers.values()})
            hidden = [h for h in headers.values() if h.get("hidden")]
            shown = [h for h in headers.values() if not h.get("hidden")]
            if not hidden or not shown:
                # all-visible is the wide unreadable table we are avoiding;
                # all-hidden is an empty one
                print(f"{path}: {len(shown)} visible / {len(hidden)} hidden columns")
                bad = True
                continue
            print(f"{path}: valid table, {len(doc['data'])} rows, "
                  f"namespaces {spaces}, {len(shown)} visible / {len(hidden)} hidden")

        else:
            print(f"{path}: unexpected plot_type {kind}")
            bad = True

    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
