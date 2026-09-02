#!/usr/bin/env python3
"""Parse the MultiQC custom content the pipeline writes.

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

        if doc["plot_type"] != "html" or "<table" not in doc["data"]:
            print(f"{path}: expected an html table in data")
            bad = True
            continue

        print(f"{path}: valid, {doc['data'].count('<tr>')} rows")

    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
