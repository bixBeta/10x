#!/usr/bin/env python3
"""Guard the workflow file against the two ways it has actually broken.

1. A column 0 heredoc inside a "run: |" block. The block is a YAML scalar held
   open by indentation, so a terminator written flush left ends the scalar and
   the rest of the step becomes garbage YAML. This has broken ci.yml three
   times; put the script in test/ and call it instead.
2. The file no longer parsing at all, which means CI does not run rather than
   failing visibly.

    python3 test/check_ci_yaml.py .github/workflows/ci.yml
"""
import re
import sys

import yaml


def main(paths):
    if not paths:
        raise SystemExit(__doc__)
    bad = False

    for path in paths:
        with open(path) as fh:
            text = fh.read()

        try:
            yaml.safe_load(text)
        except yaml.YAMLError as exc:
            print(f"{path}: does not parse: {exc}")
            bad = True
            continue

        # inside a run: | block, any line at column 0 that is not blank has
        # ended the scalar. A heredoc terminator is the way that happens.
        in_run = False
        for lineno, line in enumerate(text.split("\n"), 1):
            if re.match(r"^\s*(-\s+)?run:\s*\|", line):
                in_run = True
                continue
            if in_run:
                if not line.strip():
                    continue
                if not line.startswith(" "):
                    # back to the step list, i.e. the block legitimately ended
                    in_run = False
                    continue
            if in_run and re.match(r"^[A-Za-z_]+$", line):
                print(f"{path}:{lineno}: heredoc terminator at column 0 inside run: |")
                print(f"  {line}")
                bad = True

        # the reliable signal: a heredoc opened inside a run block
        for m in re.finditer(r"run:\s*\|(.*?)(?=\n\s{0,8}-\s|\Z)", text, re.S):
            block = m.group(1)
            for hd in re.finditer(r"<<-?'?([A-Za-z_]+)'?", block):
                term = hd.group(1)
                if re.search(rf"^{term}$", block, re.M):
                    start = text[: m.start()].count("\n") + 1
                    print(f"{path}: heredoc <<{term} inside a run block near line {start}, "
                          f"terminated at column 0 - move the script to test/")
                    bad = True

    if not bad:
        print(f"ci yaml OK ({len(paths)} file)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
