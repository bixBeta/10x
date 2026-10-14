#!/usr/bin/env python3
"""Catch backslashes that Groovy eats before the task script ever runs.

A process script lives in a Groovy triple quoted string, which processes
escapes itself: a lone backslash n becomes a real newline, and a lone trailing
backslash is swallowed. So anything the shell or an embedded python block needs
to see as a backslash has to be written doubled.

That is why every line continuation in these modules is a double backslash, and
why an embedded python writing a newline has to say it doubled too. Getting it
wrong does not fail at parse time - it fails inside the task, which on a real
run means after the counting.

    python3 test/check_module_escapes.py modules/*/main.nf
"""
import re
import sys

# """ ... """ - the script and stub blocks
BLOCK = re.compile(r'"""(.*?)"""', re.S)

# Escapes groovy is MEANT to process, and which therefore survive correctly:
#   \$  stops interpolation, the task sees a plain dollar
#   \\  the task sees one backslash - how every line continuation and every
#       regex backslash in these modules is written
#   \"  \'  quoting
# Anything else is an escape the author expected the shell or an embedded
# python to see, and groovy consumed it first.
INTENDED = set('$\\"\'')


def check(path):
    with open(path) as fh:
        text = fh.read()

    problems = []
    for block in BLOCK.finditer(text):
        body = block.group(1)
        start_line = text[: block.start()].count("\n") + 1
        for offset, line in enumerate(body.split("\n")):
            i = 0
            while i < len(line):
                if line[i] != "\\":
                    i += 1
                    continue
                nxt = line[i + 1] if i + 1 < len(line) else ""
                if nxt in INTENDED:
                    i += 2          # groovy handles it, the task sees what it should
                    continue
                problems.append((start_line + offset, line.strip(), nxt or "end of line"))
                break

    if problems:
        print(f"{path}: backslashes groovy will eat before the task sees them:")
        for lineno, line, nxt in problems:
            print(f"  line {lineno}: backslash before {nxt!r} - write it doubled")
            print(f"           {line[:100]}")
        return 1
    return 0


def main(paths):
    if not paths:
        raise SystemExit(__doc__)
    bad = sum(check(p) for p in paths)
    if not bad:
        print(f"module escapes OK ({len(paths)} files)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
