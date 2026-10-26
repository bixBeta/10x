#!/usr/bin/env python3
"""Dump the structure of a Cell Ranger web_summary.html's embedded JSON.

The counterpart to make_web_summary.py, which fabricates one: this inspects a
real one. A web summary carries the whole report as a single JSON blob on a
`const data = {...}` line, which is how MultiQC's cellranger module reads it
and how modules/websummary lifts out the atac tables and figures.

Run this before writing a parser for a mode that does not have one yet, so the
key names come from the file rather than from memory:

    python3 test/dump_websummary_structure.py web_summary_htmls/HS2_web_summary.html
    python3 test/dump_websummary_structure.py web_summary_htmls/
    python3 test/dump_websummary_structure.py 'CELLRANGER_MULTI/*/outs/per_sample_outs/*/web_summary.html'

Everything the pipeline publishes is a symlink into work/, often a chain of
two, so links are resolved and a broken one is reported as such rather than as
a missing file. Gzipped input is handled too.

Prints keys, types and plot shapes only - not metric values - so the output is
short and carries nothing but structure. Plot titles are printed because they
are what the MultiQC sections get named.
"""
import glob
import gzip
import io
import json
import os
import re
import sys

MAX_DEPTH = 6


def resolve(path):
    """The real file behind a (possibly chained) symlink, or a clear error."""
    seen = []
    cur = path
    for _ in range(40):
        seen.append(cur)
        if not os.path.islink(cur):
            break
        target = os.readlink(cur)
        if not os.path.isabs(target):
            target = os.path.join(os.path.dirname(cur), target)
        cur = os.path.normpath(target)
    else:
        raise SystemExit(f"symlink loop starting at {path}")

    if not os.path.exists(cur):
        chain = "\n     -> ".join(seen + [cur])
        raise SystemExit(
            f"broken symlink:\n     {chain}\n"
            f"\nThe target is gone - work/ was probably cleaned, or the run that\n"
            f"produced it failed. Re-run, or point this at a web_summary.html\n"
            f"inside a task directory that still exists."
        )
    return cur


def expand(args):
    """Files from paths, directories or globs, with symlinks resolved."""
    out = []
    for a in args:
        hits = glob.glob(a) if any(c in a for c in "*?[") else [a]
        if not hits:
            raise SystemExit(f"nothing matches: {a}")
        for h in hits:
            real = resolve(h)
            if os.path.isdir(real):
                found = sorted(glob.glob(os.path.join(h, "*web_summary*.htm*")))
                if not found:
                    raise SystemExit(f"no *web_summary*.html under {h}")
                out.extend((f, resolve(f)) for f in found)
            else:
                out.append((h, real))
    return out


def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), errors="replace")
    return open(path, errors="replace")


def load(path):
    with open_text(path) as fh:
        for line in fh:
            s = line.strip()
            m = re.match(r"^(?:const|var|let)\s+\w+\s*=\s*(\{.*\})\s*;?$", s)
            if m:
                try:
                    return json.loads(m.group(1))
                except json.JSONDecodeError:
                    continue
    # fall back: the longest {...} anywhere in the file
    with open_text(path) as fh:
        text = fh.read()
    best = None
    for m in re.finditer(r"\{.*\}", text):
        if best is None or len(m.group(0)) > len(best):
            best = m.group(0)
    if best:
        try:
            return json.loads(best)
        except json.JSONDecodeError:
            pass
    raise SystemExit(
        f"no embedded JSON object found in {path}\n"
        f"Is this really a Cell Ranger web_summary.html? "
        f"(size {os.path.getsize(path)} bytes)"
    )


def is_plot(node):
    return isinstance(node, dict) and "data" in node and "layout" in node


def title_of(v):
    return v.get("text") if isinstance(v, dict) else v


def describe_plot(path, node):
    layout = node.get("layout") or {}
    x = (layout.get("xaxis") or {}).get("title")
    y = (layout.get("yaxis") or {}).get("title")
    traces = node.get("data") or []
    print(f"PLOT {path}")
    print(f"     title={title_of(layout.get('title'))!r} "
          f"xlab={title_of(x)!r} ylab={title_of(y)!r} traces={len(traces)}")
    for t in traces[:4]:
        if not isinstance(t, dict):
            continue
        print(f"     trace name={t.get('name')!r} type={t.get('type')!r} "
              f"mode={t.get('mode')!r} points={len(t.get('x') or [])}")
    if len(traces) > 4:
        print(f"     ... {len(traces) - 4} more traces")


def walk(node, path="", depth=0):
    if depth > MAX_DEPTH:
        return
    if is_plot(node):
        describe_plot(path or ".", node)
        return
    if isinstance(node, dict):
        if "rows" in node and isinstance(node["rows"], list):
            rows = node["rows"]
            print(f"TABLE {path}  rows={len(rows)}")
            for r in rows:
                if isinstance(r, list) and r:
                    print(f"      - {r[0]}")
            return
        for k, v in node.items():
            child = f"{path}.{k}" if path else k
            if isinstance(v, (dict, list)):
                walk(v, child, depth + 1)
            else:
                print(f"      {child}  <{type(v).__name__}>")
    elif isinstance(node, list):
        if not node:
            print(f"      {path}  <empty list>")
            return
        print(f"LIST  {path}  len={len(node)}")
        walk(node[0], f"{path}[0]", depth + 1)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)

    for given, real in expand(sys.argv[1:]):
        print("=" * 78)
        print(f"given: {given}")
        if os.path.realpath(given) != os.path.abspath(given):
            print(f"real:  {real}")
        print("=" * 78)
        doc = load(real)
        root = doc.get("summary", doc)
        print("top level keys:", sorted(doc.keys()))
        if "summary" in doc:
            print("summary keys: ", sorted(root.keys()))
        print()
        walk(root)
        print()


if __name__ == "__main__":
    main()
