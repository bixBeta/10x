#!/usr/bin/env python3
"""Build a sample-sheet.csv for a 10x delivery.

    scripts/mk-sheet-10x.py spatial /local/Illumina/DRV/<run>/Unaligned/Project_104*** \
        --imagedir visium_images

Walks the delivery, works out one row per physical library, and derives whatever
the naming convention will support - the label, the capture area, the slide
serial, which image belongs to which area. Anything it cannot derive with
confidence is left EMPTY with a note saying so, rather than guessed at: a sheet
that is wrong in a plausible way costs far more than one with a blank to fill.

Nothing is copied, renamed or moved. The paths written are the delivery
directories exactly as found.

WHAT IS DERIVED, AND FROM WHAT

    label        field 1 of the fastq prefix, i.e. SC2630_YF1_... -> YF1
                 ( --label-field moves it, --rename fixes one up )
    fastqs       the directory the R1 files were found in
    area         a field of the prefix that looks like a capture area:
                 SC2630_YF1_V-G3_A1_... -> A1
    slide        read off the CytAssist image file name, which is where 10x
                 puts it: ..._H1-DYCH77W_A1_gingiva3.tif -> H1-DYCH77W
                 ( --slide overrides, for a run whose images are not named so )
    cytaimage    an image under --imagedir carrying BOTH a slide serial and
                 this row's capture area
    image        an image under --imagedir carrying this row's label or area
    alignment    a .json under --imagedir carrying this row's label or area

    Deliveries are one row each, so a capture area delivered as eight per lane
    directories becomes eight rows with the label repeated - which is exactly
    how the pipeline pools them.

ADDING A MODE

Every mode shares the discovery above and differs only in its columns, so that
is all a new one has to say. Subclass Mode, set:

    name        the --mode value
    required    columns always written, even when empty - the ones a human has
                to fill in if this script could not
    optional    columns written only when at least one row has a value

implement fill() to put values in a row, and add_args() if the mode needs flags
of its own. Then @register it. arc and multi are not here yet: both need to know
WHICH ASSAY each library is, which means reading it off a naming token
( _MA_ for ATAC, _CSP_ for Antibody Capture, _VDJ_ ... ) - doable, but a wrong
guess there silently counts the wrong reads together, so it wants its own pass.
"""
import argparse
import csv
import os
import re
import sys
from pathlib import Path

# <prefix>_S<n>_L<nnn>_R1_001.fastq.gz - the same shape main.nf reads --sample
# off, so a delivery this script accepts is one the pipeline accepts.
R1_RE    = re.compile(r"^(?P<prefix>.+)_S\d+_L(?P<lane>\d+)_R1(?:_001)?\.f(?:ast)?q\.gz$")

# A1 .. D1: the four capture areas of a Visium slide
AREA_RE  = re.compile(r"^[A-D]1$")

# H1-DYCH77W ( CytAssist ) or V19L29-035 ( Visium v1 / v2 )
SLIDE_RE = re.compile(r"^(?:[A-Z]\d-[A-Z0-9]{7}|V\d{2}[A-Z]\d{2}-\d{3})$")

# the house project id that sits in front of the label
PROJECT_RE = re.compile(r"^SC\d+$", re.I)

IMAGE_EXT = {".tif", ".tiff", ".jpg", ".jpeg", ".png", ".btf", ".qptiff"}


class Delivery(object):
    """One physical library: a directory, the fastq prefix inside it, its lanes."""

    def __init__(self, directory, prefix, lanes):
        self.directory = directory
        self.prefix    = prefix
        self.lanes     = sorted(lanes)
        self.label     = None
        self.derived   = None    # the label before --rename, i.e. per library
        self.area      = None

    def __repr__(self):
        return "Delivery(%s, %s)" % (self.prefix, self.directory)


def tokens(name):
    """Lower case fields of a file or directory name, for matching."""
    stem = name
    for ext in (".gz", ".tif", ".tiff", ".jpg", ".jpeg", ".png", ".json",
                ".csv", ".gpr", ".btf", ".qptiff"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    return [t for t in re.split(r"[^A-Za-z0-9-]+", stem) if t]


# --------------------------------------------------------------------------
# discovery
# --------------------------------------------------------------------------

def prefixes_in(directory):
    """prefix -> lanes, for the R1 files sitting directly in this directory."""
    found = {}
    try:
        names = sorted(os.listdir(str(directory)))
    except OSError:
        return found
    for name in names:
        m = R1_RE.match(name)
        if m:
            found.setdefault(m.group("prefix"), set()).add(int(m.group("lane")))
    return found


def find_deliveries(paths, note):
    """Deliveries under the given paths, looking one level down as well.

    A path may be the Sample_* delivery directory, the Project_* directory above
    it, or any fastq inside one - the same three things the sheet itself accepts.
    """
    out = []
    seen = set()

    for given in paths:
        p = Path(given).expanduser()
        if not p.exists():
            sys.exit("no such path: %s" % given)

        # a file stands for its parent, since --fastqs takes a directory
        if p.is_file():
            p = p.parent

        candidates = [p] + sorted(c for c in p.iterdir() if c.is_dir())

        hits = 0
        for d in candidates:
            found = prefixes_in(d)
            for prefix, lanes in sorted(found.items()):
                key = (str(d.resolve()), prefix)
                if key in seen:
                    continue
                seen.add(key)
                out.append(Delivery(d.resolve(), prefix, lanes))
                hits += 1

        if not hits:
            note("nothing found under %s - no *_S<n>_L<n>_R1_001.fastq.gz "
                 "there or one level down" % p)

    if not out:
        sys.exit("no 10x deliveries found. Point this at the Project_* or "
                 "Sample_* directory the sequencing core delivered.")
    return out


def assign_labels(deliveries, field, renames, note):
    for d in deliveries:
        fields = d.prefix.split("_")

        if len(fields) > field:
            if field == 1 and not PROJECT_RE.match(fields[0]):
                note("%s does not start with an SC<number> project id, so "
                     "field %d may not be the label. Check it, or use "
                     "--label-field." % (d.prefix, field))
            label = fields[field]
        else:
            label = fields[0]
            note("%s has no field %d; using '%s' as the label."
                 % (d.prefix, field, label))

        d.derived = label
        d.label   = renames.get(label, label)

        for f in fields:
            if AREA_RE.match(f):
                d.area = f
                break


# --------------------------------------------------------------------------
# images
# --------------------------------------------------------------------------

class ImageDir(object):
    """The files sitting in --imagedir, indexed by the tokens of their names."""

    def __init__(self, directory):
        self.directory = directory
        self.files = []
        self.used = set()
        if directory is None:
            return
        if not directory.is_dir():
            sys.exit("--imagedir is not a directory: %s" % directory)
        for f in sorted(directory.iterdir()):
            if f.is_file() and not f.name.startswith("."):
                self.files.append(f)

    def _slide_in(self, f):
        for t in tokens(f.name):
            if SLIDE_RE.match(t.upper()):
                return t.upper()
        return None

    def _matches(self, f, *wanted):
        low = [t.lower() for t in tokens(f.name)]
        return any(w and w.lower() in low for w in wanted)

    def cytaimage(self, area):
        """A tif carrying a slide serial AND this capture area: high confidence."""
        hits = [f for f in self.files
                if f.suffix.lower() in IMAGE_EXT
                and self._slide_in(f)
                and area and self._matches(f, area)]
        return hits

    def by_name(self, exts, *wanted):
        return [f for f in self.files
                if f.suffix.lower() in exts and self._matches(f, *wanted)]

    def take(self, f):
        self.used.add(f.name)
        return f.name

    def slide_of(self, f):
        return self._slide_in(f)

    def unused(self):
        return [f.name for f in self.files if f.name not in self.used]


# --------------------------------------------------------------------------
# modes
# --------------------------------------------------------------------------

MODES = {}


def register(cls):
    MODES[cls.name] = cls()
    return cls


class Mode(object):
    name     = ""
    required = ()
    optional = ()

    def add_args(self, parser):
        pass

    def fill(self, row, delivery, args, images, note):
        raise NotImplementedError


@register
class Gex(Mode):
    name     = "gex"
    required = ("label", "fastqs")
    # library is for separate GEM wells of one sample. It is filled in only when
    # one label really does cover several libraries - see library_column().
    optional = ("library",)

    def fill(self, row, d, args, images, note):
        pass


@register
class Atac(Mode):
    name     = "atac"
    required = ("label", "fastqs", "library_type")
    optional = ("library",)

    def fill(self, row, d, args, images, note):
        row["library_type"] = "Chromatin Accessibility"


@register
class Spatial(Mode):
    name     = "spatial"
    # slide / area / cytaimage are required even when empty: the pipeline
    # refuses a sheet without them, so an empty cell is the prompt to fill one.
    required = ("label", "fastqs", "slide", "area", "cytaimage")
    optional = ("image", "darkimage", "colorizedimage",
                "probeset", "alignment", "slidefile")

    def add_args(self, parser):
        parser.add_argument("--imagedir", metavar="DIR", default="visium_images",
                            help="where the microscope images live; the sheet "
                                 "names them and --imagedir tells the pipeline "
                                 "where to look (default: visium_images)")
        parser.add_argument("--slide", metavar="SERIAL",
                            help="slide serial for every row, when the image "
                                 "file names do not carry it")

    def fill(self, row, d, args, images, note):
        if d.area:
            row["area"] = d.area
        else:
            note("%s: no capture area in '%s' - fill in area by hand (A1..D1)."
                 % (d.label, d.prefix))

        hits = images.cytaimage(d.area) if d.area else []
        if len(hits) == 1:
            row["cytaimage"] = images.take(hits[0])
            serial = images.slide_of(hits[0])
            if args.slide:
                row["slide"] = args.slide
                if serial and serial != args.slide.upper():
                    note("%s: --slide says %s but %s carries %s. Using --slide."
                         % (d.label, args.slide, hits[0].name, serial))
            elif serial:
                row["slide"] = serial
        elif len(hits) > 1:
            note("%s: %d images carry area %s (%s) - cannot choose, fill "
                 "cytaimage in by hand."
                 % (d.label, len(hits), d.area, ", ".join(h.name for h in hits)))
        elif d.area:
            note("%s: no image under %s carries a slide serial and area %s - "
                 "fill cytaimage in by hand."
                 % (d.label, images.directory, d.area))

        if not row.get("slide"):
            row["slide"] = args.slide or ""
            if not row["slide"]:
                note("%s: no slide serial found. It is printed on the slide, "
                     "e.g. H1-DYCH77W; pass --slide or fill the column in."
                     % d.label)

        # a brightfield / fluorescence image is often named after the section
        # rather than the slide, so it is matched on the label or the area
        extra = [f for f in images.by_name(IMAGE_EXT, d.label, d.area)
                 if f.name != row.get("cytaimage")]
        if len(extra) == 1:
            row["image"] = images.take(extra[0])
        elif len(extra) > 1:
            note("%s: %d further images match (%s) - put the right one in "
                 "image / darkimage / colorizedimage by hand."
                 % (d.label, len(extra), ", ".join(f.name for f in extra)))

        align = images.by_name({".json"}, d.label, d.area)
        if len(align) == 1:
            row["alignment"] = images.take(align[0])
        elif len(align) > 1:
            note("%s: %d alignment json files match (%s) - pick one by hand."
                 % (d.label, len(align), ", ".join(f.name for f in align)))

        gpr = images.by_name({".gpr"}, d.label, d.area)
        if len(gpr) == 1:
            row["slidefile"] = images.take(gpr[0])


# --------------------------------------------------------------------------

def library_column(deliveries, note):
    """Which labels need a library column, i.e. cover several GEM wells.

    Repeating a label means "the same library sequenced again", and the fastq
    prefix is what tells the two apart: one prefix over several directories is
    one library across lanes or flow cells, while two prefixes under one label
    are two GEM wells, whose barcodes are not comparable and which therefore
    have to be counted separately. That is exactly what the library column says.
    """
    per_label = {}
    for d in deliveries:
        per_label.setdefault(d.label, set()).add(d.derived)

    needs = set(label for label, kinds in per_label.items() if len(kinds) > 1)
    for label in sorted(needs):
        note("%s covers %d libraries (%s) - they are separate GEM wells, so a "
             "library column was added and each is counted on its own."
             % (label, len(per_label[label]), ", ".join(sorted(per_label[label]))))
    return needs


def build(mode, deliveries, args, images, note):
    columns = list(mode.required) + list(mode.optional)
    rows = []

    wells = library_column(deliveries, note) if "library" in columns else set()

    for d in sorted(deliveries, key=lambda x: (x.label, x.derived, str(x.directory))):
        row = dict((c, "") for c in columns)
        row["label"]  = d.label
        row["fastqs"] = str(d.directory)
        if d.label in wells:
            row["library"] = d.derived
        mode.fill(row, d, args, images, note)
        rows.append(row)

    # a capture area is one physical thing, so the pipeline rejects a label
    # whose rows disagree about it. Say so here instead.
    for col in columns:
        if col in ("fastqs", "library"):
            continue
        per_label = {}
        for row in rows:
            per_label.setdefault(row["label"], set()).add(row[col])
        for label, values in sorted(per_label.items()):
            if len(values) > 1:
                note("%s: rows disagree about %s (%s). The pipeline rejects "
                     "that - every row of a label describes one library."
                     % (label, col, ", ".join(sorted(v or "<empty>" for v in values))))

    # keep the required columns whatever happens; drop optional ones nobody used
    keep = [c for c in columns
            if c in mode.required or any(r[c] for r in rows)]
    return keep, [dict((c, r[c]) for c in keep) for r in rows]


def main():
    # Shared flags live on a parent parser and each mode gets a subparser, so a
    # mode's own flags are its own: gex does not inherit spatial's --imagedir,
    # and `mk-sheet-10x.py spatial --help` lists only what spatial takes.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("paths", nargs="+", metavar="PATH",
                        help="the Project_* or Sample_* delivery directories")
    common.add_argument("-o", "--out", default="sample-sheet.csv",
                        help="where to write it (default: sample-sheet.csv)")
    common.add_argument("-n", "--dry-run", action="store_true",
                        help="print the sheet and write nothing")
    common.add_argument("-f", "--force", action="store_true",
                        help="overwrite an existing --out")
    common.add_argument("--label-field", type=int, default=1, metavar="N",
                        help="which _ separated field of the fastq prefix is "
                             "the label (default: 1, i.e. SC2630_YF1_... -> YF1)")
    common.add_argument("--rename", action="append", default=[],
                        metavar="OLD=NEW",
                        help="rename one derived label; repeatable")

    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        epilog="See the module docstring for what is derived and how to add a mode.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    modes = parser.add_subparsers(dest="mode", metavar="MODE")
    modes.required = True

    # adding a mode never means touching this function
    for name in sorted(MODES):
        sub = modes.add_parser(name, parents=[common],
                               help="a sheet for --mode %s" % name)
        MODES[name].add_args(sub)

    args = parser.parse_args()
    mode = MODES[args.mode]

    renames = {}
    for pair in args.rename:
        if "=" not in pair:
            sys.exit("--rename wants OLD=NEW, got: %s" % pair)
        old, new = pair.split("=", 1)
        renames[old] = new

    # a note is about a label or a file, but fill() runs per ROW - a capture area
    # delivered over eight lanes would otherwise repeat its note eight times
    notes = []
    def note(message):
        if message not in notes:
            notes.append(message)

    deliveries = find_deliveries(args.paths, note)
    assign_labels(deliveries, args.label_field, renames, note)

    imagedir = getattr(args, "imagedir", None)
    images = ImageDir(Path(imagedir).expanduser() if imagedir else None)

    columns, rows = build(mode, deliveries, args, images, note)

    out = sys.stdout if args.dry_run else None
    if out is None:
        target = Path(args.out)
        if target.exists() and not args.force:
            sys.exit("%s exists - move it, choose --out, or pass --force." % target)

    text = []
    text.append(",".join(columns))
    for row in rows:
        text.append(",".join(row[c] for c in columns))
    sheet = "\n".join(text) + "\n"

    if args.dry_run:
        sys.stdout.write(sheet)
    else:
        with open(args.out, "w") as fh:
            fh.write(sheet)

    labels = []
    for row in rows:
        if row["label"] not in labels:
            labels.append(row["label"])

    report = sys.stderr
    report.write("\n%s: %d row(s), %d label(s): %s\n"
                 % (args.out if not args.dry_run else "(dry run)",
                    len(rows), len(labels), ", ".join(labels)))

    for label in labels:
        mine = [r for r in rows if r["label"] == label]
        report.write("  %-10s %d delivery dir(s)" % (label, len(mine)))
        for col in columns:
            # library is per row by design, and its own note explains it
            if col in ("label", "fastqs", "library"):
                continue
            value = mine[0][col]
            if value:
                report.write("  %s=%s" % (col, value))
        report.write("\n")

    if images.directory is not None:
        left = images.unused()
        if left:
            report.write("\nnot placed, under %s:\n" % images.directory)
            for name in left:
                report.write("  %s\n" % name)
            report.write("  ( a run wide probe set goes in --probeset rather "
                         "than the sheet )\n")

    if notes:
        report.write("\ncheck these:\n")
        for message in notes:
            report.write("  - %s\n" % message)

    # only the required columns: an optional one is sparse on purpose, and an
    # empty library cell just means "this label is one library", which is right
    blanks = sorted(set(c for r in rows for c in mode.required if not r[c]))
    if blanks:
        report.write("\nempty column(s) the pipeline needs filled in: %s\n"
                     % ", ".join(blanks))

    report.write("\nNothing was copied or renamed. Read the sheet before "
                 "running it.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
