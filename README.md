# Nextflow Pipeline for 10x Genomics single cell runs

[![ci](https://github.com/bixBeta/10x/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/bixBeta/10x/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-25.04.1-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![cellranger](https://img.shields.io/badge/cellranger-8.0.1%20%7C%209.0.1%20%7C%2010.1.0-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)
[![cellranger-arc 2.2.0](https://img.shields.io/badge/cellranger--arc-2.2.0-6f42c1.svg)](https://www.10xgenomics.com/support/software/cell-ranger-arc)
[![cellranger-atac 2.2.0](https://img.shields.io/badge/cellranger--atac-2.2.0-e8710a.svg)](https://www.10xgenomics.com/support/software/cell-ranger-atac)
[![spaceranger 4.0.1](https://img.shields.io/badge/spaceranger-4.0.1-0b7285.svg)](https://www.10xgenomics.com/support/software/space-ranger)
[![multiqc 1.35](https://img.shields.io/badge/multiqc-1.35-1a9e77.svg)](https://multiqc.info/)
[![modes](https://img.shields.io/badge/modes-gex%20%7C%20atac%20%7C%20arc%20%7C%20multi%20%7C%20spatial-555.svg)](#)

<hr>

| mode | tool | status |
|---|---|---|
| `gex` | `cellranger count` | working |
| `arc` | `cellranger-arc count` | working |
| `atac` | `cellranger-atac count` | working |
| `multi` | `cellranger multi` | working — GEX + CITE-seq + VDJ |
| `spatial` | `spaceranger count` | working — Visium / Visium HD |
| — | `cellranger aggr` | working, for multi-library labels (gex) |
| — | `multiqc` | working, one report per run |

Requires Nextflow >= 24.04, tested against 25.04.1.

## Usage

Grab [`params.yaml`](params.yaml) next to your sample sheet:

```bash
wget https://raw.githubusercontent.com/bixBeta/10x/main/params.yaml
```

Set `ref` and anything else you need, then:

```bash
nextflow run bixBeta/10x -r main -params-file params.yaml
```

It carries every default with comments, so a line only needs changing if you
want something else. Note the single dash: `-params-file` is a Nextflow option,
not a pipeline param.

The command line still wins over the file, which is handy for one-off tweaks:

```bash
nextflow run bixBeta/10x -r main -params-file params.yaml --crversion 10.1.0
```

Everything can be passed directly instead, if you prefer:

```bash
nextflow run bixBeta/10x -r main --id BRC_1234 --sheet sample-sheet.csv --ref Human
```

```bash
nextflow run bixBeta/10x -r main --help
```

```bash
nextflow run bixBeta/10x -r main --listRefs
```

The run directory gathers the summaries so a whole run can be skimmed at a
glance:

```
CELLRANGER/JS4/outs   -> work/../outs        the full cellranger output
filtered_counts/JS4/  filtered_feature_bc_matrix.h5
web_summary_htmls/    JS4_web_summary.html   JS5_web_summary.html    ...
summary_metrics/      JS4_metrics_summary.csv  JS5_metrics_summary.csv ...
                      ALL_metrics_summary.csv
multiqc/              BRC_1234_multiqc_report.html
pipeline_info/        software_versions.yml
```

`ALL_metrics_summary.csv` is every per library table stacked, header kept once,
with a leading `label` column so each row can be traced back — Cell Ranger's
own `metrics_summary.csv` carries no sample identifier.
Multiome and scATAC join the same directories. Both report `summary.csv` rather
than `metrics_summary.csv`, with their own columns, so each gets its own
`ALL_summary.csv`; `CELLRANGER_ARC/<label>/` also keeps the generated
`libraries.csv`, and `CELLRANGER_AGGR/<label>/` the aggregated result.

`--mode spatial` joins them too, under `SPACERANGER/<label>/`. Space Ranger
reports a wide `metrics_summary.csv` like `gex` does, so its rows stack into
`ALL_metrics_summary.csv` the same way.

`--mode multi` mirrors `outs` the same way under `CELLRANGER_MULTI/<label>/`,
keeps the generated config beside it, and links the per-sample results out of
`outs/per_sample_outs/` as `<label>_<sample>_web_summary.html` and
`<label>_<sample>_metrics_summary.csv` — multi writes one directory per sample,
which is one for a plain run and several when multiplexed.

With `--createBam true` it also gathers every sample's alignments into one
place:

```
bams/   EW1_SC2622_EW1_SC2622_alignments.bam(.bai)   -> work/../per_sample_outs/.../count/
        EW2_SC2622_EW2_SC2622_alignments.bam(.bai)
```

Links, not copies. The per-run `unassigned_alignments.bam` and the VDJ contig
BAMs are left where Cell Ranger put them, reachable through the `outs` mirror.
With `createBam` false there is nothing to collect and no `bams/` directory.

For scATAC, `filtered_counts/<label>/` holds what Signac and ArchR actually
read — the peak matrix, the fragments file with its index, and `singlecell.csv`:

```
CELLRANGER_ATAC/CK1/outs  -> work/../outs
filtered_counts/CK1/      filtered_peak_bc_matrix.h5
                          fragments.tsv.gz  fragments.tsv.gz.tbi
                          singlecell.csv
```

### MultiQC

One report per run, in `multiqc/<id>_multiqc_report.html`. MultiQC's
`cellranger` and `cellranger_arc` modules find runs by matching content inside
`web_summary.html`, so the summaries the pipeline already collects are the whole
input — `metrics_summary.csv` is not one of its search patterns.

The **Software Versions** table at the end of the report is the pipeline's own,
not MultiQC's. `disable_version_detection` is set in
[`assets/multiqc_config.yaml`](assets/multiqc_config.yaml), and `DUMP_VERSIONS`
writes `software_versions_mqc.yml` from the `versions.yml` every process emits,
so the table reports what the tools themselves said. `report_section_order`
pins it last. MultiQC's own version is not in that table — it cannot be, since
the table is an input to the report.

For multiome, MultiQC's own `cellranger_arc` module cannot parse Cell Ranger
ARC 2.2 web summaries ([MultiQC #3609](https://github.com/MultiQC/MultiQC/issues/3609),
open), so the pipeline renders the metrics itself from `summary.csv` as three
sections — **library**, **ATAC** and **GEX** — with the headline metrics visible
and the long tail behind *Configure columns*. Splitting by assay rather than
using MultiQC's column `namespace` is deliberate: a namespace only shows in
tooltips, which left ATAC and GEX columns of the same name indistinguishable.
The module is left enabled so its richer section returns the day upstream
fixes it.

For scATAC there is no `cellranger_atac` module in MultiQC at all, so the
pipeline parses the web summary itself. A Cell Ranger web summary carries the
whole report as one JSON blob on a `const data = {...}` line — tables, plot
traces and axis labels — which is how MultiQC's own `cellranger` module reads
the gex one. [`modules/websummary`](modules/websummary/main.nf) lifts out:

| section | from | how |
|---|---|---|
| Sequencing / Cells / Mapping / Targeting / Complexity | the report's own metric tables | MultiQC tables |
| Barcode Rank · Fragment Distribution · TSS Enrichment · Insert Size · Assay Efficiency | the report's own plots | Cell Ranger's Plotly figure, embedded |

Those tables are richer than `ALL_summary.csv` — the mapping one carries the
nucleosome-free and mononucleosome fractions, which the CSV omits entirely.
The CSV-derived table is kept alongside them.

The plots are **not redrawn**. Each one in a web summary is already a Plotly
figure (`data` + `layout`), and MultiQC bundles Plotly, so the figure is handed
to `Plotly.newPlot` as Cell Ranger built it — same axes, same log scales, same
styling. Only three things change: each trace is prefixed with its library, a
library gets one colour across every plot (dashes separate traces within it),
and the per-library subtitle is dropped since several are now overlaid.

Rebuilding them as MultiQC line graphs was tried first and looked wrong: the
styling is lost, and these profiles are noisy at 1 bp — which the web summary
hides by drawing into a small panel, and which reads as a hairy band at report
width. Passing the figure through avoids the question entirely.

This adds roughly 0.5–1.5 MB to the report, mostly the barcode-rank traces.
The tSNE cluster scatters are deliberately skipped: a per-library 2D embedding
is not comparable across libraries.

The web summaries themselves are published for a human to open but deliberately
not handed to MultiQC — nothing would parse them, and the `cellranger_arc`
module trying and failing would take the whole report down.

For `--mode multi` the per-sample web summaries are passed to MultiQC as they
are. MultiQC's `cellranger` module matches on `"subcommand":"count"`, so a
`multi` summary is not parsed by it today and the report carries the Software
Versions section rather than per-sample tables. The same web-summary parsing
used for scATAC would fix that; it is not wired up yet.

`--multiqc false` skips the report, `--multiqcconfig` swaps the config, and
`--multiqcversion` selects the image.

`CELLRANGER/<label>/outs` is a single symlink to the task's real `outs`, so
everything Cell Ranger wrote — matrices, `molecule_info.h5`, `analysis/`, the
BAM if asked for — is reachable from the run directory without copying
anything. Multiome and aggregation mirror theirs the same way under
`CELLRANGER_ARC/<label>/` and `CELLRANGER_AGGR/<label>/`.

Everything published is a link into `work/`, so keep `work/` until you are
done, or dereference when archiving (`cp -rL`, `rsync -L`, `tar -h`).

## Building the sheet

`scripts/make-sheet.py` writes one for you from the delivery. It derives what
the naming convention supports and leaves everything else **empty with a note**,
rather than guessing — a sheet that is wrong in a plausible way costs far more
than one with a blank to fill.

```bash
scripts/make-sheet.py spatial /local/Illumina/DRV/<run>/Unaligned/Project_104*** \
    --imagedir visium_images -n          # -n prints it instead of writing
```

```
sample-sheet.csv: 16 row(s), 2 label(s): YF1, YF2
  YF1        8 delivery dir(s)  slide=H1-DYCH77W  area=A1  cytaimage=CAVG..._A1_yf1.tif
  YF2        8 delivery dir(s)  slide=H1-DYCH77W  area=B1  cytaimage=CAVG..._B1_yf2.tif

not placed, under visium_images:
  GE_1_VisiumHD_20261105JZS.jpg
  GE_2_VisiumHD_20261105JZS.jpg
```

| it derives | from |
|---|---|
| `label` | field 1 of the fastq prefix — `SC2630_YF1_...` → `YF1`. `--label-field` moves it, `--rename OLD=NEW` fixes one up |
| `fastqs` | the directory the R1 files were found in, one row each. A capture area delivered as eight per-lane directories becomes eight rows with the label repeated — which is how the pipeline pools them |
| `library` | added only when one label really covers several **fastq prefixes**, i.e. separate GEM wells, which is the one case that column exists for |
| `area` | a field of the prefix shaped like a capture area: `..._V-G3_A1_...` → `A1` |
| `slide` | read off the CytAssist image file name, which is where 10x puts it: `..._H1-DYCH77W_A1_yf1.tif` → `H1-DYCH77W`. `--slide` overrides |
| `cytaimage` | an image under `--imagedir` carrying **both** a slide serial and this row's area |
| `image` / `alignment` | a further image, or a `.json`, carrying this row's label or area |

Give it the `Project_*` directory, the `Sample_*` directories, or any fastq
inside one — the same three things the sheet itself accepts. Nothing is copied
or renamed; the paths written are the directories exactly as found.

Modes: `gex`, `atac`, `spatial`. `arc` and `multi` are not there yet — both need
to know *which assay* each library is, and a wrong guess silently counts the
wrong reads together, so that wants its own pass. Adding a mode means a
subclass with its column names and a `fill()`; the discovery, labelling and
reporting are shared. See the docstring.

## Sample sheet

Most runs need two columns:

```csv
label,fastqs
JS4,/local/Illumina/DRV/260***_RX_0***/Unaligned/Project_104***/Sample_SC2620_JS4_G3_***_104***_253***_L2
JS5,/local/Illumina/DRV/260***_RX_0***/Unaligned/Project_104***/Sample_SC2620_JS5_G3_***_104***_253***_L2
```

**`label`** — your short name for the library. It becomes the output folder and
Cell Ranger's `--id`.

**`fastqs`** — the 10x delivery directory, used exactly as given. Nothing is
copied or renamed. Any fastq inside it works too, and resolves to its parent.

The remaining columns are optional, one per situation:

| column | add it when |
|---|---|
| `library` | one sample was run in two GEM wells |
| `sample` | the `fastqs` path holds more than one sample |
| `library_type` | multiome, to mark the GEX and ATAC halves |
| `slide` / `area` / `cytaimage` … | spatial; see [Visium](#visium---mode-spatial) |

### The same library sequenced twice

Repeat the label. Both runs are pooled into one count.

```csv
label,fastqs
JS5,/local/Illumina/DRV/run1/.../Sample_SC2620_JS5_G3_***_104***_253***_L2
JS5,/local/Illumina/DRV/run2/.../Sample_SC2620_JS5_G3_***_104***_999***_L3
```

### Two GEM wells from one sample

Give each its own `library`. They are counted separately, since barcodes are
not comparable across GEM wells, then combined with `cellranger aggr`.

```csv
label,library,fastqs
JS6,JS6_wellA,/local/.../Sample_SC2620_JS6A_G3_***_104***_253***_L2
JS6,JS6_wellB,/local/.../Sample_SC2620_JS6B_G3_***_104***_253***_L2
```

You get a count per library plus `CELLRANGER_AGGR/JS6/`. Labels with a single
library are left alone — a re-sequenced library was already pooled into one
count, and aggregating it would only depth-normalise two halves of the same
library against each other. `--aggr false` skips the step; `--normalize`
controls Cell Ranger's own depth normalisation (`mapped`, the default, or
`none`).

### A path with several samples in it

Pointing at a `Project_*` directory is fine, but name which sample the row
means. Without it the run stops rather than guessing.

```csv
label,fastqs,sample
JS4,/local/.../Project_104***,SC2620_JS4_G3_***_104***_253***
```

### Multiome (`--mode arc`)

One label, two rows: the GEX library and the ATAC library, marked with
`library_type`.

```csv
label,fastqs,sample,library_type
JS4,/local/.../Project_104***,SC2619_JS4_BC_MG3_***_104***_25F***,Gene Expression
JS4,/local/.../Project_104***,SC2619_JS4_MA_***_104***_23C***,Chromatin Accessibility
```

The pipeline writes the `libraries.csv` that `cellranger-arc` needs and keeps a
copy with the results. `library_type` also accepts shorthands (`gex`, `atac`,
`adt`, `citeseq`, `crispr`, `cmo`), and the full 10x vocabulary for the modes
still to come.

### scATAC (`--mode atac`)

One row per library, the same shape as `gex` — `cellranger-atac` takes
`--fastqs` / `--sample`, not a `libraries.csv`. The `library_type` column can be
left out, since Chromatin Accessibility is the only thing this mode counts:

```csv
label,fastqs
CK1,/local/.../Sample_SC2617_CK1_MA_***_104***_450***_L1
```

To count the ATAC half of a **multiome** run on its own, use this mode with an
arc reference and `--chemistry ARC-v1`:

```bash
nextflow run bixBeta/10x -r main --mode atac --ref Human --chemistry ARC-v1 \
  --sheet sample-sheet.csv
```

`--mode arc` cannot do this. `cellranger-arc count` needs both libraries in its
`libraries.csv` and has no single-assay mode, and no chemistry setting changes
that. `--chemistry` is otherwise left off entirely in this mode: unlike
`cellranger`, `cellranger-atac` detects the configuration itself and takes no
`auto` literal.

`cellranger-atac aggr` is not wired up, so several libraries under one label
are counted separately rather than aggregated.

<details>
<summary>How <code>--sample</code> is worked out</summary>

Cell Ranger needs the fastq **prefix**, which is not the directory name:

```
Sample_SC2620_JS4_G3_***_104***_253***_L2/
  SC2620_JS4_G3_***_104***_253***_S3_L002_R1_001.fastq.gz
  -> --sample=SC2620_JS4_G3_***_104***_253***
```

It is read from the actual filenames, so an unexpected naming variant is caught
while reading the sheet rather than inside Cell Ranger. The `sample` column
overrides it, and is only needed when one path holds more than one library — a
`Project_*` directory rather than a single `Sample_*` one. Give it the prefix,
not the file name: a value still carrying its `_S3_L002` tail is rejected, with
the prefix it meant named in the error.

The flow cell is part of that prefix, so the same library re-sequenced on a
different flow cell has a *different* prefix. Both are collected and passed as
`--sample=prefix1,prefix2` — passing only one would silently drop a run.

</details>

### Multi (`--mode multi`)

One label carries every library of one sample — Gene Expression plus any of
Antibody Capture (CITE-seq), CRISPR and VDJ — and `cellranger multi` counts
them together. **You supply only the sheet**; the multi config CSV is generated
per label from it.

```csv
label,fastqs,sample,lanes,library_type
EW1,/local/.../Project_104***,SC2622_EW1_G5_***,2,Gene Expression
EW1,/local/.../Project_104***,SC2622_EW1_CSP_***,2,Antibody Capture
EW1,/local/.../Project_104***,SC2622_EW1_VDJ_***,5,VDJ
```

becomes `EW1_multi_config.csv`:

```
[gene-expression]
reference,<--ref>
create-bam,<--createBam>

[vdj]
reference,<--vdjref>

[feature]
reference,<--featureref>

[libraries]
fastq_id,fastqs,lanes,feature_types
SC2622_EW1_G5_***,/local/.../Project_104***,2,Gene Expression
SC2622_EW1_CSP_***,/local/.../Project_104***,2,Antibody Capture
SC2622_EW1_VDJ_***,/local/.../Project_104***,5,VDJ
```

A copy is kept at `CELLRANGER_MULTI/<label>/<label>_multi_config.csv`, next to
the results, so it can be diffed against a hand-written one.

`lanes` is **optional and multi-only**. It restricts the library to those flow
cell lanes; leave it out and `any` is used. Several lanes are separated with a
**pipe**, since the comma is the CSV delimiter — `1|2|3` — and a range works
too, `1-4`.

Libraries are written in assay order: Gene Expression, then the feature barcode
assays, then VDJ. That is how a hand-written config reads, and keeping to it
means the generated file diffs cleanly against one.

`[vdj]` and `[feature]` are written **only when a library needs them**, since an
empty section is an error in Cell Ranger rather than a no-op. So `--vdjref` and
`--featureref` are required exactly when the sheet carries a library that needs
one, and the run stops at launch naming the label and the assay if it is
missing.

<details>
<summary>Which options reach the generated config</summary>

Only these. Everything else in the `multi` config spec is **not** expressible
yet, and a param the pipeline does not declare is silently ignored:

| section | options |
|---|---|
| `[gene-expression]` | `reference`, `create-bam`, `chemistry`, `expect-cells`, `force-cells`, `include-introns`, `r1-length`, `r2-length` |
| `[vdj]` | `reference` |
| `[feature]` | `reference` |
| `[libraries]` | `fastq_id,fastqs,lanes,feature_types` |

Not supported: `probe-set`, `target-panel`, `no-secondary`,
`check-library-compatibility`, `cmo-set`, `min-assignment-confidence`,
`inner-enrichment-primers`, and the whole `[samples]` section — so CMO/hashtag
multiplexing cannot be expressed.

Watch `--r1length`: it defaults to **28**, so leaving it out writes
`r1-length,28` into `[gene-expression]`. Set it to `0` to omit the line, which
is usually what you want for 5' GEX.

</details>

### Visium (`--mode spatial`)

One row per **capture area**. `spaceranger count` takes `--fastqs` / `--sample`
exactly like the `gex` path; what makes spatial different is that it needs two
things the reads cannot supply — **which slide and area** the section sat on,
and an **image** of the tissue to align the spots to.

```csv
label,fastqs,slide,area,cytaimage,image
YF1,/local/.../Sample_SC2630_YF1_V-G3_A1_104***_23N***_L1,H1-DYCH77W,A1,CAVG10243R1_..._A1_yf1.tif,GE_1_20260212JZS.jpg
YF2,/local/.../Sample_SC2630_YF2_V-G3_B1_104***_23N***_L1,H1-DYCH77W,B1,CAVG10243R1_..._B1_yf2.tif,GE_2_20260212JZS.jpg
```

```bash
scripts/make-sheet.py spatial /local/.../Project_104*** --imagedir visium_images
nextflow run bixBeta/10x -r dev --mode spatial --ref Canine \
  --sheet sample-sheet.csv --slide H1-DYCH77W
```

`scripts/make-sheet.py` derives the label, the area, the slide serial and the
CytAssist image from the delivery — see [Building the sheet](#building-the-sheet).

**Images live in one folder.** A bare file name in the sheet is looked up under
`--imagedir`, which defaults to `visium_images` next to where you launch — so
drop the TIFFs the lab sent in there and the sheet only has to name them. A full
path works too.

```
<project>/
  sample-sheet.csv
  params.yaml
  visium_images/
    CAVG10243R1_2026-02-12_..._H1-DYCH77W_A1_yf1.tif
    GE_1_20260212JZS.jpg
```

| column | what it is |
|---|---|
| `slide` | the serial printed on the slide, e.g. `H1-DYCH77W`. One slide normally carries a whole run, so `--slide` sets it for every row and only `area` has to be in the sheet |
| `area` | the capture area the section sits in: `A1`, `B1`, `C1`, `D1` |
| `cytaimage` | the TIFF the CytAssist instrument wrote — CytAssist assays, Visium HD included |
| `image` | a brightfield H&E image |
| `darkimage` | a fluorescence image |
| `colorizedimage` | a composite colour fluorescence image |
| `probeset` | `--probe-set`, for a probe-based (FFPE / CytAssist) assay needing a panel other than the reference's |
| `alignment` | `--loupe-alignment`, the JSON Loupe Browser exports after a manual fiducial alignment |
| `slidefile` | `--slidefile`, the slide layout `.gpr`. Space Ranger normally fetches the layout from 10x, so this is what an offline node needs |

At least one image column is **required** and the sheet is rejected without one,
along with a slide that has no area, or an image name `--imagedir` does not hold
— all before anything is submitted, since an image typo should not surface an
hour into a run.

If the slide serial was lost, put one of Space Ranger's keywords in `slide`
instead — `visium-1`, `visium-2`, `visium-2-large` — and leave `area` out. That
becomes `--unknown-slide`, which is mutually exclusive with `--slide` / `--area`.

A capture area sequenced over several lanes or runs: **repeat the label**, exactly
as in `gex` mode. Spatial deliveries often arrive one directory per lane, which
is what that looks like:

```csv
label,fastqs,area,cytaimage
YF1,/local/.../Sample_SC2630_YF1_V-G3_A1_104***_23N***_L1,A1,..._A1_yf1.tif
YF1,/local/.../Sample_SC2630_YF1_V-G3_A1_104***_23N***_L2,A1,..._A1_yf1.tif
```

The `slide`, `area` and image columns must **match** on every row of a label —
they describe one physical section — and a label whose rows disagree is rejected
rather than having one of them picked at random.

> **`--r1length` does not apply here.** A Visium read is not a gex read: a
> CytAssist library carries a 43 nt R1, not 28. Reusing the gex flag would mean
> its default silently trimmed every spatial read, so this mode has its own
> `--spatialR1length` / `--spatialR2length` and passes **neither** unless asked,
> leaving Space Ranger to read the length itself. Passing `--r1length` in
> `spatial` mode is an error rather than a no-op.

`filtered_counts/<label>/` holds what Seurat's `Load10X_Spatial` and squidpy
read — the matrix is no use without the `spatial/` folder beside it:

```
SPACERANGER/YF1/outs    -> work/../outs
filtered_counts/YF1/    filtered_feature_bc_matrix.h5
                        spatial/   tissue_positions.csv
                                   scalefactors_json.json
                                   tissue_hires_image.png ...
                        cloupe.cloupe
```

Visium HD writes no top-level matrix — every bin size gets its own directory
with its own `spatial/` folder — so `binned_outputs/` is linked there instead.

MultiQC ships no `spaceranger` module, so the web summaries are published for a
human to open but deliberately not handed over; the pipeline renders the metrics
itself as the **Space Ranger (Visium)** section, the same arrangement scATAC
uses. The knee plots and spatial figures are not lifted out of the web summary
yet — the treatment `modules/websummary` gives scATAC would work here too.

## References

`--listRefs` prints the full map. Keys point at the build rather than the
species directory, and the same key picks the right build for the mode.

`atac` uses the `arc` column: `cellranger-atac --reference` takes an ATAC or an
ARC reference. `multi` uses the `gex` column for `--ref`, plus its own
`--vdjref` (`Human`, `Mouse`) for the `[vdj]` section — VDJ references are a
different build, so they are a separate key rather than derived from `--ref`.

| key | gex | arc / atac |
|---|---|---|
| `Apoculata_stonyCoral` | `240514_fromSarahArnold/apoculata` | — |
| `Canine` | `CanFam3_1` | `arc/CanFam3_Ensembl101annot` |
| `Combo_Human_Mouse` | `refdata-gex-GRCh38_and_GRCm39-2024-A` | — |
| `Combo_Human_Mouse_2020` | `refdata-gex-GRCh38-and-mm10-2020-A` | — |
| `Feline` | `Fca126` | — |
| `Horse` | `EquCab3/ENSEMBL_annot116/spaceranger/EquCab3_ENS116` | — |
| `Horse_ENS112` | `EquCab3/ENSEMBL_annot112/spaceranger/EquCab3_ENS112` | — |
| `Human` | `refdata-gex-GRCh38-2024-A` | `refdata-cellranger-arc-GRCh38-2024-A` |
| `MAIZE` | `MAIZE_CellRanger` | — |
| `Mouse` | `refdata-gex-GRCm39-2024-A` | `refdata-cellranger-arc-GRCm39-2024-A` |
| `Nematostella` | `Nematostella2_2/Nematostella2_2_standard` | `Nematostella2_2/Nematostella2_2` |
| `Nematostella_jaNemVect1` | — | `jaNemVect1.1/NCBI/jaNemVect1_1` |

Rooted at `/local/workdir/10x_analysis/REFS`. `--ref` also takes a full path.

A reference is a directory holding `reference.json`. If `--ref` points at a
parent, the run stops and names the builds it found so you can pick one.

## 10x software

Two engines, `singularity` by default:

| engine | source | resolves to |
|---|---|---|
| `singularity` | a locally built image | `<sifdir>/cellranger-<crversion>.sif` |
| `local` | a native install | `/programs/cellranger-<crversion>/cellranger` |

Images currently built, in `--sifdir`
(`/local/workdir/10x_analysis/singularity-sifs`):

| image | select with |
|---|---|
| `cellranger-8.0.1.sif` | `--crversion 8.0.1` |
| `cellranger-9.0.1.sif` | `--crversion 9.0.1` *(default)* |
| `cellranger-10.1.0.sif` | `--crversion 10.1.0` |
| `cellranger-arc-2.2.0.sif` | `--arcversion 2.2.0` *(default)* |

`--mode spatial` wants `spaceranger-<srversion>.sif`, built the same way from
the Space Ranger tarball:

```bash
bash containers/build-sif.sh spaceranger 4.0.1 spaceranger-4.0.1.tar.gz /local/workdir/10x_analysis/singularity-sifs
```

`--mode atac` wants `cellranger-atac-<atacversion>.sif`, which has to be built
first — see below.

`--crsif` / `--arcsif` / `--atacsif` / `--srsif` name an image directly.
`--crpath` / `--arcpath` / `--atacpath` / `--srpath` name a binary and force
`local` for that tool.
`--listPrograms` shows what is installed under `--programs`.

The image or binary is checked before any work starts, and the error lists what
is actually present.

Results are not comparable across Cell Ranger majors, so a project started on
one version stays on it — `crversion` in `params.yaml` is the only thing to
change, and the guard below catches a mismatch before any counting.

### Which version am I really running?

A filename is only a claim. Before counting anything the tool is asked its own
version inside the image, and the run stops if it disagrees:

```
ERROR: asked for cellranger 9.0.1 but cellranger reports 7.0.0
       Set --crversion to match, point at another image with --crsif, or pass --checkversion false.
```

This applies to `--engine local` too. The version actually reported is recorded
in `pipeline_info/software_versions.yml`.

## Building an image

`containers/build-sif.sh` builds on the machine that will run it. Nothing is
pushed anywhere.

From a 10x tarball:

```bash
bash containers/build-sif.sh cellranger 10.1.0 cellranger-10.1.0.tar.gz /local/workdir/10x_analysis/singularity-sifs
```

```bash
bash containers/build-sif.sh cellranger-atac 2.2.0 cellranger-atac-2.2.0.tar.gz /local/workdir/10x_analysis/singularity-sifs
```

An existing image can also be converted rather than rebuilt, if you have one:

```bash
bash containers/build-sif.sh cellranger 9.0.1 docker://<registry>/cellranger:9.0.1 /local/workdir/10x_analysis/singularity-sifs
```

It writes `<tool>-<version>.sif`, exactly the name the pipeline looks for, then
runs `--version` inside the image so you see what you got. Each version is its
own image: the install is baked in and cannot be upgraded in place.

Notes for the build host:

- Needs root or `--fakeroot`; set `SIF_FAKEROOT=1`.
- Needs outbound network to bootstrap `ubuntu:22.04`. The finished image does not.
- Needs several GB of scratch. Set `APPTAINER_TMPDIR` if `/tmp` is small.
- Sites that auto-bind paths need those mount points inside the image. The
  default `/workdir /local /programs` is created during the build; adjust with
  `SIF_BINDS="/workdir /local /fs"`.
- `cellranger`, `cellranger-atac` and `spaceranger` keep their executable at the top of the
  install dir, `cellranger-arc` in `bin/`. All are put on `PATH`.

## Params

| param | default | notes |
|---|---|---|
| `--id` | `BRC_ID` | BRC Project ID |
| `--sheet` | `sample-sheet.csv` | |
| `--mode` | `gex` | `gex`, `atac`, `arc`, `multi`, `spatial` |
| `--ref` | — | key from `--listRefs`, or a full path |
| `--chemistry` | `auto` | atac: omitted unless set; `ARC-v1` for a multiome ATAC library |
| `--expectCells` / `--forceCells` | — | |
| `--createBam` | `false` | multi also collects the per-sample BAMs into `bams/` |
| `--introns` | tool default | `--include-introns` |
| `--r1length` | `28` | `--r1-length`; pass `0` to omit. multi writes it into the config |
| `--r2length` | — | `--r2-length` |
| `--vdjref` | — | multi only; VDJ reference, key or path. Required for a VDJ library |
| `--featureref` | — | multi only; the feature barcode panel csv. Required for Antibody Capture / CRISPR |
| `--peaks` | — | atac only; a BED to use instead of calling peaks |
| `--dimReduce` | tool default | atac only; `lsa`, `plsa` or `pca` |
| `--imagedir` | `visium_images` | spatial only; where a bare image name in the sheet is looked up, relative to the launch dir |
| `--slide` | — | spatial only; slide serial for every row that gives none |
| `--probeset` | — | spatial only; `--probe-set` panel csv for every row that gives none |
| `--nosecondary` | `false` | spatial only; skip clustering / t-SNE / UMAP |
| `--spatialR1length` / `--spatialR2length` | — | spatial only; `--r1-length` / `--r2-length`. **Not** `--r1length` — see below |
| `--localcores` | `32` | `--localcores` **and** the CPUs reserved |
| `--localmem` | `128` | `--localmem` in GB **and** the memory reserved |
| `--maxforks` | `1` | concurrent tasks **per process** |
| `--aggr` | `true` | run `cellranger aggr` for multi-library labels |
| `--normalize` | `mapped` | aggr depth normalisation; or `none` |
| `--engine` | `singularity` | or `local` |
| `--crversion` / `--arcversion` / `--atacversion` / `--srversion` | `9.0.1` / `2.2.0` / `2.2.0` / `4.0.1` | selects the image or install. `--crversion` drives `gex` **and** `multi` |
| `--multiqc` | `true` | build the MultiQC report |
| `--multiqcversion` | `1.35` | selects `<sifdir>/multiqc-<v>.sif` |
| `--multiqcsif` / `--multiqcpath` | — | a specific image, or a binary |
| `--multiqcconfig` | shipped | a custom MultiQC config |
| `--sifdir` | see above | where the `.sif` files live |
| `--crsif` / `--arcsif` / `--atacsif` / `--srsif` | — | a specific image |
| `--programs` | `/programs` | where native installs live |
| `--crpath` / `--arcpath` / `--atacpath` / `--srpath` | — | a specific binary; forces `--engine local` |
| `--checkversion` | `true` | verify the tool's own version |

The command built per library:

```bash
cellranger count --id=<label> \
  --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam> --r1-length=<r1length> \
  --transcriptome=<ref> \
  --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]>
```

`--localcores` and `--localmem` set the Nextflow reservation *and* the Cell
Ranger flags from one value, so the allocation and what Cell Ranger believes it
has cannot drift apart. Every process retries a failed task twice before the run
stops, and no wall time is imposed.

## Development

No Cell Ranger, no containers and no real data needed:

```bash
bash test/make_test_data.sh
nextflow run . -stub-run -c test/ci.config --sheet test/sample-sheet.csv --ref test/ref
```

CI runs the same on every push, and additionally drives the real command
against a stand-in binary to check the flags, the pooling of re-sequenced
libraries, the generated `libraries.csv`, engine selection, and every rejection
case.

## Layout

```
main.nf                            params, help, reference maps, workflows
modules/cellranger/main.nf         cellranger count
modules/cellrangeraggr/main.nf     cellranger aggr
modules/cellrangerarc/main.nf      cellranger-arc count
modules/cellrangeratac/main.nf     cellranger-atac count
modules/cellrangermulti/main.nf    cellranger multi, incl. the config it generates
modules/spaceranger/main.nf        spaceranger count
modules/websummary/main.nf         web_summary.html -> multiqc custom content
modules/versions/main.nf           software_versions.yml
nextflow.config                    engines, resources, retries
containers/build-sif.sh            build a .sif locally
scripts/make-sheet.py              generate a sample-sheet.csv from a delivery
test/                              stub inputs and a stand-in binary
```

## Licence

Pipeline code: MIT. Cell Ranger and the 10x reference packages are covered by
the 10x Genomics End User Licence Agreement and are not distributed here.
