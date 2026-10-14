# Nextflow Pipeline for 10x Genomics single cell runs

[![ci](https://github.com/bixBeta/10x/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/bixBeta/10x/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-25.04.1-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![cellranger 10.1.0](https://img.shields.io/badge/cellranger-10.1.0-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)
[![cellranger 9.0.1](https://img.shields.io/badge/cellranger-9.0.1-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)
[![cellranger-arc 2.2.0](https://img.shields.io/badge/cellranger--arc-2.2.0-6f42c1.svg)](https://www.10xgenomics.com/support/software/cell-ranger-arc)
[![cellranger-atac 2.2.0](https://img.shields.io/badge/cellranger--atac-2.2.0-e8710a.svg)](https://www.10xgenomics.com/support/software/cell-ranger-atac)
[![multiqc 1.35](https://img.shields.io/badge/multiqc-1.35-1a9e77.svg)](https://multiqc.info/)

<hr>

| mode | tool | status |
|---|---|---|
| `gex` | `cellranger count` | working |
| `arc` | `cellranger-arc count` | working |
| `atac` | `cellranger-atac count` | working |
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

`--multiqc false` skips the report, `--multiqcconfig` swaps the config, and
`--multiqcversion` selects the image.

`CELLRANGER/<label>/outs` is a single symlink to the task's real `outs`, so
everything Cell Ranger wrote — matrices, `molecule_info.h5`, `analysis/`, the
BAM if asked for — is reachable from the run directory without copying
anything. Multiome and aggregation mirror theirs the same way under
`CELLRANGER_ARC/<label>/` and `CELLRANGER_AGGR/<label>/`.

Everything published is a link into `work/`, so keep `work/` until you are
done, or dereference when archiving (`cp -rL`, `rsync -L`, `tar -h`).

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

## References

`--listRefs` prints the full map. Keys point at the build rather than the
species directory, and the same key picks the right build for the mode.

`atac` uses the `arc` column: `cellranger-atac --reference` takes an ATAC or an
ARC reference.

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
| `cellranger-9.0.1.sif` | `--crversion 9.0.1` *(default)* |
| `cellranger-10.1.0.sif` | `--crversion 10.1.0` |
| `cellranger-arc-2.2.0.sif` | `--arcversion 2.2.0` *(default)* |

`--mode atac` wants `cellranger-atac-<atacversion>.sif`, which has to be built
first — see below.

`--crsif` / `--arcsif` / `--atacsif` name an image directly. `--crpath` /
`--arcpath` / `--atacpath` name a binary and force `local` for that tool.
`--listPrograms` shows what is installed under `--programs`.

The image or binary is checked before any work starts, and the error lists what
is actually present.

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
- `cellranger` and `cellranger-atac` keep their executable at the top of the
  install dir, `cellranger-arc` in `bin/`. All are put on `PATH`.

## Params

| param | default | notes |
|---|---|---|
| `--id` | `BRC_ID` | BRC Project ID |
| `--sheet` | `sample-sheet.csv` | |
| `--mode` | `gex` | `gex`, `atac`, `arc` |
| `--ref` | — | key from `--listRefs`, or a full path |
| `--chemistry` | `auto` | atac: omitted unless set; `ARC-v1` for a multiome ATAC library |
| `--expectCells` / `--forceCells` | — | |
| `--createBam` | `false` | |
| `--introns` | tool default | `--include-introns` |
| `--r1length` | `28` | `--r1-length`; pass `0` to omit |
| `--r2length` | — | `--r2-length` |
| `--peaks` | — | atac only; a BED to use instead of calling peaks |
| `--dimReduce` | tool default | atac only; `lsa`, `plsa` or `pca` |
| `--localcores` | `32` | `--localcores` **and** the CPUs reserved |
| `--localmem` | `128` | `--localmem` in GB **and** the memory reserved |
| `--maxforks` | `1` | concurrent tasks **per process** |
| `--aggr` | `true` | run `cellranger aggr` for multi-library labels |
| `--normalize` | `mapped` | aggr depth normalisation; or `none` |
| `--engine` | `singularity` | or `local` |
| `--crversion` / `--arcversion` / `--atacversion` | `9.0.1` / `2.2.0` / `2.2.0` | selects the image or install |
| `--multiqc` | `true` | build the MultiQC report |
| `--multiqcversion` | `1.35` | selects `<sifdir>/multiqc-<v>.sif` |
| `--multiqcsif` / `--multiqcpath` | — | a specific image, or a binary |
| `--multiqcconfig` | shipped | a custom MultiQC config |
| `--mqcpoints` | `200` | bins per series in the rendered curves; `0` plots every raw point |
| `--sifdir` | see above | where the `.sif` files live |
| `--crsif` / `--arcsif` / `--atacsif` | — | a specific image |
| `--programs` | `/programs` | where native installs live |
| `--crpath` / `--arcpath` / `--atacpath` | — | a specific binary; forces `--engine local` |
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
modules/websummary/main.nf         web_summary.html -> multiqc custom content
modules/versions/main.nf           software_versions.yml
nextflow.config                    engines, resources, retries
containers/build-sif.sh            build a .sif locally
test/                              stub inputs and a stand-in binary
```

## Licence

Pipeline code: MIT. Cell Ranger and the 10x reference packages are covered by
the 10x Genomics End User Licence Agreement and are not distributed here.
