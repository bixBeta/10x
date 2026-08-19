# Nextflow Pipeline for 10x Genomics single cell runs

[![ci](https://github.com/bixBeta/10x/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/bixBeta/10x/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-25.04.1-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![cellranger 10.1.0](https://img.shields.io/badge/cellranger-10.1.0-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)
[![cellranger 9.0.1](https://img.shields.io/badge/cellranger-9.0.1-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)
[![cellranger-arc 2.2.0](https://img.shields.io/badge/cellranger--arc-2.2.0-6f42c1.svg)](https://www.10xgenomics.com/support/software/cell-ranger-arc)

<hr>

| mode | tool | status |
|---|---|---|
| `gex` | `cellranger count` | working |
| `arc` | `cellranger-arc count` | working |
| — | `cellranger aggr` | working, for multi-library labels |
| `atac` | `cellranger-atac count` | planned |

Requires Nextflow >= 24.04, tested against 25.04.1.

## Usage

```bash
nextflow run bixBeta/10x -r main --id BRC_1234 --sheet sample-sheet.csv --ref Human
```

```bash
nextflow run bixBeta/10x -r main --help
```

```bash
nextflow run bixBeta/10x -r main --listRefs
```

Results land in `CELLRANGER/<label>/outs`, or `CELLRANGER_ARC/<label>/` for
multiome, plus `pipeline_info/`. Everything is published as symlinks into the
work directory, so nothing is duplicated — keep `work/` until you are done.

Two flat collections gather the per library summaries so a whole run can be
skimmed without walking each folder:

```
web_summary_htmls/    JS4_web_summary.html  JS5_web_summary.html  ...
summary_metrics/      JS4_metrics_summary.csv  JS5_metrics_summary.csv  ...
```

Multiome joins the same directories; its metrics file is `summary.csv`, so it
lands as `<label>_summary.csv`.

## Sample sheet

Most runs need two columns:

```csv
label,fastqs
JS4,/local/Illumina/DRV/260810_RX_0556_253GGLLT4/Unaligned/Project_10488923/Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2
JS5,/local/Illumina/DRV/260810_RX_0556_253GGLLT4/Unaligned/Project_10488923/Sample_SC2620_JS5_G3_Sofia_10488923_253GGLLT4_L2
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
JS5,/local/Illumina/DRV/run1/.../Sample_SC2620_JS5_G3_Sofia_10488923_253GGLLT4_L2
JS5,/local/Illumina/DRV/run2/.../Sample_SC2620_JS5_G3_Sofia_10488923_999XYZAB2_L3
```

### Two GEM wells from one sample

Give each its own `library`. They are counted separately, since barcodes are
not comparable across GEM wells, then combined with `cellranger aggr`.

```csv
label,library,fastqs
JS6,JS6_wellA,/local/.../Sample_SC2620_JS6A_G3_Scooter_10488923_253GGLLT4_L2
JS6,JS6_wellB,/local/.../Sample_SC2620_JS6B_G3_Scooter_10488923_253GGLLT4_L2
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
JS4,/local/.../Project_10488923,SC2620_JS4_G3_Reign_10488923_253GGLLT4
```

### Multiome (`--mode arc`)

One label, two rows: the GEX library and the ATAC library, marked with
`library_type`.

```csv
label,fastqs,sample,library_type
JS4,/local/.../Project_10488522,SC2619_JS4_BC_MG3_8Healthy_10488522_25FWVCLT4,Gene Expression
JS4,/local/.../Project_10488522,SC2619_JS4_MA_8Healthy_10488522_23C52HLT4,Chromatin Accessibility
```

The pipeline writes the `libraries.csv` that `cellranger-arc` needs and keeps a
copy with the results. `library_type` also accepts shorthands (`gex`, `atac`,
`adt`, `citeseq`, `crispr`, `cmo`), and the full 10x vocabulary for the modes
still to come.

<details>
<summary>How <code>--sample</code> is worked out</summary>

Cell Ranger needs the fastq **prefix**, which is not the directory name:

```
Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2/
  SC2620_JS4_G3_Reign_10488923_253GGLLT4_S3_L002_R1_001.fastq.gz
  -> --sample=SC2620_JS4_G3_Reign_10488923_253GGLLT4
```

It is read from the actual filenames, so an unexpected naming variant is caught
while reading the sheet rather than inside Cell Ranger. The `sample` column
overrides it.

The flow cell is part of that prefix, so the same library re-sequenced on a
different flow cell has a *different* prefix. Both are collected and passed as
`--sample=prefix1,prefix2` — passing only one would silently drop a run.

</details>

## References

`--listRefs` prints the full map. Keys point at the build rather than the
species directory, and the same key picks the right build for the mode.

| key | gex | arc |
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

`--crsif` / `--arcsif` name an image directly. `--crpath` / `--arcpath` name a
binary and force `local` for that tool. `--listPrograms` shows what is installed
under `--programs`.

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
- `cellranger` keeps its executable at the top of the install dir,
  `cellranger-arc` in `bin/`. Both are put on `PATH`.

## Params

| param | default | notes |
|---|---|---|
| `--id` | `BRC_ID` | BRC Project ID |
| `--sheet` | `sample-sheet.csv` | |
| `--mode` | `gex` | `gex`, `arc` |
| `--ref` | — | key from `--listRefs`, or a full path |
| `--chemistry` | `auto` | |
| `--expectCells` / `--forceCells` | — | |
| `--createBam` | `false` | |
| `--introns` | tool default | `--include-introns` |
| `--r1length` | `28` | `--r1-length`; pass `0` to omit |
| `--r2length` | — | `--r2-length` |
| `--localcores` | `32` | `--localcores` **and** the CPUs reserved |
| `--localmem` | `180` | `--localmem` in GB **and** the memory reserved |
| `--maxforks` | `2` | concurrent tasks **per process** |
| `--aggr` | `true` | run `cellranger aggr` for multi-library labels |
| `--normalize` | `mapped` | aggr depth normalisation; or `none` |
| `--engine` | `singularity` | or `local` |
| `--crversion` / `--arcversion` | `9.0.1` / `2.2.0` | selects the image or install |
| `--sifdir` | see above | where the `.sif` files live |
| `--crsif` / `--arcsif` | — | a specific image |
| `--programs` | `/programs` | where native installs live |
| `--crpath` / `--arcpath` | — | a specific binary; forces `--engine local` |
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
modules/versions/main.nf           software_versions.yml
nextflow.config                    engines, resources, retries
containers/build-sif.sh            build a .sif locally
test/                              stub inputs and a stand-in binary
```

## Licence

Pipeline code: MIT. Cell Ranger and the 10x reference packages are covered by
the 10x Genomics End User Licence Agreement and are not distributed here.
