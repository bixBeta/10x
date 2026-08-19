# Nextflow Pipeline for 10x Genomics single cell runs

[![ci](https://github.com/bixBeta/10x/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/bixBeta/10x/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-25.04.1-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![Cell Ranger](https://img.shields.io/badge/cellranger-9.0.1-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)

<hr>

| mode | tool | status |
|---|---|---|
| `gex` | `cellranger count` | working |
| `arc` | `cellranger-arc count` | working |
| `atac` | `cellranger-atac count` | v0.2 |

## Usage

```bash
nextflow run https://github.com/bixBeta/10x -r main --id TREx_1234 --sheet sample-sheet.csv --ref GRCh38
```

```bash
nextflow run https://github.com/bixBeta/10x -r main --help
```

```bash
nextflow run https://github.com/bixBeta/10x -r main --listRefs
```

## Sample sheet

```csv
label,fastqs
JS4,/local/Illumina/DRV/260810_RX_0556_253GGLLT4/Unaligned/Project_10488923/Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2
JS5,/local/Illumina/DRV/260810_RX_0556_253GGLLT4/Unaligned/Project_10488923/Sample_SC2620_JS5_G3_Sofia_10488923_253GGLLT4_L2
```

`fastqs` is the 10x delivery directory — or any fastq inside it, whichever is
easier to paste. A file resolves to its parent, since `--fastqs` takes a
directory. **Nothing is copied, staged or renamed**; the path goes to Cell
Ranger exactly as given.

`label` is yours: short, user defined, used for `--id` and the output folder.
It is never derived from the path.

`--sample` is read off the files in the directory, because Cell Ranger needs
the fastq prefix rather than the directory name:

```
Sample_SC2620_JS4_G3_Reign_10488923_253GGLLT4_L2/
  SC2620_JS4_G3_Reign_10488923_253GGLLT4_S3_L002_R1_001.fastq.gz
  -> --sample=SC2620_JS4_G3_Reign_10488923_253GGLLT4
```

It is read from the actual filenames rather than parsed out of the directory
name, so an unexpected naming variant fails in the sheet check with a clear
message instead of inside Cell Ranger.

### Same library sequenced more than once

Top-ups and extra flow cells: repeat the label. The directories are handed to
Cell Ranger as one comma-separated `--fastqs`, which is how it pools them.

```csv
label,fastqs
JS5,/local/Illumina/DRV/run1/.../Sample_SC2620_JS5_G3_Sofia_10488923_253GGLLT4_L2
JS5,/local/Illumina/DRV/run2/.../Sample_SC2620_JS5_G3_Sofia_10488923_999XYZAB2_L3
```

Note the flow cell is part of the fastq prefix, so a top-up on a different flow
cell has a *different* `--sample` value. Both are collected and passed as
`--sample=prefix1,prefix2`; missing that would silently drop the second run.

### Several libraries from one sample

Separate GEM wells: add the optional `library` column. Barcodes are not
comparable across GEM wells, so each library is counted separately while
staying tied to its `label`.

```csv
label,library,fastqs
JS6,JS6_wellA,/local/.../Sample_SC2620_JS6A_G3_Scooter_10488923_253GGLLT4_L2
JS6,JS6_wellB,/local/.../Sample_SC2620_JS6B_G3_Scooter_10488923_253GGLLT4_L2
```

### Columns per mode

| column | required | notes |
|---|---|---|
| `label` | yes | your short name; becomes `--id` and the output folder |
| `fastqs` | yes | delivery dir, or any fastq in it |
| `library` | no | one GEM well; defaults to `label` |
| `sample` | only if ambiguous | the fastq prefix; required when the path holds more than one |
| `library_type` | no | defaults to `Gene Expression` |

`library_type` uses 10x's own vocabulary, so it lands in `libraries.csv`
verbatim: `Gene Expression`, `Chromatin Accessibility`, `Antibody Capture`,
`CRISPR Guide Capture`, `Multiplexing Capture`, `VDJ-T`, `VDJ-B`. Shorthands
are accepted (`gex`, `rna`, `atac`, `adt`, `citeseq`, `hto`, `crispr`, `cmo`,
`cellplex`).

R1/R2/R3 all live inside the delivery directory, so there is no `fastq3`
column — that existed only back when the sheet listed individual files.

### Multiome (`--mode arc`)

One label carries a Gene Expression library and a Chromatin Accessibility
library. They are not counted separately: `cellranger-arc` takes a
`libraries.csv` naming both, which the pipeline writes for you.

```csv
label,fastqs,sample,library_type
JS4,/local/.../Project_10488522,SC2619_JS4_BC_MG3_8Healthy_10488522_25FWVCLT4,Gene Expression
JS4,/local/.../Project_10488522,SC2619_JS4_MA_8Healthy_10488522_23C52HLT4,Chromatin Accessibility
```

becomes `JS4_libraries.csv`, published next to the results:

```csv
fastqs,sample,library_type
/local/.../Project_10488522,SC2619_JS4_BC_MG3_8Healthy_10488522_25FWVCLT4,Gene Expression
/local/.../Project_10488522,SC2619_JS4_MA_8Healthy_10488522_23C52HLT4,Chromatin Accessibility
```

then:

```bash
cellranger-arc count --id=JS4 --reference=<ref> --libraries=JS4_libraries.csv \
  --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam>
```

Note both rows point at the same Project directory, which holds *both*
libraries — so the `sample` column is required here. One row is one library, and
a path yielding several prefixes is rejected rather than guessed at. Outputs
land in `CELLRANGER_ARC/<label>/`.

## Params

| param | default | notes |
|---|---|---|
| `--id` | `TREx_ID` | TREx Project ID |
| `--sheet` | `sample-sheet.csv` | |
| `--mode` | `gex` | `gex`, `atac`, `arc` |
| `--ref` | — | key from `--listRefs`, or a path to a transcriptome dir |
| `--chemistry` | `auto` | |
| `--expectCells` / `--forceCells` | — | |
| `--createBam` | `false` | |
| `--introns` | tool default | |
| `--r1length` | `28` | `--r1-length`; pass `0` to omit |
| `--r2length` | — | `--r2-length` |
| `--localcores` | `32` | `--localcores` **and** the CPUs reserved |
| `--localmem` | `180` | `--localmem` in GB **and** the memory reserved |
| `--maxforks` | `2` | concurrent tasks **per process** — 2 means 2 Cell Ranger runs at once |
| `--crversion` | `9.0.1` | → `/programs/cellranger-<v>/cellranger` |
| `--arcversion` | `2.2.0` | → `/programs/cellranger-arc-<v>/bin/cellranger-arc` |
| `--programs` | `/programs` | where the installs live |
| `--crpath` | — | full path to the binary; overrides `--crversion` |
| `--arcpath` | — | full path to the arc binary; overrides `--arcversion` |
| `--container` / `--arccontainer` | — | opt in to running inside an image |

Outputs land in `CELLRANGER/<label>/outs` and `pipeline_info/`.

Every flag in the usual invocation is a param, and the command built per
library is:

```bash
<crpath|cellranger> count --id=<label> \
  --localcores=<localcores> --localmem=<localmem> --create-bam=<createBam> --r1-length=<r1length> \
  --transcriptome=<ref> \
  --fastqs=<dir[,dir2]> --sample=<prefix[,prefix2]>
```

`--localcores` and `--localmem` set the Nextflow reservation *and* the Cell
Ranger flags from one value, so what the scheduler holds and what Cell Ranger
believes it has cannot drift apart.

Every process retries a failed task twice before the run stops
(`errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }`, `maxRetries 2`).
`finish` lets tasks already running complete rather than killing them, so a
late failure does not throw away hours of work on the other libraries.

## 10x software

The pipeline runs **your local install** — nothing is pulled and no container is
involved unless you ask for one. Give it a version and the path is built by
convention:

| param | default | resolves to |
|---|---|---|
| `--crversion` | `9.0.1` | `/programs/cellranger-9.0.1/cellranger` |
| `--arcversion` | `2.2.0` | `/programs/cellranger-arc-2.2.0/bin/cellranger-arc` |

Note the two layouts differ: `cellranger` keeps its executable at the top of the
install directory, `cellranger-arc` keeps it in `bin/`. The pipeline knows this.

See what is installed:

```bash
nextflow run https://github.com/bixBeta/10x -r main --listPrograms
```

Pick a version, or point at a binary directly:

```bash
nextflow run https://github.com/bixBeta/10x -r main --crversion 8.0.1 --sheet sample-sheet.csv --ref CanFam3_1
```

```bash
nextflow run https://github.com/bixBeta/10x -r main --crpath /home/me/builds/cellranger --sheet sample-sheet.csv --ref CanFam3_1
```

A different site layout is `--programs /opt/apps`.

The resolved binary is checked before any work starts, so an uninstalled version
fails immediately and tells you what *is* installed, rather than dying inside the
first task.

### Containers

Images are built **on the machine that runs them** and never pulled by the
pipeline. `containers/build-sif.sh` makes a `.sif` either by converting an image
you already have:

```bash
bash containers/build-sif.sh cellranger 9.0.1 docker://ghcr.io/bixbeta/cellranger:9.0.1 /local/workdir/singularity
```

or by building from a 10x tarball, which needs no registry at all:

```bash
bash containers/build-sif.sh cellranger-arc 2.2.0 ~/cellranger-arc-2.2.0.tar.gz /local/workdir/singularity
```

Either way it writes `<sifdir>/<tool>-<version>.sif`, which is exactly the name
the pipeline looks for. Building usually needs root or `--fakeroot`; set
`SIF_FAKEROOT=1` to add the flag.

The script knows the two install layouts — `cellranger` keeps its executable at
the top of the install dir, `cellranger-arc` in `bin/` — and puts both on `PATH`.

Set the SIF location once for everyone in a site config rather than asking
users to pass `--sifdir`, or pin an image outright:

```groovy
params.sifdir = '/local/workdir/singularity'
// or, for one specific image
params.crsif  = 'file:///workdir/TREx_shared/projects/CELLRANGER_9.0.1.sif'
```

### Which version am I actually running?

The filename is only a claim — `cellranger-9.0.1.sif` can contain anything. So
before counting anything, the pipeline asks the tool its own version inside the
image and fails if it disagrees with what you asked for:

```
ERROR: asked for cellranger 9.0.1 but cellranger reports 7.0.0
       Set --crversion to match, point at another image with --crsif, or pass --checkversion false.
```

This applies to `--engine local` too, where `--crpath` can point at any binary.
It costs milliseconds and turns a silent wrong-version result into an immediate
stop. `--checkversion false` disables it. The version Cell Ranger reported is
also recorded in `pipeline_info/software_versions.yml` for every run.

## Development

No Cell Ranger or real data needed:

```bash
bash test/make_test_data.sh
nextflow run . -stub-run -c test/ci.config --sheet test/sample-sheet.csv --ref test/ref
```

CI runs the same thing on every push, against Nextflow 25.04.1 and latest.
