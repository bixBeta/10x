# Nextflow Pipeline for 10x Genomics single cell runs

[![ci](https://github.com/bixBeta/10x/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/bixBeta/10x/actions/workflows/ci.yml)
[![build-container](https://github.com/bixBeta/10x/actions/workflows/build-containers.yml/badge.svg)](https://github.com/bixBeta/10x/actions/workflows/build-containers.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-25.04.1-23aa62.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-1d355c.svg)](https://sylabs.io/singularity/)
[![Cell Ranger](https://img.shields.io/badge/cellranger-9.0.1-blue.svg)](https://www.10xgenomics.com/support/software/cell-ranger)
[![GHCR](https://img.shields.io/badge/ghcr.io-private-lightgrey.svg)](https://github.com/bixBeta/10x/pkgs/container/cellranger)

<hr>

| mode | tool | status |
|---|---|---|
| `gex` | `cellranger count` | working |
| `atac` | `cellranger-atac count` | v0.2 |
| `arc` | `cellranger-arc count` | v0.3 |

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
label,fastq1,fastq2
SS1,SS1_R1.fastq.gz,SS1_R2.fastq.gz
SS2,SS2_R1.fastq.gz,SS2_R2.fastq.gz
```

One row per fastq pair. Fastq names are arbitrary — they get symlinked into the
bcl2fastq convention Cell Ranger requires. Use `--fastqs` if the files sit in a
`fastqs` folder in the project dir, otherwise give absolute paths.

### Same library sequenced more than once

Top-ups and extra flow cells: just repeat the label. The rows are pooled into
**one** `cellranger count` as consecutive lanes, which is what you want — same
GEM well, same cells.

```csv
label,fastq1,fastq2
SS1,run1/SS1_R1.fastq.gz,run1/SS1_R2.fastq.gz
SS1,run2/SS1_R1.fastq.gz,run2/SS1_R2.fastq.gz
```

### Several libraries from one sample

Separate GEM wells: add the optional `library` column. Barcodes are not
comparable across GEM wells, so each library is counted **separately** while
staying tied to its `label` for later aggregation.

```csv
label,library,fastq1,fastq2
SS1,SS1_wellA,A_R1.fastq.gz,A_R2.fastq.gz
SS1,SS1_wellB,B_R1.fastq.gz,B_R2.fastq.gz
```

Rows sharing `library` are pooled; rows differing in `library` are not. The two
rules compose, so a named library can still have several run rows.

### fastq3

Not used by `gex`, and rejected with a clear error if present. It returns with
`--mode atac` and `--mode arc`, where ATAC reads are R1 + R2 (16 bp barcode) +
R3 and all three have to reach `cellranger-atac`. A stray `fastq3` in a GEX
sheet almost always means the wrong `--mode`, which is why it errors rather
than being ignored.

## Params

| param | default | notes |
|---|---|---|
| `--id` | `TREx_ID` | TREx Project ID |
| `--sheet` | `sample-sheet.csv` | |
| `--mode` | `gex` | `gex`, `atac`, `arc` |
| `--ref` | — | key from `--listRefs`, or a path to a transcriptome dir |
| `--fastqs` | — | fastqs live in `fastqs/` in the project dir |
| `--chemistry` | `auto` | |
| `--expectCells` / `--forceCells` | — | |
| `--createBam` | `false` | |
| `--introns` | tool default | |
| `--crversion` | `9.0.1` | Cell Ranger version, selects the container tag |
| `--container` | — | full override, e.g. a local `.sif` |

Outputs land in `CELLRANGER/<library>/outs` and `pipeline_info/`.

## Cell Ranger version

The container tag follows `--crversion`, so a run can be pinned to any built
version without touching the code:

```bash
nextflow run https://github.com/bixBeta/10x -r main --crversion 7.2.0 --sheet sample-sheet.csv --ref GRCh38
```

## Containers

Cell Ranger is licensed software and must not be redistributed, so the images
are **private** on GHCR and this repo never holds the tarball.

Built in two layers so the Cell Ranger install is isolated:

| image | holds | rebuilt |
|---|---|---|
| `cellranger-base:<version>` | the Cell Ranger install, nothing else | only when missing |
| `cellranger:<version>` | anything on top, built `FROM` the base | every run, in seconds |

To build a version:

1. Accept the EULA on the [10x downloads page](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
   to reveal the signed link. It expires, so use it promptly.
2. Store it as the repo secret `TENX_DOWNLOAD_URL` — only the URL, not the whole
   `curl` command.
3. Run the **build-container** workflow from the Actions tab.

Additions go in the marked extension block of
`containers/cellranger/Dockerfile`; leave `Dockerfile.base` alone.

### Pulling on the server

The packages are private, so either export credentials:

```bash
export SINGULARITY_DOCKER_USERNAME=bixBeta
export SINGULARITY_DOCKER_PASSWORD=<ghcr-pat-with-read:packages>
```

or, better for a shared server, pull the SIF once into a shared cache so nobody
else needs credentials:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/singularity_cache
singularity pull docker://ghcr.io/bixbeta/cellranger:9.0.1
```

## Development

No Cell Ranger or real data needed:

```bash
bash test/make_test_data.sh
nextflow run . -stub-run -c test/ci.config --sheet test/sample-sheet.csv --ref test/ref
```

CI runs the same thing on every push, against Nextflow 25.04.1 and latest.
