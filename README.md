# 10x

A modular Nextflow (DSL2) pipeline for 10x Genomics single-cell data.

| Modality | Tool | Status |
|---|---|---|
| Gene expression (`gex`) | `cellranger count` | implemented |
| ATAC (`atac`) | `cellranger-atac count` | planned, v0.2 |
| Multiome (`arc_gex` + `arc_atac`) | `cellranger-arc count` | planned, v0.3 |

The samplesheet schema already accepts all three, so adding a modality means
adding a subworkflow — existing samplesheets keep working.

Requires Nextflow >= 24.04 (developed and tested against 25.04.1).

## Quick start

```bash
nextflow run bixBeta/10x -profile singularity \
    --input samplesheet.csv \
    --reference /refs/refdata-gex-GRCh38-2024-A \
    --outdir results
```

Pin a specific Cell Ranger version for a run:

```bash
nextflow run bixBeta/10x -profile singularity --cellranger_version 7.2.0 --input samplesheet.csv --reference /refs/refdata-gex-GRCh38-2024-A
```

Validate a samplesheet without running anything:

```bash
nextflow run bixBeta/10x --input samplesheet.csv --validate_only
```

## Samplesheet

| Column | Required | Notes |
|---|---|---|
| `sample` | yes | Biological sample; the aggregation key |
| `library` | no | One physical 10x library (one GEM well). Defaults to `sample` |
| `fastq_1` | yes | |
| `fastq_2` | yes | |
| `fastq_3` | ATAC only | Required for `atac`/`arc_atac`, forbidden otherwise |
| `modality` | no | `gex` (default), `atac`, `arc_gex`, `arc_atac` |

FASTQs can be named anything — the pipeline symlinks them into the
bcl2fastq-style names Cell Ranger insists on.

### The two "multiple runs" cases

**Same library, sequenced more than once** (top-up runs, extra flowcells) — repeat
the `sample`, leave `library` alone. The rows are merged into a single
`cellranger count` as consecutive pseudo-lanes:

```csv
sample,library,fastq_1,fastq_2,modality
pbmc,,/runs/jul/pbmc_R1.fq.gz,/runs/jul/pbmc_R2.fq.gz,gex
pbmc,,/runs/aug/pbmc_R1.fq.gz,/runs/aug/pbmc_R2.fq.gz,gex
```

**Several libraries from one sample** (separate GEM wells) — give each its own
`library`. Barcodes are not comparable across GEM wells, so each library is
counted separately; they stay linked through `sample` for later aggregation:

```csv
sample,library,fastq_1,fastq_2,modality
pbmc,pbmc_wellA,/runs/jul/A_R1.fq.gz,/runs/jul/A_R2.fq.gz,gex
pbmc,pbmc_wellB,/runs/jul/B_R1.fq.gz,/runs/jul/B_R2.fq.gz,gex
```

The two compose: a library with its own `library` name can still have several
run rows.

### Multiome (ARC)

One sample, one `arc_gex` library and one `arc_atac` library. The pipeline
generates the `libraries.csv` that `cellranger-arc count` needs.

```csv
sample,library,fastq_1,fastq_2,fastq_3,modality
brain_1,brain_1_gex,/runs/b1gex_R1.fq.gz,/runs/b1gex_R2.fq.gz,,arc_gex
brain_1,brain_1_atac,/runs/b1atac_R1.fq.gz,/runs/b1atac_R2.fq.gz,/runs/b1atac_R3.fq.gz,arc_atac
```

`arc_*` is deliberately distinct from plain `gex`/`atac`: a standalone scRNA
library and a standalone scATAC library from the same tissue are *not* multiome
and cannot go through `cellranger-arc`. Stating intent in the samplesheet lets
validation catch that instead of guessing.

## Parameters

| Param | Default | Description |
|---|---|---|
| `--input` | — | Samplesheet CSV |
| `--reference` | — | Cell Ranger transcriptome directory (GEX) |
| `--outdir` | `./results` | |
| `--cellranger_version` | `9.0.1` | Selects the GHCR image tag |
| `--cellranger_container` | — | Full image or `.sif` path; overrides the version |
| `--chemistry` | `auto` | |
| `--expect_cells` / `--force_cells` | — | |
| `--create_bam` | `false` | |
| `--include_introns` | tool default | |
| `--max_cpus` / `--max_memory` / `--max_time` | `32` / `256.GB` / `96.h` | Ceilings |
| `--validate_only` | `false` | Parse and validate the samplesheet, then stop |

Any flag not exposed above can go through `task.ext.args`:

```groovy
process {
    withName: CELLRANGER_COUNT {
        ext.args = '--nosecondary --r1-length 26'
    }
}
```

## Containers

Cell Ranger is licensed software and **must not be redistributed**, so there is
no public image and the tarball never touches this (public) repository. Build
your own private image:

1. Accept the EULA on the [10x downloads page](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
   to reveal the signed download link. It expires, so grab a fresh one per build.
2. Store it as a repository secret:

   ```bash
   gh secret set TENX_DOWNLOAD_URL --body '<signed-url>'
   ```

3. Run the **build-container** workflow from the Actions tab, choosing the tool
   and version. It pushes `ghcr.io/bixbeta/cellranger:<version>`.
4. Confirm the package is **private** in your GHCR package settings. The repo
   being public does not make the package public, but check it once.

Prefer to keep binaries off GitHub entirely? Build locally and push by hand:

```bash
cp ~/Downloads/cellranger-9.0.1.tar.gz containers/cellranger/
docker build containers/cellranger --build-arg TOOL=cellranger --build-arg VERSION=9.0.1 -t ghcr.io/bixbeta/cellranger:9.0.1
docker push ghcr.io/bixbeta/cellranger:9.0.1
```

Multiple versions can coexist as tags; `--cellranger_version` picks one per run.

### Singularity on the server

GHCR needs credentials for private images. Use a classic PAT with `read:packages`:

```bash
export SINGULARITY_DOCKER_USERNAME=bixBeta
export SINGULARITY_DOCKER_PASSWORD=<ghcr-pat>
export NXF_SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
nextflow run bixBeta/10x -profile singularity,slurm --input samplesheet.csv --reference /refs/refdata-gex-GRCh38-2024-A
```

To avoid pulling on every node, pre-build the SIF once and point at it directly:

```bash
singularity pull /shared/sif/cellranger-9.0.1.sif docker://ghcr.io/bixbeta/cellranger:9.0.1
nextflow run bixBeta/10x -profile singularity --cellranger_container /shared/sif/cellranger-9.0.1.sif --input samplesheet.csv --reference /refs/refdata-gex-GRCh38-2024-A
```

## Profiles

`docker`, `singularity`, `apptainer`, `slurm`, `debug`, and the stub profiles
`test`, `test_multi`, `test_arc`. Combine them: `-profile singularity,slurm`.

## Development

Everything is testable without Cell Ranger or real data:

```bash
bash assets/test_data/make_test_data.sh
nextflow run . -profile test_multi -stub-run
```

CI runs exactly this on every push.

## Layout

```
main.nf                                  entrypoint
workflows/tenx.nf                        dispatch by modality
subworkflows/local/samplesheet.nf        parse, validate, group
subworkflows/local/gex.nf                GEX path
modules/local/cellranger/count/main.nf   cellranger count
conf/                                    resources, publishing, test profiles
containers/cellranger/Dockerfile         private GHCR image
```

## Licence

Pipeline code: MIT. Cell Ranger and the 10x reference packages are covered by
the 10x Genomics EULA and are not distributed here.
