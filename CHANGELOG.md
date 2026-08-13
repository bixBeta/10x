# Changelog

## v0.1.0 — 2026-08-13

Initial release.

- GEX quantification with `cellranger count`, one run per library.
- Samplesheet schema with `sample` / `library` / `modality`, so re-sequencing
  runs merge as pseudo-lanes while separate GEM wells stay separate.
- `atac`, `arc_gex` and `arc_atac` rows are parsed and validated but not yet
  quantified.
- Cell Ranger version selectable per run via `--cellranger_version`.
- Docker / Singularity / Apptainer profiles, SLURM executor profile.
- Private GHCR container build workflow; stub-run CI.

## Planned

- v0.2 — ATAC via `cellranger-atac count`
- v0.3 — Multiome via `cellranger-arc count` (generated `libraries.csv`)
- v0.4 — `cellranger aggr` across libraries of a sample, MultiQC report
