# syntax=docker/dockerfile:1.7
# =============================================================================
#  OVERLAY layer: everything the pipeline needs on top of the 10x install.
# =============================================================================
#  Builds FROM the base image, which holds Cell Ranger and nothing else. Adding
#  or changing anything here rebuilds in seconds and never touches the Cell
#  Ranger install - no re-download, no re-extract, and no need for the 10x
#  download secret.
#
#  Pushed as:  ghcr.io/<owner>/<tool>:<version>
#
#  Add future tooling in the block marked below. Keep the base untouched.
# =============================================================================

ARG BASE_IMAGE
FROM ${BASE_IMAGE}

ARG TOOL=cellranger
ARG VERSION=9.0.1

# ---- extension point --------------------------------------------------------
# Extra utilities go here, e.g. samtools for post-processing or python3 for
# format conversion. One RUN per concern keeps the layers cacheable.
#
#   RUN apt-get update && apt-get install -y --no-install-recommends \
#           samtools python3-minimal \
#       && rm -rf /var/lib/apt/lists/*
# -----------------------------------------------------------------------------

LABEL org.opencontainers.image.title="${TOOL}" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.description="10x Genomics ${TOOL} ${VERSION} for the 10x Nextflow pipeline - private, redistribution prohibited by the 10x EULA" \
      org.opencontainers.image.source="https://github.com/bixBeta/10x" \
      org.opencontainers.image.licenses="LicenseRef-10x-EULA"

# Cell Ranger writes temp files relative to cwd; Nextflow sets that per task.
WORKDIR /tmp
CMD ["bash"]
