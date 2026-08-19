#!/usr/bin/env bash
# =============================================================================
#  Build a Singularity image for a 10x tool.
#
#  Each Cell Ranger version needs its own image: the install is baked in, so an
#  existing <tool>-<version>.sif cannot be upgraded in place.
#
#    ./containers/build-sif.sh <tool> <version> <source> [outdir]
#
#  <source> is a 10x tarball, which is built from scratch:
#    ./containers/build-sif.sh cellranger-arc 2.2.0 ~/cellranger-arc-2.2.0.tar.gz /local/workdir/singularity
#
#  or an existing image, which is just converted:
#    ./containers/build-sif.sh cellranger 9.0.1 docker://<registry>/cellranger:9.0.1 /local/workdir/singularity
#
#  Produces <outdir>/<tool>-<version>.sif, which is exactly what the pipeline
#  looks for, i.e. sifdir + crversion / arcversion in params.yaml
#
#  Nothing is pushed anywhere; the image stays on this machine.
#
#  Building usually needs either root or --fakeroot. Set SIF_FAKEROOT=1 to add
#  --fakeroot, or run the script under sudo.
#
#  Sites that auto-bind paths ( apptainer.conf "bind path = /workdir" ) fail the
#  build with
#      destination /workdir doesn't exist in container
#  because build has no overlay to create the mount point with. The mount points
#  are therefore created in %setup, which runs on the host before %post. Adjust
#  the list with SIF_BINDS, e.g.
#      SIF_BINDS="/workdir /local /programs /fs" bash build-sif.sh ...
#  As a last resort apptainer build also takes --no-mount bind-paths.
#
#  For the same reason the tarball is staged at /opt, not /tmp: the host /tmp is
#  bind mounted over the container's during %post, hiding anything %files put
#  there.
# =============================================================================

set -euo pipefail

TOOL="${1:-}"
VERSION="${2:-}"
SOURCE="${3:-}"
OUTDIR="${4:-$PWD}"

if [ -z "$TOOL" ] || [ -z "$VERSION" ] || [ -z "$SOURCE" ] ; then
    sed -n '2,22p' "$0"
    exit 1
fi

case "$TOOL" in
    cellranger|cellranger-arc|cellranger-atac) : ;;
    *) echo "unknown tool '$TOOL' (expected cellranger, cellranger-arc or cellranger-atac)" >&2 ; exit 1 ;;
esac

# resolved only once the inputs are known good, so a bad tarball reports itself
# rather than hiding behind a missing runtime
RUNNER=""
need_runner () {
    [ -n "$RUNNER" ] && return
    command -v singularity > /dev/null 2>&1 || command -v apptainer > /dev/null 2>&1 || {
        echo "neither singularity nor apptainer is on PATH" >&2 ; exit 1 ; }
    RUNNER=$( command -v singularity 2>/dev/null || command -v apptainer )
}

mkdir -p "$OUTDIR"
SIF="${OUTDIR}/${TOOL}-${VERSION}.sif"

FAKEROOT=""
[ "${SIF_FAKEROOT:-0}" = "1" ] && FAKEROOT="--fakeroot"

# Mount points that must exist inside the image, because the site auto-binds
# them. They are also what the pipeline binds at run time.
SIF_BINDS="${SIF_BINDS:-/workdir /local /programs}"

# ---- source is an existing image: convert it straight to a .sif -------------
# Private packages need credentials first:
#   export SINGULARITY_DOCKER_USERNAME=<user>
#   export SINGULARITY_DOCKER_PASSWORD=<token with read:packages>
case "$SOURCE" in
    docker://*|oras://*|library://*|shub://*|*.sif)
        need_runner
        echo "building ${SIF} from ${SOURCE}"
        # shellcheck disable=SC2086
        "$RUNNER" build $FAKEROOT "$SIF" "$SOURCE"

        echo
        echo "built: ${SIF}"
        "$RUNNER" exec "$SIF" "$TOOL" --version || true
        exit 0
        ;;
esac

# ---- source is a tarball: build from scratch --------------------------------
TARBALL="$SOURCE"
[ -f "$TARBALL" ] || { echo "no such tarball, and not an image uri: $TARBALL" >&2 ; exit 1 ; }

# the tarball must match the version, or the image would be mislabelled.
# Read the top level component from the first entries, tolerating a leading ./
TAR_ERR=$( mktemp )
set +o pipefail
TOP=$( tar -tzf "$TARBALL" 2>"$TAR_ERR" | head -20 | sed 's#^\./##' | awk -F/ 'NF>0 && $1!="" {print $1}' | sort -u | head -n1 )
set -o pipefail

if [ -z "$TOP" ] ; then
    echo "could not read '$TARBALL' as a gzipped tar." >&2
    echo >&2
    echo "  size : $( du -h "$TARBALL" 2>/dev/null | cut -f1 )" >&2
    echo "  type : $( file -b "$TARBALL" 2>/dev/null )" >&2
    if [ -s "$TAR_ERR" ] ; then
        echo "  tar  : $( head -n2 "$TAR_ERR" | tr '\n' ' ' )" >&2
    fi
    echo >&2
    echo "  A truncated download, or an expired 10x link that returned an error" >&2
    echo "  page instead of the archive, both land here. Check the size against" >&2
    echo "  the download page and fetch it again if it is short." >&2
    rm -f "$TAR_ERR"
    exit 1
fi
rm -f "$TAR_ERR"

if [ "$TOP" != "${TOOL}-${VERSION}" ] ; then
    echo "tarball contains '${TOP}/', expected '${TOOL}-${VERSION}/'" >&2
    echo "  build it as:  bash $0 ${TOP%-*} ${TOP##*-} $TARBALL $OUTDIR" >&2
    exit 1
fi

DEF=$( mktemp /tmp/${TOOL}-${VERSION}.XXXXXX.def )
trap 'rm -f "$DEF"' EXIT

# %files copies from the host at build time, so the tarball never has to be
# reachable afterwards. cellranger keeps its executable at the top of the
# install dir, cellranger-arc / -atac keep it in bin/ - both go on PATH.
cat > "$DEF" <<EOF
Bootstrap: docker
From: ubuntu:22.04

# Runs on the host, before the %post container is started. Auto-bound paths
# need a destination inside the image or the build cannot enter it at all.
%setup
    ROOT="\${APPTAINER_ROOTFS:-\${SINGULARITY_ROOTFS}}"
    for d in ${SIF_BINDS} ; do
        mkdir -p "\${ROOT}\${d}"
    done

# NOT /tmp: apptainer bind mounts the host /tmp during %post, which would hide
# whatever %files put there. /opt is a plain directory in the image.
%files
    ${TARBALL} /opt/tool.tar.gz

%post
    # keep the mount points in the finished image, for run time binds
    for d in ${SIF_BINDS} ; do
        mkdir -p "\${d}"
    done

    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \\
        bash ca-certificates procps zlib1g libbz2-1.0 liblzma5
    rm -rf /var/lib/apt/lists/*

    tar -xzf /opt/tool.tar.gz -C /opt
    rm -f /opt/tool.tar.gz

    test -x "/opt/${TOOL}-${VERSION}/${TOOL}" \\
      || test -x "/opt/${TOOL}-${VERSION}/bin/${TOOL}" \\
      || { echo "no ${TOOL} executable in the extracted install" >&2 ; exit 1 ; }

%environment
    export LANG=C.UTF-8
    export LC_ALL=C.UTF-8
    export PATH="/opt/${TOOL}-${VERSION}:/opt/${TOOL}-${VERSION}/bin:\$PATH"

%runscript
    exec "\$@"

%labels
    Tool ${TOOL}
    Version ${VERSION}
EOF

need_runner
echo "building ${SIF}"

# shellcheck disable=SC2086
"$RUNNER" build $FAKEROOT "$SIF" "$DEF"

echo
echo "built: ${SIF}"
"$RUNNER" exec "$SIF" "$TOOL" --version || true
echo
case "$TOOL" in
    cellranger)      VERFLAG="--crversion"   ;;
    cellranger-arc)  VERFLAG="--arcversion"  ;;
    *)               VERFLAG="--crversion"   ;;
esac
echo "use it with, in params.yaml:"
echo "  sifdir:     ${OUTDIR}"
echo "  ${VERFLAG#--}: ${VERSION}"
