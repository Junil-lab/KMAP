#!/usr/bin/env bash
# Build missing images; use --rebuild to replace images after a successful build/test.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE=all
REBUILD=0
for arg in "$@"; do
    case "$arg" in all|main|sqanti3|discovery) MODE="$arg" ;; --rebuild) REBUILD=1 ;;
        *) echo 'Usage: bash def/build_containers.sh [all|main|sqanti3|discovery] [--rebuild]' >&2; exit 2 ;; esac
done
mkdir -p "$PIPELINE_DIR/containers"
cd "$PIPELINE_DIR"
for name in main sqanti3 discovery; do
    [[ "$MODE" == all || "$MODE" == "$name" ]] || continue
    image="$PIPELINE_DIR/containers/kmap_isoseq_${name}.sif"
    if [[ -f "$image" ]] && (( ! REBUILD )); then
        singularity test "$image"
        continue
    fi
    BUILD_TMP="$(mktemp -d "$PIPELINE_DIR/containers/.build-${name}.XXXXXX")"
    trap 'rm -rf "$BUILD_TMP"' EXIT
    singularity build --fakeroot "$BUILD_TMP/image.sif" "def/kmap_isoseq_${name}.def"
    singularity test "$BUILD_TMP/image.sif"
    if [[ -e "$image" ]]; then
        backup="$PIPELINE_DIR/backups/containers_$(date +%Y%m%d_%H%M%S)"
        mkdir -p "$backup"
        mv "$image" "$backup/"
    fi
    mv "$BUILD_TMP/image.sif" "$image"
    rmdir "$BUILD_TMP"
    trap - EXIT
done
bash "$PIPELINE_DIR/bin/check_containers.sh" \
    "$PIPELINE_DIR/containers/kmap_isoseq_main.sif" \
    "$PIPELINE_DIR/containers/kmap_isoseq_sqanti3.sif" \
    "$PIPELINE_DIR/containers/kmap_isoseq_discovery.sif"
