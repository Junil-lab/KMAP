#!/bin/bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — container build script
#
# Usage:
#   bash build_containers.sh               # build everything
#   bash build_containers.sh main          # build the main container only
#   bash build_containers.sh deepvariant   # pull the DeepVariant container only
#   bash build_containers.sh methbat       # pull the MethBat container only
#
# The prebuilt main container is also available from the GitHub Release:
#   https://github.com/Junil-lab/KMAP/releases/tag/v1.0.0
#
# BUILD_DIR: where the image is built (default: $TMPDIR or /tmp).
#   --fakeroot builds cannot write to most NFS mounts, so use a local disk.
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/../containers"
DEF_DIR="${SCRIPT_DIR}"
BUILD_DIR="${BUILD_DIR:-${TMPDIR:-/tmp}}"
DV_IMAGE="docker://google/deepvariant:1.9.0-gpu"
METHBAT_IMAGE="docker://quay.io/biocontainers/methbat:1.1.0--h9ee0642_0"

MODE="${1:-all}"
mkdir -p "${OUT_DIR}"

build_main() {
    echo "[$(date)] building main container (about 15-30 min)..."
    echo "[$(date)] kmap_pbwgs_v2.sif — pbmm2/samtools/bcftools/hiphase/trgt/sawfish/CpG/jasmine"

    LOCAL_SIF="${BUILD_DIR}/kmap_pbwgs_v2_build.$$.sif"

    singularity build \
        --fakeroot \
        "${LOCAL_SIF}" \
        "${DEF_DIR}/kmap_pbwgs.def" \
        2>&1 | tee "${BUILD_DIR}/build_main_v2.log"

    echo "[$(date)] main container built: $(du -sh "${LOCAL_SIF}" | cut -f1)"

    echo "[$(date)] testing container..."
    singularity test "${LOCAL_SIF}"

    mv "${LOCAL_SIF}" "${OUT_DIR}/kmap_pbwgs_v2.sif"
    echo "[$(date)] saved: ${OUT_DIR}/kmap_pbwgs_v2.sif"
}

pull_deepvariant() {
    # DeepVariant uses the official Google image as is
    echo "[$(date)] pulling ${DV_IMAGE}..."
    singularity pull --force "${OUT_DIR}/deepvariant_1.9.0-gpu.sif" "${DV_IMAGE}"
    echo "[$(date)] saved: ${OUT_DIR}/deepvariant_1.9.0-gpu.sif"
}

pull_methbat() {
    # MethBat uses the official bioconda (biocontainers) image as is
    echo "[$(date)] pulling ${METHBAT_IMAGE}..."
    singularity pull --force "${OUT_DIR}/methbat_1.1.0.sif" "${METHBAT_IMAGE}"
    echo "[$(date)] saved: ${OUT_DIR}/methbat_1.1.0.sif"
}

case "$MODE" in
    main)        build_main ;;
    deepvariant) pull_deepvariant ;;
    methbat)     pull_methbat ;;
    all)
        build_main
        pull_deepvariant
        pull_methbat
        echo ""
        echo "=== build finished ==="
        ls -lh "${OUT_DIR}/"*.sif
        ;;
    *)
        echo "Usage: bash build_containers.sh [main|deepvariant|methbat|all]"
        exit 1
        ;;
esac
