#!/bin/bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — container build script
#
# Usage:
#   bash build_containers.sh          # build everything
#   bash build_containers.sh main      # main container only
#   bash build_containers.sh deepvariant  # copy the DV container only
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/../containers"
DEF_DIR="${SCRIPT_DIR}"

MODE="${1:-all}"

build_main() {
    echo "[$(date)] building main container (about 15-30 min)..."
    echo "[$(date)] kmap_pbwgs_v2.sif — pbmm2/samtools/hiphase/trgt/sawfish/CpG/jasmine/AnnotSV"

    # --fakeroot cannot write directly to NFS (NAS_2) → build on local SSD, then copy to NAS_2
    LOCAL_SIF="/home/Data_Drive_8TB_3/kmap_pbwgs_v2_build.sif"

    singularity build \
        --fakeroot \
        "${LOCAL_SIF}" \
        "${DEF_DIR}/kmap_pbwgs.def" \
        2>&1 | tee /tmp/build_main_v2.log

    echo "[$(date)] main container built: $(du -sh ${LOCAL_SIF} | cut -f1)"

    echo "[$(date)] testing container..."
    singularity test "${LOCAL_SIF}"

    echo "[$(date)] copying to NAS..."
    cp "${LOCAL_SIF}" "${OUT_DIR}/kmap_pbwgs_v2.sif"
    rm -f "${LOCAL_SIF}"
    echo "[$(date)] copied: $(du -sh ${OUT_DIR}/kmap_pbwgs_v2.sif | cut -f1)"
}

copy_deepvariant() {
    # DeepVariant uses the official Google SIF as is
    echo "[$(date)] copying DeepVariant SIF..."
    DV_SRC="/home/NAS/senc0711/K_map/deepvariant/deepvariant_1.9.0-gpu.sif"
    if [ -f "$DV_SRC" ]; then
        cp "$DV_SRC" "${OUT_DIR}/deepvariant_1.9.0-gpu.sif"
        echo "[$(date)] DeepVariant SIF copied: $(du -sh ${OUT_DIR}/deepvariant_1.9.0-gpu.sif | cut -f1)"
    else
        echo "[WARN] DeepVariant SIF not found — copy it manually"
        echo "  Source: https://github.com/google/deepvariant/releases"
        echo "  Destination: ${OUT_DIR}/deepvariant_1.9.0-gpu.sif"
    fi
}

copy_vep() {
    # copy the VEP SIF
    echo "[$(date)] copying VEP SIF..."
    VEP_SRC="/home/NAS/K_MAP/Pacbio/WGS/Pipeline/ver1/output/hiphase/VEP/container/ensembl-vep-release_115.0.sif"
    if [ -f "$VEP_SRC" ]; then
        cp "$VEP_SRC" "${OUT_DIR}/ensembl-vep-release_115.0.sif"
        echo "[$(date)] VEP SIF copied"
    else
        echo "[WARN] VEP SIF not found"
    fi
}

case "$MODE" in
    main)        build_main ;;
    deepvariant) copy_deepvariant ;;
    vep)         copy_vep ;;
    all)
        build_main
        copy_deepvariant
        copy_vep
        echo ""
        echo "=== build finished ==="
        ls -lh "${OUT_DIR}/"*.sif
        ;;
    *)
        echo "Usage: bash build_containers.sh [main|deepvariant|vep|all]"
        exit 1
        ;;
esac
