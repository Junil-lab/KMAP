#!/bin/bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — 컨테이너 빌드 스크립트
#
# 사용법:
#   bash build_containers.sh          # 전체 빌드
#   bash build_containers.sh main      # 메인 컨테이너만
#   bash build_containers.sh deepvariant  # DV 컨테이너만 복사
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/../containers"
DEF_DIR="${SCRIPT_DIR}"

MODE="${1:-all}"

build_main() {
    echo "[$(date)] 메인 컨테이너 빌드 시작 (약 15~30분)..."
    echo "[$(date)] kmap_pbwgs_v2.sif — pbmm2/samtools/hiphase/trgt/sawfish/CpG/jasmine/AnnotSV"

    # NFS(NAS_2)에 --fakeroot로 직접 쓰기 불가 → 로컬 SSD에 빌드 후 NAS_2로 복사
    LOCAL_SIF="/home/Data_Drive_8TB_3/kmap_pbwgs_v2_build.sif"

    singularity build \
        --fakeroot \
        "${LOCAL_SIF}" \
        "${DEF_DIR}/kmap_pbwgs.def" \
        2>&1 | tee /tmp/build_main_v2.log

    echo "[$(date)] 메인 컨테이너 빌드 완료: $(du -sh ${LOCAL_SIF} | cut -f1)"

    echo "[$(date)] 컨테이너 테스트..."
    singularity test "${LOCAL_SIF}"

    echo "[$(date)] NAS로 복사 중..."
    cp "${LOCAL_SIF}" "${OUT_DIR}/kmap_pbwgs_v2.sif"
    rm -f "${LOCAL_SIF}"
    echo "[$(date)] 복사 완료: $(du -sh ${OUT_DIR}/kmap_pbwgs_v2.sif | cut -f1)"
}

copy_deepvariant() {
    # DeepVariant는 Google 공식 SIF를 그대로 사용
    echo "[$(date)] DeepVariant SIF 복사..."
    DV_SRC="/home/NAS/senc0711/K_map/deepvariant/deepvariant_1.9.0-gpu.sif"
    if [ -f "$DV_SRC" ]; then
        cp "$DV_SRC" "${OUT_DIR}/deepvariant_1.9.0-gpu.sif"
        echo "[$(date)] DeepVariant SIF 복사 완료: $(du -sh ${OUT_DIR}/deepvariant_1.9.0-gpu.sif | cut -f1)"
    else
        echo "[WARN] DeepVariant SIF 없음 — 수동으로 복사 필요"
        echo "  출처: https://github.com/google/deepvariant/releases"
        echo "  목적: ${OUT_DIR}/deepvariant_1.9.0-gpu.sif"
    fi
}

copy_vep() {
    # VEP SIF 복사
    echo "[$(date)] VEP SIF 복사..."
    VEP_SRC="/home/NAS/K_MAP/Pacbio/WGS/Pipeline/ver1/output/hiphase/VEP/container/ensembl-vep-release_115.0.sif"
    if [ -f "$VEP_SRC" ]; then
        cp "$VEP_SRC" "${OUT_DIR}/ensembl-vep-release_115.0.sif"
        echo "[$(date)] VEP SIF 복사 완료"
    else
        echo "[WARN] VEP SIF 없음"
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
        echo "=== 빌드 완료 ==="
        ls -lh "${OUT_DIR}/"*.sif
        ;;
    *)
        echo "사용법: bash build_containers.sh [main|deepvariant|vep|all]"
        exit 1
        ;;
esac
