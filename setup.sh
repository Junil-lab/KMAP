#!/usr/bin/env bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — 초기 설정 확인 스크립트
# git clone 후 처음 한 번만 실행
#
# 이 스크립트가 하는 일:
#   1. Java 11+ 존재 확인
#   2. Singularity 존재 확인
#   3. nextflow 바이너리 존재 확인 (없으면 자동 다운로드)
#   4. containers/ 폴더 확인 (SIF 파일 위치 안내)
#   5. resources/ref/ 폴더 확인 (레퍼런스 파일 위치 안내)
# ============================================================

set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OK=0; WARN=0

check() { echo "  [OK]   $1"; }
warn()  { echo "  [WARN] $1"; WARN=$((WARN+1)); }
fail()  { echo "  [FAIL] $1"; OK=1; }

echo ""
echo "============================================================"
echo " K_map PacBio WGS Pipeline v2 — 환경 확인"
echo "============================================================"
echo ""

# ── 1. Java ──────────────────────────────────────────────────
echo "[1] Java 11+"
if java -version 2>&1 | grep -qE 'version "1[1-9]|version "[2-9][0-9]'; then
    VER=$(java -version 2>&1 | head -1)
    check "$VER"
else
    fail "Java 11+가 필요합니다. conda install -c conda-forge openjdk=17"
fi

# ── 2. Singularity ───────────────────────────────────────────
echo ""
echo "[2] Singularity"
if command -v singularity &>/dev/null; then
    check "$(singularity --version)"
else
    fail "Singularity가 없습니다. https://docs.sylabs.io/guides/latest/admin-guide/installation.html"
fi

# ── 3. Nextflow ──────────────────────────────────────────────
echo ""
echo "[3] Nextflow"
if [[ -x "${PIPELINE_DIR}/nextflow" ]]; then
    check "repo 내 nextflow 래퍼 존재"
elif command -v nextflow &>/dev/null; then
    check "시스템 nextflow: $(nextflow -version 2>&1 | head -1)"
else
    warn "nextflow를 찾을 수 없습니다 — 자동 다운로드 시도"
    curl -fsSL https://get.nextflow.io | bash
    mv nextflow "${PIPELINE_DIR}/nextflow"
    chmod +x "${PIPELINE_DIR}/nextflow"
    check "nextflow 다운로드 완료: ${PIPELINE_DIR}/nextflow"
fi

# ── 4. Containers ────────────────────────────────────────────
echo ""
echo "[4] Singularity 컨테이너 (.sif)"
MISSING_SIF=0
for sif in kmap_pbwgs_v2.sif deepvariant_1.9.0-gpu.sif ensembl-vep-release_115.0.sif; do
    if [[ -f "${PIPELINE_DIR}/containers/${sif}" ]]; then
        SIZE=$(du -sh "${PIPELINE_DIR}/containers/${sif}" | cut -f1)
        check "containers/${sif} (${SIZE})"
    else
        fail "containers/${sif} 없음"
        MISSING_SIF=$((MISSING_SIF+1))
    fi
done
if [[ $MISSING_SIF -gt 0 ]]; then
    echo ""
    echo "  SIF 파일을 containers/ 폴더에 복사하세요:"
    echo "    cp /path/to/*.sif ${PIPELINE_DIR}/containers/"
fi

# ── 5. Reference files ───────────────────────────────────────
echo ""
echo "[5] 레퍼런스 파일 (resources/ref/)"
check_ref() {
    local f="${PIPELINE_DIR}/resources/ref/$1"
    if [[ -e "$f" ]]; then
        SIZE=$(du -sh "$f" 2>/dev/null | cut -f1)
        check "resources/ref/$1 (${SIZE})"
    else
        fail "resources/ref/$1 없음"
    fi
}
check_ref "GRCh38_no_alt_analysis_set.fna"
check_ref "GRCh38_no_alt_analysis_set.fna.mmi"
check_ref "repeat_catalog.hg38.bed"
check_ref "TE_ref/teref.human.fa"
check_ref "sTELLeR/steller/steller.py"
check_ref "cache"       # VEP cache dir
check_ref "AnnotSV"     # AnnotSV DB dir

echo ""
echo "============================================================"
if [[ $OK -eq 0 && $WARN -eq 0 ]]; then
    echo " 모든 항목 확인 완료. 파이프라인 실행 준비 완료."
    echo ""
    echo " 실행 예시:"
    echo "   bash run.sh --samplesheet input/samples.csv"
elif [[ $OK -eq 0 ]]; then
    echo " 경고 ${WARN}개 — 실행은 가능하나 위 항목을 확인하세요."
else
    echo " 오류가 있습니다. 위 [FAIL] 항목을 해결 후 재실행하세요."
fi
echo "============================================================"
echo ""
exit $OK
