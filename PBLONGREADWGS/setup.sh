#!/usr/bin/env bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — environment check script
# run once after git clone
#
# What this script does:
#   1. check Java 11+
#   2. check Singularity
#   3. check the nextflow binary (downloaded automatically if missing)
#   4. check containers/ (where the SIF files go)
#   5. check resources/ref/ (where the reference files go)
# ============================================================

set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OK=0; WARN=0

check() { echo "  [OK]   $1"; }
warn()  { echo "  [WARN] $1"; WARN=$((WARN+1)); }
fail()  { echo "  [FAIL] $1"; OK=1; }

echo ""
echo "============================================================"
echo " K_map PacBio WGS Pipeline v2 — environment check"
echo "============================================================"
echo ""

# ── 1. Java ──────────────────────────────────────────────────
echo "[1] Java 11+"
if java -version 2>&1 | grep -qE 'version "1[1-9]|version "[2-9][0-9]'; then
    VER=$(java -version 2>&1 | head -1)
    check "$VER"
else
    fail "Java 11+ is required. conda install -c conda-forge openjdk=17"
fi

# ── 2. Singularity ───────────────────────────────────────────
echo ""
echo "[2] Singularity"
if command -v singularity &>/dev/null; then
    check "$(singularity --version)"
else
    fail "Singularity not found. https://docs.sylabs.io/guides/latest/admin-guide/installation.html"
fi

# ── 3. Nextflow ──────────────────────────────────────────────
echo ""
echo "[3] Nextflow"
if [[ -x "${PIPELINE_DIR}/nextflow" ]]; then
    check "nextflow wrapper found in repo"
elif command -v nextflow &>/dev/null; then
    check "system nextflow: $(nextflow -version 2>&1 | head -1)"
else
    warn "nextflow not found — trying automatic download"
    curl -fsSL https://get.nextflow.io | bash
    mv nextflow "${PIPELINE_DIR}/nextflow"
    chmod +x "${PIPELINE_DIR}/nextflow"
    check "nextflow downloaded: ${PIPELINE_DIR}/nextflow"
fi

# ── 4. Containers ────────────────────────────────────────────
echo ""
echo "[4] Singularity containers (.sif)"
MISSING_SIF=0
for sif in kmap_pbwgs_v2.sif deepvariant_1.9.0-gpu.sif; do
    if [[ -f "${PIPELINE_DIR}/containers/${sif}" ]]; then
        SIZE=$(du -sh "${PIPELINE_DIR}/containers/${sif}" | cut -f1)
        check "containers/${sif} (${SIZE})"
    else
        fail "containers/${sif} missing"
        MISSING_SIF=$((MISSING_SIF+1))
    fi
done
if [[ $MISSING_SIF -gt 0 ]]; then
    echo ""
    echo "  Copy the SIF files into containers/:"
    echo "    cp /path/to/*.sif ${PIPELINE_DIR}/containers/"
fi

# ── 5. Reference files ───────────────────────────────────────
echo ""
echo "[5] Reference files (resources/ref/)"
check_ref() {
    local f="${PIPELINE_DIR}/resources/ref/$1"
    if [[ -e "$f" ]]; then
        SIZE=$(du -sh "$f" 2>/dev/null | cut -f1)
        check "resources/ref/$1 (${SIZE})"
    else
        fail "resources/ref/$1 missing"
    fi
}
check_ref "GRCh38_no_alt_analysis_set.fna"
check_ref "GRCh38_no_alt_analysis_set.fna.mmi"
check_ref "trgt/adotto_strchive_20250827.hg38.noFailTag.bed.gz"
check_ref "sex/GRCh38.PAR.bed"            # DeepVariant haploid (XY)
check_ref "sex/expected_cn.hg38.XX.bed"   # Sawfish
check_ref "sex/expected_cn.hg38.XY.bed"   # Sawfish

echo ""
echo "============================================================"
if [[ $OK -eq 0 && $WARN -eq 0 ]]; then
    echo " All checks passed. The pipeline is ready to run."
    echo ""
    echo " Example:"
    echo "   bash run.sh --samplesheet input/samples.csv"
elif [[ $OK -eq 0 ]]; then
    echo " ${WARN} warning(s) — the pipeline can run, but check the items above."
else
    echo " There are errors. Fix the [FAIL] items above and run again."
fi
echo "============================================================"
echo ""
exit $OK
