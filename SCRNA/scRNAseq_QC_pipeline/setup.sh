#!/usr/bin/env bash
# ============================================================
# K_map scRNA-seq QC Pipeline - setup check
# ============================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OK=0
WARN=0

check() { echo "  [OK]   $1"; }
warn() { echo "  [WARN] $1"; WARN=$((WARN + 1)); }
fail() { echo "  [FAIL] $1"; OK=1; }

echo ""
echo "============================================================"
echo " K_map scRNA-seq QC Pipeline - environment check"
echo "============================================================"
echo ""

if [[ -n "${JAVA_HOME:-}" && ! -x "${JAVA_HOME}/bin/java" ]]; then
    unset JAVA_HOME
fi
unset JAVA_CMD

echo "[1] Java 17+"
java17_ok() {
    "$1" -version 2>&1 | grep -qE 'version "(1[7-9]|2[0-6])'
}

if java17_ok java; then
    check "$(java -version 2>&1 | head -1)"
else
    JAVA_CANDIDATE=""
    for candidate in \
        /home/Data_Drive_8TB_2/CluVar/tool/jdk-17.0.12/bin/java \
        /home/Data_Drive_8TB_2/CluVar/tool/*/bin/java \
        /usr/lib/jvm/*/bin/java \
        /opt/*/bin/java \
        /usr/local/*/bin/java \
        "${JAVA_HOME:-}/bin/java"; do
        if [[ -x "${candidate}" ]] && java17_ok "${candidate}"; then
            JAVA_CANDIDATE="${candidate}"
            break
        fi
    done

    if [[ -n "${JAVA_CANDIDATE}" ]]; then
        export JAVA_HOME="$(dirname "$(dirname "${JAVA_CANDIDATE}")")"
        export PATH="${JAVA_HOME}/bin:${PATH}"
        check "$("${JAVA_CANDIDATE}" -version 2>&1 | head -1)"
    else
        fail "Java 17+ is required by the bundled/current Nextflow binary."
    fi
fi

echo ""
echo "[2] Nextflow"
if [[ -x "${PIPELINE_DIR}/nextflow" ]]; then
    check "local nextflow: ${PIPELINE_DIR}/nextflow"
elif [[ -x "${PIPELINE_DIR}/../nextflow" ]]; then
    check "parent nextflow: ${PIPELINE_DIR}/../nextflow"
elif command -v nextflow >/dev/null 2>&1; then
    check "system nextflow: $(nextflow -version 2>&1 | head -1)"
else
    warn "nextflow not found. Install it or place the binary in this directory."
fi

echo ""
echo "[3] Singularity"
if command -v singularity >/dev/null 2>&1; then
    check "$(singularity --version)"
else
    fail "Singularity is required."
fi

echo ""
echo "[4] Cell Ranger"
CELLRANGER_SIF="${CELLRANGER_SIF:-${PIPELINE_DIR}/containers/cellranger.sif}"
CELLRANGER_BIN="${CELLRANGER_BIN:-}"
if [[ -f "${CELLRANGER_SIF}" ]]; then
    check "cellranger.sif ($(du -sh "${CELLRANGER_SIF}" | cut -f1))"
    singularity exec "${CELLRANGER_SIF}" cellranger --version 2>/dev/null | head -1 || true
elif [[ -n "${CELLRANGER_BIN}" && -x "${CELLRANGER_BIN}" ]]; then
    warn "cellranger.sif missing; host fallback available with -profile local,host_cellranger"
    check "${CELLRANGER_BIN}"
    "${CELLRANGER_BIN}" --version 2>/dev/null | head -1 || true
else
    fail "Missing cellranger.sif. Build it with: bash ${PIPELINE_DIR}/def/build_containers.sh cellranger"
fi

echo ""
echo "[5] Containers"
if [[ -f "${PIPELINE_DIR}/containers/cellranger.sif" ]]; then
    check "cellranger.sif ($(du -sh "${PIPELINE_DIR}/containers/cellranger.sif" | cut -f1))"
else
    fail "Missing cellranger.sif"
fi

if [[ -f "${PIPELINE_DIR}/containers/cellbender.sif" ]]; then
    check "cellbender.sif ($(du -sh "${PIPELINE_DIR}/containers/cellbender.sif" | cut -f1))"
else
    fail "Missing cellbender.sif"
fi

QC_CONDA_ENV="${QC_CONDA_ENV:-/home/Data_Drive_8TB/cwonkim99/.conda/envs/solo}"
if [[ -f "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" ]]; then
    check "scrnaseq_qc.sif ($(du -sh "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" | cut -f1))"
elif [[ -x "${QC_CONDA_ENV}/bin/solo" && -x "${QC_CONDA_ENV}/bin/python" ]]; then
    check "scrnaseq_qc.sif missing, using local_conda fallback: ${QC_CONDA_ENV}"
else
    fail "Missing scrnaseq_qc.sif and no usable QC conda env at ${QC_CONDA_ENV}"
fi

if [[ ! -f "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" ]]; then
    echo ""
    echo "  For fully containerized execution, build/pull containers with:"
    echo "    bash ${PIPELINE_DIR}/def/build_containers.sh all"
    echo "  Current fallback execution profile:"
    echo "    bash run.sh --samplesheet input/samples.example.csv -profile local_conda"
fi

echo ""
echo "[6] Default files"
for file in \
    "${PIPELINE_DIR}/main.nf" \
    "${PIPELINE_DIR}/nextflow.config" \
    "${PIPELINE_DIR}/configs/solo_params.json" \
    "${PIPELINE_DIR}/scripts/qc_for_solo.py" \
    "${PIPELINE_DIR}/scripts/filter_solo_singlets.py" \
    "${PIPELINE_DIR}/input/samples.example.csv"; do
    if [[ -f "${file}" ]]; then
        check "${file#${PIPELINE_DIR}/}"
    else
        fail "Missing ${file#${PIPELINE_DIR}/}"
    fi
done

echo ""
echo "============================================================"
if [[ "${OK}" -eq 0 && "${WARN}" -eq 0 ]]; then
    echo "All checks passed."
elif [[ "${OK}" -eq 0 ]]; then
    echo "Warnings: ${WARN}. Pipeline may still run after review."
else
    echo "Failures found. Fix [FAIL] items before running."
fi
echo "============================================================"
echo ""

exit "${OK}"
