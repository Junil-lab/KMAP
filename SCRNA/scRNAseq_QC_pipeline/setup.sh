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

# Nextflow 26 needs Java 17+. If none is found, OpenJDK 17 (Eclipse Temurin)
# is downloaded into <pipeline_dir>/.java17 — run.sh uses it automatically.
echo "[1] Java 17+"
JAVA_LOCAL="${PIPELINE_DIR}/.java17"
java17_ok() { "$1" -version 2>&1 | grep -qE 'version "(1[7-9]|[2-9][0-9])'; }
if [[ -x "${JAVA_LOCAL}/bin/java" ]] && java17_ok "${JAVA_LOCAL}/bin/java"; then
    check "local $("${JAVA_LOCAL}/bin/java" -version 2>&1 | head -1)  (${JAVA_LOCAL})"
elif command -v java >/dev/null 2>&1 && java17_ok java; then
    check "system $(java -version 2>&1 | head -1)"
else
    warn "Java 17+ not found — downloading OpenJDK 17 into ${JAVA_LOCAL}"
    case "$(uname -m)" in
        x86_64)  ARCH=x64 ;;
        aarch64) ARCH=aarch64 ;;
        *)       ARCH="" ;;
    esac
    URL="https://api.adoptium.net/v3/binary/latest/17/ga/linux/${ARCH}/jdk/hotspot/normal/eclipse"
    if [[ -n "${ARCH}" ]] && mkdir -p "${JAVA_LOCAL}" \
        && curl -fsSL "${URL}" | tar -xz -C "${JAVA_LOCAL}" --strip-components=1 \
        && java17_ok "${JAVA_LOCAL}/bin/java"; then
        check "installed $("${JAVA_LOCAL}/bin/java" -version 2>&1 | head -1)  (${JAVA_LOCAL})"
    else
        rm -rf "${JAVA_LOCAL}"
        fail "Java 17 download failed. Install it manually: conda create -p ${JAVA_LOCAL} -c conda-forge openjdk=17"
    fi
fi
if [[ -x "${JAVA_LOCAL}/bin/java" ]]; then
    export JAVA_HOME="${JAVA_LOCAL}"; export PATH="${JAVA_HOME}/bin:${PATH}"
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
    check "cellranger.sif ($(du -shL "${CELLRANGER_SIF}" | cut -f1))"
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
    check "cellranger.sif ($(du -shL "${PIPELINE_DIR}/containers/cellranger.sif" | cut -f1))"
else
    fail "Missing cellranger.sif"
fi

if [[ -f "${PIPELINE_DIR}/containers/cellbender.sif" ]]; then
    check "cellbender.sif ($(du -shL "${PIPELINE_DIR}/containers/cellbender.sif" | cut -f1))"
else
    fail "Missing cellbender.sif"
fi

QC_CONDA_ENV="${QC_CONDA_ENV:-}"   # optional: an existing conda env for -profile local_conda
if [[ -f "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" ]]; then
    check "scrnaseq_qc.sif ($(du -shL "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" | cut -f1))"
elif [[ -n "${QC_CONDA_ENV}" && -x "${QC_CONDA_ENV}/bin/solo" && -x "${QC_CONDA_ENV}/bin/python" ]]; then
    warn "scrnaseq_qc.sif missing; conda env usable with -profile local_conda --qc_conda_env ${QC_CONDA_ENV}"
else
    fail "Missing scrnaseq_qc.sif"
fi

if [[ ! -f "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" ]]; then
    echo ""
    echo "  Download the prebuilt QC container (GitHub Release, parts joined and checked):"
    echo "    bash ${PIPELINE_DIR}/def/fetch_qc_container.sh"
    echo "  or build it from def/scrnaseq_qc.def:"
    echo "    bash ${PIPELINE_DIR}/def/build_containers.sh qc"
    echo "  Or run Scanpy / SOLO from an existing conda environment:"
    echo "    bash run.sh --samplesheet input/samples.example.csv -profile local_conda --qc_conda_env /path/to/env"
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
