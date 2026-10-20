#!/usr/bin/env bash
# ============================================================
# K_map scRNA-seq QC Pipeline - container build script
#
# Usage:
#   bash def/build_containers.sh all
#   bash def/build_containers.sh cellranger
#   bash def/build_containers.sh qc
#   bash def/build_containers.sh cellbender
#
# Cell Ranger examples:
#   CELLRANGER_TARBALL=/path/to/cellranger-10.0.0.tar.gz bash def/build_containers.sh cellranger
#   CELLRANGER_DIR=/path/to/cellranger-10.0.0 bash def/build_containers.sh cellranger
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${PIPELINE_DIR}/containers"
MODE="${1:-all}"

mkdir -p "${OUT_DIR}"

SINGULARITY_BIN="${SINGULARITY_BIN:-singularity}"
LOCAL_BUILD_DIR="${LOCAL_BUILD_DIR:-/tmp/kmap_scrnaseq_qc_build}"
CELLBENDER_IMAGE="${CELLBENDER_IMAGE:-docker://us.gcr.io/broad-dsde-methods/cellbender:latest}"
CELLRANGER_VERSION="${CELLRANGER_VERSION:-10.0.0}"
CELLRANGER_TARBALL="${CELLRANGER_TARBALL:-}"
CELLRANGER_DIR="${CELLRANGER_DIR:-}"
CELLRANGER_URL="${CELLRANGER_URL:-https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VERSION}.tar.gz}"

mkdir -p "${LOCAL_BUILD_DIR}"

find_cellranger_source() {
    if [[ -n "${CELLRANGER_TARBALL}" && -f "${CELLRANGER_TARBALL}" ]]; then
        echo "tarball:${CELLRANGER_TARBALL}"
        return 0
    fi

    if [[ -n "${CELLRANGER_DIR}" && -x "${CELLRANGER_DIR}/cellranger" ]]; then
        echo "dir:${CELLRANGER_DIR}"
        return 0
    fi
    if [[ -n "${CELLRANGER_DIR}" && -x "${CELLRANGER_DIR}/bin/cellranger" ]]; then
        echo "dir:${CELLRANGER_DIR}"
        return 0
    fi

    local candidate
    for candidate in \
        "${PIPELINE_DIR}/downloads/cellranger-${CELLRANGER_VERSION}.tar.gz" \
        "${PIPELINE_DIR}/cellranger-${CELLRANGER_VERSION}.tar.gz" \
        "${SCRIPT_DIR}/cellranger-${CELLRANGER_VERSION}.tar.gz"; do
        if [[ -f "${candidate}" ]]; then
            echo "tarball:${candidate}"
            return 0
        fi
    done

    return 1
}

write_cellranger_def_from_tarball() {
    local source_tarball="$1"
    local def_file="$2"
    cat > "${def_file}" << EOF
Bootstrap: docker
From: ubuntu:22.04

%files
    ${source_tarball} /opt/cellranger.tar.gz

%post
    set -e
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \\
        ca-certificates wget curl tar gzip bzip2 xz-utils unzip \\
        procps util-linux libgomp1 libatomic1 zlib1g \\
        libglib2.0-0 libxrender1 libxtst6 libxi6 libxext6
    rm -rf /var/lib/apt/lists/*

    cd /opt
    tar -xzf /opt/cellranger.tar.gz
    rm -f /opt/cellranger.tar.gz
    CELLRANGER_ROOT=\$(find /opt -maxdepth 1 -type d -name 'cellranger-*' | sort -V | tail -1)
    test -n "\${CELLRANGER_ROOT}"
    ln -sfn "\${CELLRANGER_ROOT}" /opt/cellranger
    if [ -x /opt/cellranger/cellranger ]; then
        ln -sfn /opt/cellranger/cellranger /usr/local/bin/cellranger
    else
        ln -sfn /opt/cellranger/bin/cellranger /usr/local/bin/cellranger
    fi

%environment
    export PATH=/opt/cellranger:/opt/cellranger/bin:\$PATH

%runscript
    exec cellranger "\$@"

%test
    cellranger --version
EOF
}

write_cellranger_def_from_dir() {
    local source_dir="$1"
    local def_file="$2"
    cat > "${def_file}" << EOF
Bootstrap: docker
From: ubuntu:22.04

%files
    ${source_dir} /opt/cellranger

%post
    set -e
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \\
        ca-certificates wget curl tar gzip bzip2 xz-utils unzip \\
        procps util-linux libgomp1 libatomic1 zlib1g \\
        libglib2.0-0 libxrender1 libxtst6 libxi6 libxext6
    rm -rf /var/lib/apt/lists/*

    if [ -x /opt/cellranger/cellranger ]; then
        ln -sfn /opt/cellranger/cellranger /usr/local/bin/cellranger
    else
        ln -sfn /opt/cellranger/bin/cellranger /usr/local/bin/cellranger
    fi

%environment
    export PATH=/opt/cellranger:/opt/cellranger/bin:\$PATH

%runscript
    exec cellranger "\$@"

%test
    cellranger --version
EOF
}

build_cellranger() {
    local local_sif="${LOCAL_BUILD_DIR}/cellranger.sif"
    local def_file="${LOCAL_BUILD_DIR}/cellranger.def"
    local source_spec source_type source_path

    if ! source_spec="$(find_cellranger_source)"; then
        mkdir -p "${PIPELINE_DIR}/downloads"
        echo "[$(date)] Cell Ranger source was not found locally."
        echo "[$(date)] 10x may require browser/license-mediated download. Trying direct download once:"
        echo "[$(date)]   ${CELLRANGER_URL}"
        if wget -O "${PIPELINE_DIR}/downloads/cellranger-${CELLRANGER_VERSION}.tar.gz" "${CELLRANGER_URL}"; then
            source_spec="tarball:${PIPELINE_DIR}/downloads/cellranger-${CELLRANGER_VERSION}.tar.gz"
        else
            rm -f "${PIPELINE_DIR}/downloads/cellranger-${CELLRANGER_VERSION}.tar.gz"
            cat >&2 << EOF
[ERROR] Could not download Cell Ranger automatically.

Download Cell Ranger v${CELLRANGER_VERSION} from the official 10x Genomics page after accepting the license:
  https://www.10xgenomics.com/support/software/cell-ranger/latest

Then rerun one of:
  CELLRANGER_TARBALL=/path/to/cellranger-${CELLRANGER_VERSION}.tar.gz bash def/build_containers.sh cellranger
  CELLRANGER_DIR=/path/to/cellranger-${CELLRANGER_VERSION} bash def/build_containers.sh cellranger
EOF
            exit 1
        fi
    fi

    source_type="${source_spec%%:*}"
    source_path="${source_spec#*:}"
    case "${source_type}" in
        tarball)
            echo "[$(date)] Building Cell Ranger SIF from tarball: ${source_path}"
            write_cellranger_def_from_tarball "${source_path}" "${def_file}"
            ;;
        dir)
            echo "[$(date)] Building Cell Ranger SIF from installed directory: ${source_path}"
            write_cellranger_def_from_dir "${source_path}" "${def_file}"
            ;;
        *)
            echo "[ERROR] Unknown Cell Ranger source type: ${source_type}" >&2
            exit 1
            ;;
    esac

    "${SINGULARITY_BIN}" build --fakeroot "${local_sif}" "${def_file}"
    echo "[$(date)] Testing Cell Ranger container..."
    "${SINGULARITY_BIN}" test "${local_sif}"
    cp "${local_sif}" "${OUT_DIR}/cellranger.sif"
    echo "[$(date)] Wrote ${OUT_DIR}/cellranger.sif"
}

build_qc() {
    local local_sif="${LOCAL_BUILD_DIR}/scrnaseq_qc.sif"
    echo "[$(date)] Building Scanpy/SOLO container..."
    "${SINGULARITY_BIN}" build --fakeroot "${local_sif}" "${SCRIPT_DIR}/scrnaseq_qc.def"
    echo "[$(date)] Testing Scanpy/SOLO container..."
    "${SINGULARITY_BIN}" test "${local_sif}"
    cp "${local_sif}" "${OUT_DIR}/scrnaseq_qc.sif"
    echo "[$(date)] Wrote ${OUT_DIR}/scrnaseq_qc.sif"
}

pull_cellbender() {
    local local_sif="${LOCAL_BUILD_DIR}/cellbender.sif"
    echo "[$(date)] Pulling CellBender container: ${CELLBENDER_IMAGE}"
    "${SINGULARITY_BIN}" pull --force "${local_sif}" "${CELLBENDER_IMAGE}"
    cp "${local_sif}" "${OUT_DIR}/cellbender.sif"
    echo "[$(date)] Wrote ${OUT_DIR}/cellbender.sif"
}

case "${MODE}" in
    cellranger)
        build_cellranger
        ;;
    qc)
        build_qc
        ;;
    cellbender)
        pull_cellbender
        ;;
    all)
        build_cellranger
        build_qc
        pull_cellbender
        ls -lh "${OUT_DIR}"/*.sif
        ;;
    *)
        echo "Usage: bash def/build_containers.sh [all|cellranger|qc|cellbender]" >&2
        exit 1
        ;;
esac
