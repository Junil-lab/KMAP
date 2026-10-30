#!/usr/bin/env bash
# ============================================================
# K_map scRNA-seq QC Pipeline - download the prebuilt QC container
#
# scrnaseq_qc.sif (~4.5 GB) is larger than GitHub's 2 GB per-file release limit,
# so the release holds it in parts. This script downloads the parts, checks each
# against SHA256SUMS, joins them, checks the joined image and writes
# containers/scrnaseq_qc.sif. Interrupted downloads resume on the next run.
#
# Usage:
#   bash def/fetch_qc_container.sh                  # default release
#   bash def/fetch_qc_container.sh scrna-qc-v2.0.0  # a specific release tag
#   KMAP_RELEASE_REPO=owner/repo bash def/fetch_qc_container.sh
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)/containers"
REPO="${KMAP_RELEASE_REPO:-Junil-lab/KMAP}"
TAG="${1:-scrna-qc-v2.0.0}"
BASE_URL="https://github.com/${REPO}/releases/download/${TAG}"

mkdir -p "${OUT_DIR}/.download"
cd "${OUT_DIR}/.download"

echo "[$(date)] release ${REPO} ${TAG}"
curl -fsSL --retry 3 -o SHA256SUMS "${BASE_URL}/SHA256SUMS"

IMAGE=$(awk '$2 !~ /\.part[0-9]+$/ {print $2}' SHA256SUMS)
PARTS=$(awk '$2 ~ /\.part[0-9]+$/ {print $2}' SHA256SUMS | sort)
[[ -n "${IMAGE}" && -n "${PARTS}" ]] || { echo "[ERROR] unexpected SHA256SUMS content" >&2; exit 1; }

for part in ${PARTS}; do
    if [[ -f "${part}" ]] && grep " ${part}\$" SHA256SUMS | sha256sum -c --status -; then
        echo "  [OK] ${part} (already downloaded)"
        continue
    fi
    echo "  downloading ${part}..."
    curl -fL --retry 5 -C - -o "${part}" "${BASE_URL}/${part}" || \
        curl -fL --retry 5 -o "${part}" "${BASE_URL}/${part}"
    grep " ${part}\$" SHA256SUMS | sha256sum -c -
done

echo "[$(date)] joining parts into ${IMAGE}..."
cat ${PARTS} > "${IMAGE}"
grep " ${IMAGE}\$" SHA256SUMS | sha256sum -c -

mv "${IMAGE}" "${OUT_DIR}/scrnaseq_qc.sif"
rm -f ${PARTS} SHA256SUMS
cd "${OUT_DIR}" && rmdir --ignore-fail-on-non-empty .download 2>/dev/null || true
echo "[$(date)] saved: ${OUT_DIR}/scrnaseq_qc.sif"
