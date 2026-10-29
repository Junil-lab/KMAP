#!/usr/bin/env bash
# ============================================================
# K_map scRNA-seq preprocessing/QC pipeline runner
#
# Required:
#   --samplesheet input/samples.example.csv
#
# Common options:
#   --outdir DIR
#   --keep_work true|false  (default false: clean non-Cell-Ranger work after success)
#   --threads N
#   --use_gpu true|false
#   --cuda_device 0
#   --cellranger_sif containers/cellranger.sif
#   --transcriptome /path/to/refdata-gex-...
#   --mt_threshold 20  (optional; omit to keep all cells by mitochondrial percentage)
#   --hvg_n_top 2000
#   --solo_expected_doublet_rate 0.06
#   -resume
#   -profile local|local_conda|slurm
# ============================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NXF_HOME="${NXF_HOME:-${PIPELINE_DIR}/.nextflow_home}"
mkdir -p "${NXF_HOME}"

if [[ -n "${JAVA_HOME:-}" && ! -x "${JAVA_HOME}/bin/java" ]]; then
    unset JAVA_HOME
fi
unset JAVA_CMD

java17_ok() {
    "$1" -version 2>&1 | grep -qE 'version "(1[7-9]|2[0-6])'
}

if ! java17_ok java; then
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
    else
        echo "[ERROR] Java 17+ is required by this Nextflow binary." >&2
        echo "        Install OpenJDK 17 or set JAVA_HOME to a valid Java 17+ directory." >&2
        exit 1
    fi
fi

if [[ -x "${PIPELINE_DIR}/nextflow" ]]; then
    NXF="${PIPELINE_DIR}/nextflow"
elif [[ -x "${PIPELINE_DIR}/../nextflow" ]]; then
    NXF="${PIPELINE_DIR}/../nextflow"
elif command -v nextflow >/dev/null 2>&1; then
    NXF="nextflow"
else
    echo "[ERROR] nextflow was not found. Run setup.sh or install Nextflow." >&2
    exit 1
fi

MAIN="${PIPELINE_DIR}/main.nf"
CONFIG="${PIPELINE_DIR}/nextflow.config"

SAMPLESHEET=""
OUTDIR="${PIPELINE_DIR}/results"
PROFILE="local"
PROFILE_SET=false
KEEP_WORK=false
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep_work) KEEP_WORK="$2"; shift 2 ;;
        --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; EXTRA_ARGS+=("--outdir" "$2"); shift 2 ;;
        --threads) EXTRA_ARGS+=("--threads" "$2"); shift 2 ;;
        --memory_gb) EXTRA_ARGS+=("--memory_gb" "$2"); shift 2 ;;
        --use_gpu) EXTRA_ARGS+=("--use_gpu" "$2"); shift 2 ;;
        --cuda_device) EXTRA_ARGS+=("--cuda_device" "$2"); shift 2 ;;
        --cellranger_sif) EXTRA_ARGS+=("--cellranger_sif" "$2"); shift 2 ;;
        --cellranger_cmd) EXTRA_ARGS+=("--cellranger_cmd" "$2"); shift 2 ;;
        --cellranger_bin) EXTRA_ARGS+=("--cellranger_bin" "$2"); shift 2 ;;
        --transcriptome) EXTRA_ARGS+=("--transcriptome" "$2"); shift 2 ;;
        --mt_threshold) EXTRA_ARGS+=("--mt_threshold" "$2"); shift 2 ;;
        --min_genes) EXTRA_ARGS+=("--min_genes" "$2"); shift 2 ;;
        --min_counts) EXTRA_ARGS+=("--min_counts" "$2"); shift 2 ;;
        --min_cells_gene) EXTRA_ARGS+=("--min_cells_gene" "$2"); shift 2 ;;
        --hvg_n_top) EXTRA_ARGS+=("--hvg_n_top" "$2"); shift 2 ;;
        --cellbender_fpr) EXTRA_ARGS+=("--cellbender_fpr" "$2"); shift 2 ;;
        --cellbender_epochs) EXTRA_ARGS+=("--cellbender_epochs" "$2"); shift 2 ;;
        --solo_expected_doublet_rate) EXTRA_ARGS+=("--solo_expected_doublet_rate" "$2"); shift 2 ;;
        --solo_json) EXTRA_ARGS+=("--solo_json" "$2"); shift 2 ;;
        --solo_extra_args) EXTRA_ARGS+=("--solo_extra_args" "$2"); shift 2 ;;
        --qc_conda_env) EXTRA_ARGS+=("--qc_conda_env" "$2"); shift 2 ;;
        -profile) PROFILE="$2"; PROFILE_SET=true; shift 2 ;;
        -resume) EXTRA_ARGS+=("-resume"); shift ;;
        *) EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "${SAMPLESHEET}" ]]; then
    cat << EOF
Usage:
  bash run.sh --samplesheet input/samples.example.csv [options]

Required:
  --samplesheet CSV

Common:
  --outdir DIR
  --keep_work true|false  (default false: clean non-Cell-Ranger work after success)
  --threads N
  --use_gpu true|false
  --cuda_device 0
  --cellranger_sif containers/cellranger.sif
  --transcriptome /path/to/refdata-gex
  --mt_threshold 20  (optional; omitted by default)
  --hvg_n_top 2000
  --solo_expected_doublet_rate 0.06
  --qc_conda_env /path/to/solo_conda_env
  -resume
  -profile local|local_conda|slurm
EOF
    exit 1
fi

if [[ "${KEEP_WORK}" != "true" && "${KEEP_WORK}" != "false" ]]; then
    echo "[ERROR] --keep_work must be true or false" >&2
    exit 2
fi
if [[ "${KEEP_WORK}" = "false" ]]; then
    python3 -c 'import sys; sys.exit(sys.version_info < (3, 9))' || { echo "[ERROR] Python 3.9+ required for work cleanup" >&2; exit 2; }
fi

mkdir -p "${PIPELINE_DIR}/containers"

if [[ "${PROFILE_SET}" = "false" && ! -f "${PIPELINE_DIR}/containers/scrnaseq_qc.sif" ]]; then
    DEFAULT_QC_CONDA="/home/Data_Drive_8TB/cwonkim99/.conda/envs/solo"
    if [[ -x "${DEFAULT_QC_CONDA}/bin/solo" && -x "${DEFAULT_QC_CONDA}/bin/python" ]]; then
        PROFILE="local_conda"
    fi
fi

if [[ "${PROFILE}" != *"host_cellranger"* && ! -f "${PIPELINE_DIR}/containers/cellranger.sif" ]]; then
    echo "[WARN] Missing Cell Ranger SIF: ${PIPELINE_DIR}/containers/cellranger.sif" >&2
    echo "       Build it with: bash def/build_containers.sh cellranger" >&2
    echo "       If you intentionally want a host binary fallback, use -profile local,host_cellranger --cellranger_bin /path/to/cellranger" >&2
fi

mkdir -p "${OUTDIR}/pipeline_info" "${OUTDIR}/logs"

MEM_LOG="${OUTDIR}/logs/system_memory.tsv"
echo -e "timestamp\ttotal_MB\tused_MB\tfree_MB" > "${MEM_LOG}"
(
    while true; do
        TS=$(date '+%Y-%m-%d %H:%M:%S')
        MEM=$(free -m | awk '/^Mem:/{print $2"\t"$3"\t"$4}')
        echo -e "${TS}\t${MEM}" >> "${MEM_LOG}"
        sleep 5
    done
) &
MEM_PID=$!
trap 'kill "${MEM_PID}" 2>/dev/null || true' EXIT INT TERM

LOG_FILE="${OUTDIR}/pipeline_info/run_$(date '+%Y%m%d_%H%M%S').log"

echo "============================================================"
echo " K_map scRNA-seq QC Pipeline"
echo " $(date)"
echo " pipeline dir : ${PIPELINE_DIR}"
echo " samplesheet  : ${SAMPLESHEET}"
echo " outdir       : ${OUTDIR}"
echo " profile      : ${PROFILE}"
echo " mem log      : ${MEM_LOG}"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo " extra args   : ${EXTRA_ARGS[*]}"
echo "============================================================"

"${NXF}" run "${MAIN}" \
    -c "${CONFIG}" \
    --samplesheet "${SAMPLESHEET}" \
    -profile "${PROFILE}" \
    "${EXTRA_ARGS[@]}" \
    -with-trace "${LOG_FILE%.log}.trace.csv" \
    2>&1 | tee "${LOG_FILE}"

echo "============================================================"
if [[ "${KEEP_WORK}" = "false" ]]; then
    python3 "${PIPELINE_DIR}/scripts/cleanup_work.py" \
        --trace "${LOG_FILE%.log}.trace.csv" \
        --outdir "${OUTDIR}" --apply \
        2>&1 | tee "${LOG_FILE%.log}.cleanup.log"
fi

echo "Done: $(date)"
echo "Log: ${LOG_FILE}"
echo "Memory log: ${MEM_LOG}"
echo "============================================================"
