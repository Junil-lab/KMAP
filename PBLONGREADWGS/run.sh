#!/usr/bin/env bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — run wrapper
# Fully self-contained: all paths are set relative to this script
#
# Required arguments:
#   --samplesheet   path to the input CSV
#
# Optional arguments:
#   --outdir        output directory (default: <pipeline_dir>/results)
#   --threads       number of CPU threads (default: 16)
# Sex is set per sample with the samplesheet sex column (XX|XY) — used by TRGT/DeepVariant/Sawfish
#   --use_gpu       use GPU for DeepVariant true|false (default: true)
#   --vcf_dir       existing VCF directory when start_step=hiphase
#   -resume         reuse the Nextflow cache (restart an interrupted run)
#   -profile        local|slurm (default: local)
#
# Examples:
#   # default run (from pbmm2)
#   bash run.sh --samplesheet input/samples.csv
#
#   # set the output directory
#   bash run.sh --samplesheet input/samples.csv --outdir /data/results
#
#   # restart an interrupted run
#   bash run.sh --samplesheet input/samples.csv -resume
#
#   # restart from hiphase (start_step=hiphase in the samplesheet)
#   bash run.sh --samplesheet input/samples.csv --vcf_dir /data/results/HG002 -resume
#
#   # run without GPU, 32 threads
#   bash run.sh --samplesheet input/samples.csv --use_gpu false --threads 32
#
#   # run on a SLURM cluster
#   bash run.sh --samplesheet input/samples.csv -profile slurm
# ============================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Java 17+ — <pipeline_dir>/.java17 (installed by setup.sh) first, then the system java
if [[ -x "${PIPELINE_DIR}/.java17/bin/java" ]]; then
    export JAVA_HOME="${PIPELINE_DIR}/.java17"
    export PATH="${JAVA_HOME}/bin:${PATH}"
fi
if ! java -version 2>&1 | grep -qE 'version "(1[7-9]|[2-9][0-9])'; then
    echo "[ERROR] Java 17+ is required (Nextflow 26). Run: bash setup.sh   (installs it into .java17/)"
    exit 1
fi

# nextflow binary: wrapper in the repo first, otherwise PATH
if [[ -x "${PIPELINE_DIR}/nextflow" ]]; then
    NXF="${PIPELINE_DIR}/nextflow"
elif command -v nextflow &>/dev/null; then
    NXF="nextflow"
else
    echo "[ERROR] nextflow not found."
    echo "  Install: curl -s https://get.nextflow.io | bash && mv nextflow ${PIPELINE_DIR}/"
    exit 1
fi
MAIN="${PIPELINE_DIR}/main.nf"
CONFIG="${PIPELINE_DIR}/nextflow.config"

# ── Argument parsing ───────────────────────────────────────────
SAMPLESHEET=""
OUTDIR="${PIPELINE_DIR}/results"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2";      EXTRA_ARGS+=("--outdir" "$2"); shift 2 ;;
        --threads)     EXTRA_ARGS+=("--threads" "$2");    shift 2 ;;
        --use_gpu)     EXTRA_ARGS+=("--use_gpu" "$2");    shift 2 ;;
        --vcf_dir)     EXTRA_ARGS+=("--vcf_dir" "$2");    shift 2 ;;
        *)             EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "$SAMPLESHEET" ]]; then
    echo ""
    echo "Usage: bash run.sh --samplesheet <CSV> [options]"
    echo ""
    echo "Required:"
    echo "  --samplesheet   path to the input CSV"
    echo ""
    echo "Options:"
    echo "  --outdir        output directory (default: <pipeline_dir>/results)"
    echo "  --threads       number of CPU threads (default: 16)"
    echo "  --use_gpu       use GPU true|false (default: true)"
    echo "  --vcf_dir       existing VCF directory when restarting at hiphase"
    echo "  -resume         restart an interrupted run"
    echo "  -profile        local|slurm (default: local)"
    echo ""
    echo "Examples:"
    echo "  bash run.sh --samplesheet input/samples.csv"
    echo "  bash run.sh --samplesheet input/samples.csv --outdir /data/results -resume"
    echo "  bash run.sh --samplesheet input/samples.csv --use_gpu false --threads 32"
    exit 1
fi

mkdir -p "${OUTDIR}/pipeline_info"

# ── Check samplesheet BAM paths ──────────────────────────────
# BAMs must be inside the pipeline directory (e.g. <pipeline_dir>/input/bam/)
# for external paths, create a symlink: ln -s /external/path/sample.bam input/bam/
INPUT_BAM_DIR="${PIPELINE_DIR}/input/bam"
mkdir -p "${INPUT_BAM_DIR}"

BAM_ERROR=0
# find columns by header name (sample_id, bam, sex required)
HEADER="$(head -1 "${SAMPLESHEET}" | tr -d '\r')"
col_idx() { echo "$HEADER" | tr ',' '\n' | grep -nx "$1" | cut -d: -f1; }
C_ID=$(col_idx sample_id); C_BAM=$(col_idx bam); C_SEX=$(col_idx sex)
if [[ -z "$C_ID" || -z "$C_BAM" || -z "$C_SEX" ]]; then
    echo "[ERROR] samplesheet header must contain sample_id, bam and sex columns: ${HEADER}"
    exit 1
fi
while IFS= read -r line; do
    line="$(echo "$line" | tr -d '\r')"
    [[ -z "$line" ]] && continue
    sample_id=$(echo "$line" | cut -d, -f"$C_ID")
    bam=$(echo "$line" | cut -d, -f"$C_BAM")
    sex=$(echo "$line" | cut -d, -f"$C_SEX" | tr '[:lower:]' '[:upper:]')
    if [[ "$sex" != "XX" && "$sex" != "XY" ]]; then
        echo "[ERROR] ${sample_id}: sex is '${sex}' — must be XX or XY"
        BAM_ERROR=1
    fi
    if [[ ! -f "$bam" ]]; then
        echo "[ERROR] BAM file not found: ${bam}"
        BAM_ERROR=1
    elif [[ "$bam" != "${PIPELINE_DIR}"* ]]; then
        echo "[ERROR] BAM is outside the pipeline directory: ${bam}"
        echo "        → move it into input/bam/ or create a symlink:"
        echo "          ln -s ${bam} ${INPUT_BAM_DIR}/$(basename "$bam")"
        BAM_ERROR=1
    else
        echo "[OK]    ${sample_id} (${sex}): $(basename "$bam")"
    fi
done < <(tail -n +2 "${SAMPLESHEET}")

if [[ $BAM_ERROR -ne 0 ]]; then
    echo ""
    echo "[ERROR] check the BAM paths / sex values in the samplesheet."
    echo "        all BAMs must be inside ${INPUT_BAM_DIR}/."
    exit 1
fi

# ── Bind real locations of symlinked references ──────────────
# resources/ref/ files may be symlinks to another disk (e.g. a shared NAS);
# the containers only see projectDir, so bind the link targets too.
REF_BINDS=$(find "${PIPELINE_DIR}/resources/ref" -type l -exec readlink -f {} \; 2>/dev/null \
            | xargs -r -n1 dirname | grep -v "^${PIPELINE_DIR}" | sort -u | paste -sd, -)
if [[ -n "$REF_BINDS" ]]; then
    export KMAP_EXTRA_BINDS="${REF_BINDS}"   # read by nextflow.config (singularity.runOptions)
fi

# ── external memory monitor (every 5 s, for the whole run) ──
LOG_DIR="${OUTDIR}/logs"
mkdir -p "$LOG_DIR"
MEM_LOG="${LOG_DIR}/system_memory.tsv"
echo -e "timestamp\ttotal_MB\tused_MB\tfree_MB" > "$MEM_LOG"
(
    while true; do
        TS=$(date '+%Y-%m-%d %H:%M:%S')
        MEM=$(free -m | awk '/^Mem:/{print $2"\t"$3"\t"$4}')
        echo -e "${TS}\t${MEM}" >> "$MEM_LOG"
        sleep 5
    done
) &
MEM_PID=$!
trap "kill $MEM_PID 2>/dev/null || true" EXIT INT TERM

# ── Print run info ─────────────────────────────────────────────
echo "============================================================"
echo " K_map PacBio WGS Pipeline v2"
echo " $(date)"
echo " pipeline dir : ${PIPELINE_DIR}"
echo " samplesheet  : ${SAMPLESHEET}"
echo " outdir       : ${OUTDIR}"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo " Extra options: ${EXTRA_ARGS[*]}"
echo " mem log      : ${MEM_LOG}"
echo "============================================================"

# ── Run Nextflow ──────────────────────────────────────────────
LOG_FILE="${OUTDIR}/pipeline_info/run_$(date '+%Y%m%d_%H%M%S').log"

"$NXF" run "$MAIN" \
    -c "$CONFIG" \
    --samplesheet "$SAMPLESHEET" \
    -profile local \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee "$LOG_FILE"

echo "============================================================"
echo " Finished: $(date)"
echo " Log: ${LOG_FILE}"
echo " Memory log: ${MEM_LOG}"
echo "============================================================"
