#!/usr/bin/env bash
# Run from any directory. Relative BAM paths are relative to the CSV directory.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLESHEET=''
OUTDIR="${PIPELINE_DIR}/results"
THREADS=16
PROFILE=local
CHECK_ONLY=0
EXTRA_ARGS=()
GLOBAL_ARGS=()
# Override parameters are validated using the effective Nextflow configuration.
usage() {
    cat <<'HELP'
Usage: bash run.sh --samplesheet input/samples.csv [options]
  --methods LIST          Comma-separated isoquant,isocall,flair (default: all three)
  --analysis_mode MODE    auto (default), single, multi, or both
                          single: independent analysis of each CSV row
                          multi: cohort isocall/FLAIR; requires >=2 samples
                          both: sample and cohort outputs; requires >=2 samples
                          auto: sample outputs plus cohort when >=2 samples
                          IsoQuant always runs per sample when selected.
  --threads N             CPUs per task (default: 16)
  --outdir PATH           Results directory (use a new directory for each selection)
  -profile local|slurm    Execution profile (default: local)
  -resume                 Reuse matching retained task results
  --check-only            Validate inputs and containers without starting analysis
Examples:
  bash run.sh --samplesheet input/sample1.csv --analysis_mode single
  bash run.sh --samplesheet input/samples.csv --methods isocall --analysis_mode multi
  bash run.sh --samplesheet input/samples.csv --methods flair --analysis_mode single
  bash run.sh --samplesheet input/samples.csv --methods isocall,flair --analysis_mode both
HELP
}
while (( $# )); do
    case "$1" in
        --*=*) set -- "${1%%=*}" "${1#*=}" "${@:2}"; continue ;;
        -c|-C)
            (( $# >= 2 )) && [[ -f "$2" ]] || { echo "ERROR: config file required for $1" >&2; exit 2; }
            GLOBAL_ARGS+=("$1" "$(realpath -e "$2")"); shift 2 ;;
        --samplesheet|--outdir|--threads|-profile)
            (( $# >= 2 )) && [[ -n "$2" && "$2" != -* ]] || { echo "ERROR: missing value for $1" >&2; exit 2; }
            case "$1" in
                --samplesheet) SAMPLESHEET="$2" ;;
                --outdir) OUTDIR="$2" ;;
                --threads) THREADS="$2" ;;
                -profile) PROFILE="$2" ;;
            esac
            shift 2 ;;
        --methods|--analysis_mode)
            (( $# >= 2 )) && [[ -n "$2" && "$2" != -* ]] || { echo "ERROR: missing value for $1" >&2; exit 2; }
            EXTRA_ARGS+=("$1" "$2"); shift 2 ;;
        --check-only) CHECK_ONLY=1; shift ;;
        --help|-h) usage; exit 0 ;;
        *) EXTRA_ARGS+=("$1"); shift ;;
    esac
done
[[ -n "$SAMPLESHEET" ]] || { usage >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'ERROR: threads must be a positive integer' >&2; exit 2; }
command -v python3 >/dev/null || { echo 'ERROR: python3 is required' >&2; exit 1; }
command -v singularity >/dev/null || { echo 'ERROR: singularity is required' >&2; exit 1; }
source "$PIPELINE_DIR/bin/select_java.sh"
select_java || { echo 'ERROR: Java 17+ required; run bash setup.sh' >&2; exit 1; }
if [[ -x "$PIPELINE_DIR/nextflow" ]]; then NXF="$PIPELINE_DIR/nextflow";
elif command -v nextflow >/dev/null; then NXF="$(command -v nextflow)";
else echo 'ERROR: nextflow is required; install a pinned release in the pipeline directory' >&2; exit 1; fi
export NXF_VER=26.04.3
SAMPLESHEET="$(realpath -e "$SAMPLESHEET")"
OUTDIR="$(realpath -m "$OUTDIR")"
python3 "$PIPELINE_DIR/bin/validate_samples.py" "$SAMPLESHEET"
# Extract effective params including -c/-params-file and CLI overrides.
CHECK_DIR="$(mktemp -d)"
trap 'rm -rf "$CHECK_DIR"' EXIT
"$NXF" "${GLOBAL_ARGS[@]}" config "$PIPELINE_DIR" -profile "$PROFILE" -flat > "$CHECK_DIR/effective.config"
python3 "$PIPELINE_DIR/bin/check_config.py" "$CHECK_DIR/effective.config" "$CHECK_DIR/containers.txt" --samplesheet "$SAMPLESHEET" "${EXTRA_ARGS[@]}"
mapfile -t IMAGES < "$CHECK_DIR/containers.txt"
bash "$PIPELINE_DIR/bin/check_containers.sh" "${IMAGES[@]}"
if (( CHECK_ONLY )); then echo 'Preflight passed; no analysis launched'; exit 0; fi
mkdir -p "$OUTDIR/pipeline_info"
# Keep the original CSV path unchanged for stable resume; the workflow resolves BAMs.
RUN_LOG="$OUTDIR/pipeline_info/run_$(date +%Y%m%d_%H%M%S).log"
"$NXF" "${GLOBAL_ARGS[@]}" run "$PIPELINE_DIR/main.nf" -profile "$PROFILE" \
    --samplesheet "$SAMPLESHEET" --outdir "$OUTDIR" --threads "$THREADS" \
    "${EXTRA_ARGS[@]}" 2>&1 | tee "$RUN_LOG"
echo "Pipeline completed: $OUTDIR"
