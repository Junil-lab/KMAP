#!/usr/bin/env bash
# ============================================================
# K_map PacBio WGS Pipeline v2 — 실행 래퍼
# 완전 독립형: 모든 경로는 스크립트 위치 기준으로 자동 설정됨
#
# 필수 인자:
#   --samplesheet   입력 CSV 파일 경로
#
# 선택 인자:
#   --outdir        결과 출력 디렉토리 (기본: <pipeline_dir>/results)
#   --threads       CPU 스레드 수 (기본: 16)
#   --karyotype     TRGT 성별 설정 XY|XX (기본: XY)
#   --use_gpu       DeepVariant GPU 사용 true|false (기본: true)
#   --vcf_dir       start_step=hiphase 일 때 기존 VCF 경로
#   -resume         Nextflow 캐시 재사용 (중단된 파이프라인 재시작)
#   -profile        local|slurm (기본: local)
#
# 사용 예:
#   # 기본 실행 (pbmm2부터)
#   bash run.sh --samplesheet input/samples.csv
#
#   # 결과 경로 지정
#   bash run.sh --samplesheet input/samples.csv --outdir /data/results
#
#   # 중단된 파이프라인 재시작
#   bash run.sh --samplesheet input/samples.csv -resume
#
#   # hiphase 단계부터 재시작 (samplesheet의 start_step=hiphase)
#   bash run.sh --samplesheet input/samples.csv --vcf_dir /data/results/HG002 -resume
#
#   # GPU 없이 실행, 스레드 32개
#   bash run.sh --samplesheet input/samples.csv --use_gpu false --threads 32
#
#   # SLURM 클러스터 실행
#   bash run.sh --samplesheet input/samples.csv -profile slurm
# ============================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Java 17+ — 시스템 java 우선 사용, 없으면 conda 환경에서 탐색
if ! java -version 2>&1 | grep -qE 'version "1[1-9]|version "[2-9][0-9]'; then
    # 시스템 java가 없거나 버전 미달 → conda 환경 탐색
    JAVA_CANDIDATES=$(find /home /opt /usr/local -name "java" -path "*/bin/java" 2>/dev/null \
                      | xargs -I{} sh -c '{} -version 2>&1 | grep -qE "version \"1[1-9]|version \"[2-9][0-9]" && echo {}' 2>/dev/null \
                      | head -1)
    if [[ -n "$JAVA_CANDIDATES" ]]; then
        export JAVA_HOME="$(dirname "$(dirname "$JAVA_CANDIDATES")")"
        export PATH="${JAVA_HOME}/bin:${PATH}"
    else
        echo "[ERROR] Java 11+ 가 필요합니다. java를 설치하거나 JAVA_HOME을 설정하세요."
        exit 1
    fi
fi

# nextflow 바이너리: repo 내 래퍼 우선, 없으면 시스템 PATH에서 탐색
if [[ -x "${PIPELINE_DIR}/nextflow" ]]; then
    NXF="${PIPELINE_DIR}/nextflow"
elif command -v nextflow &>/dev/null; then
    NXF="nextflow"
else
    echo "[ERROR] nextflow를 찾을 수 없습니다."
    echo "  설치: curl -s https://get.nextflow.io | bash && mv nextflow ${PIPELINE_DIR}/"
    exit 1
fi
MAIN="${PIPELINE_DIR}/main.nf"
CONFIG="${PIPELINE_DIR}/nextflow.config"

# ── 인자 파싱 ──────────────────────────────────────────────────
SAMPLESHEET=""
OUTDIR="${PIPELINE_DIR}/results"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2";      EXTRA_ARGS+=("--outdir" "$2"); shift 2 ;;
        --threads)     EXTRA_ARGS+=("--threads" "$2");    shift 2 ;;
        --karyotype)   EXTRA_ARGS+=("--karyotype" "$2");  shift 2 ;;
        --use_gpu)     EXTRA_ARGS+=("--use_gpu" "$2");    shift 2 ;;
        --vcf_dir)     EXTRA_ARGS+=("--vcf_dir" "$2");    shift 2 ;;
        *)             EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "$SAMPLESHEET" ]]; then
    echo ""
    echo "사용법: bash run.sh --samplesheet <CSV> [옵션]"
    echo ""
    echo "필수:"
    echo "  --samplesheet   입력 CSV 파일 경로"
    echo ""
    echo "선택 옵션:"
    echo "  --outdir        결과 출력 디렉토리 (기본: <pipeline_dir>/results)"
    echo "  --threads       CPU 스레드 수 (기본: 16)"
    echo "  --karyotype     TRGT 성별 XY|XX (기본: XY)"
    echo "  --use_gpu       GPU 사용 true|false (기본: true)"
    echo "  --vcf_dir       hiphase 재시작 시 기존 VCF 경로"
    echo "  -resume         중단된 파이프라인 재시작"
    echo "  -profile        local|slurm (기본: local)"
    echo ""
    echo "예시:"
    echo "  bash run.sh --samplesheet input/samples.csv"
    echo "  bash run.sh --samplesheet input/samples.csv --outdir /data/results -resume"
    echo "  bash run.sh --samplesheet input/samples.csv --use_gpu false --threads 32"
    exit 1
fi

mkdir -p "${OUTDIR}/pipeline_info"

# ── samplesheet BAM 경로 검증 ────────────────────────────────
# BAM은 반드시 <pipeline_dir>/input/bam/ 안에 있어야 함
# 외부 경로 사용 시: ln -s /외부/경로/sample.bam input/bam/ 으로 심링크 생성
INPUT_BAM_DIR="${PIPELINE_DIR}/input/bam"
mkdir -p "${INPUT_BAM_DIR}"

BAM_ERROR=0
while IFS=',' read -r sample_id bam start_step; do
    [[ "$sample_id" == "sample_id" ]] && continue
    bam="$(echo "$bam" | tr -d '\r')"
    if [[ ! -f "$bam" ]]; then
        echo "[ERROR] BAM 파일 없음: ${bam}"
        BAM_ERROR=1
    elif [[ "$bam" != "${PIPELINE_DIR}"* ]]; then
        echo "[ERROR] BAM이 pipeline 폴더 밖에 있습니다: ${bam}"
        echo "        → input/bam/ 안으로 이동하거나 심링크를 생성하세요:"
        echo "          ln -s ${bam} ${INPUT_BAM_DIR}/$(basename "$bam")"
        BAM_ERROR=1
    else
        echo "[OK]    ${sample_id}: $(basename "$bam")"
    fi
done < "${SAMPLESHEET}"

if [[ $BAM_ERROR -ne 0 ]]; then
    echo ""
    echo "[ERROR] samplesheet의 BAM 경로를 확인하세요."
    echo "        모든 BAM은 ${INPUT_BAM_DIR}/ 안에 있어야 합니다."
    exit 1
fi

# ── 외부 메모리 모니터 (5초 간격, 파이프라인 실행 중 내내 기록) ──
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

# ── 실행 정보 출력 ─────────────────────────────────────────────
echo "============================================================"
echo " K_map PacBio WGS Pipeline v2"
echo " $(date)"
echo " pipeline dir : ${PIPELINE_DIR}"
echo " samplesheet  : ${SAMPLESHEET}"
echo " outdir       : ${OUTDIR}"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo " 추가 옵션    : ${EXTRA_ARGS[*]}"
echo " mem log      : ${MEM_LOG}"
echo "============================================================"

# ── Nextflow 실행 ──────────────────────────────────────────────
LOG_FILE="${OUTDIR}/pipeline_info/run_$(date '+%Y%m%d_%H%M%S').log"

"$NXF" run "$MAIN" \
    -c "$CONFIG" \
    --samplesheet "$SAMPLESHEET" \
    -profile local \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee "$LOG_FILE"

echo "============================================================"
echo " 완료: $(date)"
echo " 로그: ${LOG_FILE}"
echo " 메모리 로그: ${MEM_LOG}"
echo "============================================================"
