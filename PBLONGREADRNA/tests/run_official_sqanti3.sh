#!/usr/bin/env bash
# Compare an official image against the completed default smoke run.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIF="$(realpath -e "$1")"
BASE="$(realpath -e "$2")"
OUT="$(realpath -m "$3")"
cd "$PIPELINE_DIR"
NXF_WORK="$BASE/work" bash run.sh \
    --samplesheet "$BASE/data/samples.csv" --outdir "$OUT" \
    --sqanti3_sif "$SIF" --sqanti3_qc auto \
    --threads 2 -c tests/smoke.config -resume \
    --genome_fa "$BASE/data/genome.fa" --gencode_gtf "$BASE/data/annotation.gtf" \
    --polya_motif_list "$BASE/data/polyA_motifs.txt" \
    --polya_peak_bed "$BASE/data/polya.bed" --cage_peak_bed "$BASE/data/cage.bed"
python3 tests/check_smoke_outputs.py "$OUT"
python3 tests/compare_sqanti3_outputs.py "$BASE/results" "$OUT"
