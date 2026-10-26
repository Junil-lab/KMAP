#!/usr/bin/env bash
# Run after run_smoke.sh using its synthetic inputs and reusable task cache.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(realpath -e "$1")"
cd "$PIPELINE_DIR"
for mode in single multi both; do
    NXF_WORK="$TEST_DIR/work" bash run.sh \
        --samplesheet "$TEST_DIR/data/samples.csv" --outdir "$TEST_DIR/$mode" \
        --methods isocall,flair --analysis_mode "$mode" --threads 2 \
        -c tests/smoke.config -resume \
        --genome_fa "$TEST_DIR/data/genome.fa" --gencode_gtf "$TEST_DIR/data/annotation.gtf" \
        --polya_motif_list "$TEST_DIR/data/polyA_motifs.txt" \
        --polya_peak_bed "$TEST_DIR/data/polya.bed" --cage_peak_bed "$TEST_DIR/data/cage.bed"
    python3 tests/check_mode_outputs.py "$TEST_DIR/$mode" "$mode"
done
