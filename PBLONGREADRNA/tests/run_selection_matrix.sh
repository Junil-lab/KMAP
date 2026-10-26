#!/usr/bin/env bash
# Run after run_smoke.sh. Uses real tools or their verified matching cached tasks.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(realpath -e "${1:-$PIPELINE_DIR/test_runs/smoke}")"
cd "$PIPELINE_DIR"
mkdir -p "$TEST_DIR/selection"
for case in \
    single:isoquant single:isocall single:flair \
    single:isoquant,isocall single:isoquant,flair single:isocall,flair \
    single:isoquant,isocall,flair \
    multi:isocall multi:flair multi:isoquant,isocall,flair \
    both:isoquant,isocall,flair; do
    mode="${case%%:*}"; methods="${case#*:}"
    name="${mode}_${methods//,/_}"
    out="$TEST_DIR/selection/$name"
    echo "TEST: $name"
    NXF_WORK="$TEST_DIR/work" bash run.sh \
        --samplesheet "$TEST_DIR/data/samples.csv" --outdir "$out" \
        --methods "$methods" --analysis_mode "$mode" --threads 2 \
        -c tests/smoke.config -resume \
        --genome_fa "$TEST_DIR/data/genome.fa" --gencode_gtf "$TEST_DIR/data/annotation.gtf" \
        --polya_motif_list "$TEST_DIR/data/polyA_motifs.txt" \
        --polya_peak_bed "$TEST_DIR/data/polya.bed" --cage_peak_bed "$TEST_DIR/data/cage.bed" \
        > "$TEST_DIR/selection/$name.log" 2>&1
    python3 tests/check_selection_outputs.py "$out" "$mode" "$methods"
done
