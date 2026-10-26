#!/usr/bin/env bash
# One-row synthetic sample sheet: the same form as a user's single-sample CSV.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(realpath -e "${1:-$PIPELINE_DIR/test_runs/smoke}")"
cd "$PIPELINE_DIR"
python3 - "$TEST_DIR/data" <<'PY'
import csv,sys
from pathlib import Path
p=Path(sys.argv[1])
with (p/'samples.csv').open() as f:
 r=csv.DictReader(f);fields=r.fieldnames;row=next(r)
with (p/'sample1.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerow(row)
PY
NXF_WORK="$TEST_DIR/work" bash run.sh \
    --samplesheet "$TEST_DIR/data/sample1.csv" --outdir "$TEST_DIR/one_sample" \
    --threads 2 -c tests/smoke.config -resume \
    --genome_fa "$TEST_DIR/data/genome.fa" --gencode_gtf "$TEST_DIR/data/annotation.gtf" \
    --polya_motif_list "$TEST_DIR/data/polyA_motifs.txt" \
    --polya_peak_bed "$TEST_DIR/data/polya.bed" --cage_peak_bed "$TEST_DIR/data/cage.bed"
python3 tests/check_selection_outputs.py "$TEST_DIR/one_sample" auto isoquant,isocall,flair SmokeA
