#!/usr/bin/env bash
# Runs real pbmm2, IsoQuant, SQANTI3/ORF/HTML reporting, and merged QC on synthetic data.
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SMOKE_DIR="$(realpath -m "${1:-$PIPELINE_DIR/test_runs/smoke}")"
mkdir -p "$SMOKE_DIR/data"
singularity exec --cleanenv "$PIPELINE_DIR/containers/kmap_isoseq_main.sif" \
    python "$PIPELINE_DIR/tests/make_smoke_data.py" "$SMOKE_DIR/data"
cd "$PIPELINE_DIR"
NXF_WORK="$SMOKE_DIR/work" bash "$PIPELINE_DIR/run.sh" \
    --samplesheet "$SMOKE_DIR/data/samples.csv" --outdir "$SMOKE_DIR/results" \
    --threads 2 -c "$PIPELINE_DIR/tests/smoke.config" \
    --genome_fa "$SMOKE_DIR/data/genome.fa" --gencode_gtf "$SMOKE_DIR/data/annotation.gtf" \
    --polya_motif_list "$SMOKE_DIR/data/polyA_motifs.txt" \
    --polya_peak_bed "$SMOKE_DIR/data/polya.bed" --cage_peak_bed "$SMOKE_DIR/data/cage.bed"
python3 "$PIPELINE_DIR/tests/check_smoke_outputs.py" "$SMOKE_DIR/results"
