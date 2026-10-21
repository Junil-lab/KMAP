#!/usr/bin/env bash
set -euo pipefail

BASE="/home/Data_Drive_8TB_3/K_map/haplotype_pbcpg/method1_phased_bam"
THREADS=16

for TYPE in Normal Tumor; do
  OUT_DIR="$BASE/$TYPE/joint_segment"
  COLLECTION="$OUT_DIR/input_collection.tsv"
  mkdir -p "$OUT_DIR"

  # build the input collection TSV (header: identifier \t filename \t labels)
  if [ "$TYPE" = "Normal" ]; then
    cat > "$COLLECTION" <<EOF
identifier	filename	labels
m84039_230412_183900_s2	$BASE/Normal/m84039_230412_183900_s2	Normal
m84039_230415_005427_s4	$BASE/Normal/m84039_230415_005427_s4	Normal
EOF
  else
    cat > "$COLLECTION" <<EOF
identifier	filename	labels
m84039_230414_235240_s2	$BASE/Tumor/m84039_230414_235240_s2	Tumor
m84039_230415_002321_s3	$BASE/Tumor/m84039_230415_002321_s3	Tumor
EOF
  fi

  echo "[$TYPE joint-segment] start: $(date)"

  methbat joint-segment \
    --input-collection "$COLLECTION" \
    --output-prefix "$OUT_DIR/joint" \
    --threads "$THREADS" \
    --min-asm-abs-delta-mean 0.5 \
    --min-asm-frac 0.6

  echo "[$TYPE joint-segment] done: $(date)"
done

echo "[ALL DONE] output: $BASE/{Normal,Tumor}/joint_segment/"
