#!/usr/bin/env bash
# Build the example input: HG002 HiFi reads for chr20 + chrX + chrY as an unaligned BAM
# (with MM/ML 5mC tags), streamed from PacBio's public Revio HG002 data (needs internet).
#
# Usage (from the PBLONGREADWGS folder, using samtools inside the main container):
#   singularity exec -B $PWD containers/kmap_pbwgs_v2.sif \
#     bash example/HG002_chr20XY/input/download_input.sh input/bam
#
# Output:
#   <outdir>/HG002.chr20XY.hifi_reads.bam      (~2.2 GB, 370,018 reads, md5 36765a3750cc9f0578d43970722aa3b7)
#   <outdir>/HG002_chr20XY.samplesheet.csv     (ready to pass to run.sh)
set -euo pipefail
OUT=${1:-.}
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)
U=https://downloads.pacbcloud.com/public/revio/2022Q4/HG002-rep1/analysis/HG002.m84011_220902_175841_s1.GRCh38.bam
STEP=5000000
declare -A LEN=([chr20]=64444167 [chrX]=156040895 [chrY]=57227415)   # GRCh38 lengths
mkdir -p "$OUT/chunks"; cd "$OUT"
n=0
for chr in chr20 chrX chrY; do
  for ((s=0; s<${LEN[$chr]}; s+=STEP)); do
    e=$((s+STEP)); [ $e -gt ${LEN[$chr]} ] && e=${LEN[$chr]}
    out=chunks/$(printf "%03d" $n).${chr}_${s}.bam; n=$((n+1))
    [ -f $out.ok ] && continue
    ok=0
    for try in 1 2 3 4 5; do
      # primary reads whose start falls in this chunk only → no duplicates at chunk boundaries
      if samtools view -@4 -b -F 0x900 -e "pos>$s" -o $out $U ${chr}:$((s+1))-$e 2>chunks/err.log && samtools quickcheck $out; then ok=1; break; fi
      echo "retry $try $chr:$s $(tail -1 chunks/err.log)"; sleep 20
    done
    [ $ok -eq 1 ] || { echo "CHUNK FAILED $chr:$s"; exit 1; }
    touch $out.ok; echo "ok $out"
  done
done
samtools cat -o sub.aln.bam $(ls chunks/*.bam | sort)
# drop alignment and phasing/kinetics tags → unaligned HiFi BAM that keeps MM/ML
samtools reset -@4 -x HP,PS,ip,pw,fi,ri,fp,rp,MC,mc,rm --no-PG -o HG002.chr20XY.hifi_reads.bam sub.aln.bam
samtools quickcheck -u HG002.chr20XY.hifi_reads.bam
echo "reads $(samtools view -c -@4 HG002.chr20XY.hifi_reads.bam)"
rm -rf sub.aln.bam chunks "$(basename "$U").bai"   # .bai = remote index cached by samtools
cat > HG002_chr20XY.samplesheet.csv <<CSV
sample_id,bam,sex,start_step
HG002_XY,$OUT/HG002.chr20XY.hifi_reads.bam,XY,pbmm2
CSV
echo "samplesheet: $OUT/HG002_chr20XY.samplesheet.csv"
