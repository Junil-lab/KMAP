# Example run: HG002 (chr20 + chrX + chrY)

A small end-to-end test of the pipeline on HG002 (male, XY), used to check the sex-aware calling on chrX/chrY.

## Input

| | |
|---|---|
| Source | PacBio public Revio HG002-rep1 (2022Q4), `HG002.m84011_220902_175841_s1.GRCh38.bam` |
| Subset | primary reads on chr20, chrX, chrY (~5% of the genome), converted back to an unaligned HiFi BAM with MM/ML tags |
| Size | 2.2 GB, 370,018 reads — **not stored in git** (over GitHub's 100 MB file limit) |

Rebuild it with [input/download_input.sh](input/download_input.sh), then run:

```bash
bash input/download_input.sh /path/to/PBLONGREADWGS/input/bam
# edit the bam path in input/samplesheet.csv, then
taskset -c 0-4 bash ../../run.sh --samplesheet input/samplesheet.csv --outdir results --threads 5
```

## Output

[output/](output/) holds the results of that run, **except BAM files** (aligned, haplotagged and TRGT spanning BAMs are 0.1–2.6 GB each).

```
output/HG002_XY/
├── 02_aligned_bam/       flagstat
├── 03_small_variants/    DeepVariant VCF + gVCF
├── 04_structural_var/    Sawfish SV VCF, run stats
├── 05_tandem_repeats/    TRGT VCF
├── 06a_phased_bam/       HiPhase summary / blocks / stats
├── 06b_phased_SNP/       phased SNV VCF
├── 06c_phased_SV/        phased SV VCF
├── 06d_phased_TR/        phased TR VCF
└── 07_CpG_methylation/   combined / hap1 / hap2 CpG BED
output/qc_summary/, output/qc_summary.csv, output/pipeline_info/trace.csv
```

## Results

Run on 5 CPU cores + 1 GPU: **3 h 16 min** in total.

| Step | Time |
|---|---|
| PBMM2 | 33 min |
| DeepVariant | 2 h 29 min |
| Sawfish | 6 min |
| TRGT | 1.6 min |
| HiPhase | 3 min |
| CpG tools | 2 min |

Sex-aware calls (XY):

| Tool | chrX / chrY outside PAR | PAR |
|---|---|---|
| DeepVariant | chrX: 129,036 hom-alt, 0 het; chrY: 0 het | 6,981 het (diploid kept) |
| Sawfish | 1,101 haploid calls; the 14 diploid records are all BND | – |
| TRGT | all haploid (chrY: 3,529 called, 17 missing) | chrX PAR also haploid (TRGT `--karyotype XY` behaviour) |

Phasing (chr20): 88.5% of SNVs phased, NG50 1.02 Mb.
