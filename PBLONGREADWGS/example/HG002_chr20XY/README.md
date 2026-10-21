# Example run: HG002 (chr20 + chrX + chrY)

A small end-to-end test of the pipeline on HG002 (male, XY). Use it to

- check that your installation works before running real samples (~3 h on 5 cores instead of a day or more for a whole genome),
- see what every output file looks like ([output/](output/)),
- check the sex-aware calling on chrX/chrY.

## Input

| | |
|---|---|
| Source | PacBio public Revio HG002-rep1 (2022Q4), `HG002.m84011_220902_175841_s1.GRCh38.bam` on `downloads.pacbcloud.com` |
| Subset | primary reads on chr20, chrX, chrY (~5% of the genome), converted back to an unaligned HiFi BAM with MM/ML (5mC) tags |
| Size | 2.2 GB, 370,018 reads — **not stored in git** (over GitHub's file size limit); rebuild it with [input/download_input.sh](input/download_input.sh) |

`download_input.sh` streams only the chr20/X/Y regions of the public BAM over HTTPS with `samtools` (the whole-genome file is never downloaded), 5 Mb at a time with retries, then strips the alignment, phasing and kinetics tags so the result looks like a HiFi BAM straight from the sequencer. It needs an internet connection; `samtools` comes from the main container.

## How to run it

All commands are run from the `PBLONGREADWGS` folder of a fresh clone.

### 0. Set up the pipeline (once)

Follow [Running → 1](../../README.md#1-check-the-environment-once) in the main README: download the two containers, put the GRCh38 reference in `resources/ref/`, and run `bash setup.sh` until every line is `[OK]`.

### 1. Download the input (~2.2 GB)

```bash
singularity exec -B $PWD containers/kmap_pbwgs_v2.sif \
  bash example/HG002_chr20XY/input/download_input.sh input/bam
```

This writes

```
input/bam/HG002.chr20XY.hifi_reads.bam       # 370,018 reads
input/bam/HG002_chr20XY.samplesheet.csv      # sample_id=HG002_XY, sex=XY, start_step=pbmm2
```

It takes about 30 min (57 chunks streamed from PacBio's server; 28 min in our test). The result is byte-identical to the BAM used for the reference run (md5 `36765a3750cc9f0578d43970722aa3b7`). If it stops (network error), run the same command again: finished chunks are kept and skipped.

**On the K-MAP lab server** the same BAM is already on the NAS (md5 `36765a3750cc9f0578d43970722aa3b7`), so skip the download and link it:

```bash
mkdir -p input/bam
ln -s /home/NAS_2/K_MAP/KMAP_Pipeline_example_data/HG002_chr20XY/HG002.chr20XY.hifi_reads.bam input/bam/
sed "s#PIPELINE_DIR#$PWD#" /home/NAS_2/K_MAP/KMAP_Pipeline_example_data/HG002_chr20XY/HG002_chr20XY.samplesheet.csv \
  > input/bam/HG002_chr20XY.samplesheet.csv
```

### 2. Run the pipeline

```bash
taskset -c 0-4 bash run.sh \
  --samplesheet input/bam/HG002_chr20XY.samplesheet.csv \
  --outdir results/HG002_example --threads 5
# no GPU: add --use_gpu false (DeepVariant becomes much slower)
```

Reference run on 5 CPU cores + 1 GPU: **3 h 16 min**.

| Step | Time |
|---|---|
| PBMM2 | 33 min |
| DeepVariant | 2 h 29 min |
| Sawfish | 6 min |
| TRGT | 1.6 min |
| HiPhase | 3 min |
| CpG tools | 2 min |

### 3. Compare with the reference results

Your results should match [output/](output/) (same tool versions → same numbers; small differences in DeepVariant are possible on a different GPU/CPU). Quick check:

```bash
R=results/HG002_example/HG002_XY
cat results/HG002_example/qc_summary.csv
singularity exec -B $PWD containers/kmap_pbwgs_v2.sif bash -c "
  bcftools view -H -f PASS $R/03_small_variants/HG002_XY.deepvariant.vcf.gz chr20 | wc -l
  bcftools view -H -f PASS $R/04_structural_var/HG002_XY.sv.vcf.gz chr20 chrX chrY | wc -l
  bcftools view -H $R/05_tandem_repeats/HG002_XY.trgt.vcf.gz chr20 chrX chrY | awk '\$10!~/^\./' | wc -l
"
```

| Check | Expected |
|---|---|
| `qc_summary.csv`: total_reads / total_sv / phased_variants / phase_blocks / block_n50 | 399,105 / 3,088 / 146,382 / 428 / 1,021,679 |
| DeepVariant PASS variants on chr20 | 127,155 |
| Sawfish PASS SVs on chr20/X/Y | 2,669 |
| TRGT loci with a genotype on chr20/X/Y | 50,705 |
| chr20 phasing (`06a_phased_bam/*.summary.tsv`) | 53,814 of 60,808 het SNVs phased (88.5%), NG50 1.02 Mb |
| CpGs on chr20 with coverage ≥ 10 (`07_CpG_methylation/*.combined.bed.gz`) | 747,520 |

`total_reads` (399,105) is larger than the input (370,018) because it counts supplementary alignments.

### 4. Check the sex-aware calls (XY)

```bash
singularity exec -B $PWD containers/kmap_pbwgs_v2.sif bash -c "
  # chrX outside the PARs: expect no heterozygous (0/1, 1/2) calls
  # reference run: 129,036 x 1/1, 627 x 2/2, 6 x 3/3, 20 x 0/0
  bcftools view -H -f PASS -r chrX:2781480-155701382 $R/03_small_variants/HG002_XY.deepvariant.vcf.gz \
    | cut -f10 | cut -d: -f1 | sort | uniq -c
"
```

| Tool | chrX / chrY outside PAR | PAR |
|---|---|---|
| DeepVariant | chrX: 129,036 hom-alt, **0 het**; chrY: 0 het | 6,981 het (diploid kept) |
| Sawfish | 1,101 haploid calls (`GT=1`); the 14 diploid records are all BND | – |
| TRGT | all haploid (chrY: 3,529 called, 17 missing) | chrX PAR also haploid (TRGT `--karyotype XY` behaviour) |

If chrX shows many `0/1` calls, the `sex` column was not `XY` or the PAR BED is missing.

### 5. Look at the data

The output BAMs are not in git, but your own run has them. Load in IGV (genome: hg38):

- `02_aligned_bam/HG002_XY.aligned.sorted.bam` — colour by base modification (5mC) to see methylation,
- `06a_phased_bam/HG002_XY.haplotagged.bam` — group by tag `HP` to see the two haplotypes,
- `05_tandem_repeats/HG002_XY.trgt.spanning.sorted.bam` — reads over each repeat.

## Output in this folder

[output/](output/) holds the results of the reference run, **except BAM files** (aligned, haplotagged and TRGT spanning BAMs are 0.1–2.6 GB each). See [What each folder contains](../../README.md#what-each-folder-contains) for the file contents.

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

> **Note on `04_structural_var`**: the SV VCF contains whole-chromosome `<DEL>` CNV records with `CN=0` on every chromosome and contig other than chr20/X/Y (e.g. chr1:1-248956000). These are artefacts of the test input, which only has reads on chr20, chrX and chrY; they do not appear in whole-genome runs. That is also why the totals (`total_sv`, all-chromosome counts) are larger than the chr20/X/Y numbers above.
