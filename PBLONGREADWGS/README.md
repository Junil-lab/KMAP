# PBLONGREADWGS — K-MAP PacBio HiFi WGS pipeline

A Nextflow (DSL2) + Singularity pipeline that takes one PacBio HiFi BAM per sample and runs alignment → SNV/indel → SV → tandem repeats → phasing → CpG methylation in one go. Reference: GRCh38 (`GRCh38_no_alt_analysis_set`).

## How it works

```
HiFi BAM (one BAM per sample, whole genome)
  │
  ├─ 01 JASMINE      kinetics (ip/pw) → 5mC tags (MM/ML)   * only when start_step=jasmine
  ├─ 02 PBMM2        alignment to GRCh38
  ├─ 03 DEEPVARIANT  SNV/indel (GPU)
  ├─ 04 SAWFISH      structural variants (SV)
  ├─ 05 TRGT         tandem repeat genotyping (adotto + STRchive catalog)
  ├─ 06 HIPHASE      joint SNV + SV + TR phasing, haplotagged BAM
  └─ 07 CPG_TOOLS    CpG methylation (combined / hap1 / hap2)
       └─ QC summary (per sample + all samples)
```

### Sex handling

The `sex` value (XX/XY) in the samplesheet is passed to three tools. It is **required**; the run stops before starting if it is empty or anything else.

| Tool | XX | XY |
|---|---|---|
| DeepVariant | default (diploid) | `--haploid_contigs chrX,chrY` + PAR BED (PARs stay diploid) |
| Sawfish | `expected_cn.hg38.XX.bed` | `expected_cn.hg38.XY.bed` + `--treat-single-copy-as-haploid` |
| TRGT | `--karyotype XX` | `--karyotype XY` (chrX PAR is also treated as haploid — TRGT's own behaviour) |

For cell lines whose karyotype differs from the donor's sex (e.g. a line that lost one X), set the value **from the data** (one X copy, no Y → XY).

## Requirements

| Item | Location | Notes |
|---|---|---|
| Java 17+ | system, or **installed automatically** by `setup.sh` | Nextflow 26 needs Java 17+. If missing, `setup.sh` downloads OpenJDK 17 into `.java17/` and `run.sh` uses it |
| Singularity / Apptainer | system | |
| GPU (optional) | – | for DeepVariant; without one use `--use_gpu false` |
| `containers/kmap_pbwgs_v2.sif` | [GitHub Release v1.0.0](https://github.com/hi6167579-droid/KMAP_Pipeline/releases/tag/v1.0.0) | pbmm2, jasmine, sawfish, trgt, hiphase, pb-CpG-tools, samtools/bcftools — can be built with `def/build_containers.sh` |
| `containers/deepvariant_1.9.0-gpu.sif` | Docker Hub | `google/deepvariant:1.9.0-gpu` |
| `resources/ref/GRCh38_no_alt_analysis_set.fna` (+ `.fai`, `.mmi`) | distributed separately | |
| `resources/ref/trgt/adotto_strchive_20250827.hg38.noFailTag.bed.gz` | **in this repo** | catalog from the PacBio HiFi-human-WGS-WDL v4.0.0 reference bundle (Zenodo) with the `;INCLUDE_FAIL_READS` tag removed (same processing as the WDL) |
| `resources/ref/sex/*.bed` | **in this repo** | PAR BED, Sawfish expected copy number BEDs |

The GRCh38 reference and BAMs are not in the repo because of their size (`.gitignore`). Ask the lab admin for their location on the lab server.

## Running

### 1. Check the environment (once)

```bash
git clone https://github.com/hi6167579-droid/KMAP_Pipeline.git
cd KMAP_Pipeline/PBLONGREADWGS

# containers
wget -O containers/kmap_pbwgs_v2.sif \
  https://github.com/hi6167579-droid/KMAP_Pipeline/releases/download/v1.0.0/kmap_pbwgs_v1.0.0.sif
singularity pull containers/deepvariant_1.9.0-gpu.sif docker://google/deepvariant:1.9.0-gpu

# GRCh38 reference (+ .fai, .mmi) goes in resources/ref/ (see Requirements)
bash setup.sh          # installs Java 17 if needed; checks Singularity / nextflow / SIFs / references; all [OK] means ready
```

### 2. Write a samplesheet

Copy `input/samples.csv` and fill it in.

```csv
sample_id,bam,sex,start_step
SampleA,/absolute/path/SampleA.hifi_reads.bam,XX,pbmm2
SampleB,/absolute/path/SampleB.hifi_reads.bam,XY,pbmm2
```

| Column | Value |
|---|---|
| `sample_id` | name of the output folder |
| `bam` | absolute path. Must be **inside the pipeline folder** → if it is elsewhere, `ln -s /external/path/x.bam input/bam/` |
| `sex` | `XX` or `XY` (required) |
| `start_step` | `jasmine` — BAM with kinetics (ip/pw) but no methylation tags<br>`pbmm2` — HiFi BAM that already has methylation tags (MM/ML) (most cases)<br>`deepvariant` / `hiphase` — restart from the middle |

### 3. Run

```bash
bash run.sh --samplesheet input/my_samples.csv --outdir /path/to/results --threads 16

# common options
#   -resume              continue from where it stopped (reuses the cache)
#   --use_gpu false      DeepVariant without GPU
#   --vcf_dir <dir>      existing VCF location when start_step=hiphase
#   export NXF_WORK=/fast/local/disk/work   # move intermediate files (default: <pipeline folder>/.nf_work)
```

On a shared server, limit the cores: `taskset -c 0-4 bash run.sh ... --threads 5`

## Outputs

```
<outdir>/<sample_id>/
├── 01_methyl_bam/        (only if the jasmine step ran)
├── 02_aligned_bam/       *.aligned.sorted.bam, flagstat
├── 03_small_variants/    *.deepvariant.vcf.gz, *.g.vcf.gz
├── 04_structural_var/    *.sv.vcf.gz, run.stats.json
├── 05_tandem_repeats/    *.trgt.vcf.gz, *.trgt.spanning.sorted.bam
├── 06a_phased_bam/       *.haplotagged.bam, summary/blocks/stats
├── 06b_phased_SNP/       *.snv.phased.vcf.gz
├── 06c_phased_SV/        *.sv.phased.vcf.gz
├── 06d_phased_TR/        *.tr.phased.vcf.gz
├── 07_CpG_methylation/   *.combined / hap1 / hap2 .bed.gz
└── logs/
<outdir>/qc_summary/      QC table for all samples
<outdir>/pipeline_info/   timeline / report / trace, run logs
```

## Example

[example/HG002_chr20XY/](example/HG002_chr20XY/) has a test run on HG002 (XY, chr20 + chrX + chrY): the script to rebuild the input BAM, the samplesheet, and all outputs except BAMs.

## Notes

- **Runtime**: HG002 chr20+X+Y (~5% of the genome) on 5 cores took 3 h 16 min, 76% of it DeepVariant. Give a whole-genome 30x run plenty of cores.
- **Methylation**: the pipeline summarises the methylation tags already in the BAM. 5hmC/6mA are only available if the tags were produced at sequencing time and the kinetics were kept.
- **Check TRGT logs**: TRGT exits 0 even if some loci fail. Check the `[ERROR]` lines in the `05_tandem_repeats` task log.
- **Removed steps**: variant annotation (VEP, AnnotSV) and sTELLeR (TE detection) were removed in this version.
- **Reference point**: the tool set mostly matches PacBio's official [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) v4.0.0; versions here are DeepVariant 1.9.0 / TRGT 5.0.0 / HiPhase 1.6.0 / Sawfish 2.2.1 / pb-CpG-tools 3.0.0.
