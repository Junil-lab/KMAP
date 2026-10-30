<h1 align="center">KMAP Pipeline</h1>

<p align="center">
  Reproducible Nextflow + Singularity pipelines for the K-MAP multi-omics project<br>
  PacBio HiFi whole-genome · PacBio long-read RNA · single-cell RNA-seq
</p>

<p align="center">
  <a href="https://github.com/Junil-lab/KMAP/releases"><img alt="release" src="https://img.shields.io/github/v/release/Junil-lab/KMAP?label=release&color=2ea44f"></a>
  <img alt="Nextflow" src="https://img.shields.io/badge/Nextflow-DSL2%20%E2%89%A5%2026.04-0DC09D">
  <img alt="Singularity" src="https://img.shields.io/badge/run%20with-Singularity%20%2F%20Apptainer-1d355c">
  <img alt="Reference" src="https://img.shields.io/badge/reference-GRCh38-blue">
  <img alt="Java" src="https://img.shields.io/badge/Java-17%2B-orange">
</p>

---

## Pipelines

| Module | Data | What it does | Status |
|---|---|---|---|
| [**PBLONGREADWGS**](PBLONGREADWGS/) | PacBio HiFi WGS (BAM) | alignment → SNV/indel → SV → tandem repeats → joint phasing → CpG methylation, sex-aware on chrX/chrY | ✅ available |
| [**PBLONGREADRNA**](PBLONGREADRNA/) | PacBio Iso-Seq / Kinnex | pbmm2 → IsoQuant / isocall / FLAIR → SQANTI3 + ORF; sample/cohort modes | ✅ available |
| [**SCRNA**](SCRNA/scRNAseq_QC_pipeline/) | 10x Genomics scRNA-seq (FASTQ) | Cell Ranger → CellBender → QC → SOLO doublet removal | ✅ available |

Every module follows the same pattern: **one samplesheet in, one command to run, all tools inside Singularity containers**, so results do not depend on what is installed on the server.

### What is ready in each module

✅ ready · ⚠️ works with extra steps or limits · ❌ not available yet

| | PBLONGREADWGS | SCRNA | PBLONGREADRNA |
|---|---|---|---|
| **Input** | HiFi BAM | 10x FASTQ, raw H5 or CellBender H5 | FLNC BAM |
| **Start from the middle** | `start_step`: jasmine / pbmm2 / deepvariant / hiphase | `start_step`: cellranger / cellbender / qc | – (choose callers with `--methods`) |
| **Main container** | ✅ GitHub Release `v1.0.0` (`kmap_pbwgs_v1.0.0.sif`) | ✅ GitHub Release `scrna-qc-v2.0.0` (`scrnaseq_qc.sif`, 3 parts joined by `def/fetch_qc_container.sh`); also buildable from `def/scrnaseq_qc.def` | ✅ GitHub Release `isoseq-v1.0.0` (main + discovery) |
| **Other containers** | ✅ DeepVariant (Docker Hub) and MethBat (quay.io) are pulled with one command | ⚠️ Cell Ranger: build from your own 10x download (license)<br>✅ CellBender pulled from its public image (tag `latest`, version not pinned) | ✅ official SQANTI3 image fetched by `setup.sh` |
| **Reference** | GRCh38 from NCBI (+ index); TRGT catalog, sex BEDs and CpG islands are in the repo | 10x `refdata-gex` download | hg38 + GENCODE v45 (see `resources/ref/`) |
| **Java 17 / Nextflow** | ✅ installed by `setup.sh` | ⚠️ Java 17 must already be installed | ✅ installed by `setup.sh` |
| **Lab-server paths in the code** | ✅ none | ⚠️ fallback conda environment and Java paths of the lab server (`nextflow.config`, `run.sh`, `setup.sh`); not needed when the QC container is used | ✅ none |
| **Example data** | ✅ HG002 chr20/X/Y (download script + expected outputs) | ❌ | ❌ |

---

## PBLONGREADWGS at a glance

```mermaid
flowchart TB
    A[HiFi BAM + sex XX/XY] --> J{5mC tags?}
    J -- no --> B0[01 Jasmine<br/>5mC from kinetics]
    J -- yes --> B[02 pbmm2<br/>align to GRCh38]
    B0 --> B
    B --> C[03 DeepVariant<br/>SNV / indel]
    B --> D[04 Sawfish<br/>SV / CNV]
    B --> E[05 TRGT<br/>tandem repeats]
    C & D & E --> F[06 HiPhase<br/>joint phasing]
    F --> G[07 pb-CpG-tools<br/>CpG, hap1 / hap2]
    F --> M[07b MethBat<br/>5mC / 5hmC / 6mA]
    F --> H[(phased VCFs<br/>haplotagged BAM)]
    G --> Q[QC summary]
```

| Step | Tool | Version |
|---|---|---|
| 5mC from kinetics | Jasmine | 2.7.99 |
| Alignment | pbmm2 | 26.1.99 |
| Small variants | DeepVariant | 1.9.0 (GPU) |
| Structural variants | Sawfish | 2.2.1 |
| Tandem repeats | TRGT (adotto + STRchive catalog, 937k loci) | 5.0.0 |
| Phasing | HiPhase | 1.6.0 |
| CpG methylation | pb-CpG-tools | 3.0.0 |
| Methylation pileup + CpG island profile | MethBat | 1.1.0 |
| Utilities | samtools / bcftools | 1.23.1 |

**Sex-aware calling** — the samplesheet `sex` column switches DeepVariant haploid calling, Sawfish expected copy number and TRGT karyotype together. On HG002 (XY): 0 heterozygous DeepVariant calls on chrX outside the PARs, and haploid TRGT and Sawfish genotypes on chrX/chrY.

The tool set mirrors PacBio's official [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) v4.0.0.

---

## PBLONGREADRNA at a glance

```mermaid
flowchart TB
    A[FLNC BAM per sample] --> B[pbmm2 ISOSEQ<br/>align once per sample]
    B --> I1[IsoQuant]
    B --> I2[isocall<br/>profile + call]
    B --> I3[FLAIR<br/>transcriptome]
    I1 & I2 & I3 --> S[SQANTI3 QC<br/>+ ORF: TD2 / PSAURON]
    I2 & I3 -. 2+ samples .-> C[cohort transcriptome<br/>isocall / FLAIR]
    C --> S
    S --> R[HTML report<br/>+ QC summary]
```

| Step | Tool | Version |
|---|---|---|
| Alignment | pbmm2 (ISOSEQ preset) | 26.1.99 |
| Transcript discovery | IsoQuant | 3.13.1 |
| Transcript discovery | isocall | 1.3.0 |
| Transcript discovery | FLAIR | 3.0.0 (pinned commit + patch) |
| Structural QC | SQANTI3 | 6.0.1 (official image) |
| ORF prediction | TD2 / PSAURON | in the SQANTI3 image |

Choose callers with `--methods isoquant,isocall,flair`; with more than one sample, isocall and FLAIR also build a joint cohort transcriptome (`--run_cohort false` to skip). Details: [PBLONGREADRNA/README.md](PBLONGREADRNA/README.md).

---

## SCRNA at a glance

```mermaid
flowchart TB
    A[10x FASTQ] --> B[01 Cell Ranger<br/>count matrix]
    B --> C[02 CellBender<br/>ambient RNA removal]
    C --> D[03 Scanpy QC<br/>cell QC, all genes kept]
    D --> E[04 SOLO<br/>doublet scores]
    E --> F[05 singlet filter]
    D --> F
    F --> G[(singlet H5AD<br/>all genes + raw)]
    F --> Q[06 QC summary]
    X[raw H5] -. start_step=cellbender .-> C
    Y[CellBender H5] -. start_step=qc .-> D
```

| Step | Tool | Version |
|---|---|---|
| Count matrix | Cell Ranger | 10.0.0 (built from your own 10x download) |
| Background removal | CellBender | 0.3.0 (public image, tag `latest`) |
| Cell QC, HVGs | Scanpy | 1.11.5 |
| Doublet detection | SOLO (solo-sc) + scvi-tools | 1.2 + 1.4.0.post1 |
| Deep learning backend | PyTorch | 2.9.0 (CUDA 12.8) |

Mitochondrial filtering is off by default (`--mt_threshold 20` to enable); all input genes stay in the final H5AD and `.raw`. Details: [SCRNA/scRNAseq_QC_pipeline/README.md](SCRNA/scRNAseq_QC_pipeline/README.md).

---

## Quick start (WGS)

```bash
git clone https://github.com/Junil-lab/KMAP.git
cd KMAP/PBLONGREADWGS

# containers
wget -O containers/kmap_pbwgs_v2.sif \
  https://github.com/Junil-lab/KMAP/releases/download/v1.0.0/kmap_pbwgs_v1.0.0.sif
singularity pull containers/deepvariant_1.9.0-gpu.sif docker://google/deepvariant:1.9.0-gpu
singularity pull containers/methbat_1.1.0.sif docker://quay.io/biocontainers/methbat:1.1.0--h9ee0642_0

# reference: see PBLONGREADWGS/README.md → Reference genome
bash setup.sh            # installs Java 17 if needed, checks everything

bash run.sh --samplesheet input/samples.csv --outdir results --threads 16
```

```csv
sample_id,bam,sex,start_step
SampleA,/abs/path/SampleA.hifi_reads.bam,XX,pbmm2
```

New here? Run the **[HG002 example](PBLONGREADWGS/example/HG002_chr20XY/)** first (~3 h on 5 cores) and compare your numbers with the reference results.

---

## Repository layout

```
KMAP/
├── PBLONGREADWGS/          PacBio HiFi WGS pipeline
│   ├── main.nf, nextflow.config, run.sh, setup.sh
│   ├── def/                container recipes
│   ├── resources/ref/      sex BEDs + TRGT catalog (GRCh38 FASTA downloaded separately)
│   └── example/HG002_chr20XY/   test input script, expected outputs
├── PBLONGREADRNA/          PacBio long-read RNA pipeline
└── SCRNA/
    └── scRNAseq_QC_pipeline/   10x scRNA-seq preprocessing and QC
```

## Documentation

- [PBLONGREADWGS/README.md](PBLONGREADWGS/README.md) — requirements, samplesheet, options, every output file explained
- [PBLONGREADWGS/example/HG002_chr20XY/](PBLONGREADWGS/example/HG002_chr20XY/) — step-by-step test run with expected results
- [PBLONGREADRNA/README.md](PBLONGREADRNA/README.md) — tool selection, sample/cohort modes, SQANTI3 and ORF outputs
- [SCRNA/scRNAseq_QC_pipeline/README.md](SCRNA/scRNAseq_QC_pipeline/README.md) — scRNA-seq workflow and QC defaults
- [Releases](https://github.com/Junil-lab/KMAP/releases) — prebuilt containers

