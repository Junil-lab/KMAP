# SCRNA — K-MAP 10x scRNA-seq preprocessing / QC pipeline

A Nextflow (DSL2) + Singularity pipeline that takes 10x Genomics FASTQs per sample and runs Cell Ranger → CellBender → Scanpy QC → SOLO doublet detection in one go, ending with a QC-passed, singlet-only full-gene H5AD for each sample.

## How it works

```
10x FASTQ (one folder per sample)
  │
  ├─ 01 CELLRANGER_COUNT              FASTQ → raw / filtered feature-barcode matrix
  ├─ 02 CELLBENDER_REMOVE_BACKGROUND  ambient RNA / background removal (GPU)
  ├─ 03 SCANPY_QC_FOR_SOLO            cell QC (all genes kept), gene prefilter + HVG subset for SOLO only
  ├─ 04 SOLO_DOUBLETS                 doublet scoring (GPU), then score correction (correct_solo_scores.py)
  ├─ 05 FILTER_SOLO_SINGLETS          SOLO labels added to the full matrix, singlet-only matrix
  └─ 06 QC_MERGE                      QC summary for all samples
```

You can also start later: from an existing raw matrix (`start_step=cellbender`) or an existing CellBender output (`start_step=qc`).

### QC defaults

| Filter | Default | Option |
|---|---|---|
| mitochondrial fraction | **off** — `pct_counts_mt` is computed and recorded but not used to remove cells | `--mt_threshold 20` removes cells with `pct_counts_mt >= 20` |
| detected genes | remove cells with `< 200` genes | `--min_genes` |
| total UMI | remove cells with `< 500` counts | `--min_counts` |
| gene prefilter | genes seen in `>= 5` cells — applied **only to the SOLO training copy**, not to the final matrices | `--min_cells_gene` |
| HVGs for SOLO | `2000` (`seurat_v3`, falls back to log-normalized Seurat, then dispersion) | `--hvg_n_top` |
| CellBender | FPR `0.01`, `150` epochs | `--cellbender_fpr`, `--cellbender_epochs` |

SOLO model settings are in `configs/solo_params.json` (`n_hidden` 384, `n_latent` 64, `n_layers` 1, `cl_hidden` 128, `cl_layers` 1, `dropout_rate` 0.2, `learning_rate` 0.001, `valid_pct` 0.1, `batch_size` 128), following the [Calico SOLO](https://github.com/calico/solo) README.

## Requirements

| Item | Location | Notes |
|---|---|---|
| Java 17+ | system, or **installed automatically** by `setup.sh` | Nextflow 26 needs Java 17+. If missing, `setup.sh` downloads OpenJDK 17 into `.java17/` and `run.sh` uses it |
| Nextflow | **in this folder** (`./nextflow`, version 26.04.3) | downloads itself on first use |
| Singularity / Apptainer | system | |
| GPU (recommended) | – | CellBender and SOLO; without one use `--use_gpu false` (much slower) |
| `containers/cellranger.sif` | build from the 10x tarball (see [Containers](#containers)) | Cell Ranger needs a license agreement, so it cannot be redistributed |
| `containers/cellbender.sif` | Google Artifact Registry | `us.gcr.io/broad-dsde-methods/cellbender:latest` |
| `containers/scrnaseq_qc.sif` | GitHub Release `scrna-qc-v2.0.0` (`bash def/fetch_qc_container.sh`), or build with `def/scrnaseq_qc.def` | Scanpy 1.11.5, scvi-tools 1.4.0.post1, SOLO 1.2, PyTorch 2.9.0 (CUDA 12.8). 4.5 GB, stored as 3 parts in the release |
| 10x transcriptome reference | 10x Genomics (see [Reference](#reference)) | not downloaded by the pipeline |

> `def/scrnaseq_qc.def` installs every package from public sources (PyPI, the PyTorch CUDA 12.8 index and the official calico/solo source at commit `d94baa3`), pinned to the versions the pipeline was validated with, so it builds on any server. On a 10x PBMC 1k test it gave results identical to the earlier lab image (same QC cells, HVGs, matrix values, doublet barcodes and scores).

### Containers

Cell Ranger: download the tarball from the [10x Genomics download page](https://www.10xgenomics.com/support/software/cell-ranger/latest) (v10.0 at the time of writing) and package it:

```bash
CELLRANGER_TARBALL=/path/to/cellranger-10.0.0.tar.gz bash def/build_containers.sh cellranger
# or from an extracted directory
CELLRANGER_DIR=/path/to/cellranger-10.0.0 bash def/build_containers.sh cellranger
```

The QC container (Scanpy / SOLO) — download the prebuilt image:

```bash
bash def/fetch_qc_container.sh     # downloads 3 parts, checks SHA-256, joins them into containers/scrnaseq_qc.sif
```

The others:

```bash
bash def/build_containers.sh cellbender
bash def/build_containers.sh qc     # optional: build the QC image yourself instead of downloading (~30 min)
# or all three at once, with downloads/cellranger-10.0.0.tar.gz in place
bash def/build_containers.sh all
```

To use a Cell Ranger installed on the host instead of the SIF: `-profile local,host_cellranger --cellranger_bin /path/to/cellranger`.

### Reference

Download and extract a 10x reference, then give its path in the samplesheet `transcriptome` column or with `--transcriptome`:

```bash
mkdir -p reference && cd reference
wget -c https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz    # human GRCh38
wget -c https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz    # mouse GRCm39
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
```

## Running

### 1. Check the environment (once)

```bash
git clone https://github.com/hi6167579-droid/KMAP_Pipeline.git
cd KMAP_Pipeline/SCRNA/scRNAseq_QC_pipeline

# containers (see Containers above), then
bash setup.sh          # installs Java 17 if needed; checks nextflow / Singularity / SIFs; "All checks passed." means ready
```

### 2. Write a samplesheet

Copy `input/samples.example.csv` and fill it in. Only the columns for your `start_step` are needed; leave the rest empty.

```csv
sample_id,fastqs,sample,transcriptome,start_step,raw_h5,cellbender_h5,expected_cells,chemistry,create_bam,include_introns,localmem_gb,cuda_device,cellbender_expected_cells,cellbender_total_droplets_included,expected_doublet_rate,expected_doublets
SampleA,/abs/path/fastq/SampleA,SampleA,/abs/path/refdata-gex-GRCh38-2024-A,cellranger,,,,auto,true,false,200,0,,,,
SampleB,,,,qc,,/abs/path/SampleB/cellbender_cleaned_filtered.h5,,,,,,0,,,,
```

| Column | Value |
|---|---|
| `sample_id` | name of the output folder and barcode prefix |
| `fastqs` | Cell Ranger FASTQ directory (`start_step=cellranger`) |
| `sample` | Cell Ranger `--sample` (FASTQ file prefix) |
| `transcriptome` | extracted Cell Ranger reference directory |
| `start_step` | `cellranger` — from FASTQ<br>`cellbender` — from a raw matrix (`raw_h5`)<br>`qc` — from a CellBender output (`cellbender_h5`) |
| `raw_h5` | `raw_feature_bc_matrix.h5`, when `start_step=cellbender` |
| `cellbender_h5` | CellBender cleaned H5, when `start_step=qc` |
| `expected_cells`, `chemistry`, `create_bam`, `include_introns`, `localmem_gb` | Cell Ranger settings (empty = Cell Ranger default) |
| `cuda_device` | GPU index for this sample |
| `cellbender_expected_cells`, `cellbender_total_droplets_included` | CellBender settings (empty = automatic) |
| `expected_doublet_rate` | only if known, e.g. `0.06`; otherwise SOLO's own score threshold is used |
| `expected_doublets` | an explicit number of doublets instead of a rate |

### 3. Run

```bash
bash run.sh --samplesheet input/my_samples.csv --outdir /path/to/results --threads 24 --use_gpu true --cuda_device 0

# common options
#   -resume                  continue from where it stopped (reuses the cache)
#   --use_gpu false          CellBender / SOLO on CPU
#   --keep_work true         keep all intermediate files (see Disk use and clean-up)
#   --memory_gb 200          memory for Cell Ranger / GPU steps
#   --mt_threshold, --min_genes, --min_counts, --hvg_n_top, --solo_expected_doublet_rate
#   -profile local_conda --qc_conda_env /path/to/env   Scanpy / SOLO from conda instead of scrnaseq_qc.sif
#   -profile slurm           submit to SLURM (queue: gpu)
#   export NXF_WORK=/fast/local/disk/work   # intermediate files (default: ~/nf_work/K_map_scrnaseq_qc)
```

On a shared server, limit the cores: `taskset -c 0-4 bash run.sh ... --threads 5`

### Disk use and clean-up

Nextflow keeps every intermediate file in the work directory, and outputs are **copied** from there to `--outdir`. To avoid keeping everything twice, `run.sh` cleans up after a **successful** run:

| | |
|---|---|
| Kept | the whole Cell Ranger work directory (including `outs`), inputs, references, all of `--outdir` |
| Removed | CellBender checkpoints / posterior / temporary files, and the work directories of the QC, SOLO, singlet and summary steps |
| Safety | only successful tasks in the current run's trace are touched, and every output is first compared with its published copy (SHA-256); any mismatch stops the clean-up before anything is deleted |
| Failed run | nothing is cleaned, so you can debug and `-resume` |

- Cleaned steps are recomputed on a later `-resume`; Cell Ranger is not, because its work is kept.
- Earlier runs and other work directories are never touched. Running `nextflow run main.nf` directly skips the clean-up.
- Needs Python 3.9+ on the host (standard library only). `--keep_work true` turns it off.
- Records: `pipeline_info/run_*.trace.csv` (tasks of the run) and `pipeline_info/run_*.cleanup.log` (what was removed).

## Outputs

```
<outdir>/<sample_id>/
├── 01_cellranger/
├── 02_cellbender/
├── 03_scanpy_qc/
├── 04_solo/
└── 05_singlets/
<outdir>/scrnaseq_qc_summary.csv   QC table for all samples
<outdir>/pipeline_info/            timeline / report / trace, run and clean-up logs
```

### What each folder contains

| Folder | File | Values |
|---|---|---|
| `01_cellranger` | `<id>/outs/raw_feature_bc_matrix.h5` | all barcodes, unfiltered (input to CellBender) |
| | `<id>/outs/filtered_feature_bc_matrix.h5` | Cell Ranger's called cells |
| | `<id>/outs/metrics_summary.csv`, `web_summary.html` | sequencing / mapping / cell-calling metrics |
| | `<id>/outs/possorted_genome_bam.bam` (+ `.bai`) | aligned reads (only with `create_bam=true`) |
| `02_cellbender` | `*.cellbender_cleaned.h5` | background-corrected counts |
| | `*.cellbender_cleaned*.csv`, `*.pdf` | called cell barcodes, metrics, training / cell-probability plots |
| `03_scanpy_qc` | `*_qc_full.h5ad` | cells that passed QC, **all input genes** in `.X`, `layers["counts"]` and `.raw`; no SOLO labels yet |
| | `*_qc_hvg2000.h5ad` | SOLO training copy: `.X` = HVGs after the gene prefilter, `.raw` = all input genes (not for final analysis) |
| | `*.qc_metrics.csv` | cells / genes before and after QC, cells removed (total and by the MT filter), median counts / genes / % MT before and after |
| `04_solo` | `*.solo/` | SOLO model output: `soloed.h5ad`, doublet scores (corrected, see [SOLO score correction](#solo-score-correction)), `*_double_softmax_original` archives, `score_correction.json`, plots |
| `05_singlets` | `*_qc_solo_annotated_full.h5ad` | every QC-passed cell and all input genes (+ `.raw`), with `solo_is_doublet`, `solo_doublet_probability`, `solo_softmax_scores`, `solo_logit_scores`, `solo_double_softmax_scores_original` in `obs` |
| | `*_qc_solo_singlets_full.h5ad` | only `solo_is_doublet == False` cells, all input genes (+ `.raw`) — **the matrix for downstream analysis** |
| | `*.solo_summary.csv` | cells after QC, doublets, singlets, doublet %, gene counts (`n_genes_full`, `n_genes_raw`) |
| | `*.solo_doublet_barcodes.csv` | barcodes called doublets, with scores |
| `scrnaseq_qc_summary.csv` | | one row per sample: the QC metrics and SOLO summary combined |
| `pipeline_info/` | | run time and memory per step (`trace.csv`, `timeline.html`, `report.html`) |

The annotated file lets you compare results before and after doublet removal, or keep the doublet calls only as metadata.

More detail with example tables: [docs/scrnaseq_qc_pipeline_io_tools.html](docs/scrnaseq_qc_pipeline_io_tools.html) and [docs/scrnaseq_qc_pipeline_step_io_examples.xlsx](docs/scrnaseq_qc_pipeline_step_io_examples.xlsx).

## Notes

- **Doublet rate**: if `expected_doublet_rate` is empty, SOLO decides with its own score threshold. For 10x data, setting the rate per sample (from the loaded cell number) is better when known.
- **GPU**: CellBender and SOLO use `--nv`; pick the card per sample with `cuda_device` or for the whole run with `--cuda_device`.
- **Cell Ranger license**: the Cell Ranger SIF is built locally from your own download and is not distributed with this repository.

## Full input genes and AnnData.raw

After cell QC (default `genes >= 200`, `counts >= 500`; plus `mt < 20%` if `--mt_threshold 20` is given) and **before any gene filtering**, the pipeline stores every input gene with `adata.raw = adata.copy()`. Here raw means the CellBender background-corrected counts used as input, not the uncorrected Cell Ranger counts.

- `*_qc_full.h5ad`: cell QC only; all input genes in `.X`, `layers["counts"]` and `.raw`.
- `*_qc_hvg*.h5ad`: the gene prefilter (`min_cells_gene`) and HVG selection are applied to this training copy only. `.X` / `layers["counts"]` hold the HVGs (used to train SOLO), `.raw` holds all input genes.
- `*_qc_solo_annotated_full.h5ad` and `*_qc_solo_singlets_full.h5ad`: all QC-passed cells (or singlets only) with all input genes, `.raw` and the SOLO calls.
- If SOLO is missing any QC cell barcode, the run stops instead of silently dropping cells.

In the QC summary, `n_genes_after_qc` and `n_genes_raw` are the number of preserved genes, `n_genes_solo_prefilter` the number after the SOLO gene prefilter, and `n_hvg` the number of genes SOLO was trained on. `.var["solo_gene_filter_pass"]` and `.var["highly_variable"]` mark the genes used for training. Older gene-filtered full files cannot recover the dropped genes; rerun QC from the original input.

## SOLO score correction

The SOLO CLI on our server saves the output of `predict(soft=True)` — already a probability — as `logit_scores` and then applies softmax again. After SOLO finishes, `scripts/correct_solo_scores.py` checks this relation and restores `softmax_scores` (and `solo_softmax_scores`) to the original doublet probability; `solo_doublet_probability` holds the same value.

- `solo_logit_scores` keeps its old values for compatibility, but they are probabilities, not logits.
- The old double-softmax values are kept in the `*_double_softmax_original` files and the `solo_double_softmax_scores_original` column.
- The singlet / doublet calls themselves are unchanged.
- The PDF plots written by the SOLO CLI show the uncorrected scores; this is recorded in `score_correction.json`.

## QC consistency checks (2026-10-24)

Sample IDs must be unique and contain only letters, digits, dots, underscores or
hyphens, starting with a letter or digit. QC and SOLO must contain exactly the same
unique cell barcodes (order may differ). Missing/extra barcodes or missing/invalid
doublet labels stop the run before final H5AD files are written.

The combined summary requires one QC and one SOLO summary for every input sample.
It rejects missing files, duplicate samples and inconsistent cell counts; doublets
plus singlets must equal the QC-passed cell count. The successful-run work cleanup
therefore also remains blocked when these checks fail.

Run the tests inside the QC environment/container with
`python3 -m unittest discover -s tests -v` (AnnData, NumPy and pandas required).
