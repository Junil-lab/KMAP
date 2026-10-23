# K_map scRNA-seq preprocessing/QC pipeline

A Nextflow + Singularity pipeline for 10x Genomics scRNA-seq preprocessing and QC, built on the same structure as the PacBio WGS pipeline in this repository.

## Workflow

1. `CELLRANGER_COUNT`: FASTQ -> 10x raw/filtered matrix
2. `CELLBENDER_REMOVE_BACKGROUND`: raw matrix -> background-corrected H5
3. `SCANPY_QC_FOR_SOLO`: basic QC, mitochondrial filter, HVG subset
4. `SOLO_DOUBLETS`: SOLO doublet scoring
5. `FILTER_SOLO_SINGLETS`: adds SOLO labels/scores to the full-gene QC matrix and writes a singlet-only matrix
6. `QC_MERGE`: merges the per-sample QC summaries

## Default QC

- remove cells with `pct_counts_mt >= 20`
- remove cells with `n_genes_by_counts < 200`
- remove cells with `total_counts < 500`
- gene prefilter: `min_cells = 5`
- HVG: `2000` genes
- HVGs are computed with `seurat_v3` first; if that fails, the log-normalized Seurat flavour is used, and dispersion as the last fallback.

## SOLO defaults

`configs/solo_params.json`:

- `n_hidden`: 384
- `n_latent`: 64
- `n_layers`: 1
- `cl_hidden`: 128
- `cl_layers`: 1
- `dropout_rate`: 0.2
- `learning_rate`: 0.001
- `valid_pct`: 0.1
- `batch_size`: 128

The SOLO CLI and example parameters follow the Calico SOLO README: https://github.com/calico/solo

## Setup

```bash
cd KMAP_Pipeline/SCRNA/scRNAseq_QC_pipeline
```

### Cell Ranger SIF

By default Cell Ranger runs from its own Singularity image, not from a host installation:

```text
containers/cellranger.sif
```

10x Genomics Cell Ranger requires accepting a license, so a fully automatic download may be blocked depending on the server or account. Download the latest Cell Ranger tarball from the official page and package it as a SIF:

```text
https://www.10xgenomics.com/support/software/cell-ranger/latest
```

The latest Cell Ranger in the 10x documentation is v10.0. From the tarball:

```bash
CELLRANGER_TARBALL=/path/to/cellranger-10.0.0.tar.gz \
  bash def/build_containers.sh cellranger
```

From an already extracted Cell Ranger directory:

```bash
CELLRANGER_DIR=/path/to/cellranger-10.0.0 \
  bash def/build_containers.sh cellranger
```

Build the other containers:

```bash
bash def/build_containers.sh qc
bash def/build_containers.sh cellbender
bash setup.sh
```

If `downloads/cellranger-10.0.0.tar.gz` is already in place, everything can be built at once:

```bash
bash def/build_containers.sh all
```

Use the fallback profile only if you must use a host Cell Ranger installation:

```bash
bash run.sh \
  --samplesheet input/samples.example.csv \
  -profile local,host_cellranger \
  --cellranger_bin /path/to/cellranger
```

### Transcriptome reference

The pipeline does not download references. Download and extract a 10x reference first, then give its path in the samplesheet `transcriptome` column or with `--transcriptome`.

Common 10x Cell Ranger reference tarballs:

```text
mouse GRCm39:
  https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz

human GRCh38:
  https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
```

Example:

```bash
mkdir -p reference
cd reference
wget -c https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz
tar -xzf refdata-gex-GRCm39-2024-A.tar.gz
```

If a SOLO conda environment already exists on the server, the pipeline can run before `scrnaseq_qc.sif` is built:

```bash
bash run.sh \
  --samplesheet input/samples.example.csv \
  -profile local_conda \
  --qc_conda_env /path/to/solo_conda_env
```

## Samplesheet

Example file:

```bash
input/samples.example.csv
```

Columns:

- `sample_id`: sample name, used for the output directory and the barcode prefix
- `fastqs`: Cell Ranger input FASTQ directory
- `sample`: Cell Ranger `--sample`
- `transcriptome`: extracted Cell Ranger reference directory
- `start_step`: `cellranger`, `cellbender` or `qc`
- `raw_h5`: raw feature matrix H5, when `start_step=cellbender`
- `cellbender_h5`: CellBender cleaned H5, when `start_step=qc`
- `expected_doublet_rate`: only if known, e.g. `0.06`
- `expected_doublets`: an explicit number of doublets instead of a rate

## Run

From FASTQ:

```bash
bash run.sh \
  --samplesheet input/samples.example.csv \
  --outdir results \
  --threads 24 \
  --use_gpu true \
  --cuda_device 0
```

To run only QC/SOLO on existing CellBender output, set `start_step=qc` and `cellbender_h5=/path/to/cellbender_cleaned.h5` in the samplesheet.

## Outputs

Main outputs per sample:

- `01_cellranger/`: Cell Ranger output
- `02_cellbender/`: CellBender cleaned H5
- `03_scanpy_qc/`: `*_qc_full.h5ad`, `*_qc_hvg2000.h5ad`, `*.qc_metrics.csv`
- `04_solo/`: `soloed.h5ad`, SOLO score arrays/plots
- `05_singlets/`: `*_qc_solo_annotated_full.h5ad`, `*_qc_solo_singlets_full.h5ad`, doublet barcode CSV
- `scrnaseq_qc_summary.csv`: QC summary for all samples

`*_qc_solo_annotated_full.h5ad` keeps every cell that passed QC and adds the `solo_is_doublet`, `solo_softmax_scores` and `solo_logit_scores` columns to `obs`. `*_qc_solo_singlets_full.h5ad` keeps only the cells with `solo_is_doublet == False`; this is the final matrix for downstream analysis.

## Work directory retention

By default, `bash run.sh` cleans up the **work directories of every step except Cell Ranger** after the whole run succeeds.
During the run, the intermediate files needed by later steps are created and kept. Only after Nextflow exits successfully and
the copies to `results` are finished does the clean-up start, and it only touches successful tasks recorded in the trace of the current run.
Every output is first checked against its published copy (SHA-256) before anything is deleted.

- Kept: the whole Cell Ranger work directory (including `outs`), the original inputs, references, and all of `results`.
- Removed: CellBender checkpoints/posterior/temporary files, and the whole work directory of the QC/SOLO/singlet/summary steps.
- Failed runs: nothing is cleaned, so the work directory stays available for debugging and `-resume`.
- Copy check fails: the clean-up stops before deleting anything. Results are kept; see the clean-up log for the error.
- Cleaned steps may be recomputed on a later `-resume`. Cell Ranger work is kept.
- Earlier finished runs and other work directories are never cleaned automatically.
- Running `nextflow run main.nf` directly skips this post-processing. Use `run.sh` for automatic clean-up.
- Needs Python 3.9+ on the host (standard library only).

To keep all work directories:

```bash
bash run.sh --samplesheet input/samples.example.csv --keep_work true
```

The per-run trace is written to `results/pipeline_info/run_*.trace.csv` and the clean-up record to `run_*.cleanup.log`.
