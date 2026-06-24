# K_map scRNA-seq preprocessing/QC pipeline

이 파이프라인은 기존 `/home/NAS_2/K_MAP/PacBio_WGS_pipeline_v2/`의 Nextflow/Singularity 구조를 기준으로 만든 scRNA-seq 전처리 및 QC 파이프라인입니다.

## Workflow

1. `CELLRANGER_COUNT`: FASTQ -> 10x raw/filtered matrix
2. `CELLBENDER_REMOVE_BACKGROUND`: raw matrix -> background-corrected H5
3. `SCANPY_QC_FOR_SOLO`: basic QC, mitochondrial filter, HVG subset 생성
4. `SOLO_DOUBLETS`: SOLO doublet scoring
5. `FILTER_SOLO_SINGLETS`: full-gene QC matrix에 SOLO label/score를 붙이고 singlet-only matrix 생성
6. `QC_MERGE`: 샘플별 QC 요약 병합

## Default QC

- `pct_counts_mt >= 20` cell 제거
- `n_genes_by_counts < 200` cell 제거
- `total_counts < 500` cell 제거
- gene prefilter: `min_cells = 5`
- HVG: `2000` genes
- HVG 계산은 `seurat_v3`를 우선 사용하고 실패하면 log-normalized Seurat 방식, 마지막으로 dispersion fallback을 사용합니다.

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

SOLO CLI와 예시 파라미터는 Calico SOLO README 기준입니다: https://github.com/calico/solo

## Setup

```bash
cd /home/NAS_2/K_MAP/PacBio_WGS_pipeline_v2/scRNAseq_QC_pipeline
```

### Cell Ranger SIF

Cell Ranger는 기본적으로 host 설치본이 아니라 전용 Singularity image를 사용합니다.

```text
containers/cellranger.sif
```

10x Genomics Cell Ranger는 라이선스 동의가 필요한 소프트웨어라서, 완전 자동 다운로드가 서버/계정에 따라 막힐 수 있습니다. 공식 페이지에서 최신 Cell Ranger tarball을 받은 뒤 SIF로 패키징하는 방식을 권장합니다.

공식 페이지:

```text
https://www.10xgenomics.com/support/software/cell-ranger/latest
```

현재 10x 공식 문서의 최신 Cell Ranger는 v10.0입니다. tarball을 받은 경우:

```bash
CELLRANGER_TARBALL=/path/to/cellranger-10.0.0.tar.gz \
  bash def/build_containers.sh cellranger
```

이미 압축 해제된 Cell Ranger 디렉터리가 있으면 그 디렉터리로도 SIF를 만들 수 있습니다.

```bash
CELLRANGER_DIR=/path/to/cellranger-10.0.0 \
  bash def/build_containers.sh cellranger
```

나머지 컨테이너는 다음처럼 준비합니다.

```bash
bash def/build_containers.sh qc
bash def/build_containers.sh cellbender
bash setup.sh
```

`downloads/cellranger-10.0.0.tar.gz`를 미리 넣어둔 상태라면 한 번에 빌드할 수도 있습니다.

```bash
bash def/build_containers.sh all
```

기존 host Cell Ranger를 꼭 써야 하는 경우에만 fallback profile을 사용합니다.

```bash
bash run.sh \
  --samplesheet input/samples.example.csv \
  -profile local,host_cellranger \
  --cellranger_bin /path/to/cellranger
```

### Transcriptome reference

파이프라인은 reference를 자동 다운로드하지 않습니다. 사용자가 먼저 10x reference를 다운로드하고 압축 해제한 뒤, samplesheet의 `transcriptome` 열이나 실행 옵션 `--transcriptome`으로 경로를 지정합니다.

자주 쓰는 10x Cell Ranger reference tarball:

```text
mouse GRCm39:
  https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz

human GRCh38:
  https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
```

예시:

```bash
mkdir -p reference
cd reference
wget -c https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz
tar -xzf refdata-gex-GRCm39-2024-A.tar.gz
```

현재 서버에 기존 SOLO conda 환경이 있으면 `scrnaseq_qc.sif`를 만들기 전에도 실행할 수 있습니다.

```bash
bash run.sh \
  --samplesheet input/samples.example.csv \
  -profile local_conda
```

기본 conda 환경 경로:

```bash
/home/Data_Drive_8TB/cwonkim99/.conda/envs/solo
```

## Samplesheet

예시 파일:

```bash
input/samples.example.csv
```

주요 열:

- `sample_id`: 결과 디렉터리와 barcode prefix로 사용할 샘플명
- `fastqs`: Cell Ranger 입력 FASTQ 디렉터리
- `sample`: Cell Ranger `--sample`
- `transcriptome`: 압축 해제된 Cell Ranger reference 디렉터리
- `start_step`: `cellranger`, `cellbender`, `qc`
- `raw_h5`: `start_step=cellbender`일 때 raw feature matrix H5
- `cellbender_h5`: `start_step=qc`일 때 CellBender cleaned H5
- `expected_doublet_rate`: 알 때만 입력. 예: `0.06`
- `expected_doublets`: rate 대신 명시적인 doublet 수를 줄 때 입력

## Run

FASTQ부터 실행:

```bash
bash run.sh \
  --samplesheet input/samples.example.csv \
  --outdir results \
  --threads 24 \
  --use_gpu true \
  --cuda_device 0
```

기존 CellBender 결과부터 QC/SOLO만 실행하려면 samplesheet에서 `start_step=qc`, `cellbender_h5=/path/to/cellbender_cleaned.h5`를 지정합니다.

## Outputs

샘플별 주요 출력:

- `01_cellranger/`: Cell Ranger output
- `02_cellbender/`: CellBender cleaned H5
- `03_scanpy_qc/`: `*_qc_full.h5ad`, `*_qc_hvg2000.h5ad`, `*.qc_metrics.csv`
- `04_solo/`: `soloed.h5ad`, SOLO score arrays/plots
- `05_singlets/`: `*_qc_solo_annotated_full.h5ad`, `*_qc_solo_singlets_full.h5ad`, doublet barcode CSV
- `scrnaseq_qc_summary.csv`: 전체 샘플 QC 요약

`*_qc_solo_annotated_full.h5ad`는 QC를 통과한 모든 cell을 유지하며 `obs`에 `solo_is_doublet`, `solo_softmax_scores`, `solo_logit_scores` 컬럼을 추가합니다. `*_qc_solo_singlets_full.h5ad`는 이 중 `solo_is_doublet == False`인 cell만 남긴 최종 분석용 matrix입니다.
