# PBLONGREADWGS — K-MAP PacBio HiFi WGS 파이프라인

PacBio HiFi 리드(BAM) 하나로 정렬 → SNV/indel → SV → tandem repeat → phasing → CpG 메틸화까지 한 번에 돌리는 Nextflow (DSL2) + Singularity 파이프라인입니다. 레퍼런스는 GRCh38 (`GRCh38_no_alt_analysis_set`) 입니다.

## 어떻게 굴러가나

```
HiFi BAM (샘플 1개 = BAM 1개, 전장)
  │
  ├─ 01 JASMINE      kinetics(ip/pw) → 5mC 태그(MM/ML)   ※ start_step=jasmine 일 때만
  ├─ 02 PBMM2        GRCh38 정렬
  ├─ 03 DEEPVARIANT  SNV/indel (GPU)
  ├─ 04 SAWFISH      구조변이(SV)
  ├─ 05 TRGT         tandem repeat 유전형 (adotto + STRchive 카탈로그)
  ├─ 06 HIPHASE      SNV + SV + TR 동시 phasing, haplotag BAM
  └─ 07 CPG_TOOLS    CpG 메틸화 (combined / hap1 / hap2)
       └─ QC 요약 (샘플별 + 전체)
```

### 성별(sex) 처리

샘플시트의 `sex` 값(XX/XY)은 아래 세 도구에 반영됩니다. **필수 값**이며, 비어 있거나 다른 값이면 시작 전에 멈춥니다.

| 도구 | XX | XY |
|---|---|---|
| DeepVariant | 기본(이배체) | `--haploid_contigs chrX,chrY` + PAR BED (PAR 은 이배체 유지) |
| Sawfish | `expected_cn.hg38.XX.bed` | `expected_cn.hg38.XY.bed` + `--treat-single-copy-as-haploid` |
| TRGT | `--karyotype XX` | `--karyotype XY` (chrX PAR 도 haploid 로 처리됨 — TRGT 자체 동작) |

공여자 성별과 실제 염색체 수가 다른 세포주(예: X 를 하나 잃은 세포주)는 **데이터 기준**으로 넣는 것을 권장합니다 (X 1 copy, Y 없음 → XY).

## 준비물

| 항목 | 위치 | 비고 |
|---|---|---|
| Java 17+ | 시스템 또는 conda | Nextflow 26 은 Java 17 이상 필요 (`conda create -p .java17 -c conda-forge openjdk=17`) |
| Singularity / Apptainer | 시스템 | |
| GPU (선택) | – | DeepVariant 용. 없으면 `--use_gpu false` |
| `containers/kmap_pbwgs_v2.sif` | 별도 배포 | pbmm2, jasmine, sawfish, trgt, hiphase, pb-CpG-tools, samtools/bcftools — `def/build_containers.sh` 로 빌드 가능 |
| `containers/deepvariant_1.9.0-gpu.sif` | 별도 배포 | `google/deepvariant:1.9.0-gpu` |
| `resources/ref/GRCh38_no_alt_analysis_set.fna` (+ `.fai`, `.mmi`) | 별도 배포 | |
| `resources/ref/trgt/adotto_strchive_20250827.hg38.noFailTag.bed.gz` | 별도 배포 | PacBio HiFi-human-WGS-WDL v4.0.0 레퍼런스 번들(Zenodo)의 카탈로그에서 `;INCLUDE_FAIL_READS` 태그를 제거한 것 (WDL 과 동일 처리) |
| `resources/ref/sex/*.bed` | **저장소에 포함** | PAR BED, Sawfish expected copy number BED |

SIF·레퍼런스·BAM 은 용량 때문에 저장소에 없습니다 (`.gitignore`). 연구실 서버의 위치는 관리자에게 문의하세요.

## 실행 방법

### 1. 환경 확인 (처음 한 번)

```bash
git clone <this repo>
cd PBLONGREADWGS
bash setup.sh          # Java / Singularity / nextflow / SIF / 레퍼런스를 점검, 모두 [OK] 이면 준비 완료
```

### 2. 샘플시트 작성

`input/samples.csv` 를 복사해서 작성합니다.

```csv
sample_id,bam,sex,start_step
SampleA,/절대경로/SampleA.hifi_reads.bam,XX,pbmm2
SampleB,/절대경로/SampleB.hifi_reads.bam,XY,pbmm2
```

| 열 | 값 |
|---|---|
| `sample_id` | 결과 폴더 이름 |
| `bam` | 절대경로. **파이프라인 폴더 안**에 있어야 함 → 밖에 있으면 `ln -s /외부/경로/x.bam input/bam/` |
| `sex` | `XX` 또는 `XY` (필수) |
| `start_step` | `jasmine` — kinetics(ip/pw)만 있고 메틸화 태그가 없는 BAM<br>`pbmm2` — 메틸화 태그(MM/ML)가 이미 있는 HiFi BAM (대부분 이것)<br>`deepvariant` / `hiphase` — 중간부터 재시작 |

### 3. 실행

```bash
bash run.sh --samplesheet input/my_samples.csv --outdir /결과/경로 --threads 16

# 자주 쓰는 옵션
#   -resume              중단된 지점부터 이어서 (캐시 재사용)
#   --use_gpu false      GPU 없이 DeepVariant
#   --vcf_dir <dir>      start_step=hiphase 일 때 기존 VCF 위치
#   export NXF_WORK=/빠른/로컬/디스크/work   # 중간 파일 위치 변경 (기본: 파이프라인 폴더/.nf_work)
```

공유 서버에서는 코어 수를 제한해서 돌리세요: `taskset -c 0-4 bash run.sh ... --threads 5`

## 결과

```
<outdir>/<sample_id>/
├── 01_methyl_bam/        (jasmine 단계를 돌린 경우)
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
<outdir>/qc_summary/      전 샘플 QC 표
<outdir>/pipeline_info/   timeline / report / trace, 실행 로그
```

## 참고 / 주의

- **실행 시간 감**: HG002 chr20+X+Y (전장의 약 5%) 를 5 코어로 돌렸을 때 3시간 16분, 그중 DeepVariant 가 76%. 전장 30x 는 코어를 충분히 주는 것을 권장합니다.
- **메틸화**: BAM 에 들어 있는 메틸화 태그를 집계합니다. 5hmC/6mA 는 시퀀싱 단계에서 태그가 만들어졌고 kinetics 가 남아 있어야만 얻을 수 있습니다.
- **TRGT 로그 확인**: TRGT 는 일부 locus 가 실패해도 exit 0 으로 끝납니다. `05_tandem_repeats` 작업 로그의 `[ERROR]` 줄을 확인하세요.
- **제거된 단계**: 변이 어노테이션(VEP, AnnotSV)과 sTELLeR(TE 검출)는 이 버전에서 제거되었습니다.
- **비교 기준**: 도구 구성은 PacBio 공식 [HiFi-human-WGS-WDL](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) v4.0.0 과 대부분 같고, 버전은 DeepVariant 1.9.0 / TRGT 5.0.0 / HiPhase 1.6.0 / Sawfish 2.2.1 / pb-CpG-tools 3.0.0 입니다.
