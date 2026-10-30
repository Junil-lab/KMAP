#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
// K_map scRNA-seq preprocessing and QC pipeline
//
// Steps:
//   01 cellranger_count      FASTQ -> raw/filtered 10x matrix
//   02 cellbender_background raw 10x H5 -> background-corrected H5
//   03 scanpy_qc_for_solo   mito/count/gene QC -> full H5AD + HVG H5AD
//   04 solo_doublets        SOLO doublet scoring
//   05 filter_singlets      annotate full QC H5AD with SOLO labels and write singlet-only H5AD
//   06 qc_merge             merged per-sample QC summary
//
// samplesheet CSV columns:
//   sample_id,fastqs,sample,transcriptome,start_step,raw_h5,cellbender_h5,
//   expected_cells,chemistry,create_bam,include_introns,cuda_device,
//   expected_doublet_rate,expected_doublets
//
// start_step:
//   cellranger  - requires fastqs and transcriptome
//   cellbender  - requires raw_h5
//   qc          - requires cellbender_h5
// ============================================================

def rowValue(row, key, fallback = null) {
    if (!row.containsKey(key)) {
        return fallback
    }
    def value = row[key]
    if (value == null) {
        return fallback
    }
    def text = value.toString().trim()
    return text ? text : fallback
}

def requireText(value, label, sampleId) {
    if (value == null || !value.toString().trim()) {
        throw new IllegalArgumentException("${sampleId}: '${label}' is required for this start_step")
    }
    return value.toString().trim()
}

// ============================================================
// STEP 01: Cell Ranger count
// ============================================================
process CELLRANGER_COUNT {
    tag "${meta.id}"
    label 'cellranger'
    publishDir { "${params.outdir}/${meta.id}/01_cellranger" }, mode: 'copy'

    input:
    tuple val(meta), path(fastqs, stageAs: "input_fastqs"), path(transcriptome, stageAs: "transcriptome")

    output:
    tuple val(meta), path("${meta.id}/outs/raw_feature_bc_matrix.h5"), emit: raw_h5
    path "${meta.id}/outs/filtered_feature_bc_matrix.h5", optional: true
    path "${meta.id}/outs/metrics_summary.csv", optional: true
    path "${meta.id}/outs/web_summary.html", optional: true
    path "${meta.id}/outs/possorted_genome_bam.bam", optional: true
    path "${meta.id}/outs/possorted_genome_bam.bam.bai", optional: true
    path "${meta.id}.cellranger.log"

    script:
    """
    set -euo pipefail

    CELLRANGER_CMD="${params.cellranger_cmd}"
    if ! command -v "\${CELLRANGER_CMD}" >/dev/null 2>&1; then
        if [ -n "${params.cellranger_bin}" ] && [ -x "${params.cellranger_bin}" ]; then
            CELLRANGER_CMD="${params.cellranger_bin}"
        else
            echo "[ERROR] cellranger command was not found in PATH." >&2
            echo "        Expected SIF: ${params.cellranger_sif}" >&2
            echo "        Build it with: bash def/build_containers.sh cellranger" >&2
            exit 2
        fi
    fi

    if [ ! -d "${transcriptome}" ]; then
        echo "[ERROR] transcriptome reference directory does not exist in task workdir: ${transcriptome}" >&2
        exit 2
    fi

    ARGS=(
        count
        --id="${meta.id}"
        --fastqs="${fastqs}"
        --sample="${meta.cr_sample}"
        --transcriptome="${transcriptome}"
        --localcores="${task.cpus}"
        --localmem="${meta.localmem_gb}"
        --create-bam="${meta.create_bam}"
    )

    if [ -n "${meta.chemistry}" ] && [ "${meta.chemistry}" != "auto" ]; then
        ARGS+=(--chemistry="${meta.chemistry}")
    fi
    if [ -n "${meta.expected_cells}" ] && [ "${meta.expected_cells}" != "null" ]; then
        ARGS+=(--expect-cells="${meta.expected_cells}")
    fi
    if [ "${meta.include_introns}" = "true" ]; then
        ARGS+=(--include-introns=true)
    fi

    "\${CELLRANGER_CMD}" "\${ARGS[@]}" 2>&1 | tee "${meta.id}.cellranger.log"
    """
}

// ============================================================
// STEP 02: CellBender remove-background
// ============================================================
process CELLBENDER_REMOVE_BACKGROUND {
    tag "${meta.id}"
    label 'gpu'
    publishDir { "${params.outdir}/${meta.id}/02_cellbender" }, mode: 'copy'

    input:
    tuple val(meta), path(raw_h5)

    output:
    tuple val(meta), path("${meta.id}.cellbender_cleaned.h5"), emit: cleaned_h5
    path "${meta.id}.cellbender.log"
    path "${meta.id}.cellbender_cleaned*.csv", optional: true
    path "${meta.id}.cellbender_cleaned*.pdf", optional: true

    script:
    """
    set -euo pipefail

    mkdir -p tmp_cellbender
    export TMPDIR="\$PWD/tmp_cellbender"
    export TEMPDIR="\$PWD/tmp_cellbender"
    export CUDA_VISIBLE_DEVICES="${meta.cuda_device}"

    CUDA_ARG=""
    if [ "${params.use_gpu}" = "true" ]; then
        CUDA_ARG="--cuda"
    fi

    EXTRA_ARGS=()
    if [ -n "${meta.cellbender_expected_cells}" ] && [ "${meta.cellbender_expected_cells}" != "null" ]; then
        EXTRA_ARGS+=(--expected-cells "${meta.cellbender_expected_cells}")
    fi
    if [ -n "${meta.cellbender_total_droplets_included}" ] && [ "${meta.cellbender_total_droplets_included}" != "null" ]; then
        EXTRA_ARGS+=(--total-droplets-included "${meta.cellbender_total_droplets_included}")
    fi
    if [ -n "${meta.cellbender_empty_drops}" ] && [ "${meta.cellbender_empty_drops}" != "null" ]; then
        EXTRA_ARGS+=(--empty-drop-training-fraction "${meta.cellbender_empty_drops}")
    fi

    cellbender remove-background \\
        --input "${raw_h5}" \\
        --output "${meta.id}.cellbender_cleaned.h5" \\
        --fpr "${params.cellbender_fpr}" \\
        --epochs "${params.cellbender_epochs}" \\
        --checkpoint-mins "${params.cellbender_checkpoint_mins}" \\
        \${CUDA_ARG} \\
        "\${EXTRA_ARGS[@]}" \\
        2>&1 | tee "${meta.id}.cellbender.log"
    """
}

// ============================================================
// STEP 03: Scanpy QC and HVG subset for SOLO
// ============================================================
process SCANPY_QC_FOR_SOLO {
    tag "${meta.id}"
    label 'standard'
    publishDir { "${params.outdir}/${meta.id}/03_scanpy_qc" }, mode: 'copy'

    input:
    tuple val(meta), path(cleaned_h5)
    path qc_script

    output:
    tuple val(meta),
          path("${meta.id}_qc_full.h5ad"),
          path("${meta.id}_qc_hvg${params.hvg_n_top}.h5ad"),
          path("${meta.id}.qc_metrics.csv"),
          emit: qc

    script:
    """
    set -euo pipefail
    mkdir -p numba_cache mpl_cache
    export NUMBA_CACHE_DIR="\$PWD/numba_cache"
    export MPLCONFIGDIR="\$PWD/mpl_cache"

    MT_ARGS=()
    if [ -n "${params.mt_threshold == null ? '' : params.mt_threshold}" ]; then
        MT_ARGS=(--mt-threshold "${params.mt_threshold}")
    fi

    python3 "${qc_script}" \\
        --input "${cleaned_h5}" \\
        --sample-id "${meta.id}" \\
        --output-full "${meta.id}_qc_full.h5ad" \\
        --output-hvg "${meta.id}_qc_hvg${params.hvg_n_top}.h5ad" \\
        --metrics "${meta.id}.qc_metrics.csv" \\
        "\${MT_ARGS[@]}" \\
        --min-genes "${params.min_genes}" \\
        --min-counts "${params.min_counts}" \\
        --min-cells-gene "${params.min_cells_gene}" \\
        --hvg-n-top "${params.hvg_n_top}"
    """
}

// ============================================================
// STEP 04: SOLO doublet detection
// ============================================================
process SOLO_DOUBLETS {
    tag "${meta.id}"
    label 'gpu'
    publishDir { "${params.outdir}/${meta.id}/04_solo" }, mode: 'copy'

    input:
    tuple val(meta), path(hvg_h5ad)
    path solo_json
    path score_script

    output:
    tuple val(meta), path("${meta.id}.solo"), emit: solo_dir

    script:
    """
    set -euo pipefail

    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export CUDA_VISIBLE_DEVICES="${meta.cuda_device}"
    mkdir -p numba_cache mpl_cache
    export NUMBA_CACHE_DIR="\$PWD/numba_cache"
    export MPLCONFIGDIR="\$PWD/mpl_cache"

    EXPECTED_ARGS=()
    if [ -n "${meta.expected_doublets}" ] && [ "${meta.expected_doublets}" != "null" ]; then
        EXPECTED_ARGS=(-e "${meta.expected_doublets}")
    elif [ -n "${meta.expected_doublet_rate}" ] && [ "${meta.expected_doublet_rate}" != "null" ]; then
        EXPECTED_N=\$(python3 - "${hvg_h5ad}" "${meta.expected_doublet_rate}" << 'PYEOF'
import math
import sys
import anndata as ad

adata = ad.read_h5ad(sys.argv[1], backed="r")
rate = float(sys.argv[2])
print(max(1, int(round(adata.n_obs * rate))))
PYEOF
)
        EXPECTED_ARGS=(-e "\${EXPECTED_N}")
    fi

    PLOT_ARG=""
    if [ "${params.solo_plot}" = "true" ]; then
        PLOT_ARG="-p"
    fi

    GPU_ARG=""
    if [ "${params.use_gpu}" = "true" ]; then
        GPU_ARG="-g"
    fi

    solo \\
        -j "${solo_json}" \\
        -d "${hvg_h5ad}" \\
        -o "${meta.id}.solo" \\
        -a \\
        \${GPU_ARG} \\
        \${PLOT_ARG} \\
        -r "${params.solo_doublet_ratio}" \\
        --set-reproducible-seed "${params.solo_seed}" \\
        "\${EXPECTED_ARGS[@]}" \\
        ${params.solo_extra_args}

    test -s "${meta.id}.solo/soloed.h5ad"
    python3 "${score_script}" --solo-dir "${meta.id}.solo"
    """
}

// ============================================================
// STEP 05: Remove SOLO doublets from full QC matrix
// ============================================================
process FILTER_SOLO_SINGLETS {
    tag "${meta.id}"
    label 'standard'
    publishDir { "${params.outdir}/${meta.id}/05_singlets" }, mode: 'copy'

    input:
    tuple val(meta), path(full_h5ad), path(solo_dir)
    path filter_script

    output:
    tuple val(meta),
          path("${meta.id}_qc_solo_annotated_full.h5ad"),
          path("${meta.id}_qc_solo_singlets_full.h5ad"),
          path("${meta.id}.solo_summary.csv"),
          path("${meta.id}.solo_doublet_barcodes.csv"),
          emit: singlets

    script:
    """
    set -euo pipefail

    python3 "${filter_script}" \\
        --full-h5ad "${full_h5ad}" \\
        --solo-h5ad "${solo_dir}/soloed.h5ad" \\
        --sample-id "${meta.id}" \\
        --annotated-output "${meta.id}_qc_solo_annotated_full.h5ad" \\
        --output "${meta.id}_qc_solo_singlets_full.h5ad" \\
        --summary "${meta.id}.solo_summary.csv" \\
        --doublet-barcodes "${meta.id}.solo_doublet_barcodes.csv"
    """
}

// ============================================================
// STEP 06: Merge QC summaries
// ============================================================
process QC_MERGE {
    label 'standard'
    publishDir params.outdir, mode: 'copy'

    input:
    path qc_metrics
    path solo_metrics
    path merge_script
    val expected_samples

    output:
    path "scrnaseq_qc_summary.csv"

    script:
    def expectedArgs = expected_samples.collect { "'${it}'" }.join(' ')
    """
    set -euo pipefail

    python3 "${merge_script}" \\
        --qc-files ${qc_metrics} \\
        --solo-files ${solo_metrics} \\
        --expected-samples ${expectedArgs} --output scrnaseq_qc_summary.csv
    """
}

// ============================================================
// WORKFLOW
// ============================================================
workflow {
    if (!params.samplesheet) {
        error "Missing required parameter: --samplesheet"
    }

    ch_samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sampleId = rowValue(row, 'sample_id')
            if (!sampleId) {
                throw new IllegalArgumentException("samplesheet row is missing 'sample_id'")
            }

            if (!(sampleId ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error "Invalid sample_id: ${sampleId}"
            }
            def startStep = rowValue(row, 'start_step', 'cellranger').toLowerCase()
            if (!(startStep in ['cellranger', 'cellbender', 'qc'])) {
                throw new IllegalArgumentException("${sampleId}: start_step must be one of cellranger, cellbender, qc")
            }

            def meta = [
                id: sampleId,
                cr_sample: rowValue(row, 'sample', sampleId),
                transcriptome: rowValue(row, 'transcriptome', params.transcriptome),
                localmem_gb: rowValue(row, 'localmem_gb', params.memory_gb),
                expected_cells: rowValue(row, 'expected_cells', params.cellranger_expect_cells),
                chemistry: rowValue(row, 'chemistry', params.chemistry),
                create_bam: rowValue(row, 'create_bam', params.create_bam.toString()),
                include_introns: rowValue(row, 'include_introns', params.include_introns.toString()),
                cuda_device: rowValue(row, 'cuda_device', params.cuda_device.toString()),
                cellbender_expected_cells: rowValue(row, 'cellbender_expected_cells', rowValue(row, 'expected_cells', null)),
                cellbender_total_droplets_included: rowValue(row, 'cellbender_total_droplets_included', params.cellbender_total_droplets_included),
                cellbender_empty_drops: rowValue(row, 'cellbender_empty_drops', params.cellbender_empty_drops),
                expected_doublet_rate: rowValue(row, 'expected_doublet_rate', params.solo_expected_doublet_rate),
                expected_doublets: rowValue(row, 'expected_doublets', null)
            ]

            tuple(meta,
                  rowValue(row, 'fastqs', null),
                  rowValue(row, 'raw_h5', null),
                  rowValue(row, 'cellbender_h5', null),
                  startStep)
        }
        .toList()
        .flatMap { entries ->
            def ids = entries.collect { it[0].id }
            if (!ids || ids.unique(false).size() != ids.size()) {
                error "Samplesheet must contain nonempty, unique sample IDs"
            }
            entries
        }

    ch_qc_script = Channel.value(file("${projectDir}/scripts/qc_for_solo.py", checkIfExists: true))
    ch_filter_script = Channel.value(file("${projectDir}/scripts/filter_solo_singlets.py", checkIfExists: true))
    ch_merge_script = Channel.value(file("${projectDir}/scripts/merge_qc.py", checkIfExists: true))
    ch_solo_json = Channel.value(file(params.solo_json, checkIfExists: true))
    ch_score_script = Channel.value(file("${projectDir}/scripts/correct_solo_scores.py", checkIfExists: true))

    ch_samples.branch {
        from_cellranger: it[4] == 'cellranger'
        from_cellbender: it[4] == 'cellbender'
        from_qc:         it[4] == 'qc'
    }.set { by_start }

    ch_cellranger_in = by_start.from_cellranger.map { meta, fastqs, raw_h5, cellbender_h5, startStep ->
        def fastqPath = requireText(fastqs, 'fastqs', meta.id)
        def transcriptomePath = requireText(meta.transcriptome, 'transcriptome', meta.id)
        tuple(meta, file(fastqPath, checkIfExists: true), file(transcriptomePath, checkIfExists: true))
    }

    ch_raw_direct = by_start.from_cellbender.map { meta, fastqs, raw_h5, cellbender_h5, startStep ->
        def rawPath = requireText(raw_h5, 'raw_h5', meta.id)
        tuple(meta, file(rawPath, checkIfExists: true))
    }

    ch_clean_direct = by_start.from_qc.map { meta, fastqs, raw_h5, cellbender_h5, startStep ->
        def cleanPath = requireText(cellbender_h5, 'cellbender_h5', meta.id)
        tuple(meta, file(cleanPath, checkIfExists: true))
    }

    CELLRANGER_COUNT(ch_cellranger_in)
    ch_raw_for_cellbender = CELLRANGER_COUNT.out.raw_h5.mix(ch_raw_direct)

    CELLBENDER_REMOVE_BACKGROUND(ch_raw_for_cellbender)
    ch_clean_for_qc = CELLBENDER_REMOVE_BACKGROUND.out.cleaned_h5.mix(ch_clean_direct)

    SCANPY_QC_FOR_SOLO(ch_clean_for_qc, ch_qc_script)

    ch_solo_in = SCANPY_QC_FOR_SOLO.out.qc.map { meta, full_h5ad, hvg_h5ad, metrics ->
        tuple(meta, hvg_h5ad)
    }
    SOLO_DOUBLETS(ch_solo_in, ch_solo_json, ch_score_script)

    ch_full_by_id = SCANPY_QC_FOR_SOLO.out.qc.map { meta, full_h5ad, hvg_h5ad, metrics ->
        tuple(meta.id, meta, full_h5ad, metrics)
    }
    ch_solo_by_id = SOLO_DOUBLETS.out.solo_dir.map { meta, solo_dir ->
        tuple(meta.id, meta, solo_dir)
    }
    ch_filter_in = ch_full_by_id
        .join(ch_solo_by_id)
        .map { id, metaA, full_h5ad, metrics, metaB, solo_dir ->
            tuple(metaA, full_h5ad, solo_dir)
        }
    FILTER_SOLO_SINGLETS(ch_filter_in, ch_filter_script)

    ch_qc_metrics = SCANPY_QC_FOR_SOLO.out.qc.map { meta, full_h5ad, hvg_h5ad, metrics -> metrics }.toList()
    ch_solo_metrics = FILTER_SOLO_SINGLETS.out.singlets.map { meta, annotated_h5ad, singlets_h5ad, summary, doublet_barcodes -> summary }.toList()
    ch_expected_samples = ch_samples.map { meta, fastqs, raw, clean, step -> meta.id }.toList()
    QC_MERGE(ch_qc_metrics, ch_solo_metrics, ch_merge_script, ch_expected_samples)
}
