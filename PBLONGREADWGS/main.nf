#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
// K_map PacBio HiFi WGS Pipeline v2
// Fully self-contained: every tool runs inside a Singularity container
//
// Steps:
//   01 methyl_bam        ← JASMINE  (kinetics BAMs only)
//   02 aligned_bam       ← PBMM2
//   03 small_variants    ← DEEPVARIANT (GPU)
//   04 structural_var    ← SAWFISH
//   05 tandem_repeats    ← TRGT
//   06 phasing           ← HIPHASE (joint SNV+SV+TR)
//   07 CpG_methylation   ← pb-CpG-tools
//   (VEP / AnnotSV annotation and sTELLeR TE detection steps have been removed)
//
// samplesheet CSV: sample_id, bam, sex, start_step
//   sex: XX | XY  (required; used for TRGT karyotype, DeepVariant haploid calling, Sawfish expected CN)
//   start_step: jasmine | pbmm2 | deepvariant | hiphase
// ============================================================

// ============================================================
// STEP 01: JASMINE — kinetics tags → MM/ML methylation tags
// ============================================================
process JASMINE {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/01_methyl_bam" }, mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.methyl.bam"), path("${meta.id}.methyl.bam.bai"), emit: methyl_bam
    path "${meta.id}.jasmine.log"

    script:
    """
    jasmine -j ${params.threads} ${bam} ${meta.id}.methyl.bam 2>&1 | tee ${meta.id}.jasmine.log
    samtools index ${meta.id}.methyl.bam
    """
}

// ============================================================
// STEP 02: PBMM2 — alignment to GRCh38_no_alt_analysis_set
// ============================================================
process PBMM2 {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/02_aligned_bam" }, mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.aligned.sorted.bam"), path("${meta.id}.aligned.sorted.bam.bai"), emit: aligned_bam
    path "${meta.id}.pbmm2.log"
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat

    script:
    """
    pbmm2 align \\
        ${params.ref_mmi} \\
        ${bam} \\
        ${meta.id}.aligned.sorted.bam \\
        --preset HIFI \\
        --sort \\
        --sample ${meta.id} \\
        -j ${params.threads} \\
        --log-file ${meta.id}.pbmm2.log \\
        --bam-index BAI

    samtools flagstat ${meta.id}.aligned.sorted.bam > ${meta.id}.flagstat.txt
    """
}

// ============================================================
// STEP 03: DEEPVARIANT — SNV/Indel calling (GPU)
// Reference: GRCh38_no_alt_analysis_set (no unplaced contigs → avoids exit 245)
// ============================================================
process DEEPVARIANT {
    tag "${meta.id}"
    label 'high_mem'
    container params.dv_sif
    containerOptions params.use_gpu ? "--nv" : ""
    publishDir { "${params.outdir}/${meta.id}/03_small_variants" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.deepvariant.vcf.gz"), path("${meta.id}.deepvariant.vcf.gz.tbi"), emit: snv_vcf
    path "${meta.id}.deepvariant.g.vcf.gz"
    path "${meta.id}.deepvariant.g.vcf.gz.tbi"

    script:
    // XY: call chrX/chrY outside the PARs as haploid
    def haploid_args = meta.sex == 'XY' ? "--haploid_contigs=chrX,chrY --par_regions_bed=${params.par_bed}" : ""
    """
    mkdir -p tmp_dv
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=PACBIO \\
        --ref=${params.ref} \\
        --reads=${bam} \\
        --output_vcf=${meta.id}.deepvariant.vcf.gz \\
        --output_gvcf=${meta.id}.deepvariant.g.vcf.gz \\
        --intermediate_results_dir=tmp_dv \\
        --num_shards=${params.threads} \\
        --vcf_stats_report=false \\
        ${haploid_args}
    """
}

// ============================================================
// STEP 04: SAWFISH — SV + CNV calling
// ============================================================
process SAWFISH {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/04_structural_var" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.sv.vcf.gz"), path("${meta.id}.sv.vcf.gz.tbi"), emit: sv_vcf
    tuple val(meta), path("${meta.id}.sv.run.stats.json"), emit: stats

    script:
    def expected_cn = meta.sex == 'XY' ? params.expected_cn_xy : params.expected_cn_xx
    """
    sawfish discover \\
        --bam ${bam} \\
        --ref ${params.ref} \\
        --expected-cn ${expected_cn} \\
        --output-dir discover \\
        --threads ${params.threads}

    # regions with expected CN 1 (chrX/chrY in XY) are genotyped as haploid SVs too
    sawfish joint-call \\
        --sample discover \\
        --treat-single-copy-as-haploid \\
        --output-dir call \\
        --threads ${params.threads}

    cp call/genotyped.sv.vcf.gz     ${meta.id}.sv.vcf.gz
    cp call/genotyped.sv.vcf.gz.tbi ${meta.id}.sv.vcf.gz.tbi
    cp call/run.stats.json          ${meta.id}.sv.run.stats.json
    """
}

// ============================================================
// STEP 05: TRGT — Tandem Repeat Genotyping
// ============================================================
process TRGT {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/05_tandem_repeats" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.trgt.vcf.gz"), path("${meta.id}.trgt.vcf.gz.tbi"), emit: trgt_vcf
    path "${meta.id}.trgt.spanning.sorted.bam"
    path "${meta.id}.trgt.spanning.sorted.bam.bai"

    script:
    """
    trgt genotype \\
        --threads ${params.threads} \\
        --genome ${params.ref} \\
        --reads ${bam} \\
        --repeats ${params.trgt_catalog} \\
        --output-prefix ${meta.id}.trgt \\
        --karyotype ${meta.sex} \\
        --preset wgs

    samtools sort -@ ${params.threads} \\
        -o ${meta.id}.trgt.spanning.sorted.bam \\
        ${meta.id}.trgt.spanning.bam

    samtools index ${meta.id}.trgt.spanning.sorted.bam
    bcftools sort -O z -o ${meta.id}.trgt.sorted.vcf.gz ${meta.id}.trgt.vcf.gz
    mv ${meta.id}.trgt.sorted.vcf.gz ${meta.id}.trgt.vcf.gz
    bcftools index -f --tbi ${meta.id}.trgt.vcf.gz
    """
}

// ============================================================
// STEP 06: HIPHASE — joint SNV + SV + TR phasing
// ============================================================
process HIPHASE {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/06a_phased_bam" }, mode: 'copy', pattern: "*.haplotagged.bam*"
    publishDir { "${params.outdir}/${meta.id}/06a_phased_bam" }, mode: 'copy', pattern: "*.blocks.tsv"
    publishDir { "${params.outdir}/${meta.id}/06a_phased_bam" }, mode: 'copy', pattern: "*.stats.csv"
    publishDir { "${params.outdir}/${meta.id}/06a_phased_bam" }, mode: 'copy', pattern: "*.summary.tsv"
    publishDir { "${params.outdir}/${meta.id}/06b_phased_SNP" }, mode: 'copy', pattern: "*.snv.phased*"
    publishDir { "${params.outdir}/${meta.id}/06c_phased_SV" },  mode: 'copy', pattern: "*.sv.phased*"
    publishDir { "${params.outdir}/${meta.id}/06d_phased_TR" },  mode: 'copy', pattern: "*.tr.phased*"

    input:
    tuple val(meta), path(bam), path(bai),
          path(sv_vcf), path(sv_tbi),
          path(tr_vcf), path(tr_tbi),
          path(snv_vcf), path(snv_tbi)

    output:
    tuple val(meta), path("${meta.id}.haplotagged.bam"), path("${meta.id}.haplotagged.bam.bai"), emit: haplotagged_bam
    tuple val(meta), path("${meta.id}.snv.phased.vcf.gz"), path("${meta.id}.snv.phased.vcf.gz.tbi"), emit: snv_phased
    tuple val(meta), path("${meta.id}.sv.phased.vcf.gz"),  path("${meta.id}.sv.phased.vcf.gz.tbi"),  emit: sv_phased
    tuple val(meta), path("${meta.id}.tr.phased.vcf.gz"),  path("${meta.id}.tr.phased.vcf.gz.tbi"),  emit: tr_phased
    tuple val(meta), path("${meta.id}.summary.tsv"), emit: summary
    path "${meta.id}.blocks.tsv"
    path "${meta.id}.stats.csv"

    script:
    """
    bcftools view -f PASS -Oz -o snv.pass.vcf.gz ${snv_vcf}
    bcftools index --tbi snv.pass.vcf.gz
    bcftools view -f PASS -Oz -o sv.pass.vcf.gz ${sv_vcf}
    bcftools index --tbi sv.pass.vcf.gz

    hiphase \\
        --threads ${params.threads} \\
        --reference ${params.ref} \\
        --bam ${bam} \\
        --vcf snv.pass.vcf.gz \\
        --output-vcf ${meta.id}.snv.phased.vcf.gz \\
        --vcf sv.pass.vcf.gz \\
        --output-vcf ${meta.id}.sv.phased.vcf.gz \\
        --vcf ${tr_vcf} \\
        --output-vcf ${meta.id}.tr.phased.vcf.gz \\
        --output-bam ${meta.id}.haplotagged.bam \\
        --summary-file ${meta.id}.summary.tsv \\
        --blocks-file  ${meta.id}.blocks.tsv \\
        --stats-file   ${meta.id}.stats.csv

    samtools index ${meta.id}.haplotagged.bam
    """
}

// ============================================================
// STEP 07: CPG_TOOLS — CpG methylation (haplotype-aware)
// ============================================================
process CPG_TOOLS {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/07_CpG_methylation" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.combined.bed.gz"), emit: cpg_bed
    path "${meta.id}.combined.bed.gz.tbi"
    path "${meta.id}.hap1.bed.gz",     optional: true
    path "${meta.id}.hap1.bed.gz.tbi", optional: true
    path "${meta.id}.hap2.bed.gz",     optional: true
    path "${meta.id}.hap2.bed.gz.tbi", optional: true

    script:
    """
    aligned_bam_to_cpg_scores \\
        --bam ${bam} \\
        --ref ${params.ref} \\
        --output-prefix ${meta.id} \\
        --pileup-mode model \\
        --modsites-mode reference \\
        --threads ${params.threads}

    tabix -f -p bed ${meta.id}.combined.bed.gz
    [ -f ${meta.id}.hap1.bed.gz ] && tabix -f -p bed ${meta.id}.hap1.bed.gz || true
    [ -f ${meta.id}.hap2.bed.gz ] && tabix -f -p bed ${meta.id}.hap2.bed.gz || true
    """
}

// ============================================================
// QC_COLLECT / QC_MERGE
// ============================================================
process QC_COLLECT {
    tag "${meta.id}"
    container params.main_sif
    errorStrategy 'terminate'
    publishDir "${params.outdir}/qc_summary", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(sv_vcf), path(sv_tbi), path(hiphase_summary), path(cpg_bed)
    path qc_script

    output:
    path "${meta.id}.qc.csv", emit: qc_csv

    script:
    """
    set -euo pipefail
    samtools flagstat "${bam}" > alignment.flagstat.txt
    bcftools view -H "${sv_vcf}" | wc -l > sv_count.txt
    python3 "${qc_script}" collect --sample "${meta.id}" \
        --flagstat alignment.flagstat.txt --sv-count sv_count.txt \
        --summary "${hiphase_summary}" --output "${meta.id}.qc.csv"
    """
}

process QC_MERGE {
    container params.main_sif
    errorStrategy 'terminate'
    publishDir params.outdir, mode: 'copy'

    input:
    path csv_files
    val expected_samples
    path qc_script

    output:
    path "qc_summary.csv"

    script:
    def expectedArgs = expected_samples.collect { "'${it}'" }.join(' ')
    """
    python3 "${qc_script}" merge --files ${csv_files} \
        --expected-samples ${expectedArgs} --output qc_summary.csv
    """
}

// ============================================================
// MEMORY_MONITOR
// ============================================================
process MEMORY_MONITOR {
    // Runs after QC_MERGE finishes (= end of the whole pipeline).
    // Memory during the run is recorded by the background loop in run.sh (system_memory.tsv).
    // This process records a per-sample snapshot at completion + the list of processes that ran.
    tag "${meta.id}"
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/logs" }, mode: 'copy'

    input:
    val meta
    val done_signal   // QC_MERGE output file name (completion signal)

    output:
    path "pipeline_done.tsv"

    script:
    """
    TS=\$(date '+%Y-%m-%d %H:%M:%S')
    MEM=\$(free -m | awk '/^Mem:/{print \$2"\\t"\$3"\\t"\$4}')
    echo -e "sample_id\ttimestamp\ttotal_MB\tused_MB\tfree_MB\tqc_file" > pipeline_done.tsv
    echo -e "${meta.id}\\t\${TS}\\t\${MEM}\\t${done_signal}" >> pipeline_done.tsv
    """
}

// ============================================================
// WORKFLOW
// ============================================================
workflow {

    ch_input = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sex        = (row.sex ?: '').trim().toUpperCase()
            if (!(sex in ['XX', 'XY'])) {
                error "samplesheet: sex of ${row.sample_id} is '${row.sex}' — must be XX or XY"
            }
            if (!(row.sample_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error "Invalid sample_id: ${row.sample_id}"
            }
            def meta       = [id: row.sample_id, sex: sex]
            def bam        = file(row.bam)
            def start_step = row.start_step ?: 'pbmm2'
            if (!(start_step in ['jasmine', 'pbmm2', 'deepvariant', 'hiphase'])) {
                error "Invalid start_step for ${meta.id}: ${start_step}"
            }
            tuple(meta, bam, start_step)
        }
        .toList()
        .flatMap { entries ->
            def ids = entries.collect { it[0].id }
            if (!ids || ids.unique(false).size() != ids.size()) {
                error "Samplesheet must contain nonempty, unique sample IDs"
            }
            entries
        }

    ch_input.branch {
        jasmine:     it[2] == 'jasmine'
        pbmm2:       it[2] == 'pbmm2'
        deepvariant: it[2] == 'deepvariant'
        hiphase:     it[2] == 'hiphase'
    }.set { ch_branched }

    // Memory monitor runs after QC_MERGE (captures the end of the whole pipeline)

    // STEP 01: JASMINE
    ch_for_jasmine = ch_branched.jasmine.map { meta, bam, step -> tuple(meta, bam) }
    jasmine_out    = JASMINE(ch_for_jasmine)
    ch_methyl_for_pbmm2 = jasmine_out.methyl_bam.map { meta, bam, bai -> tuple(meta, bam) }

    // STEP 02: PBMM2
    ch_for_pbmm2 = ch_branched.pbmm2
        .map { meta, bam, step -> tuple(meta, bam) }
        .mix(ch_methyl_for_pbmm2)
    pbmm2_out = PBMM2(ch_for_pbmm2)

    // ch_aligned: PBMM2 output → reused by several steps
    pbmm2_out.aligned_bam
        .multiMap { meta, bam, bai ->
            for_dv:       tuple(meta, bam, bai)
            for_sawfish:  tuple(meta, bam, bai)
            for_trgt:     tuple(meta, bam, bai)
            for_hiphase:  tuple(meta, bam, bai)
        }
        .set { ch_aligned }

    // supply the BAM when start_step=deepvariant
    ch_dv_bam = ch_branched.deepvariant
        .map { meta, bam, step -> tuple(meta, bam, file("${bam}.bai")) }
    ch_dv_bam
        .multiMap { meta, bam, bai ->
            for_dv:       tuple(meta, bam, bai)
            for_sawfish:  tuple(meta, bam, bai)
            for_trgt:     tuple(meta, bam, bai)
            for_hiphase:  tuple(meta, bam, bai)
        }
        .set { ch_dv_entry }

    // STEP 03: DEEPVARIANT
    ch_for_dv = ch_dv_entry.for_dv.mix(ch_aligned.for_dv)
    dv_out    = DEEPVARIANT(ch_for_dv)

    // STEP 04: SAWFISH — runs after DV (ordering enforced by joining on the DV VCF)
    ch_for_sawfish = ch_dv_entry.for_sawfish
        .mix(ch_aligned.for_sawfish)
        .join(dv_out.snv_vcf, by: 0)
        .map { meta, bam, bai, snv_vcf, snv_tbi -> tuple(meta, bam, bai) }
    sawfish_out = SAWFISH(ch_for_sawfish)

    // STEP 05: TRGT — runs after SAWFISH (ordering enforced by joining on the SAWFISH stats)
    ch_for_trgt = ch_dv_entry.for_trgt
        .mix(ch_aligned.for_trgt)
        .join(sawfish_out.stats, by: 0)
        .map { meta, bam, bai, stats -> tuple(meta, bam, bai) }
    trgt_out = TRGT(ch_for_trgt)

    // STEP 06: HIPHASE — runs after DV + SAWFISH + TRGT all finish
    ch_for_hiphase = ch_dv_entry.for_hiphase
        .mix(ch_aligned.for_hiphase)
        .join(sawfish_out.sv_vcf, by: 0)
        .join(trgt_out.trgt_vcf,  by: 0)
        .join(dv_out.snv_vcf,     by: 0)
        .map { meta, bam, bai, sv_vcf, sv_tbi, tr_vcf, tr_tbi, snv_vcf, snv_tbi ->
            tuple(meta, bam, bai, sv_vcf, sv_tbi, tr_vcf, tr_tbi, snv_vcf, snv_tbi)
        }

    ch_hiphase_extra = ch_branched.hiphase
        .map { meta, bam, step ->
            def bai      = file("${bam}.bai")
            def vcf_base = params.vcf_dir ? "${params.vcf_dir}/${meta.id}" : "${params.outdir}/${meta.id}"
            def snv_vcf  = file("${vcf_base}/03_small_variants/${meta.id}.deepvariant.vcf.gz")
            def snv_tbi  = file("${vcf_base}/03_small_variants/${meta.id}.deepvariant.vcf.gz.tbi")
            def sv_vcf   = file("${vcf_base}/04_structural_var/${meta.id}.sv.vcf.gz")
            def sv_tbi   = file("${vcf_base}/04_structural_var/${meta.id}.sv.vcf.gz.tbi")
            def tr_vcf   = file("${vcf_base}/05_tandem_repeats/${meta.id}.trgt.vcf.gz")
            def tr_tbi   = file("${vcf_base}/05_tandem_repeats/${meta.id}.trgt.vcf.gz.tbi")
            tuple(meta, bam, bai, sv_vcf, sv_tbi, tr_vcf, tr_tbi, snv_vcf, snv_tbi)
        }

    ch_all_hiphase = ch_for_hiphase.mix(ch_hiphase_extra)
    hiphase_out = HIPHASE(ch_all_hiphase)

    // STEP 07: CPG_TOOLS — haplotagged BAM → CpG methylation
    cpg_out = CPG_TOOLS(hiphase_out.haplotagged_bam)

    // Recompute alignment statistics for every entry route. Wait for CpG too.
    ch_qc_script = Channel.value(file("${projectDir}/scripts/qc_summary.py", checkIfExists: true))
    ch_expected_samples = ch_input.map { meta, bam, step -> meta.id }.toList()
    ch_qc_input = ch_all_hiphase
        .map { meta, bam, bai, sv_vcf, sv_tbi, tr_vcf, tr_tbi, snv_vcf, snv_tbi ->
            tuple(meta, bam, bai, sv_vcf, sv_tbi)
        }
        .join(hiphase_out.summary, by: 0)
        .join(cpg_out.cpg_bed, by: 0)
    qc_per_sample = QC_COLLECT(ch_qc_input, ch_qc_script)
    qc_merge_out = QC_MERGE(qc_per_sample.qc_csv.toList(), ch_expected_samples, ch_qc_script)

    // MEMORY_MONITOR: runs after QC_MERGE (= end of the whole pipeline)
    ch_meta      = ch_input.map { meta, bam, step -> meta }
    ch_done_sig  = qc_merge_out.map { f -> f.name }
    MEMORY_MONITOR(ch_meta.combine(ch_done_sig).map { meta, sig -> meta },
                   ch_meta.combine(ch_done_sig).map { meta, sig -> sig }.first())
}
