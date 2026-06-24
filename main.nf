#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
// K_map PacBio HiFi WGS Pipeline v2
// 완전 독립형 — 모든 툴이 Singularity 컨테이너 안에 포함
//
// 단계:
//   01 methyl_bam        ← JASMINE  (kinetics BAM만)
//   02 aligned_bam       ← PBMM2
//   03 small_variants    ← DEEPVARIANT (GPU)
//   04 structural_var    ← SAWFISH
//   05 tandem_repeats    ← TRGT
//   06 phasing           ← HIPHASE (SNV+SV+TR 동시)
//   07 CpG_methylation   ← pb-CpG-tools
//   08 SNP_annotation    ← VEP 115
//   09 SV_annotation     ← AnnotSV
//   10 TR_annotation     ← AnnotSV
//   11 TE_detection      ← sTELLeR
//   12 TE_annotation     ← AnnotSV
//
// samplesheet CSV: sample_id, bam, start_step
//   start_step: jasmine | pbmm2 | deepvariant | hiphase
// ============================================================

// ============================================================
// STEP 01: JASMINE — kinetics 태그 → MM/ML 메틸화 태그 변환
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
// STEP 02: PBMM2 — GRCh38_no_alt_analysis_set 정렬
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
// 레퍼런스: GRCh38_no_alt_analysis_set (unplaced contig 제거 → exit 245 방지)
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
        --vcf_stats_report=false
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
    """
    sawfish discover \\
        --bam ${bam} \\
        --ref ${params.ref} \\
        --output-dir discover \\
        --threads ${params.threads}

    sawfish joint-call \\
        --sample discover \\
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
        --karyotype ${params.karyotype ?: 'XY'} \\
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
// STEP 06: HIPHASE — SNV + SV + TR 동시 phasing
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
// STEP 07: CPG_TOOLS — CpG 메틸화 (haplotype-aware)
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
// STEP 08: VEP — SNV annotation (VEP SIF 사용)
// ============================================================
process SNP_ANNOTATION {
    tag "${meta.id}"
    label 'high_mem'
    container params.vep_sif
    publishDir { "${params.outdir}/${meta.id}/10_SNP_annotation" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.snv.phased.vep115.merged.vcf.gz"), emit: vep_vcf
    path "${meta.id}.snv.phased.vep115.merged.vcf.gz.stdout.log"

    script:
    """
    vep \\
        --offline \\
        --cache \\
        --dir_cache ${params.vep_cache} \\
        --species homo_sapiens \\
        --assembly GRCh38 \\
        --fasta ${params.vep_ref_1kg} \\
        --merged \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --input_file ${vcf} \\
        --output_file ${meta.id}.snv.phased.vep115.merged.vcf.gz \\
        --force_overwrite \\
        --fork ${params.threads} \\
        --symbol --mane --canonical --biotype --numbers \\
        --variant_class --hgvs --check_existing \\
        --af --af_gnomade --af_gnomadg --max_af \\
        --sift b --polyphen b \\
        --no_stats --no_progress \\
        > ${meta.id}.snv.phased.vep115.merged.vcf.gz.stdout.log 2>&1
    """
}

// ============================================================
// STEP 09: ANNOTSV_SV — SV annotation
// ============================================================
process ANNOTSV_SV {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/11_SV_annotation" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.sv.annotsv.tsv"), emit: sv_annot
    path "${meta.id}.sv.annotsv.unannotated.tsv"

    script:
    """
    AnnotSV \\
        -SVinputFile    ${vcf} \\
        -outputFile     ${meta.id}.sv.annotsv \\
        -annotationsDir ${params.annotsv_db} \\
        -genomeBuild    GRCh38 \\
        -annotationMode full
    find . -name "${meta.id}.sv.annotsv.tsv" ! -path "./${meta.id}.sv.annotsv.tsv" -exec mv {} . \\;
    find . -name "${meta.id}.sv.annotsv.unannotated.tsv" ! -path "./${meta.id}.sv.annotsv.unannotated.tsv" -exec mv {} . \\;
    """
}

// ============================================================
// STEP 10: ANNOTSV_TR — TR annotation
// ============================================================
process ANNOTSV_TR {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/12_TR_annotation" }, mode: 'copy'
    // errorStrategy: nextflow.config withName:'STELLER|ANNOTSV_TR' 에서 관리 (retry 1 → ignore)

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.tr.annotsv.tsv"), emit: tr_annot, optional: true
    path "${meta.id}.tr.annotsv.unannotated.tsv", optional: true

    script:
    """
    AnnotSV \\
        -SVinputFile    ${vcf} \\
        -outputFile     ${meta.id}.tr.annotsv \\
        -annotationsDir ${params.annotsv_db} \\
        -genomeBuild    GRCh38 \\
        -annotationMode full
    find . -name "${meta.id}.tr.annotsv.tsv" ! -path "./${meta.id}.tr.annotsv.tsv" -exec mv {} . \\; 2>/dev/null || true
    find . -name "${meta.id}.tr.annotsv.unannotated.tsv" ! -path "./${meta.id}.tr.annotsv.unannotated.tsv" -exec mv {} . \\; 2>/dev/null || true
    """
}

// ============================================================
// STEP 11: STELLER — TE detection
// ============================================================
process STELLER {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/08_transposable_elements" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.steller.sorted.vcf.gz"), path("${meta.id}.steller.sorted.vcf.gz.tbi"), emit: te_vcf
    path "${meta.id}.steller.log"

    script:
    """
    # minimap2는 컨테이너 PATH에 없을 수 있으므로 명시
    export PATH=/opt/conda/envs/kmap/bin:\$PATH

    # steller -o 는 sample prefix — 실제 출력 파일은 <prefix>_repeats.vcf
    STELLER_EXIT=0
    python3 ${params.steller_py} \\
        -R ${params.ref} \\
        -tf ${params.te_fasta} \\
        -b ${bam} \\
        -s pb -r 3 \\
        -o ${meta.id} \\
        2>&1 | tee ${meta.id}.steller.log || STELLER_EXIT=\$?

    # steller 출력: ${meta.id}_repeats.vcf
    RAW_VCF="${meta.id}_repeats.vcf"

    # sTELLeR 실패 시: 빈 VCF 생성 후 경고 (downstream은 계속, retry는 config에서 관리)
    if [ \${STELLER_EXIT} -ne 0 ] || [ ! -s "\${RAW_VCF}" ]; then
        echo "WARNING: sTELLeR failed or produced empty output for ${meta.id} (exit=\${STELLER_EXIT}). Creating empty VCF." >&2
        printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${meta.id}\n' > "\${RAW_VCF}"
    fi

    # TE 이름 콜론 fix + FORMAT 헤더 정의 추가 (bcftools 파싱 에러 방지)
    # HT 필드에 콤마(1,2,0) 포함 케이스 때문에 String 타입으로 헤더 정의 필요
    awk 'BEGIN{OFS="\t"}
    /^##fileformat/{print; next}
    /^##FORMAT/{print; next}
    /^##/{print; next}
    /^#CHROM/{
        print "##FORMAT=<ID=TE,Number=1,Type=String,Description=\"Transposable element\">"
        print "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Transposon start\">"
        print "##FORMAT=<ID=HT,Number=.,Type=String,Description=\"Haplotag\">"
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        print "##INFO=<ID=MEINFO,Number=.,Type=String,Description=\"Mobile element info\">"
        print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">"
        print "##INFO=<ID=SVTPE,Number=1,Type=String,Description=\"SV type (steller)\">"
        print; next
    }
    {
        n=split(\$10,v,":");
        if(n==5){ \$10=v[1]"_"v[2]":"v[3]":"v[4]":"v[5] };
        print
    }' "\${RAW_VCF}" > ${meta.id}.steller.fixed.vcf

    bcftools sort -Oz -o ${meta.id}.steller.sorted.vcf.gz ${meta.id}.steller.fixed.vcf
    bcftools index --tbi ${meta.id}.steller.sorted.vcf.gz
    """
}

// ============================================================
// STEP 12: ANNOTSV_TE — TE annotation
// ============================================================
process ANNOTSV_TE {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/09_TE_annotation" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.steller.annotated.tsv"), emit: te_annot
    path "${meta.id}.steller.unannotated.tsv", optional: true

    script:
    """
    AnnotSV \\
        -SVinputFile    ${vcf} \\
        -outputFile     ${meta.id}.steller.annotated \\
        -annotationsDir ${params.annotsv_db} \\
        -genomeBuild    GRCh38 \\
        -annotationMode full
    find . -name "${meta.id}.steller.annotated.tsv" ! -path "./${meta.id}.steller.annotated.tsv" -exec mv {} . \\;
    find . -name "${meta.id}.steller.unannotated.tsv" ! -path "./${meta.id}.steller.unannotated.tsv" -exec mv {} . \\; 2>/dev/null || true
    """
}

// ============================================================
// QC_COLLECT / QC_MERGE
// ============================================================
process QC_COLLECT {
    tag "${meta.id}"
    container params.main_sif
    publishDir "${params.outdir}/qc_summary", mode: 'copy'

    input:
    tuple val(meta), path(flagstat), path(sv_stats), path(hiphase_summary)

    output:
    path "${meta.id}.qc.csv", emit: qc_csv

    script:
    """
    python3 << 'PYEOF'
import csv, json, os
sample = "${meta.id}"
qc = {"sample_id": sample, "status": "success", "error": ""}

try:
    lines = open("${flagstat}").readlines()
    total  = int(next(l.split()[0] for l in lines if "in total" in l))
    mapped = int(next(l.split()[0] for l in lines if "mapped (" in l))
    qc["mapping_rate"] = f"{mapped/total*100:.2f}%"
    qc["total_reads"]  = total
    qc["mapped_reads"] = mapped
except Exception as e:
    qc["mapping_rate"] = qc["total_reads"] = qc["mapped_reads"] = "NA"
    qc["error"] += f"flagstat:{e};"

try:
    d = json.load(open("${sv_stats}"))
    qc["total_sv"] = d["score_stats"]["vcf_output_record_count"]
except Exception as e:
    qc["total_sv"] = "NA"; qc["error"] += f"sv_stats:{e};"

try:
    phased, total_v, blocks, ng50 = 0, 0, 0, 0
    with open("${hiphase_summary}") as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            phased  += int(row.get("num_phased", 0))
            total_v += int(row.get("num_variants", 0))
            blocks  += int(row.get("num_blocks", 0))
            ng50     = max(ng50, int(row.get("block_ng50", 0)))
    qc["phased_variants"] = phased
    qc["total_variants"]  = total_v
    qc["phase_blocks"]    = blocks
    qc["block_n50"]       = ng50
except Exception as e:
    qc["phased_variants"] = qc["total_variants"] = qc["phase_blocks"] = qc["block_n50"] = "NA"
    qc["error"] += f"hiphase:{e};"

fields = ["sample_id","status","error",
          "mapping_rate","total_reads","mapped_reads",
          "total_sv","phased_variants","total_variants","phase_blocks","block_n50"]
with open("${meta.id}.qc.csv","w",newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerow({k: qc.get(k,"NA") for k in fields})
PYEOF
    """
}

process QC_MERGE {
    container params.main_sif
    publishDir params.outdir, mode: 'copy'

    input:
    path csv_files

    output:
    path "qc_summary.csv"

    script:
    """
    python3 << 'PYEOF'
import csv, glob
files = sorted(glob.glob("*.qc.csv"))
if not files:
    open("qc_summary.csv","w").write("no data\\n")
else:
    rows, header = [], None
    for f in files:
        with open(f) as fh:
            reader = csv.DictReader(fh)
            if header is None:
                header = reader.fieldnames
            rows.extend(reader)
    with open("qc_summary.csv","w",newline="") as out:
        w = csv.DictWriter(out, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    print(f"QC summary: {len(rows)} samples")
PYEOF
    """
}

// ============================================================
// MEMORY_MONITOR
// ============================================================
process MEMORY_MONITOR {
    // QC_MERGE 완료(= 파이프라인 전체 종료) 후 실행됨.
    // 실행 중 메모리 기록은 run.sh의 외부 백그라운드 루프(system_memory.tsv)가 담당.
    // 이 프로세스는 각 샘플별 완료 시점 스냅샷 + 실행된 프로세스 목록을 기록한다.
    tag "${meta.id}"
    container params.main_sif
    publishDir { "${params.outdir}/${meta.id}/logs" }, mode: 'copy'

    input:
    val meta
    val done_signal   // QC_MERGE 출력 파일명 (완료 신호)

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
            def meta       = [id: row.sample_id]
            def bam        = file(row.bam)
            def start_step = row.start_step ?: 'pbmm2'
            tuple(meta, bam, start_step)
        }

    ch_input.branch {
        jasmine:     it[2] == 'jasmine'
        pbmm2:       it[2] == 'pbmm2'
        deepvariant: it[2] == 'deepvariant'
        hiphase:     it[2] == 'hiphase'
    }.set { ch_branched }

    // 메모리 모니터는 QC_MERGE 완료 후 실행 (파이프라인 전체 종료 시점 포착)

    // STEP 01: JASMINE
    ch_for_jasmine = ch_branched.jasmine.map { meta, bam, step -> tuple(meta, bam) }
    jasmine_out    = JASMINE(ch_for_jasmine)
    ch_methyl_for_pbmm2 = jasmine_out.methyl_bam.map { meta, bam, bai -> tuple(meta, bam) }

    // STEP 02: PBMM2
    ch_for_pbmm2 = ch_branched.pbmm2
        .map { meta, bam, step -> tuple(meta, bam) }
        .mix(ch_methyl_for_pbmm2)
    pbmm2_out = PBMM2(ch_for_pbmm2)

    // ch_aligned: PBMM2 출력 → 여러 단계 재사용
    pbmm2_out.aligned_bam
        .multiMap { meta, bam, bai ->
            for_dv:       tuple(meta, bam, bai)
            for_sawfish:  tuple(meta, bam, bai)
            for_trgt:     tuple(meta, bam, bai)
            for_hiphase:  tuple(meta, bam, bai)
        }
        .set { ch_aligned }

    // start_step=deepvariant 일 때 BAM 공급
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

    // STEP 04: SAWFISH — DV 완료 후 실행 (DV VCF join으로 순서 강제)
    ch_for_sawfish = ch_dv_entry.for_sawfish
        .mix(ch_aligned.for_sawfish)
        .join(dv_out.snv_vcf, by: 0)
        .map { meta, bam, bai, snv_vcf, snv_tbi -> tuple(meta, bam, bai) }
    sawfish_out = SAWFISH(ch_for_sawfish)

    // STEP 05: TRGT — SAWFISH 완료 후 실행 (SAWFISH stats join으로 순서 강제)
    ch_for_trgt = ch_dv_entry.for_trgt
        .mix(ch_aligned.for_trgt)
        .join(sawfish_out.stats, by: 0)
        .map { meta, bam, bai, stats -> tuple(meta, bam, bai) }
    trgt_out = TRGT(ch_for_trgt)

    // STEP 06: HIPHASE — DV + SAWFISH + TRGT 모두 완료 후 실행
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

    hiphase_out = HIPHASE(ch_for_hiphase.mix(ch_hiphase_extra))

    // STEP 07: CPG_TOOLS — haplotagged BAM → CpG methylation (가장 먼저)
    cpg_out = CPG_TOOLS(hiphase_out.haplotagged_bam)

    // STEP 08: STELLER — CPG_TOOLS 완료 후 실행 (cpg_bed join으로 순서 강제)
    ch_for_steller = hiphase_out.haplotagged_bam
        .join(cpg_out.cpg_bed, by: 0)
        .map { meta, bam, bai, cpg -> tuple(meta, bam, bai) }
    steller_out = STELLER(ch_for_steller)

    // STEP 09: ANNOTSV_TE — STELLER 완료 후 실행
    annotsv_te_out = ANNOTSV_TE(steller_out.te_vcf)

    // STEP 10: SNP_ANNOTATION (VEP) — ANNOTSV_TE 완료 후 실행 (te_annot join으로 순서 강제)
    ch_for_snp = hiphase_out.snv_phased
        .join(annotsv_te_out.te_annot, by: 0)
        .map { meta, vcf, tbi, te -> tuple(meta, vcf, tbi) }
    SNP_ANNOTATION(ch_for_snp)

    // STEP 11: ANNOTSV_SV — SNP_ANNOTATION 과 병렬 (ANNOTSV_TE 완료 후)
    ch_for_annotsv_sv = hiphase_out.sv_phased
        .join(annotsv_te_out.te_annot, by: 0)
        .map { meta, vcf, tbi, te -> tuple(meta, vcf, tbi) }
    ANNOTSV_SV(ch_for_annotsv_sv)

    // STEP 12: ANNOTSV_TR — SNP_ANNOTATION 과 병렬 (ANNOTSV_TE 완료 후)
    ch_for_annotsv_tr = hiphase_out.tr_phased
        .join(annotsv_te_out.te_annot, by: 0)
        .map { meta, vcf, tbi, te -> tuple(meta, vcf, tbi) }
    ANNOTSV_TR(ch_for_annotsv_tr)

    // QC
    ch_qc_input = pbmm2_out.flagstat
        .join(sawfish_out.stats,      by: 0)
        .join(hiphase_out.summary,    by: 0)
        .map { meta, flagstat, sv_stats, hiphase_sum ->
            tuple(meta, flagstat, sv_stats, hiphase_sum)
        }
    qc_per_sample = QC_COLLECT(ch_qc_input)
    qc_merge_out  = QC_MERGE(qc_per_sample.qc_csv.collect())

    // MEMORY_MONITOR: QC_MERGE 완료(= 파이프라인 전체 종료) 후 실행
    ch_meta      = ch_input.map { meta, bam, step -> meta }
    ch_done_sig  = qc_merge_out.map { f -> f.name }
    MEMORY_MONITOR(ch_meta.combine(ch_done_sig).map { meta, sig -> meta },
                   ch_meta.combine(ch_done_sig).map { meta, sig -> sig }.first())
}
