def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process PBMM2 {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/samples/${meta.id}/alignment/pbmm2" }, mode: 'copy'
    input:
    tuple val(meta), path(flnc_bam)
    path genome_fa
    path genome_fai
    output:
    tuple val(meta), path("${meta.id}.aligned.bam"), path("${meta.id}.aligned.bam.bai"), emit: aligned_bam
    path "${meta.id}.pbmm2.log"
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat
    script:
    """
    pbmm2 align --preset ISOSEQ --sort -j ${task.cpus} \
        ${quote(genome_fa)} ${quote(flnc_bam)} ${meta.id}.aligned.bam \
        2>&1 | tee ${meta.id}.pbmm2.log
    samtools index -@ ${task.cpus} ${meta.id}.aligned.bam
    samtools flagstat ${meta.id}.aligned.bam > ${meta.id}.flagstat.txt
    """
}

// start_step=discovery: use an existing aligned BAM instead of running pbmm2
process ALIGNED_INPUT {
    tag "${meta.id}"
    label 'low_mem'
    container params.main_sif
    publishDir { "${params.outdir}/samples/${meta.id}/alignment/input" }, mode: 'copy', pattern: '*.flagstat.txt'
    input:
    tuple val(meta), path(aligned_bam), path(aligned_bai)
    output:
    tuple val(meta), path(aligned_bam), path(aligned_bai), emit: aligned_bam
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat
    script:
    """
    samtools quickcheck ${quote(aligned_bam)}
    if ! samtools view -H ${quote(aligned_bam)} | grep -q '^@HD.*SO:coordinate'; then
        echo "${meta.id}: aligned_bam must be coordinate-sorted (pbmm2 align --sort or samtools sort)" >&2
        exit 1
    fi
    samtools flagstat -@ ${task.cpus} ${quote(aligned_bam)} > ${meta.id}.flagstat.txt
    """
}

process ISOQUANT {
    tag "${meta.id}"
    label 'high_mem'
    container params.main_sif
    publishDir { "${params.outdir}/samples/${meta.id}/isoquant/discovery" }, mode: 'copy'
    input:
    tuple val(meta), path(aligned_bam), path(aligned_bai)
    path genome_fa
    path genome_fai
    path gencode_gtf
    output:
    tuple val(meta), path("${meta.id}/${meta.id}.transcript_models.gtf"), emit: transcript_gtf
    // Counts corresponding to transcript_models.gtf (known + novel models).
    tuple val(meta), path("${meta.id}/${meta.id}.discovered_transcript_counts.tsv"),
        path("${meta.id}/${meta.id}.discovered_transcript_tpm.tsv"), emit: counts
    // Preserve reference-only quantification separately.
    path "${meta.id}/${meta.id}.transcript_counts.tsv"
    path "${meta.id}/${meta.id}.transcript_tpm.tsv"
    path "${meta.id}/${meta.id}.gene_counts.tsv"
    path "${meta.id}/${meta.id}.gene_tpm.tsv"
    path "${meta.id}/${meta.id}.extended_annotation.gtf"
    path "isoquant.log"
    script:
    """
    export HOME="\$PWD/.home"
    mkdir -p "\$HOME"
    export NUMBA_CACHE_DIR="\$PWD/.numba_cache"
    mkdir -p "\$NUMBA_CACHE_DIR"
    isoquant --reference ${quote(genome_fa)} --genedb ${quote(gencode_gtf)} \
        --complete_genedb --bam ${quote(aligned_bam)} --data_type pacbio_ccs \
        --fl_data --model_construction_strategy fl_pacbio \
        --output . --prefix ${meta.id} --labels ${meta.id} \
        --threads ${task.cpus} --sqanti_output --check_canonical --count_exons
    """
}
