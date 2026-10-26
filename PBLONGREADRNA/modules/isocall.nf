def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }


process ISOCALL_PREP {
    label 'low_mem'
    container params.discovery_sif
    publishDir "${params.outdir}/reference/isocall", mode: 'copy'
    input:
    path gencode_gtf
    output:
    path "known_isoforms.gz", emit: known
    script:
    """
    isocall prep-isoforms --gtf ${quote(gencode_gtf)} --output known_isoforms.gz
    """
}

process ISOCALL_PROFILE {
    tag "${meta.id}"
    label 'high_mem'
    container params.discovery_sif
    publishDir { "${params.outdir}/samples/${meta.id}/isocall/profile" }, mode: 'copy'
    input:
    tuple val(meta), path(aligned_bam), path(aligned_bai)
    output:
    tuple val(meta), path("${meta.id}.profile.gz"), emit: profile
    path "${meta.id}.profile.log"
    script:
    """
    isocall -v profile --reads ${quote(aligned_bam)} --sample ${quote(meta.id)} \
        --io-threads ${task.cpus} --output ${meta.id}.profile.gz \
        2>&1 | tee ${meta.id}.profile.log
    """
}

process ISOCALL_MERGE {
    tag "${meta.scope}:${meta.id}"
    label 'high_mem'
    container params.discovery_sif
    publishDir { "${params.outdir}/${meta.scope == 'cohort' ? 'cohort' : 'samples/' + meta.id}/isocall/discovery" }, mode: 'copy'
    input:
    tuple val(meta), path(profiles, stageAs: 'profiles/*')
    output:
    tuple val(meta), path("${meta.prefix}.merged.gz"), emit: merged
    script:
    def profile_args = [profiles].flatten().collect { quote(it) }.join(' ')
    """
    isocall -v merge --profiles ${profile_args} --output ${meta.prefix}.merged.gz
    """
}

process ISOCALL_CALL {
    tag "${meta.scope}:${meta.id}"
    label 'high_mem'
    container params.discovery_sif
    publishDir { "${params.outdir}/${meta.scope == 'cohort' ? 'cohort' : 'samples/' + meta.id}/isocall/discovery" }, mode: 'copy'
    input:
    tuple val(meta), path(merged_profile)
    path known_isoforms
    path genome_fa
    path genome_fai
    output:
    tuple val(meta), path("${meta.prefix}.isoforms.gtf"), emit: gtf
    path "${meta.prefix}.isoforms.gtf.gz"
    path "${meta.prefix}.count_matrix.txt"
    path "${meta.prefix}.closest_known.txt"
    path "${meta.prefix}.isocall.log"
    script:
    """
    isocall -v call --merged-profile ${quote(merged_profile)} \
        --known-isoforms ${quote(known_isoforms)} --reference ${quote(genome_fa)} \
        --threads ${task.cpus} --output-prefix ${meta.prefix} \
        2>&1 | tee ${meta.prefix}.isocall.log
    gzip -dc ${meta.prefix}.isoforms.gtf.gz > ${meta.prefix}.isoforms.gtf
    """
}
