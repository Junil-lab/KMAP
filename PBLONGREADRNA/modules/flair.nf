def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }


process FLAIR_TRANSCRIPTOME {
    tag "${meta.id}"
    label 'high_mem'
    container params.discovery_sif
    publishDir { "${params.outdir}/samples/${meta.id}/flair/discovery" }, mode: 'copy'
    input:
    tuple val(meta), path(aligned_bam), path(aligned_bai)
    path genome_fa
    path genome_fai
    path gencode_gtf
    output:
    tuple val(meta), path("${meta.id}.isoforms.gtf"), path("${meta.id}.isoforms.bed"),
        path("${meta.id}.isoforms.fa"), path("${meta.id}.isoform.read.map.txt"), emit: models
    path "${meta.id}.isoform.counts.txt"
    path "${meta.id}.flair_transcriptome.log"
    script:
    """
    export HOME="\$PWD/.home"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$HOME" "\$TMPDIR"
    flair transcriptome -b ${quote(aligned_bam)} -g ${quote(genome_fa)} \
        -f ${quote(gencode_gtf)} -o ${meta.id} -t ${task.cpus} \
        -s 3 --stringent --check_splice 2>&1 | tee ${meta.id}.flair_transcriptome.log
    """
}

process FLAIR_COMBINE {
    tag 'cohort'
    label 'high_mem'
    container params.discovery_sif
    publishDir "${params.outdir}/cohort/flair/discovery", mode: 'copy'
    input:
    tuple val(meta), path(model_files, stageAs: 'models/*')
    output:
    tuple val(meta), path("${meta.prefix}.gtf"), path("${meta.prefix}.bed"), path("${meta.prefix}.fa"), emit: models
    path "${meta.prefix}.counts.tsv"
    path "combine_manifest.tsv"
    path "flair_combine.log"
    script:
    def manifest = meta.samples.collect { id -> "${id}\tisoform\tmodels/${id}.isoforms.bed\tmodels/${id}.isoforms.fa\tmodels/${id}.isoform.read.map.txt" }.join('\n')
    """
    cat > combine_manifest.tsv <<'MANIFEST'
${manifest}
MANIFEST
    flair combine -m combine_manifest.tsv -o ${meta.prefix} --convert_gtf -s \
        2>&1 | tee flair_combine.log
    """
}
