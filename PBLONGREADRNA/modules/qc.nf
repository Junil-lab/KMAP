def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process SQANTI3 {
    tag "${meta.prefix}"
    label 'high_mem'
    container params.sqanti3_sif
    publishDir { "${params.outdir}/${meta.scope == 'cohort' ? 'cohort' : 'samples/' + meta.id}/${meta.method}/sqanti3" }, mode: 'copy'
    input:
    tuple val(meta), path(transcript_gtf)
    path genome_fa
    path genome_fai
    path gencode_gtf
    path annotations, stageAs: 'annotations/*'
    output:
    tuple val(meta), path("${meta.prefix}_classification.txt"), emit: classification
    path "${meta.prefix}_junctions.txt"
    path "${meta.prefix}_corrected.gtf"
    path "${meta.prefix}_corrected.fasta"
    path "${meta.prefix}_SQANTI3_report.html"
    path "${meta.prefix}_corrected.faa", optional: true
    path "${meta.prefix}.sqanti3.log"
    script:
    def opts = []
    if (params.polya_motif_list) opts << "--polyA_motif_list ${quote('annotations/' + file(params.polya_motif_list).name)}"
    if (params.polya_peak_bed) opts << "--polyA_peak ${quote('annotations/' + file(params.polya_peak_bed).name)}"
    if (params.cage_peak_bed) opts << "--CAGE_peak ${quote('annotations/' + file(params.cage_peak_bed).name)}"
    """
    export HOME="\$PWD/.home"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$HOME" "\$XDG_CACHE_HOME"
    CUDA_VISIBLE_DEVICES="" bash ${quote(projectDir.toString() + "/bin/run_sqanti3.sh")} ${quote(params.sqanti3_qc)} \
        --isoforms ${quote(transcript_gtf)} --refGTF ${quote(gencode_gtf)} \
        --refFasta ${quote(genome_fa)} --output ${meta.prefix} --dir . \
        --report html --cpus ${task.cpus} --include_ORF ${opts.join(' ')} \
        2>&1 | tee ${meta.prefix}.sqanti3.log
    """
}


process QC_COLLECT {
    tag "${meta.prefix}"
    label 'low_mem'
    container params.main_sif
    publishDir "${params.outdir}/qc_summary", mode: 'copy'
    input:
    tuple val(meta), path(flagstats, stageAs: 'flagstats/*'), path(classification)
    output:
    path "${meta.prefix}.qc.csv"
    script:
    def flag_args = [flagstats].flatten().collect { quote(it) }.join(' ')
    """
    python3 ${quote(projectDir.toString() + '/bin/qc_summary.py')} collect \
        --sample ${quote(meta.id)} --scope ${quote(meta.scope)} --method ${quote(meta.method)} \
        --tissue ${quote(meta.tissue)} --flagstat ${flag_args} \
        --classification ${quote(classification)} --output ${meta.prefix}.qc.csv
    """
}

process QC_MERGE {
    label 'low_mem'
    container params.main_sif
    publishDir "${params.outdir}/qc_summary", mode: 'copy'
    input:
    path qc_csvs, stageAs: 'qc_csvs/*'
    val expected_keys
    output:
    path 'qc_summary.csv'
    script:
    """
    python3 ${quote(projectDir.toString() + '/bin/qc_summary.py')} merge --input-dir qc_csvs \
        --expected ${quote(expected_keys.join(','))} --output qc_summary.csv
    """
}
