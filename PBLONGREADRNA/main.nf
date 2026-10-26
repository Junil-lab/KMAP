#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PBMM2; ISOQUANT } from './modules/core'
include { FLAIR_TRANSCRIPTOME; FLAIR_COMBINE } from './modules/flair'
include { ISOCALL_PREP; ISOCALL_PROFILE; ISOCALL_MERGE as ISOCALL_SAMPLE_MERGE; ISOCALL_MERGE as ISOCALL_COHORT_MERGE; ISOCALL_CALL as ISOCALL_SAMPLE_CALL; ISOCALL_CALL as ISOCALL_COHORT_CALL } from './modules/isocall'
include { SQANTI3; QC_COLLECT; QC_MERGE } from './modules/qc'

def annotate(meta, method, scope = 'sample', samples = null) {
    meta + [method: method, scope: scope, prefix: "${scope}.${meta.id}.${method}".toString(),
            samples: samples ?: [meta.id]]
}
def qcKey(meta) { "${meta.scope}|${meta.id}|${meta.method}".toString() }

workflow {
    if (!params.samplesheet) error '--samplesheet is required'
    if (!(params.threads.toString() ==~ /[1-9][0-9]*/)) error '--threads must be a positive integer'
    def methods = params.methods.toString().split(',').collect { it.trim().toLowerCase() }.unique()
    if (!methods || methods.any { !(it in ['isoquant', 'isocall', 'flair']) })
        error '--methods must be a comma-separated subset of isoquant,isocall,flair'
    def sheet = file(params.samplesheet, checkIfExists: true)
    def genome = file(params.genome_fa, checkIfExists: true)
    def genome_index = file(params.genome_fa + '.fai', checkIfExists: true)
    def annotation = file(params.gencode_gtf, checkIfExists: true)
    file(params.main_sif, checkIfExists: true)
    file(params.sqanti3_sif, checkIfExists: true)
    if (methods.any { it in ['isocall', 'flair'] }) file(params.discovery_sif, checkIfExists: true)
    def optional_files = [params.polya_motif_list, params.polya_peak_bed, params.cage_peak_bed]
        .findAll { it }.collect { file(it, checkIfExists: true) }
    if (optional_files.collect { it.name }.unique().size() != optional_files.size())
        error 'Optional annotation filenames must be unique'

    def validation = ['python3', "${projectDir}/bin/validate_samples.py", sheet.toString(), '--json'].execute()
    def validated_json = new StringBuffer()
    def validation_errors = new StringBuffer()
    validation.waitForProcessOutput(validated_json, validation_errors)
    if (validation.exitValue() != 0) error validation_errors.toString().trim()
    def rows = new groovy.json.JsonSlurper().parseText(validated_json.toString())
    def samples = rows.collect { row ->
        [[id: row.sample_id.toString(), tissue: row.tissue.toString()], file(row.flnc_bam, checkIfExists: true)]
    }
    def sample_ids = samples.collect { it[0].id }.sort()
    def mode = params.analysis_mode.toString().toLowerCase()
    if (!(mode in ['auto','single','multi','both']))
        error '--analysis_mode must be auto, single, multi, or both'
    if (mode in ['multi','both'] && samples.size() < 2)
        error '--analysis_mode multi/both requires at least two samples'
    if (mode in ['multi','both'] && !methods.any { it in ['isocall','flair'] })
        error '--analysis_mode multi/both requires isocall or flair in --methods'
    def do_sample = mode != 'multi'
    def do_cohort = mode in ['multi','both'] || (mode == 'auto' && params.run_cohort && samples.size() > 1)
    log.info "Analysis mode: ${mode}; sample-level FLAIR/isocall QC: ${do_sample}; cohort: ${do_cohort}"

    def cohort_base = [id: 'cohort', tissue: 'COHORT']
    def expected = samples.collectMany { row -> methods.findAll { it == 'isoquant' || do_sample }.collect { method -> qcKey(annotate(row[0], method)) } }
    if (do_cohort) expected.addAll(methods.findAll { it in ['isocall','flair'] }.collect { qcKey(annotate(cohort_base,it,'cohort',sample_ids)) })
    ch_input = Channel.fromList(samples)

    def iq_gtfs = Channel.empty()
    def iq_flags = Channel.empty()
    PBMM2(ch_input, genome, genome_index)
    if ('isoquant' in methods) {
        ISOQUANT(PBMM2.out.aligned_bam, genome, genome_index, annotation)
        iq_gtfs = ISOQUANT.out.transcript_gtf.map { meta, gtf -> [annotate(meta,'isoquant'),gtf] }
        iq_flags = PBMM2.out.flagstat.map { meta, flag ->
            def m = annotate(meta,'isoquant'); [qcKey(m),m,[flag]]
        }
    }

    def flair_gtfs = Channel.empty()
    def flair_cohort_gtfs = Channel.empty()
    def flair_flags = Channel.empty()
    def flair_cohort_flags = Channel.empty()
    def isocall_gtfs = Channel.empty()
    def isocall_cohort_gtfs = Channel.empty()
    def isocall_flags = Channel.empty()
    def isocall_cohort_flags = Channel.empty()
    if (methods.any { it in ['isocall','flair'] }) {
        if ('flair' in methods) {
            FLAIR_TRANSCRIPTOME(PBMM2.out.aligned_bam, genome, genome_index, annotation)
            if (do_sample) {
                flair_gtfs = FLAIR_TRANSCRIPTOME.out.models.map { meta, gtf, bed, fa, map -> [annotate(meta,'flair'),gtf] }
                flair_flags = PBMM2.out.flagstat.map { meta, flag ->
                    def m = annotate(meta,'flair'); [qcKey(m),m,[flag]]
                }
            }
            if (do_cohort) {
                ch_flair_cohort = FLAIR_TRANSCRIPTOME.out.models.collect(flat:false).map { models ->
                    def ordered = models.sort { a,b -> a[0].id <=> b[0].id }
                    [annotate(cohort_base,'flair','cohort',sample_ids), ordered.collectMany { [it[2],it[3],it[4]] }]
                }
                FLAIR_COMBINE(ch_flair_cohort)
                flair_cohort_gtfs = FLAIR_COMBINE.out.models.map { meta,gtf,bed,fa -> [meta,gtf] }
                flair_cohort_flags = PBMM2.out.flagstat.collect(flat:false).map { entries ->
                    def m = annotate(cohort_base,'flair','cohort',sample_ids); [qcKey(m),m,entries.collect { it[1] }]
                }

            }
        }
        if ('isocall' in methods) {
            ISOCALL_PREP(annotation)
            ISOCALL_PROFILE(PBMM2.out.aligned_bam)
            if (do_sample) {
                ch_sample_profiles = ISOCALL_PROFILE.out.profile.map { meta, profile -> [annotate(meta,'isocall'),[profile]] }
                ISOCALL_SAMPLE_MERGE(ch_sample_profiles)
                ISOCALL_SAMPLE_CALL(ISOCALL_SAMPLE_MERGE.out.merged,ISOCALL_PREP.out.known,genome,genome_index)
                isocall_gtfs = ISOCALL_SAMPLE_CALL.out.gtf
                isocall_flags = PBMM2.out.flagstat.map { meta,flag ->
                    def m = annotate(meta,'isocall'); [qcKey(m),m,[flag]]
                }
            }
            if (do_cohort) {
                ch_cohort_profiles = ISOCALL_PROFILE.out.profile.collect(flat:false).map { entries ->
                    [annotate(cohort_base,'isocall','cohort',sample_ids),entries.sort { a,b -> a[0].id <=> b[0].id }.collect { it[1] }]
                }
                ISOCALL_COHORT_MERGE(ch_cohort_profiles)
                ISOCALL_COHORT_CALL(ISOCALL_COHORT_MERGE.out.merged,ISOCALL_PREP.out.known,genome,genome_index)
                isocall_cohort_gtfs = ISOCALL_COHORT_CALL.out.gtf
                isocall_cohort_flags = PBMM2.out.flagstat.collect(flat:false).map { entries ->
                    def m = annotate(cohort_base,'isocall','cohort',sample_ids); [qcKey(m),m,entries.collect { it[1] }]
                }
            }
        }
    }
    all_gtfs = iq_gtfs.mix(flair_gtfs,flair_cohort_gtfs,isocall_gtfs,isocall_cohort_gtfs)
    all_flags = iq_flags.mix(flair_flags,flair_cohort_flags,isocall_flags,isocall_cohort_flags)
    SQANTI3(all_gtfs,genome,genome_index,annotation,Channel.value(optional_files))
    qc_inputs = all_flags.join(SQANTI3.out.classification.map { meta,cls -> [qcKey(meta),cls] },
        by:0,failOnDuplicate:true,failOnMismatch:true).map { key,meta,flags,cls -> [meta,flags,cls] }
    QC_COLLECT(qc_inputs)
    QC_MERGE(QC_COLLECT.out.collect(),Channel.value(expected))
}
