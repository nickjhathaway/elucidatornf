#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMPLICON_READ_QC } from '../../../modules/local/amplicon_read_qc/main.nf'
include { AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR } from "../../../modules/local/amplicon_cluster_auto_SeekDeep_flag_generator/main.nf"
include { GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA } from '../../../modules/local/gen_target_info_from_genomes/main.nf'


workflow ILLUMINA_READ_QC {
    take:
    input_fastq_dir_raw //a directory with all fastqs to be analyzed
    primers_fnp_raw  // a file with primer info, 3 columsn, target,fwd_primer,rev_primer
    genome_dir_raw  // a directory with genomes
    gff_dir_raw // a directory with annotation gff files, only used when the target info is generated
    target_info_dir_raw // the output directory of an already run genTargetInfoFromGenomes, leave blank to generate it
    primers_errors_allowed // errors to allow in primers when generating the target info
    paired_end_length //the paired end length of the input, only used when the target info is generated
    output_dir_raw // the output directory
    main:

    def input_fastq_dir = file(input_fastq_dir_raw)
    def primers_fnp = file(primers_fnp_raw)
    def genome_dir = file(genome_dir_raw)
    def output_dir = file(output_dir_raw)
    def errors = []

    if (!input_fastq_dir.isDirectory()) {
        errors << "Directory does not exist for --input_fastq_dir: '${input_fastq_dir_raw}'"
    }
    if (!primers_fnp.isFile()) {
        errors << "File does not exist for --primers_fnp: '${primers_fnp_raw}'"
    }
    if (!genome_dir.isDirectory()) {
        errors << "Directory does not exist for --genome_dir: '${genome_dir_raw}'"
    }

    // either an already run genTargetInfoFromGenomes or what's needed to generate it
    def target_info_dir = target_info_dir_raw ? file(target_info_dir_raw) : null
    if (target_info_dir != null) {
        ["forSeekDeep/overlapStatuses.txt", "forSeekDeep/refSeqs", "locationsByGenome"].each { expected ->
            if (!file("${target_info_dir}/${expected}").exists()) {
                errors << "--illumina_read_qc_target_info_dir should be the output directory of genTargetInfoFromGenomes but it has no ${expected}: '${target_info_dir_raw}'"
            }
        }
    } else {
        if (!gff_dir_raw || !file(gff_dir_raw).isDirectory()) {
            errors << "--gff_dir has to be a directory unless --illumina_read_qc_target_info_dir is given: '${gff_dir_raw}'"
        }
        if (!paired_end_length) {
            errors << "--illumina_clustering_paired_end_length must be supplied unless --illumina_read_qc_target_info_dir is given"
        }
    }

    def custom_ref_seqs_dir = params.illumina_clustering_ref_seqs_dir ? file(params.illumina_clustering_ref_seqs_dir) : null
    if (custom_ref_seqs_dir != null && !custom_ref_seqs_dir.isDirectory()) {
        errors << "Directory does not exist for --illumina_clustering_ref_seqs_dir: '${custom_ref_seqs_dir}'"
    }

    def sample_key = SampleKey.load(SampleKey.keyFnp(params))
    errors += ReadQc.validate(params, genome_dir, sample_key)

    if (!errors.isEmpty()) {
        def joined = errors.collect { it -> "  - ${it}" }.join("\n")
        error """
Input validation failed in subworkflow ILLUMINA_READ_QC:

${joined}

Please correct the above issues and re-run.
"""
    }

    def read_qc_dir = "${output_dir}/readQC"

    def for_seek_deep_info_ch
    def locations_by_genome_ch
    if (target_info_dir != null) {
        for_seek_deep_info_ch = channel.value(file("${target_info_dir}/forSeekDeep"))
        locations_by_genome_ch = channel.value(file("${target_info_dir}/locationsByGenome"))
    } else {
        // the output directory as a string, as a file its modified time ends up in the hash and publishing changes that
        GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA(primers_fnp, "${output_dir}/primerInfo", genome_dir, file(gff_dir_raw), paired_end_length, 0, primers_errors_allowed, params.max_cpus)
        for_seek_deep_info_ch = GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA.out.for_seek_deep_info
        locations_by_genome_ch = GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA.out.locations_by_genome
    }
    def ref_seqs_dir_ch = custom_ref_seqs_dir != null ?
        channel.value(custom_ref_seqs_dir) :
        for_seek_deep_info_ch.map { seek_deep_info_dir -> file("${seek_deep_info_dir}/refSeqs") }

    // the same flags the extraction would get, so the primers are looked for in the same way
    AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR(input_fastq_dir, primers_fnp, "illumina", params.max_cpus, read_qc_dir, "", false)

    def fastq_input_ch = channel.fromFilePairs(
        ["${input_fastq_dir}/*_{1,2}{_001,}.fastq{.gz,}",
        "${input_fastq_dir}/*_{R1,R2}{_001,}.fastq{.gz,}"]
    ).ifEmpty { error "No paired fastqs found in --input_fastq_dir: '${input_fastq_dir}'" }

    def host_genome = params.illumina_read_qc_host_genome_fnp ? file(params.illumina_read_qc_host_genome_fnp) : null

    AMPLICON_READ_QC(
        fastq_input_ch.toSortedList { a, b -> a[0] <=> b[0] }.map { fastq_pairs -> ReadQc.samplesTuple(fastq_pairs) },
        primers_fnp,
        for_seek_deep_info_ch,
        locations_by_genome_ch,
        ref_seqs_dir_ch,
        genome_dir,
        params.illumina_read_qc_primary_genome,
        host_genome != null ? host_genome.getParent() : file(params.empty_file_fnp),
        host_genome != null ? host_genome.getName() : "",
        AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR.out.out_seekdeep_extractor_flags,
        params.max_cpus,
        read_qc_dir
    )
}
