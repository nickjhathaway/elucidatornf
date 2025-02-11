#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE } from '../../../modules/local/gen_target_info_from_genomes'
include { EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/extractor_by_kmer_matching'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'



/**
process EXTRACTOR_BY_KMER_MATCHING {
    label 'process_single'


    input:
    path fastq_fnp
    val sample_name
    path primers_fnp
    path kmers_sets
    path length_cut_offs_per_target
    val min_len_cut_off
*/

workflow NANOPORE_AMPLICON_CLUSTERING {
    take:
    input_fastq_dir //a directory with all fastqs to be analyzed
    primers_fnp  // a file with primer info, 3 columsn, target,fwd_primer,rev_primer
    genome_dir  // a directory with genomes
    gff_dir // a directory with annotation gff files
    primers_errors_allowed // errors to allow in primers
    meta_fnp //a sample meta file

    main:

    def primer_info_dir = file("${params.outdir}/primerInfo")
    primer_info_dir.mkdirs()
    GEN_TARGET_INFO_FROM_GENOMES_NANOPORE(primers_fnp, primer_info_dir, genome_dir, gff_dir, primers_errors_allowed, params.resources.max_cpus)

    def fastq_input_ch = Channel.fromPath(file("${input_fastq_dir}/*.fastq.gz"))

    // def kmers_sets_ch = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info.map{ for_seekdeep_dir ->
    //     file("${for_seekdeep_dir}/allKmers.tab.txt.gz")
    // }
    // def len_cut_offs_ch = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info.map{ for_seekdeep_dir ->
    //     file("${for_seekdeep_dir}/lenCutOffs.txt")
    // }

    def input_to_extractor = fastq_input_ch.combine(GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info)
            .map{fastq_fnp, info_dir ->
                tuple(
                    fastq_fnp,
                    file(fastq_fnp).name.replaceAll(".fastq.gz"),
                    primers_fnp,
                    file("${info_dir}/allKmers.tab.txt.gz"),
                    file("${info_dir}/lenCutOffs.txt",
                    params.nanopore_clustering_min_len)
                )
            }

    EXTRACTOR_BY_KMER_MATCHING(input_to_extractor)

}  //PATHWEAVER_EXTRACT_REGIONS_FULL

