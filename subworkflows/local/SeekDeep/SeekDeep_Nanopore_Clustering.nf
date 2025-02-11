#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE } from '../../../modules/local/gen_target_info_from_genomes'
include { EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/extractor_by_kmer_matching'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { CONCATENATE_EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/concatenate_extractor_by_kmer_matching/main.nf'
include { AMPLICON_CLUSTER_BY_KMER_SIMILARITY } from '../../../modules/local/amplicon_cluster_by_kmer_similarity/main.nf'
include { AMPLICON_POPULATION_CLUSTERING } from '../../../modules/local/amplicon_population_cluster/main.nf'
include { CONCATENATE_AMPLICON_POPULATION_CLUSTERING } from '../../../modules/local/concatenate_amplicon_population_clustering/main.nf'
include { GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA } from '../../../modules/local/gen_target_info_from_genomes/main.nf'


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
    def extraction_reports_dir = file("${params.outdir}/extractionReports")
    extraction_reports_dir.mkdirs()
    def final_results_dir = file("${params.outdir}/finalResults")
    final_results_dir.mkdirs()
    GEN_TARGET_INFO_FROM_GENOMES_NANOPORE(primers_fnp, primer_info_dir, genome_dir, gff_dir, primers_errors_allowed, params.resources.max_cpus)

    def fastq_input_ch = Channel.fromPath(file("${input_fastq_dir}/*.fastq.gz"))

    def input_to_extractor = fastq_input_ch.combine(GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info)
            .map{fastq_fnp, info_dir ->
                tuple(
                    fastq_fnp,
                    file(fastq_fnp).name.replaceAll(".fastq.gz", ""),
                    primers_fnp,
                    file("${info_dir}/allKmers.tab.txt.gz"),
                    file("${info_dir}/lenCutOffs.txt"),
                    params.nanopore_clustering_min_len
                )
            }

    //extract by kmer sets
    EXTRACTOR_BY_KMER_MATCHING(input_to_extractor)

    //concatenate extraction files
    CONCATENATE_EXTRACTOR_BY_KMER_MATCHING(
        EXTRACTOR_BY_KMER_MATCHING.out.extraction_profile_per_target | collect,
        EXTRACTOR_BY_KMER_MATCHING.out.extraction_stats | collect,
        extraction_reports_dir
    )

    //get the targets that have coverage
    // def targets = CONCATENATE_EXTRACTOR_BY_KMER_MATCHING.out.targets_with_passing_reads
    //     .splitText()
    //     .map{samp ->
    //         samp.trim() // Remove any whitespace
    //         }
    def targets = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.targets_with_extractins
        .splitText()
        .map{samp ->
            samp.trim() // Remove any whitespace
            }

    // create a list of input to cluster for the input that actually exists
    def input_to_clustering = targets.combine(EXTRACTOR_BY_KMER_MATCHING.out.sample_dir)
        .map{tar, samp, sampdir ->
                def extracted_tar_sample_fnp = file("${sampdir}/${tar}${samp}.fastq.gz")
                if(extracted_tar_sample_fnp.exists()){
                    tuple(extracted_tar_sample_fnp,  samp, tar, params.nanopore_clustering_ncpus)
                }
            }

    //cluster per target per sample
    AMPLICON_CLUSTER_BY_KMER_SIMILARITY(input_to_clustering)


    // Use Groovy to count the files in the directory matching the glob pattern
    // def fastq_count = file("${input_fastq_dir}/*.fastq.gz").size()
    // log.info "fastq_count is ${fastq_count}"
    // Get the count of fastq files
    // def fastq_count = fastq_files.size()
    def ref_seqs_dir_ch = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info.map{ seek_deep_info_dir ->
            file("${seek_deep_info_dir}/refSeqs")}
    def input_to_population_clustering = AMPLICON_CLUSTER_BY_KMER_SIMILARITY.out.output_results
        //.groupTuple(by : 1, size: fastq_input_ch.count())
        // .groupTuple(by : 1, size: fastq_count)
        .groupTuple(by : 1)
        .map { samples, target, target_fastqs ->
            tuple(
                target_fastqs, // List of fastq files for this target
                target,
                params.nanopore_clustering_ncpus,
                meta_fnp,
                ref_seqs_dir_ch,
                params.nanopore_clustering_min_sample_read_count
            )
        }
    // population clustering
    AMPLICON_POPULATION_CLUSTERING(input_to_population_clustering)

    //gather all clustered results
    CONCATENATE_AMPLICON_POPULATION_CLUSTERING(AMPLICON_POPULATION_CLUSTERING.out.clustering_results | collect, final_results_dir)

    if (params.do_variant_calling){
        def variant_call_dir = file("${final_results_dir}/variantCalls")
        def known_amino_acid_changes_fnp = params.empty_file_fnp
        genome_dir_top = genome_dir.getParent()
        if(file("${genome_dir_top}/info/drug_resistant_aaPositions.tsv").exists()){
            known_amino_acid_changes_fnp = file("${genome_dir_top}/info/drug_resistant_aaPositions.tsv")
        }
        def bed_ch = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.locations_by_genome.map{loc_dir ->
                file("${loc_dir}/${params.vc_primary_genome}_inner.bed")}
        VARIANT_CALL_ON_HAP_TABLE (
            bed_ch,
            CONCATENATE_AMPLICON_POPULATION_CLUSTERING.out.all_selected_clusters_info,
            genome_dir,
            params.vc_primary_genome,
            file("${gff_dir}/${params.vc_primary_genome}.gff"),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp,
            params.vc_getting_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
    }
}  //PATHWEAVER_EXTRACT_REGIONS_FULL

