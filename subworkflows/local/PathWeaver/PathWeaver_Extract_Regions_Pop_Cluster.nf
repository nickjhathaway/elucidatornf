#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXTRACT_REGION_ASSEMBLIES } from '../../../modules/local/extract_region_assemblies'
include { PATHWEAVER_POP_CLUSTERING } from '../../../modules/local/pathweaver_pop_clustering'
include { GET_SUB_SEGMENTS_FROM_FASTA } from '../../../modules/local/get_sub_segments_from_PW_pop_res/main.nf'
include { CONCATENATE_SUB_SEGMENT_LOCS } from '../../../modules/local/concatenate_sub_segment_locs/main.nf'
include { CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS } from '../../../modules/local/concatenate_extract_region_assemblies_results/main.nf'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { VARIANT_CALL_ON_HAP_TABLE as VARIANT_CALL_ON_HAP_TABLE_2 } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { CREATE_TWOBIT_FILE_FOR_GENOME } from '../../../modules/local/create_twobit_file_for_genome/main.nf'
include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE } from '../../../modules/local/gen_target_info_from_genomes'
include { GEN_TARGET_INFO_FROM_GENOME_NANOPORE as GEN_TARGET_INFO_FROM_GENOME } from '../../../modules/local/gen_target_info_from_genomes/main.nf'
include { DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE } from '../../../modules/local/determine_genomic_location_from_seqs_table/main.nf'
include { GET_INTERSECTING_GENE_INFO_FOR_REGIONS } from '../../../modules/local/get_intersecting_gene_info_for_regions/main.nf'
include { GET_GENE_RECORDS_FOR_GENE_IDS } from '../../../modules/local/get_gene_records_for_gene_ids/main.nf'
include { REMOVE_TANDEM_REPEATS_FROM_REGIONS } from '../../../modules/local/remove_tandem_repeats_from_regions/main.nf'




workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp_ch  // Path to the BED file
    genome_fnp // Path to the genome file
    // genome_twobit_fnp // Path to the genome twobit file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:


    // Create output directories
    pop_clus_dir = file("${results_dir}/reports/populationClustering/")
    reports_dir = file("${results_dir}/reports/rawResults/")
    pop_clus_dir.mkdirs()
    reports_dir.mkdirs()
    // Load samples from file and create a channel
    if ("EMPTY_FILE.txt" == file("${samples_file}").name){
        samples = Channel.fromPath("${bams_dir}/*${params.bams_file_ending}")
            .map{ samp_file ->
                file(samp_file).name.replaceAll("${params.bams_file_ending}", "")
        }
    } else {
        samples = Channel
            .fromPath("${samples_file}")
            .splitText()
            .map{samp ->
                samp.trim() // Remove any whitespace
                }
    }


    // Construct a channel of tuples (BAM, BAI, bed, genome_dir, primary_genome, results_dir)
    input_ch = bed_fnp_ch.combine(samples).map {bed_fnp, samp ->
        tuple(
            samp,
            file("${bams_dir}/${samp}${params.bams_file_ending}"),
            file("${bams_dir}/${samp}${params.bams_file_ending}.bai"),
            file("${bed_fnp}"),
            file("${genome_fnp}").getParent(),
            file("${genome_fnp}").getBaseName(),
            file("${results_dir}").baseName
        )
    }
    // run PathWeaver on each
    EXTRACT_REGION_ASSEMBLIES(input_ch)

    // concatenate the results files
    CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(EXTRACT_REGION_ASSEMBLIES.out | collect, reports_dir.toString())

    // run population clustering
    PATHWEAVER_POP_CLUSTERING(file("${results_dir}").baseName, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, pop_clus_dir.toString())

    emit:
    samples = samples
    assemblies = EXTRACT_REGION_ASSEMBLIES.out
    pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir
    pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info
    pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING.out.targets_with_results

}


workflow EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES {
    take:
    pop_clustering_res_pop_clustering_dir //pop_clus_results
    pop_clustering_res_targets_with_results //the targets with results
    bed_fnp_ch  // Path to the BED file
    genome_fnp // Path to the genome file
    sub_var_regions //Path to results output directory

    main:
    targets = pop_clustering_res_targets_with_results.splitText()
                .map { it.trim() }

    targets_input_ch = targets.combine(pop_clustering_res_pop_clustering_dir)
        .map{ tar, pop_clustering_dir ->
            def sub_path_to_pop_seqs = "${tar}/population/popSeqsWithMetaWtihSampleName.fasta.gz"
            def full_path_to_pop_seqs = "${pop_clustering_dir}/${sub_path_to_pop_seqs}"
            tuple(tar, full_path_to_pop_seqs)
        }.combine(bed_fnp_ch).map{tar_name, pop_seqs_fnp, bed_fnp ->
            tuple(
                file("${pop_seqs_fnp}"),
                file("${bed_fnp}"),
                file("${genome_fnp}").getParent(),
                file("${genome_fnp}").getBaseName(),
                tar_name,
                params.correction_occurence_cut_off,
                params.low_freq_cut_off
            )
        }

    GET_SUB_SEGMENTS_FROM_FASTA(targets_input_ch)
    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")

    CONCATENATE_SUB_SEGMENT_LOCS(
        GET_SUB_SEGMENTS_FROM_FASTA.out.ref_variable_expanded_genomic_0_bed | collect,
        GET_SUB_SEGMENTS_FROM_FASTA.out.ref_sharedLocs_genomic_0_bed | collect,
        gff_fnp,
        sub_var_regions.toString()
    )
    emit:
    variable_regions = CONCATENATE_SUB_SEGMENT_LOCS.out.all_variable_regions
    all_conserved_regions = CONCATENATE_SUB_SEGMENT_LOCS.out.all_conserved_regions
    combined_all_sub_regions = CONCATENATE_SUB_SEGMENT_LOCS.out.combined_all_sub_regions
}

workflow PATHWEAVER_EXTRACT_REGIONS_WITH_GENE_IDS_FULL {
    take:
    samples_file // Path to the list of sample names (one per line)
    bams_dir     // Directory containing BAM files
    gene_ids_fnp  // Path to the BED file
    genome_fnp   // Path to the genome file
    results_dir  // Directory to store results
    meta_fnp     // meta file

    main:

    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_dir_info = file("${genome_fnp}").getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_twobit_fnp = file("${genome_dir_info}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }

    GET_GENE_RECORDS_FOR_GENE_IDS(file(gene_ids_fnp), genome_twobit_fnp, gff_fnp)
    REMOVE_TANDEM_REPEATS_FROM_REGIONS(GET_GENE_RECORDS_FOR_GENE_IDS.out.out_allTranscripts_bed, genome_twobit_fnp)



    regionsInfoDir = file("${results_dir}/regions/")
    regionsInfoDir.mkdirs()
    GET_INTERSECTING_GENE_INFO_FOR_REGIONS(REMOVE_TANDEM_REPEATS_FROM_REGIONS.out.out_with_tandems_removed, gff_fnp, genome_twobit_fnp, regionsInfoDir)

    full_results_dir = file("${results_dir}/PathWeaverResults")
    PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER(samples_file, bams_dir, REMOVE_TANDEM_REPEATS_FROM_REGIONS.out.out_with_tandems_removed, genome_fnp, full_results_dir, meta_fnp)


    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/PathWeaverResults/variantCalls")
        meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
        }
        VARIANT_CALL_ON_HAP_TABLE (
            REMOVE_TANDEM_REPEATS_FROM_REGIONS.out.out_with_tandems_removed,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            params.vc_getting_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
    }

    sub_var_regions = file("${full_results_dir}/subVarRegions/")
    sub_var_regions_raw = file("${sub_var_regions}/rawResults/")
    sub_var_regions_pop_clus_dir = file("${sub_var_regions}/populationClustering/")

    if (params.run_sub_segments_determination){
        sub_var_regions_pop_clus_dir.mkdirs()
        sub_var_regions.mkdirs()

        // extract out regions
        EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES(
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            REMOVE_TANDEM_REPEATS_FROM_REGIONS.out.out_with_tandems_removed,
            genome_fnp,
            sub_var_regions)

        // Load samples from file and create a channel
        var_dir_name =  file("${results_dir}").baseName + "_variableRegions"
        input_ch_with_var_regions = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.samples.combine(EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions)
            .map{ samp, var_bed_fnp ->
                tuple(samp, var_bed_fnp)
            }.map{
                tuple(it[0],
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}"),
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}.bai"),
                    file("${it[1]}"),
                    file("${genome_fnp}").getParent(),
                    file("${genome_fnp}").getBaseName(),
                    var_dir_name
                )
            }
        // run PathWeaver on the small variable regions
        EXTRACT_REGION_ASSEMBLIES(input_ch_with_var_regions)

        // concatenate the results files
        CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_raw.toString())

        // run population clustering
        PATHWEAVER_POP_CLUSTERING(var_dir_name, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_pop_clus_dir.toString())

        if (params.do_variant_calling){
            subvars_variant_call_dir = file("${sub_var_regions}/variantCalls")
            meta_fnp_for_variant_calling2_ch = Channel.fromPath(params.meta_fnp)
            if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
                //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
                meta_fnp_for_variant_calling2_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
            }
            VARIANT_CALL_ON_HAP_TABLE_2 (
                EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions,
                PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
                file("${genome_fnp}").getParent(),
                genome_base_name,
                file(gff_fnp),
                file(known_amino_acid_changes_fnp),
                params.vc_variant_frequency_cut_off,
                params.vc_variant_occurrence_cut_off,
                meta_fnp_for_variant_calling2_ch,
                params.vc_getting_pairwise_comps,
                params.meta_fields_to_calc_pop_diffs,
                file(subvars_variant_call_dir),
                params.vc_extra_args
            )
        }
    }

    workflow.onComplete {
        def outputDir = file("${results_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params()
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime()
    }
}

workflow PATHWEAVER_EXTRACT_REGIONS_WITH_SEQS_TABLE_FULL {
    take:
    samples_file // Path to the list of sample names (one per line)
    bams_dir     // Directory containing BAM files
    seqs_table_fnp  // Path to the BED file
    seqs_table_seqs_col //
    seqs_table_name_col //
    seqs_table_target_col //
    genome_fnp   // Path to the genome file
    results_dir  // Directory to store results
    meta_fnp     // meta file

    main:

    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_dir_info = file("${genome_fnp}").getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_twobit_fnp = file("${genome_dir_info}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }
    seqsInfoDir = file("${results_dir}/seqsTableLocationInfo/")
    seqsInfoDir.mkdirs()
    DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE(file(seqs_table_fnp), file(genome_fnp), seqs_table_seqs_col, seqs_table_name_col, seqs_table_target_col, params.resources.max_cpus, seqsInfoDir)

    regionsInfoDir = file("${results_dir}/regions/")
    regionsInfoDir.mkdirs()
    GET_INTERSECTING_GENE_INFO_FOR_REGIONS(DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE.out.targets_bed, gff_fnp, genome_twobit_fnp, regionsInfoDir)

    full_results_dir = file("${results_dir}/PathWeaverResults")
    PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER(samples_file, bams_dir, DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE.out.targets_bed, genome_fnp, full_results_dir, meta_fnp)


    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/PathWeaverResults/variantCalls")
        meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
        }
        VARIANT_CALL_ON_HAP_TABLE (
            DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE.out.targets_bed,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            params.vc_getting_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
    }

    sub_var_regions = file("${full_results_dir}/subVarRegions/")
    sub_var_regions_raw = file("${sub_var_regions}/rawResults/")
    sub_var_regions_pop_clus_dir = file("${sub_var_regions}/populationClustering/")

    if (params.run_sub_segments_determination){
        sub_var_regions_pop_clus_dir.mkdirs()
        sub_var_regions.mkdirs()

        // extract out regions
        EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES(
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE.out.targets_bed,
            genome_fnp,
            sub_var_regions)

        // Load samples from file and create a channel
        var_dir_name =  file("${results_dir}").baseName + "_variableRegions"
        input_ch_with_var_regions = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.samples.combine(EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions)
            .map{ samp, var_bed_fnp ->
                tuple(samp, var_bed_fnp)
            }.map{
                tuple(it[0],
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}"),
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}.bai"),
                    file("${it[1]}"),
                    file("${genome_fnp}").getParent(),
                    file("${genome_fnp}").getBaseName(),
                    var_dir_name
                )
            }
        // run PathWeaver on the small variable regions
        EXTRACT_REGION_ASSEMBLIES(input_ch_with_var_regions)

        // concatenate the results files
        CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_raw.toString())

        // run population clustering
        PATHWEAVER_POP_CLUSTERING(var_dir_name, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_pop_clus_dir.toString())

        if (params.do_variant_calling){
            subvars_variant_call_dir = file("${sub_var_regions}/variantCalls")
            meta_fnp_for_variant_calling2_ch = Channel.fromPath(params.meta_fnp)
            if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
                //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
                meta_fnp_for_variant_calling2_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
            }
            VARIANT_CALL_ON_HAP_TABLE_2 (
                EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions,
                PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
                file("${genome_fnp}").getParent(),
                genome_base_name,
                file(gff_fnp),
                file(known_amino_acid_changes_fnp),
                params.vc_variant_frequency_cut_off,
                params.vc_variant_occurrence_cut_off,
                meta_fnp_for_variant_calling2_ch,
                params.vc_getting_pairwise_comps,
                params.meta_fields_to_calc_pop_diffs,
                file(subvars_variant_call_dir),
                params.vc_extra_args
            )
        }
    }

    workflow.onComplete {
        def outputDir = file("${results_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params()
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime()
    }
}

workflow PATHWEAVER_EXTRACT_REGIONS_WITH_PRIMERS_FULL {
    take:
    samples_file // Path to the list of sample names (one per line)
    bams_dir     // Directory containing BAM files
    primers_fnp  // Path to the BED file
    genome_fnp   // Path to the genome file
    results_dir  // Directory to store results
    meta_fnp     // meta file

    main:

    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_dir_info = file("${genome_fnp}").getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_twobit_fnp = file("${genome_dir_info}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }
    primerInfoDir = file("${results_dir}/primerInfo/")
    gff_dir = file("${top_genome_info}/info/gff/")
    genome_dir = file("${genome_fnp}").getParent()
    primerInfoDir.mkdirs()
    GEN_TARGET_INFO_FROM_GENOME(primers_fnp, primerInfoDir, genome_dir, gff_dir, genome_base_name, params.nanopore_primers_errors_allowed, params.resources.max_cpus)

    regionsInfoDir = file("${results_dir}/regions/")
    regionsInfoDir.mkdirs()
    GET_INTERSECTING_GENE_INFO_FOR_REGIONS(GEN_TARGET_INFO_FROM_GENOME.out.inner_bed, gff_fnp, genome_twobit_fnp, regionsInfoDir)

    full_results_dir = file("${results_dir}/PathWeaverResults")
    PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER(samples_file, bams_dir, GEN_TARGET_INFO_FROM_GENOME.out.inner_bed, genome_fnp, full_results_dir, meta_fnp)


    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/PathWeaverResults/variantCalls")
        meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
        }
        VARIANT_CALL_ON_HAP_TABLE (
            GEN_TARGET_INFO_FROM_GENOME.out.inner_bed,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            params.vc_getting_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
    }

    sub_var_regions = file("${full_results_dir}/subVarRegions/")
    sub_var_regions_raw = file("${sub_var_regions}/rawResults/")
    sub_var_regions_pop_clus_dir = file("${sub_var_regions}/populationClustering/")

    if (params.run_sub_segments_determination){
        sub_var_regions_pop_clus_dir.mkdirs()
        sub_var_regions.mkdirs()

        // extract out regions
        EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES(
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            GEN_TARGET_INFO_FROM_GENOME.out.inner_bed,
            genome_fnp,
            sub_var_regions)

        // Load samples from file and create a channel
        var_dir_name =  file("${results_dir}").baseName + "_variableRegions"
        input_ch_with_var_regions = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.samples.combine(EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions)
            .map{ samp, var_bed_fnp ->
                tuple(samp, var_bed_fnp)
            }.map{
                tuple(it[0],
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}"),
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}.bai"),
                    file("${it[1]}"),
                    file("${genome_fnp}").getParent(),
                    file("${genome_fnp}").getBaseName(),
                    var_dir_name
                )
            }
        // run PathWeaver on the small variable regions
        EXTRACT_REGION_ASSEMBLIES(input_ch_with_var_regions)

        // concatenate the results files
        CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_raw.toString())

        // run population clustering
        PATHWEAVER_POP_CLUSTERING(var_dir_name, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_pop_clus_dir.toString())

        if (params.do_variant_calling){
            subvars_variant_call_dir = file("${sub_var_regions}/variantCalls")
            meta_fnp_for_variant_calling2_ch = Channel.fromPath(params.meta_fnp)
            if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
                //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
                meta_fnp_for_variant_calling2_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
            }
            VARIANT_CALL_ON_HAP_TABLE_2 (
                EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions,
                PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
                file("${genome_fnp}").getParent(),
                genome_base_name,
                file(gff_fnp),
                file(known_amino_acid_changes_fnp),
                params.vc_variant_frequency_cut_off,
                params.vc_variant_occurrence_cut_off,
                meta_fnp_for_variant_calling2_ch,
                params.vc_getting_pairwise_comps,
                params.meta_fields_to_calc_pop_diffs,
                file(subvars_variant_call_dir),
                params.vc_extra_args
            )
        }
    }

    workflow.onComplete {
        def outputDir = file("${results_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params()
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime()
    }
}

workflow PATHWEAVER_EXTRACT_REGIONS_FULL {
    take:
    samples_file // Path to the list of sample names (one per line)
    bams_dir     // Directory containing BAM files
    input_bed_fnp// Path to the BED file
    genome_fnp   // Path to the genome file
    results_dir  // Directory to store results
    meta_fnp     // meta file

    main:

    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_dir_info = file("${genome_fnp}").getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_twobit_fnp = file("${genome_dir_info}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }

    full_results_dir = file("${results_dir}/PathWeaverResults")
    bed_fnp_ch = Channel.fromPath(input_bed_fnp)

    regionsInfoDir = file("${results_dir}/regions/")
    regionsInfoDir.mkdirs()
    GET_INTERSECTING_GENE_INFO_FOR_REGIONS(bed_fnp_ch, gff_fnp, genome_twobit_fnp, regionsInfoDir)

    PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER(samples_file, bams_dir, bed_fnp_ch, genome_fnp, full_results_dir, meta_fnp)


    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/PathWeaverResults/variantCalls")
        meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
        }
        VARIANT_CALL_ON_HAP_TABLE (
            bed_fnp_ch,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            params.vc_getting_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
    }

    sub_var_regions = file("${full_results_dir}/subVarRegions/")
    sub_var_regions_raw = file("${sub_var_regions}/rawResults/")
    sub_var_regions_pop_clus_dir = file("${sub_var_regions}/populationClustering/")

    if (params.run_sub_segments_determination){
        sub_var_regions_pop_clus_dir.mkdirs()
        sub_var_regions.mkdirs()

        // extract out regions
        EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES(
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            bed_fnp_ch,
            genome_fnp,
            sub_var_regions)

        // Load samples from file and create a channel
        var_dir_name =  file("${results_dir}").baseName + "_variableRegions"
        input_ch_with_var_regions = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.samples.combine(EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions)
            .map{ samp, var_bed_fnp ->
                tuple(samp, var_bed_fnp)
            }.map{
                tuple(it[0],
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}"),
                    file("${bams_dir}/${it[0]}${params.bams_file_ending}.bai"),
                    file("${it[1]}"),
                    file("${genome_fnp}").getParent(),
                    file("${genome_fnp}").getBaseName(),
                    var_dir_name
                )
            }
        // run PathWeaver on the small variable regions
        EXTRACT_REGION_ASSEMBLIES(input_ch_with_var_regions)

        // concatenate the results files
        CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_raw.toString())

        // run population clustering
        PATHWEAVER_POP_CLUSTERING(var_dir_name, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_pop_clus_dir.toString())

        if (params.do_variant_calling){
            subvars_variant_call_dir = file("${sub_var_regions}/variantCalls")
            meta_fnp_for_variant_calling2_ch = Channel.fromPath(params.meta_fnp)
            if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
                //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
                meta_fnp_for_variant_calling2_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
            }
            VARIANT_CALL_ON_HAP_TABLE_2 (
                EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions,
                PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
                file("${genome_fnp}").getParent(),
                genome_base_name,
                file(gff_fnp),
                file(known_amino_acid_changes_fnp),
                params.vc_variant_frequency_cut_off,
                params.vc_variant_occurrence_cut_off,
                meta_fnp_for_variant_calling2_ch,
                params.vc_getting_pairwise_comps,
                params.meta_fields_to_calc_pop_diffs,
                file(subvars_variant_call_dir),
                params.vc_extra_args
            )
        }
    }

    workflow.onComplete {
        def outputDir = file("${results_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params()
        record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime()
    }
}



def record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_params() {
    def output = file("${params.pw_results_dir}/run/parameters.tsv")
    output.withWriter { writer ->
        params.each { k, v ->
            writer.println("${k}\t${v}")}
    }
}

/* Record runtime information
 *
 * Records runtime and environment information and writes summary to a tabulated (tsv) file
 *
 */
def record_PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_runtime() {

    def output = file("${params.pw_results_dir}/run/runtime.tsv")
    output.withWriter { writer ->
        writer.println("PipelineVersion\t${workflow.manifest.version}")
        writer.println("ContainerEngine\t${workflow.containerEngine}")
        writer.println("Duration\t${workflow.duration}")
        writer.println("CommandLine\t${workflow.commandLine}")
        writer.println("CommitId\t${workflow.commitId}")
        writer.println("Complete\t${workflow.complete}")
        writer.println("ConfigFiles\t${workflow.configFiles.join(', ')}")
        writer.println("Container\t${workflow.container}")
        writer.println("ErrorMessage\t${workflow.errorMessage}")
        writer.println("ErrorReport\t${workflow.errorReport}")
        writer.println("ExitStatus\t${workflow.exitStatus}")
        writer.println("HomeDir\t${workflow.homeDir}")
        writer.println("LaunchDir\t${workflow.launchDir}")
        writer.println("Manifest\t${workflow.manifest}")
        writer.println("Profile\t${workflow.profile}")
        writer.println("ProjectDir\t${workflow.projectDir}")
        writer.println("Repository\t${workflow.repository}")
        writer.println("Resume\t${workflow.resume}")
        writer.println("Revision\t${workflow.revision}")
        writer.println("RunName\t${workflow.runName}")
        writer.println("ScriptFile\t${workflow.scriptFile}")
        writer.println("ScriptId\t${workflow.scriptId}")
        writer.println("ScriptName\t${workflow.scriptName}")
        writer.println("SessionId\t${workflow.sessionId}")
        writer.println("Start\t${workflow.start}")
        writer.println("StubRun\t${workflow.stubRun}")
        writer.println("Success\t${workflow.success}")
        writer.println("UserName\t${workflow.userName}")
        writer.println("WorkDir\t${workflow.workDir}")
        writer.println("NextflowBuild\t${nextflow.build}")
        writer.println("NextflowTimestamp\t${nextflow.timestamp}")
        writer.println("NextflowVersion\t${nextflow.version}")
    }
}
