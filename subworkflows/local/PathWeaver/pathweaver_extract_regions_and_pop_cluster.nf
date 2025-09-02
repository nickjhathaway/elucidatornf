#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { EXTRACT_REGION_ASSEMBLIES } from '../../../modules/local/extract_region_assemblies'
include { CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS } from '../../../modules/local/concatenate_extract_region_assemblies_results/main.nf'
include { PATHWEAVER_POP_CLUSTERING } from '../../../modules/local/pathweaver_pop_clustering'
include { PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED } from '../../../modules/local/pathweaver_pop_clustering'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { VARIANT_CALL_ON_HAP_TABLE as VARIANT_CALL_ON_HAP_TABLE_ON_TRIMMED } from '../../../modules/local/variant_call_on_hap_table/main.nf'



workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp_ch  // Path to the BED file

    genome_fnp // Path to the genome file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:

    // get possible genome
    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }

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

    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/variantCalls")
        meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        log.info(file(params.meta_fnp).baseName)
        log.info(file(params.meta_fnp).name)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).name ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{file("${it}/info/sampleMetaData.tab.txt")}
        }

        VARIANT_CALL_ON_HAP_TABLE (
            bed_fnp_ch,
            PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            !params.vc_no_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
    }
    emit:
    samples = samples
    assemblies = EXTRACT_REGION_ASSEMBLIES.out
    pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir
    pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info
    pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING.out.targets_with_results
    variant_calling_reports_dir = VARIANT_CALL_ON_HAP_TABLE.out.reports
}
workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_WITH_TRIM_BED {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp_ch  // Path to the BED file
    trim_bed_fnp_ch  // Path to the BED file
    genome_fnp // Path to the genome file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:

    // get possible genome
    top_genome_info = file("${genome_fnp}").getParent().getParent()
    genome_base_name = file("${genome_fnp}").getBaseName()
    genome_dir_info = file("${genome_fnp}").getParent()
    genome_twobit_fnp = file("${genome_dir_info}/${genome_base_name}.2bit")
    gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${top_genome_info}/info/drug_resistant_aaPositions.tsv").exists()){
        known_amino_acid_changes_fnp = file("${top_genome_info}/info/drug_resistant_aaPositions.tsv")
    }

    // Create output directories
    full_pop_clus_dir = file("${results_dir}/reports/full/populationClustering/")
    trimmed_pop_clus_dir = file("${results_dir}/reports/populationClustering/")
    reports_dir = file("${results_dir}/reports/rawResults/")
    full_pop_clus_dir.mkdirs()
    trimmed_pop_clus_dir.mkdirs()
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

    //run population clustering
    //run with full region
    PATHWEAVER_POP_CLUSTERING(              file("${results_dir}").baseName, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, full_pop_clus_dir.toString())
    //run with inner trimmed region
    PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED(file("${results_dir}").baseName, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, trimmed_pop_clus_dir.toString(),
                                            trim_bed_fnp_ch, genome_twobit_fnp)

    if (params.do_variant_calling){
        variant_call_dir = file("${results_dir}/reports/full/variantCalls")
        trimmed_variant_call_dir = file("${results_dir}/variantCalls")
        meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        trimmed_meta_fnp_for_variant_calling_ch =  Channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).name ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir.map{file("${it}/info/sampleMetaData.tab.txt")}
            trimmed_meta_fnp_for_variant_calling_ch = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.pop_clustering_dir.map{file("${it}/info/sampleMetaData.tab.txt")}
        }
        VARIANT_CALL_ON_HAP_TABLE (
            bed_fnp_ch,
            PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp_for_variant_calling_ch,
            false,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
        VARIANT_CALL_ON_HAP_TABLE_ON_TRIMMED (
            trim_bed_fnp_ch,
            PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.all_selected_clusters_info,
            file("${genome_fnp}").getParent(),
            genome_base_name,
            file(gff_fnp),
            file(known_amino_acid_changes_fnp),
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            trimmed_meta_fnp_for_variant_calling_ch,
            !params.vc_no_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(trimmed_variant_call_dir),
            params.vc_extra_args
        )
    }
    emit:
    samples = samples
    assemblies = EXTRACT_REGION_ASSEMBLIES.out
    pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING.out.pop_clustering_dir
    pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING.out.all_selected_clusters_info
    pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING.out.targets_with_results
    variant_calling_reports_dir = VARIANT_CALL_ON_HAP_TABLE.out.reports

    trimmed_pop_clustering_res_pop_clustering_dir = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.pop_clustering_dir
    trimmed_pop_clustering_res_all_selected_clusters_info = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.all_selected_clusters_info
    trimmed_pop_clustering_res_targets_with_results = PATHWEAVER_POP_CLUSTERING_WITH_TRIM_BED.out.targets_with_results
    trimmed_variant_calling_reports_dir = VARIANT_CALL_ON_HAP_TABLE_ON_TRIMMED.out.reports
}
