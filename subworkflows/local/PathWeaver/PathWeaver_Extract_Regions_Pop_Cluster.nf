#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EXTRACT_REGION_ASSEMBLIES } from '../../../modules/local/extract_region_assemblies'
include { PATHWEAVER_POP_CLUSTERING } from '../../../modules/local/pathweaver_pop_clustering'
include { GET_SUB_SEGMENTS_FROM_FASTA } from '../../../modules/local/get_sub_segments_from_PW_pop_res/main.nf'
include { CONCATENATE_SUB_SEGMENT_LOCS } from '../../../modules/local/concatenate_sub_segment_locs/main.nf'
include { CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS } from '../../../modules/local/concatenate_extract_region_assemblies_results/main.nf'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'


workflow PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp  // Path to the BED file
    genome_fnp // Path to the genome file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:


    // Create output directories
    def pop_clus_dir = file("${results_dir}/reports/populationClustering/")
    def reports_dir = file("${results_dir}/reports/")
    pop_clus_dir.mkdirs()

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
    input_ch = samples.map { samp ->
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
    bed_fnp  // Path to the BED file
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
        }.map{
            tuple(
            file("${it[1]}"),
            file("${bed_fnp}"),
            file("${genome_fnp}").getParent(),
            file("${genome_fnp}").getBaseName(),
            it[0],
            params.correction_occurence_cut_off,
            params.low_freq_cut_off
            )
        }

    GET_SUB_SEGMENTS_FROM_FASTA(targets_input_ch)
    def top_genome_info = file("${genome_fnp}").getParent().getParent()
    def genome_base_name = file("${genome_fnp}").getBaseName()
    def gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")

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

workflow PATHWEAVER_EXTRACT_REGIONS_FULL {
    take:
    samples_file  // Path to the list of sample names (one per line)
    bams_dir  // Directory containing BAM files
    bed_fnp  // Path to the BED file
    genome_fnp // Path to the genome file
    results_dir // Directory to store results
    meta_fnp //meta file

    main:

    def top_genome_info = file("${genome_fnp}").getParent().getParent()
    def genome_base_name = file("${genome_fnp}").getBaseName()
    def gff_fnp = file("${top_genome_info}/info/gff/${genome_base_name}.gff")
    def known_amino_acid_changes_fnp = file("${top_genome_info}/info/pf_drug_resistant_aaPositions.tsv")

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(results_dir)
    if (results_dir_obj.exists()){
        results_dir_obj.deleteDir()
    }
    results_dir_obj.mkdirs()
    def full_results_dir = file("${results_dir}/PathWeaverResults")
    PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER(samples_file, bams_dir, bed_fnp, genome_fnp, full_results_dir, meta_fnp)


    if (params.do_variant_calling){
        def variant_call_dir = file("${results_dir}/PathWeaverResults/variantCalls")
        def meta_fnp_for_variant_calling_ch = Channel.fromPath(params.meta_fnp)
        if ("EMPTY_FILE.txt" != file(params.meta_fnp).baseName ){
            //meta data was supplied, should use the meta data from the population clustering because it will sometimes filter and collapse samples
            meta_fnp_for_variant_calling_ch = PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir.map{file("${it[0]}/info/sampleMetaData.tab.txt")}
        }
        VARIANT_CALL_ON_HAP_TABLE (
            file(bed_fnp),
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

    def sub_var_regions = file("${full_results_dir}/subVarRegions/")
    def sub_var_regions_pop_clus_dir = file("${sub_var_regions}/populationClustering/")

    if (params.run_sub_segments_determination){
        sub_var_regions_pop_clus_dir.mkdirs()

        // extract out regions
        EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES(
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            bed_fnp,
            genome_fnp,
            sub_var_regions)

        // Load samples from file and create a channel
        def var_dir_name =  file("${results_dir}").baseName + "_variableRegions"
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
        CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS(EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions.toString())

        // run population clustering
        PATHWEAVER_POP_CLUSTERING(var_dir_name, file("${meta_fnp}"), EXTRACT_REGION_ASSEMBLIES.out | collect, sub_var_regions_pop_clus_dir.toString())
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
