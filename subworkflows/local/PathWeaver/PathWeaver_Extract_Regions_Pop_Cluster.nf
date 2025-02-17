#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// include { paramsSummaryMap       } from 'plugin/nf-schema'


include { EXTRACT_REGION_ASSEMBLIES } from '../../../modules/local/extract_region_assemblies'
include { PATHWEAVER_POP_CLUSTERING } from '../../../modules/local/pathweaver_pop_clustering'
include { GET_SUB_SEGMENTS_FROM_FASTA } from '../../../modules/local/get_sub_segments_from_PW_pop_res/main.nf'
include { CONCATENATE_SUB_SEGMENT_LOCS } from '../../../modules/local/concatenate_sub_segment_locs/main.nf'
include { CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS } from '../../../modules/local/concatenate_extract_region_assemblies_results/main.nf'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { VARIANT_CALL_ON_HAP_TABLE as VARIANT_CALL_ON_HAP_TABLE_SUBREGIONS } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { CREATE_TWOBIT_FILE_FOR_GENOME } from '../../../modules/local/create_twobit_file_for_genome/main.nf'
include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE } from '../../../modules/local/gen_target_info_from_genomes'
include { GEN_TARGET_INFO_FROM_GENOME_NANOPORE as GEN_TARGET_INFO_FROM_GENOME } from '../../../modules/local/gen_target_info_from_genomes/main.nf'
include { DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE } from '../../../modules/local/determine_genomic_location_from_seqs_table/main.nf'
include { GET_INTERSECTING_GENE_INFO_FOR_REGIONS } from '../../../modules/local/get_intersecting_gene_info_for_regions/main.nf'
include { GET_INTERSECTING_GENE_INFO_FOR_REGIONS as GET_INTERSECTING_GENE_INFO_FOR_REGIONS_SUBREGIONS } from '../../../modules/local/get_intersecting_gene_info_for_regions/main.nf'

include { GET_GENE_RECORDS_FOR_GENE_IDS } from '../../../modules/local/get_gene_records_for_gene_ids/main.nf'
include { REMOVE_TANDEM_REPEATS_FROM_REGIONS } from '../../../modules/local/remove_tandem_repeats_from_regions/main.nf'
include { PW_CREATE_DASHBOARDS } from '../../../modules/local/pw_create_dashboards/main.nf'
include { PW_CREATE_DASHBOARDS as PW_CREATE_DASHBOARDS_SUBREGIONS } from '../../../modules/local/pw_create_dashboards/main.nf'
include { PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER } from "./pathweaver_extract_regions_and_pop_cluster.nf"
include { PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER as PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_SUBREGIONS } from "./pathweaver_extract_regions_and_pop_cluster.nf"



workflow PATHWEAVER_EXTRACT_REGIONS_WITH_BED {
    take:
    samples_file // Path to the list of sample names (one per line)
    bams_dir     // Directory containing BAM files
    bed_ch       // a channel that emits a bed file
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
    full_results_dir.mkdirs()
    regionsInfoDir = file("${results_dir}/regions/")
    regionsInfoDir.mkdirs()

    GET_INTERSECTING_GENE_INFO_FOR_REGIONS(bed_ch, gff_fnp, genome_twobit_fnp, regionsInfoDir)

    PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER(samples_file, bams_dir, bed_ch, genome_fnp, full_results_dir, meta_fnp)

    if (params.do_variant_calling){
        //copy over dashboard quarto document
        PW_CREATE_DASHBOARDS(file("${results_dir}"),  file("${projectDir}/etc/pw_basic_report.qmd"),
            params.render_pw_report,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            GET_INTERSECTING_GENE_INFO_FOR_REGIONS.out.bed_withGeneInfo_tsv,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.variant_calling_reports_dir
        )
    }

    if (params.run_sub_segments_determination){
        sub_var_results_dir = file("${full_results_dir}/subVarRegions/")
        sub_var_full_results_dir = file("${sub_var_results_dir}/PathWeaverResults/")
        sub_var_full_results_dir.mkdirs()
        sub_var_regionsInfoDir = file("${sub_var_results_dir}/regions/")
        sub_var_regionsInfoDir.mkdirs()

        // extract out regions
        EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES(
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_pop_clustering_dir,
            PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER.out.pop_clustering_res_targets_with_results,
            bed_ch,
            genome_fnp,
            sub_var_results_dir)

        GET_INTERSECTING_GENE_INFO_FOR_REGIONS_SUBREGIONS(EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions, gff_fnp, genome_twobit_fnp, sub_var_regionsInfoDir)

        PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_SUBREGIONS(samples_file, bams_dir, EXTRACT_VARIABLE_REGIONS_FROM_PATHWEAVER_ASSEMBLIES.out.variable_regions, genome_fnp, sub_var_full_results_dir, meta_fnp)

        if (params.do_variant_calling){
            //copy over dashboard quarto document
            PW_CREATE_DASHBOARDS_SUBREGIONS(sub_var_results_dir,
                file("${projectDir}/etc/pw_basic_report.qmd"),
                params.render_pw_report,
                PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_SUBREGIONS.out.pop_clustering_res_targets_with_results,
                GET_INTERSECTING_GENE_INFO_FOR_REGIONS_SUBREGIONS.out.bed_withGeneInfo_tsv,
                PATHWEAVER_EXTRACT_REGIONS_AND_POP_CLUSTER_SUBREGIONS.out.variant_calling_reports_dir
            )
        }
    }



    workflow.onComplete {
        def outputDir = file("${results_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_PATHWEAVER_EXTRACT_REGIONS_WITH_BED_params()
        record_PATHWEAVER_EXTRACT_REGIONS_WITH_BED_runtime()
    }
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


    GET_GENE_RECORDS_FOR_GENE_IDS(file(gene_ids_fnp), genome_twobit_fnp, gff_fnp)
    REMOVE_TANDEM_REPEATS_FROM_REGIONS(GET_GENE_RECORDS_FOR_GENE_IDS.out.out_allTranscripts_bed, genome_twobit_fnp)

    PATHWEAVER_EXTRACT_REGIONS_WITH_BED(samples_file, bams_dir, GET_GENE_RECORDS_FOR_GENE_IDS.out.out_allTranscripts_bed, genome_fnp, results_dir, meta_fnp)
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

    seqsInfoDir = file("${results_dir}/seqsTableLocationInfo/")
    seqsInfoDir.mkdirs()
    DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE(file(seqs_table_fnp), file(genome_fnp), seqs_table_seqs_col, seqs_table_name_col, seqs_table_target_col, params.resources.max_cpus, seqsInfoDir)

    PATHWEAVER_EXTRACT_REGIONS_WITH_BED(samples_file, bams_dir, DETERMINE_GENOMIC_LOCATION_FROM_SEQS_TABLE.out.targets_bed, genome_fnp, results_dir, meta_fnp)
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
    genome_base_name = file("${genome_fnp}").getBaseName()
    primerInfoDir = file("${results_dir}/primerInfo/")
    gff_dir = file("${top_genome_info}/info/gff/")
    genome_dir = file("${genome_fnp}").getParent()
    primerInfoDir.mkdirs()
    GEN_TARGET_INFO_FROM_GENOME(primers_fnp, primerInfoDir, genome_dir, gff_dir, genome_base_name, params.nanopore_primers_errors_allowed, params.resources.max_cpus)

    PATHWEAVER_EXTRACT_REGIONS_WITH_BED(samples_file, bams_dir, GEN_TARGET_INFO_FROM_GENOME.out.inner_bed, genome_fnp, results_dir, meta_fnp)
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

    bed_fnp_ch = Channel.fromPath(input_bed_fnp)

    PATHWEAVER_EXTRACT_REGIONS_WITH_BED(samples_file, bams_dir, bed_fnp_ch, genome_fnp, results_dir, meta_fnp)
}

def record_PATHWEAVER_EXTRACT_REGIONS_WITH_BED_params() {
    def output = file("${params.outdir}/run/parameters.tsv")
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
def record_PATHWEAVER_EXTRACT_REGIONS_WITH_BED_runtime() {

    def output = file("${params.outdir}/run/runtime.tsv")
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
