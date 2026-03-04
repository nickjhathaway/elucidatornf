#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// include { paramsSummaryMap       } from 'plugin/nf-schema'


include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE as GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL } from '../../../modules/local/gen_target_info_from_genomes'
include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE } from '../../../modules/local/gen_target_info_from_genomes'
include { EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/extractor_by_kmer_matching'
include { VARIANT_CALL_ON_HAP_TABLE } from '../../../modules/local/variant_call_on_hap_table/main.nf'
include { CONCATENATE_EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/concatenate_extractor_by_kmer_matching/main.nf'
include { AMPLICON_CLUSTER_BY_KMER_SIMILARITY } from '../../../modules/local/amplicon_cluster_by_kmer_similarity/main.nf'
include { AMPLICON_POPULATION_CLUSTERING } from '../../../modules/local/amplicon_population_cluster/main.nf'
include { CONCATENATE_AMPLICON_POPULATION_CLUSTERING } from '../../../modules/local/concatenate_amplicon_population_clustering/main.nf'
include { GEN_TARGET_INFO_FROM_GENOMES_ILLUMINA } from '../../../modules/local/gen_target_info_from_genomes/main.nf'
include { AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR } from "../../../modules/local/amplicon_cluster_auto_SeekDeep_flag_generator/main.nf"
include { AMPLICON_NANOPORE_CREATE_DASHBOARDS } from "../../../modules/local/amplicon_nanopore_create_dashboards/main.nf"

include { completionSummary } from '../../../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'


process WRITE_SAMPLE_NAME_KEY {
  tag "sample_name_key"

  publishDir { meta_dir }, mode: 'copy', overwrite: true

  input:
    val lines
    val meta_dir

  output:
    path "sample_name_key.tsv"

  script:
    def header = "old_sample_name\tnew_sample_name\n"

    // sort by old_sample_name (first column)
    def sorted = (lines ?: []).sort { a, b ->
      def aOld = a.split('\t', 2)[0]
      def bOld = b.split('\t', 2)[0]
      aOld <=> bOld
    }

    def body = sorted.join('\n')
    def text = header + (body ? body : "")

    """
    cat > sample_name_key.tsv <<'EOF'
${text}
EOF
    """
}




workflow NANOPORE_AMPLICON_CLUSTERING {
    take:
    input_fastq_dir_raw //a directory with all fastqs to be analyzed
    primers_fnp_raw  // a file with primer info, 3 columsn, target,fwd_primer,rev_primer
    genome_dir_raw  // a directory with genomes
    gff_dir_raw // a directory with annotation gff files
    primers_errors_allowed // errors to allow in primers
    meta_fnp_raw //a sample meta file
    output_dir_raw // the output directory
    main:

    // convert into file type objects
    genome_dir = file(genome_dir_raw)
    primers_fnp = file(primers_fnp_raw)
    input_fastq_dir = file(input_fastq_dir_raw)
    gff_dir = file(gff_dir_raw)
    meta_fnp = file(meta_fnp_raw)
    output_dir = file(output_dir_raw)
    //
    // Set up for validation
    //
    errors = []


    // Map from internal var -> CLI-style name (for messages)
    required_params = [
        input_fastq_dir : '--input_fastq_dir',
        primers_fnp     : '--primers_fnp',
        genome_dir      : '--genome_dir',
        gff_dir         : '--gff_dir',
        output_dir      : '--outdir'
    ]

    // Map from internal var -> actual value passed to subworkflow
    inputs = [
        input_fastq_dir : input_fastq_dir,
        primers_fnp     : primers_fnp,
        genome_dir      : genome_dir_raw,
        gff_dir         : gff_dir,
        output_dir      : output_dir
    ]
    //
    // 2) Existence checks for paths (only if non-blank)
    //
    path_checks = [
        input_fastq_dir : "directory",
        primers_fnp     : "file",
        genome_dir      : "directory",
        gff_dir         : "directory"
    ]

    path_checks.each { key, type ->
        def raw = inputs[key]
        def obj = file(raw)
        if (!obj.exists()) {
            errors << "Path does not exist for ${required_params[key]}: '${raw}'"
        } else if (type == "directory" && !obj.isDirectory()) {
            errors << "Expected a directory for ${required_params[key]} but found a file: '${raw}'"
        } else if (type == "file" && !obj.isFile()) {
            errors << "Expected a file for ${required_params[key]} but found a directory: '${raw}'"
        }
    }

    //
    // 3) Fail fast if validation failed
    //
    if (!errors.isEmpty()) {
        def joined = errors.collect {it -> "  - ${it}" }.join("\n")
        error """
Input validation failed in subworkflow NANOPORE_AMPLICON_CLUSTERING:

${joined}

Please correct the above issues and re-run.
"""
    }

    def outdir = output_dir.toString().trim()

    //check if inputs exists


    //
    // Create output directory (overwrite if exists)
    //
    def results_dir_obj = file(outdir)
    if (results_dir_obj.exists()) {
        results_dir_obj.deleteDir()
    }
    results_dir_obj.mkdirs()

    //
    // Inner primers logic
    //
    if (params.nanopore_extractor_use_inner_primers) {
        params.nanopore_clustering_lower_base = "upper"
    }

    def primer_info_dir = file("${output_dir}/primerInfo")
    primer_info_dir.mkdirs()

    def meta_info_dir = file("${output_dir}/meta")
    meta_info_dir.mkdirs()

    def extraction_reports_dir = file("${output_dir}/extractionReports")
    extraction_reports_dir.mkdirs()
    def final_results_dir = file("${output_dir}/finalResults")
    final_results_dir.mkdirs()

    // optional sample renaming map
    def rename_key = SampleRename.loadRenameKey(params.nanopore_rename_key_fnp)

    primers_fnp_ch = primers_fnp

    if (params.nanopore_extractor_use_inner_primers){
        def initial_primer_info_dir = file("${output_dir}/primerInfo/initialPrimerInfo")
        initial_primer_info_dir.mkdirs()
        GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL(primers_fnp, initial_primer_info_dir, genome_dir, gff_dir, primers_errors_allowed, params.resources.max_cpus)
        primers_fnp_ch = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL.out.inner_primers
    }
    GEN_TARGET_INFO_FROM_GENOMES_NANOPORE(primers_fnp_ch, primer_info_dir, genome_dir, gff_dir, primers_errors_allowed, params.resources.max_cpus)

    AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR(file("${input_fastq_dir}"), primers_fnp_ch, "nanopore", params.resources.max_cpus, extraction_reports_dir)

    fastq_input_ch = channel.fromPath(file("${input_fastq_dir}/*.fastq.gz"))

    input_to_extractor = fastq_input_ch.combine(GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info)
            .combine(AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR.out.out_seekdeep_extractor_flags)
            .map{fastq_fnp, info_dir, auto_flags_fnp ->
                tuple(
                    fastq_fnp,
                    file(fastq_fnp).name.replaceAll(".fastq.gz", ""),
                    primers_fnp,
                    file("${info_dir}/allKmers.tab.txt.gz"),
                    file("${info_dir}/lenCutOffs.txt"),
                    auto_flags_fnp,
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
    def targets = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.targets_with_extractions
        .splitText()
        .map{samp ->
            samp.trim() // Remove any whitespace
            }


    // create a list of input to cluster for the input that actually exists
    def clustered_inputs_with_key = targets.combine(EXTRACTOR_BY_KMER_MATCHING.out.sample_dir)
        .map { tar, samp, sampdir ->

            def extracted_tar_sample_fnp = file("${sampdir}/${tar}${samp}.fastq.gz")
            if( !extracted_tar_sample_fnp.exists() ) return null

            // apply optional renaming (fallback to original if not present)
            def samp_renamed = rename_key.containsKey(samp) ? rename_key[samp] : samp

            // emit BOTH: the clustering tuple + the mapping tuple
            tuple(
                tuple(extracted_tar_sample_fnp, samp_renamed, tar, params.nanopore_clustering_ncpus),
                tuple(samp, samp_renamed)
            )
        }
        .filter { entry -> entry != null }

    // split into the two channels we want
    def input_to_clustering = clustered_inputs_with_key.map { entry -> entry[0] }
    def sample_key_ch       = clustered_inputs_with_key.map { entry -> entry[1] }.distinct()

    // write sample rename key (old -> new) for samples that had extracted reads
    // sample_key_ch is a channel of tuples: (oldName, newName)
    def sample_key_lines_ch = sample_key_ch
        .map { oldName, newName -> "${oldName}\t${newName}" }
        .distinct()
        .collect()

    WRITE_SAMPLE_NAME_KEY(sample_key_lines_ch, meta_info_dir.toString())
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
        .combine(ref_seqs_dir_ch)
        .map { _samples, target, target_fastqs, ref_seqs_dir ->
                tuple(
                    target_fastqs, // List of fastq files for this target
                    target,
                    params.nanopore_clustering_ncpus,
                    meta_fnp,
                    // GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info,
                    ref_seqs_dir,
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
            !params.vc_no_pairwise_comps,
            params.meta_fields_to_calc_pop_diffs,
            file(variant_call_dir),
            params.vc_extra_args
        )
        AMPLICON_NANOPORE_CREATE_DASHBOARDS(file("${output_dir}"),
            file("${projectDir}/etc/nanopore_clustering_report.qmd"),
            params.nanopore_render_clustering_report,
            VARIANT_CALL_ON_HAP_TABLE.out.reports,
            CONCATENATE_EXTRACTOR_BY_KMER_MATCHING.out.all_extraction_stats
        )
    } else {
        AMPLICON_NANOPORE_CREATE_DASHBOARDS(file("${output_dir}"),
            file("${projectDir}/etc/nanopore_clustering_report.qmd"),
            params.nanopore_render_clustering_report,
            CONCATENATE_AMPLICON_POPULATION_CLUSTERING.out.all_selected_clusters_info,
            CONCATENATE_EXTRACTOR_BY_KMER_MATCHING.out.all_extraction_stats
        )
    }



    workflow.onComplete {
        def outputDir = file("${output_dir}/run")
        if (!outputDir.exists()) {
            outputDir.mkdirs()
        }
        record_NANOPORE_AMPLICON_CLUSTERING_params()
        record_NANOPORE_AMPLICON_CLUSTERING_runtime()
        completionSummary(false)
    }
}  //NANOPORE_AMPLICON_CLUSTERING



def record_NANOPORE_AMPLICON_CLUSTERING_params() {
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
def record_NANOPORE_AMPLICON_CLUSTERING_runtime() {

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
