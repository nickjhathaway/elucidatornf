#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// include { paramsSummaryMap       } from 'plugin/nf-schema'


include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE as GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL } from '../../../modules/local/gen_target_info_from_genomes'
include { GEN_TARGET_INFO_FROM_GENOMES_NANOPORE } from '../../../modules/local/gen_target_info_from_genomes'
include { EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/extractor_by_kmer_matching'
include { COMBINE_VARIANT_CALLS } from '../../../modules/local/combine_variant_calls/main.nf'
include { popClusteringVariantCallArgs; noPopClusteringVariantCalls; pmoPreferredGenomes; variantCallGenomes; readTargetGenomes; variantCallGenomeFor } from '../variant_call_setup.nf'
include { ASSIGN_TARGET_GENOMES } from '../../../modules/local/assign_target_genomes/main.nf'
include { CONCATENATE_EXTRACTOR_BY_KMER_MATCHING } from '../../../modules/local/concatenate_extractor_by_kmer_matching/main.nf'
include { AMPLICON_CLUSTER_BY_KMER_SIMILARITY } from '../../../modules/local/amplicon_cluster_by_kmer_similarity/main.nf'
include { AMPLICON_POPULATION_CLUSTERING } from '../../../modules/local/amplicon_population_cluster/main.nf'
include { CONCATENATE_AMPLICON_POPULATION_CLUSTERING } from '../../../modules/local/concatenate_amplicon_population_clustering/main.nf'
include { MERGE_AMPLICON_PMOS } from '../../../modules/local/merge_amplicon_pmos/main.nf'
include { AMPLICON_CREATE_DASHBOARDS } from "../../../modules/local/amplicon_create_dashboards/main.nf"
include { WRITE_AMPLICON_SAMPLE_KEY } from "../../../modules/local/write_amplicon_sample_key/main.nf"
include { SCAN_RAW_INPUT_SAMPLE_NAMES } from "../../../modules/local/scan_raw_input_sample_names/main.nf"

include { completionSummary } from '../../../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'


workflow NANOPORE_AMPLICON_CLUSTERING {
    take:
    input_fastq_dir_raw //a directory with all fastqs (.fastq(.gz)) and/or unaligned bams (.bam, e.g. from dorado demux) to be analyzed, searched recursively with params.nanopore_input_recursive
    primers_fnp_raw  // a file with primer info, 3 columsn, target,fwd_primer,rev_primer
    genome_dir_raw  // a directory with genomes
    gff_dir_raw // a directory with annotation gff files
    genomes_ready // value channel, the genome reading waits on it (for --download_plasmodium_genomes indexing)
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
    // primers of regions within the targets to trim the reads to, primers_fnp are then the amplified primers (only used
    // to work out which target a read is) and everything from the extracted reads on is by region
    def trim_to_regions = params.nanopore_trim_primers_fnp ? true : false
    def trim_primers_fnp = trim_to_regions ? file(params.nanopore_trim_primers_fnp) : null
    // with --nanopore_extractor_use_inner_primers the reads are trimmed the same way, to inner primers made from the ends of each target's sequences
    // (given on the command line it's the string "false", which groovy would take as true)
    def use_inner_primers = params.nanopore_extractor_use_inner_primers.toString().toBoolean()
    def trim_reads = trim_to_regions || use_inner_primers
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
    // 2b) Variant calling requires a primary genome. If it is left unset the
    // default value (false) is interpolated into the file names, so the
    // workflow tries to load "false.gff" / "false_inner.bed" and fails with a
    // confusing missing-file error far downstream.
    //
    if (params.do_variant_calling) {
        def vc_genome = (params.vc_primary_genome == null) ? "" : params.vc_primary_genome.toString().trim()
        if (vc_genome == "" || vc_genome.toLowerCase() == "false") {
            errors << "--vc_primary_genome must be supplied when --do_variant_calling is set"
        } else {
            def gff_obj = file("${gff_dir}/${vc_genome}.gff")
            if (!gff_obj.exists()) {
                errors << "Variant calling GFF not found for --vc_primary_genome '${vc_genome}': '${gff_obj}' (expected in --gff_dir)"
            }
        }
        variantCallGenomes().each { genome ->
            ["${genome_dir_raw}/${genome}.fasta", "${gff_dir}/${genome}.gff"].each { fnp ->
                if (!file(fnp).exists()) {
                    errors << "Variant calling genome file not found for '${genome}': '${fnp}' (--vc_primary_genome and --vc_other_genomes need a fasta in --genome_dir and a gff in --gff_dir)"
                }
            }
        }
    }

    if (trim_to_regions) {
        if (!trim_primers_fnp.exists()) {
            errors << "Path does not exist for --nanopore_trim_primers_fnp: '${params.nanopore_trim_primers_fnp}'"
        }
        if (use_inner_primers) {
            errors << "--nanopore_trim_primers_fnp and --nanopore_extractor_use_inner_primers can't be used together"
        }
    }

    errors += SampleMeta.validate(meta_fnp)

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
    // Create output directory (overwrite if exists), keeping genomes downloaded into it
    //
    def results_dir_obj = file(outdir)
    OutputDir.clear(results_dir_obj, PlasmodiumGenomes.enabled(params) ? PlasmodiumGenomes.rootDir(params) : null, params.trace_timestamp, [workflow.launchDir, workflow.workDir])
    results_dir_obj.mkdirs()

    //
    // Inner primers logic
    //
    if (use_inner_primers) {
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

    // optional sample key, input name -> rows of sample, replicate and barcodes (see lib/SampleKey.groovy),
    // narrowed to the selected samples when there's a samples keep file (--amplicon_samples_keep_file)
    def full_sample_key = SampleKey.load(SampleKey.keyFnp(params))
    def samples_keep = SampleSelection.load(params.amplicon_samples_keep_file)
    def sample_key = SampleSelection.filterKey(samples_keep, full_sample_key)
    def keyed_names = SampleKey.keyedNames(full_sample_key)
    // barcoded inputs get their own id file with only their MIDs, named by replicate so the extractor
    // writes <target><replicate>.fastq.gz for them
    def mid_id_dir = file("${meta_info_dir}/midIdFiles")
    def barcoded_inputs = sample_key.findAll { _input, rows -> rows[0].fbarcode }
    barcoded_inputs.each { input, rows ->
        SampleKey.writeIdFile(primers_fnp, rows, "${mid_id_dir}/${input}_ids.tsv")
    }
    if (trim_reads && !barcoded_inputs.isEmpty()) {
        error "${use_inner_primers ? '--nanopore_extractor_use_inner_primers' : '--nanopore_trim_primers_fnp'} can't be used with barcodes in the sample key, they are next to the amplified primers which aren't looked for when trimming to regions"
    }
    SampleKey.writeCombinedIdFile(primers_fnp, sample_key, "${mid_id_dir}/allMids_ids.tsv")

    primers_fnp_ch = genomes_ready.map { _ready -> primers_fnp }

    // when trimming, the amplified targets' k-mers and lengths, for working out which target a read is, and the primers to trim to,
    // either the given ones or the inner primers made from the amplified targets' sequences
    def trim_primers_ch = channel.empty()
    if (trim_reads) {
        def amplified_primer_info_dir = file("${output_dir}/primerInfo/amplifiedPrimerInfo")
        amplified_primer_info_dir.mkdirs()
        GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL(primers_fnp_ch, amplified_primer_info_dir.toString(), genome_dir, gff_dir, primers_errors_allowed, params.max_cpus)
        trim_primers_ch = use_inner_primers
            ? GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL.out.inner_primers
            : genomes_ready.map { _ready -> trim_primers_fnp }
        primers_fnp_ch = trim_primers_ch
    }
    // the primers the results are for
    def result_primers_ch = trim_reads ? trim_primers_ch : genomes_ready.map { _ready -> primers_fnp }
    GEN_TARGET_INFO_FROM_GENOMES_NANOPORE(primers_fnp_ch, primer_info_dir.toString(), genome_dir, gff_dir, primers_errors_allowed, params.max_cpus)

    // given on the command line these are the string "false", which groovy would take as true
    def input_recursive = params.nanopore_input_recursive.toString().toBoolean()
    def prefix_with_top_dir = params.nanopore_sample_name_prefix_top_dir.toString().toBoolean()
    def bam_sample_name_from_sm = params.nanopore_bam_sample_name_from_sm.toString().toBoolean()

    // name every input file up front (only bam headers are read) so sample name collisions,
    // e.g. the same barcode from different runs, stop the run before the long extraction step
    SCAN_RAW_INPUT_SAMPLE_NAMES(input_fastq_dir, "nanopore", input_recursive, prefix_with_top_dir, bam_sample_name_from_sm, meta_info_dir.toString())

    def scanned_samples_ch = SCAN_RAW_INPUT_SAMPLE_NAMES.out.sample_names
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.status == "ok" }

    // sample key inputs are matched against the scanned names, which carry the top dir prefix when it's on,
    // so warn about entries that won't match anything rather than silently leaving those samples as is
    if (!sample_key.isEmpty()) {
        scanned_samples_ch
            .toList()
            .subscribe { rows ->
                def scanned_names = rows.collect { row -> row.sample } as Set
                def unmatched = sample_key.keySet().findAll { input -> !scanned_names.contains(input) }.sort()
                if (!unmatched.isEmpty()) {
                    def base_name_hits = unmatched.findAll { input -> rows.any { row -> row.baseName == input } }
                    def hint = base_name_hits.isEmpty() ? "" :
                        ". ${base_name_hits.join(", ")} match sample names before the top directory prefix, " +
                        "the sample key needs the prefixed names (see ${meta_info_dir}/rawInputSampleNames.tsv)"
                    log.warn "sample key inputs that don't match any input sample and won't be used: ${unmatched.join(", ")}${hint}"
                }
            }
    }

    if (samples_keep != null) {
        scanned_samples_ch
            .toList()
            .subscribe { rows ->
                def unmatched = SampleSelection.unmatched(samples_keep, full_sample_key, rows.collect { row -> row.sample })
                if (!unmatched.isEmpty()) {
                    def base_name_hits = unmatched.findAll { name -> rows.any { row -> row.baseName == name } }
                    def hint = base_name_hits.isEmpty() ? "" :
                        ". ${base_name_hits.join(", ")} match sample names before the top directory prefix, " +
                        "the keep file needs the prefixed names (see ${meta_info_dir}/rawInputSampleNames.tsv)"
                    log.warn "--amplicon_samples_keep_file names that don't match any input sample and won't be processed: ${unmatched.join(", ")}${hint}"
                }
            }
    }

    // file paths in the scan are relative to the input directory, ignored files (e.g. unclassified) are skipped
    def selected_samples_ch = scanned_samples_ch
        .filter { row -> SampleSelection.selectsInput(samples_keep, full_sample_key, row.sample) }
        .toList()
        .flatMap { rows ->
            // population clustering fails on samples that aren't in the meta, so nothing is extracted until they're checked
            def sample_names = rows.collectMany { row -> SampleKey.rowsForInput(sample_key, row.sample).collect { key_row -> key_row.sample } }
            def missing_meta = SampleMeta.missingMessage(SampleMeta.load(meta_fnp), sample_names)
            if (missing_meta) {
                error "${missing_meta}"
            }
            rows
        }
    def reads_input_ch = selected_samples_ch
        .map { row -> tuple(file("${input_fastq_dir}/${row.file}"), row.sample) }
    if (samples_keep != null) {
        reads_input_ch = reads_input_ch.ifEmpty { error "None of the samples in --amplicon_samples_keep_file match an input sample in --input_fastq_dir (see ${meta_info_dir}/rawInputSampleNames.tsv)" }
    }

    // the info on the targets the reads are assigned to and, when trimming to regions, on the regions
    def extractor_info_ch = trim_reads
        ? GEN_TARGET_INFO_FROM_GENOMES_NANOPORE_INITIAL.out.for_seek_deep_info
            .combine(GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info)
            .combine(trim_primers_ch)
        : GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info.map { info_dir -> tuple(info_dir, [], []) }

    input_to_extractor = reads_input_ch.combine(extractor_info_ch)
            .map{reads_fnp, sample_name, info_dir, trim_info_dir, trim_primers ->
                // an input left out of the key is a sample named after itself, which can't also be a keyed name
                if (!sample_key.containsKey(sample_name) && keyed_names.contains(sample_name)) {
                    error "input ${sample_name} isn't in the sample key but is the name of a sample or replicate in it, add it to the key"
                }
                tuple(
                    reads_fnp,
                    sample_name,
                    barcoded_inputs.containsKey(sample_name) ? file("${mid_id_dir}/${sample_name}_ids.tsv") : primers_fnp,
                    file("${info_dir}/allKmers.tab.txt.gz"),
                    file("${info_dir}/lenCutOffs.txt"),
                    trim_reads ? trim_primers : [],
                    trim_reads ? file("${trim_info_dir}/allKmers.tab.txt.gz") : [],
                    trim_reads ? file("${trim_info_dir}/lenCutOffs.txt") : [],
                    params.nanopore_clustering_min_len
                )
            }

    //extract by kmer sets
    EXTRACTOR_BY_KMER_MATCHING(input_to_extractor)

    //concatenate extraction files
    CONCATENATE_EXTRACTOR_BY_KMER_MATCHING(
        EXTRACTOR_BY_KMER_MATCHING.out.extraction_profile_per_target | collect,
        EXTRACTOR_BY_KMER_MATCHING.out.extraction_stats | collect,
        EXTRACTOR_BY_KMER_MATCHING.out.extraction_profile_per_mid.collect().ifEmpty([]),
        extraction_reports_dir.toString()
    )

    //get the targets that have coverage
    def targets = GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.targets_with_extractions
        .splitText()
        .map{samp ->
            samp.trim() // Remove any whitespace
            }


    // create a list of input to cluster for the input that actually exists. An input holds a single sample,
    // or several MID tagged ones when it has barcodes in the sample key. Each sample/replicate is clustered
    // on its own and population clustering combines the replicates of a sample
    def clustered_inputs_with_key = targets.combine(EXTRACTOR_BY_KMER_MATCHING.out.sample_dir)
        .flatMap { tar, input, sampdir ->
            SampleKey.rowsForInput(sample_key, input).findResults { row ->
                def extracted_tar_sample_fnp = file("${sampdir}/${tar}${row.extraction_name}.fastq.gz")
                if( !extracted_tar_sample_fnp.exists() ) return null

                // emit BOTH: the clustering tuple + the sample key line
                tuple(
                    tuple(extracted_tar_sample_fnp, row.sample, row.replicate, tar, params.nanopore_clustering_ncpus),
                    SampleKey.keyLine(row)
                )
            }
        }

    // split into the two channels we want
    def input_to_clustering = clustered_inputs_with_key.map { entry -> entry[0] }

    // write the sample key for the samples that had extracted reads
    def sample_key_lines_ch = clustered_inputs_with_key.map { entry -> entry[1] }
        .unique()
        .collect()

    WRITE_AMPLICON_SAMPLE_KEY(sample_key_lines_ch, meta_info_dir.toString())
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
        .combine(result_primers_ch)
        .map { _samples, target, target_fastqs, ref_seqs_dir, result_primers_fnp ->
                tuple(
                    target_fastqs, // List of fastq files for this target
                    target,
                    params.nanopore_clustering_ncpus,
                    meta_fnp,
                    // GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.for_seek_deep_info,
                    ref_seqs_dir,
                    params.nanopore_clustering_min_sample_read_count,
                    result_primers_fnp,
                    params.pmo_panel_name
                )
        }
    // the genome each target is variant called against, the first of vc_primary_genome and vc_other_genomes it's in
    def target_genomes_ch = channel.value([:])
    def assigned_genomes_ch = channel.value(file(params.empty_file_fnp))
    if (params.do_variant_calling) {
        ASSIGN_TARGET_GENOMES(GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.locations_by_genome, variantCallGenomes().join(","), primer_info_dir.toString())
        assigned_genomes_ch = ASSIGN_TARGET_GENOMES.out.assigned
        target_genomes_ch = assigned_genomes_ch.map { assigned_dir -> readTargetGenomes(assigned_dir) }
    }
    def input_to_population_clustering_with_genome = input_to_population_clustering
        .combine(target_genomes_ch)
        .map { entry ->
            def target_genomes = entry[-1]
            def pop_input = entry[0..-2]
            pop_input + [variantCallGenomeFor(target_genomes, pop_input[1])]
        }
    // population clustering
    // population clustering calls variants on the final haplotypes when variant calling is on
    def vc_known_amino_acid_changes_fnp = params.empty_file_fnp
    if(file("${genome_dir.getParent()}/info/drug_resistant_aaPositions.tsv").exists()){
        vc_known_amino_acid_changes_fnp = file("${genome_dir.getParent()}/info/drug_resistant_aaPositions.tsv")
    }
    def vc_setup_ch = noPopClusteringVariantCalls()
    if (params.do_variant_calling) {
        def vc_args = popClusteringVariantCallArgs(!params.vc_no_pairwise_comps)
        vc_setup_ch = assigned_genomes_ch.map{ assigned_dir ->
            tuple(genome_dir, params.vc_primary_genome, gff_dir, file(vc_known_amino_acid_changes_fnp), assigned_dir, vc_args)}
    }
    AMPLICON_POPULATION_CLUSTERING(input_to_population_clustering_with_genome, vc_setup_ch)

    //gather all clustered results
    CONCATENATE_AMPLICON_POPULATION_CLUSTERING(
        AMPLICON_POPULATION_CLUSTERING.out.clustering_results | collect,
        AMPLICON_POPULATION_CLUSTERING.out.sample_clusters_info | collect,
        final_results_dir.toString())

    //one PMO for all the targets with each target's genome location and genes, and the sample meta
    MERGE_AMPLICON_PMOS(
        AMPLICON_POPULATION_CLUSTERING.out.pmo | collect,
        meta_fnp,
        GEN_TARGET_INFO_FROM_GENOMES_NANOPORE.out.locations_by_genome,
        genome_dir,
        gff_dir,
        pmoPreferredGenomes(),
        final_results_dir.toString())

    if (params.do_variant_calling){

        def variant_call_dir = file("${final_results_dir}/variantCalls")

        //the reports across targets from the calls population clustering made, one set per genome, the primary genome's in variantCalls
        def combine_input = AMPLICON_POPULATION_CLUSTERING.out.variant_calls
            .groupTuple()
            .combine(assigned_genomes_ch)
            .multiMap { genome, target_calls, assigned_dir ->
                bed: file("${assigned_dir}/${genome}_variantCall.bed")
                genome: genome
                gff: file("${gff_dir}/${genome}.gff")
                known: genome == params.vc_primary_genome ? file(vc_known_amino_acid_changes_fnp) : file(params.empty_file_fnp)
                calls: target_calls
                pub_dir: genome == params.vc_primary_genome ? variant_call_dir.toString() : "${variant_call_dir}/otherGenomes/${genome}".toString()
            }
        COMBINE_VARIANT_CALLS (
            combine_input.bed,
            CONCATENATE_AMPLICON_POPULATION_CLUSTERING.out.all_selected_clusters_info,
            genome_dir,
            combine_input.genome,
            combine_input.gff,
            combine_input.known,
            params.vc_variant_frequency_cut_off,
            params.vc_variant_occurrence_cut_off,
            meta_fnp,
            combine_input.calls,
            "",
            "",
            combine_input.pub_dir,
            params.vc_extra_args
        )
        AMPLICON_CREATE_DASHBOARDS(file("${output_dir}").toString(),
            file("${projectDir}/etc/nanopore_clustering_report.qmd"),
            params.nanopore_render_clustering_report,
            COMBINE_VARIANT_CALLS.out.reports.collect().map { reports -> reports[0] },
            CONCATENATE_EXTRACTOR_BY_KMER_MATCHING.out.all_extraction_stats
        )
    } else {
        AMPLICON_CREATE_DASHBOARDS(file("${output_dir}").toString(),
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
