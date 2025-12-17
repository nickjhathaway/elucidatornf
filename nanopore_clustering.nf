#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {NANOPORE_AMPLICON_CLUSTERING} from './subworkflows/local/SeekDeep/SeekDeep_Nanopore_Clustering.nf'


workflow {

main:
    // Validate inputs
    //
    // Required param definitions
    //
    def required_params = [
        input_fastq_dir : '--input_fastq_dir',
        primers_fnp     : '--primers_fnp',
        genome_dir      : '--genome_dir',
        gff_dir         : '--gff_dir',
        outdir          : '--outdir'
    ]

    //
    // Accumulate errors
    //
    def errors = []

    //
    // Check for null or empty-string inputs
    //
    required_params.each { key, flag ->
        def v = params[key]
        if (v == null || (v instanceof String && v.trim() == "")) {
            errors << "Missing required parameter: ${flag}"
        }
    }

    //
    // Existence checks (only if param is non-empty)
    //
    def path_checks = [
        input_fastq_dir : "directory",
        primers_fnp     : "file",
        genome_dir      : "directory",
        gff_dir         : "directory"
    ]

    path_checks.each { key, type ->
        def v = params[key]
        if (v != null && v.toString().trim() != "") {
            def obj = file(v)
            if (!obj.exists()) {
                errors << "Path does not exist for ${required_params[key]}: '${v}'"
            } else if (type == "directory" && !obj.isDirectory()) {
                errors << "Expected a directory for ${required_params[key]} but found a file: '${v}'"
            } else if (type == "file" && !obj.isFile()) {
                errors << "Expected a file for ${required_params[key]} but found a directory: '${v}'"
            }
        }
    }

    //
    // If any validation failed → print full error listing
    //
    if (!errors.isEmpty()) {
        def joined = errors.collect {it -> "  - ${it}" }.join("\n")
        error """
Input validation failed:

${joined}

Please correct the above issues and re-run.
"""
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (results_dir_obj.exists()){
        results_dir_obj.deleteDir()
    }
    results_dir_obj.mkdirs()

    if (params.nanopore_extractor_use_inner_primers){
        params.nanopore_clustering_lower_base = "upper"
    }
    NANOPORE_AMPLICON_CLUSTERING(file(params.input_fastq_dir), file(params.primers_fnp), file(params.genome_dir), file(params.gff_dir), params.nanopore_primers_errors_allowed, file(params.meta_fnp), file(params.outdir))
}

