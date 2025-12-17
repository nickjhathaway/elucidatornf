#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {NANOPORE_AMPLICON_CLUSTERING} from './subworkflows/local/SeekDeep/SeekDeep_Nanopore_Clustering.nf'


workflow {

main:
    // Validate inputs
    def required_params = [
        input_fastq_dir : '--input_fastq_dir',
        primers_fnp     : '--primers_fnp',
        genome_dir      : '--genome_dir',
        gff_dir         : '--gff_dir',
        outdir          : '--outdir'
    ]

    def missing = required_params.findAll { key, _flag -> params[key] == null }

    if (!missing.isEmpty()) {
        def missing_list = missing.collect {it -> it.value }.join('\n  ')
        error """
Missing required parameters:
  ${missing_list}
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

