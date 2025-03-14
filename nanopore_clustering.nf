#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {NANOPORE_AMPLICON_CLUSTERING} from './subworkflows/local/SeekDeep/SeekDeep_Nanopore_Clustering.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null  || params.primers_fnp == null || params.genome_dir == null || params.gff_dir == null || params.outdir == null) {
        error "flags '--input_fastq_dir', '--primers_fnp', '--genome_dir', '--outdir', and '--gff_dir' must be specified!"
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

