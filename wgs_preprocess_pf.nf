#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {RUN_WGS_PREPROCESS_PF} from './subworkflows/local/WGS/run_wgs_preprocess_pf.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null || params.outdir == null) {
        error "flags '--input_fastq_dir' and '--outdir' must be specified!"
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (!results_dir_obj.exists()){
        results_dir_obj.mkdirs()
    }

    RUN_WGS_PREPROCESS_PF(file(params.input_fastq_dir), file(params.outdir))
}

