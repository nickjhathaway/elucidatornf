#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {RUN_WGS_PREPROCESS_PF} from './subworkflows/local/WGS/run_wgs_preprocess_pf.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null ||
        params.outdir == null ||
        params.wgs_host1_filter_genome_fasta_fnp == null ||
        params.wgs_host1_filter_contigs_fnp == null ||
        params.wgs_host2_filter_genome_fasta_fnp == null ||
        params.wgs_host2_filter_contigs_fnp == null||
        params.wgs_final_genome_fasta_fnp == null) {
        error "flags '--input_fastq_dir', '--outdir', '--wgs_host1_filter_genome_fasta_fnp', '--wgs_host1_filter_contigs_fnp', '--wgs_host2_filter_genome_fasta_fnp', '--wgs_host2_filter_contigs_fnp', and '--wgs_final_genome_fasta_fnp' must be specified!"
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (!results_dir_obj.exists()){
        results_dir_obj.mkdirs()
    }

    RUN_WGS_PREPROCESS_PF(file(params.input_fastq_dir),
     file(params.outdir),
     file(params.wgs_host1_filter_genome_fasta_fnp),
     file(params.wgs_host1_filter_contigs_fnp),
     file(params.wgs_host2_filter_genome_fasta_fnp),
     file(params.wgs_host2_filter_contigs_fnp),
     file(params.wgs_final_genome_fasta_fnp)
     )
}

