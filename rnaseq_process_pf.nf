#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {RUN_RNASEQ_PROCESS_PF} from './subworkflows/local/rna_seq/run_rna_seq_pf_process.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null ||
        params.outdir == null ||
        params.rnaseq_host1_filter_genome_fasta_fnp == null ||
        params.rnaseq_host1_filter_contigs_fnp == null ||
        params.rnaseq_host2_filter_genome_fasta_fnp == null ||
        params.rnaseq_host2_filter_contigs_fnp == null||
        params.rnaseq_salmon_index == null) {
        error "flags '--input_fastq_dir', '--outdir', '--rnaseq_host1_filter_genome_fasta_fnp', '--rnaseq_host1_filter_contigs_fnp', '--rnaseq_host2_filter_genome_fasta_fnp', '--rnaseq_host2_filter_contigs_fnp', and '--rnaseq_salmon_index' must be specified!"
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (!results_dir_obj.exists()){
        results_dir_obj.mkdirs()
    }

    RUN_RNASEQ_PROCESS_PF(file(params.input_fastq_dir),
                        file(params.outdir),
                        file(params.rnaseq_host1_filter_genome_fasta_fnp),
                        file(params.rnaseq_host1_filter_contigs_fnp),
                        file(params.rnaseq_host2_filter_genome_fasta_fnp),
                        file(params.rnaseq_host2_filter_contigs_fnp),
                        file(params.rnaseq_salmon_index)
     )
}

