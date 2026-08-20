#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
    Host read depletion only.

    Runs the host-filter chain from wgs_preprocess_pf.nf (minimap2 + bwa against each of
    two host genomes) without fastp trimming and without the final Pf mapping, so the
    output is the input data with host reads removed and nothing else changed. Any read
    whose mate got filtered off is discarded so the output stays strictly paired end.
*/

include { RUN_WGS_HOST_DEPLETION_PF } from './subworkflows/local/WGS/run_wgs_host_depletion_pf.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null ||
        params.outdir == null ||
        params.wgs_host1_filter_genome_fasta_fnp == null ||
        params.wgs_host1_filter_contigs_fnp == null ||
        params.wgs_host2_filter_genome_fasta_fnp == null ||
        params.wgs_host2_filter_contigs_fnp == null) {
        error "flags '--input_fastq_dir', '--outdir', '--wgs_host1_filter_genome_fasta_fnp', '--wgs_host1_filter_contigs_fnp', '--wgs_host2_filter_genome_fasta_fnp', and '--wgs_host2_filter_contigs_fnp' must be specified!"
    }

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (!results_dir_obj.exists()){
        results_dir_obj.mkdirs()
    }

    RUN_WGS_HOST_DEPLETION_PF(file(params.input_fastq_dir),
     file(params.outdir),
     file(params.wgs_host1_filter_genome_fasta_fnp),
     file(params.wgs_host1_filter_contigs_fnp),
     file(params.wgs_host2_filter_genome_fasta_fnp),
     file(params.wgs_host2_filter_contigs_fnp)
     )
}
