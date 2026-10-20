#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {RUN_RNASEQ_PROCESS_PF} from './subworkflows/local/rna_seq/run_rna_seq_pf_process.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null || params.outdir == null) {
        error "flags '--input_fastq_dir' and '--outdir' must be specified!"
    }

    // Reference files: each from its own flag if set, otherwise its expected name under
    // --rnaseq_ref_files_dir (see lib/RnaseqRefFiles.groovy). Everything is checked up
    // front so a wrong or missing file fails here rather than hours into the run.
    def refs = RnaseqRefFiles.resolve(params) { p -> file(p) }
    def ref_errors = refs.errors + RnaseqRefFiles.validate(refs.files, refs.sources)
    if (ref_errors) {
        error "Problems with the RNA-seq reference files:\n  - " + ref_errors.join("\n  - ") +
            "\nSee the README in the reference directory for the expected names, or set the individual --rnaseq_*/--single_cell_* flags."
    }
    log.info "RNA-seq reference files:\n" + RnaseqRefFiles.describe(refs.files, refs.sources)

    // Create output directory if not exists and overwrite if it does
    def results_dir_obj = file(params.outdir)
    if (!results_dir_obj.exists()){
        results_dir_obj.mkdirs()
    }

    RUN_RNASEQ_PROCESS_PF(file(params.input_fastq_dir),
                        file(params.outdir),
                        refs.files.rnaseq_host1_filter_genome_fasta_fnp,
                        refs.files.rnaseq_host1_filter_contigs_fnp,
                        refs.files.rnaseq_host2_filter_genome_fasta_fnp,
                        refs.files.rnaseq_host2_filter_contigs_fnp,
                        refs.files.rnaseq_salmon_index,
                        refs.files.single_cell_seurat_r_object_fnp,
                        refs.files.rnaseq_kmer_filter_table_fnp,
                        refs.files.rnaseq_kmer_sort_table_fnp
     )
}

