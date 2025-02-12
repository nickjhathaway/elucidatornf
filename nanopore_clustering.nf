#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {NANOPORE_AMPLICON_CLUSTERING} from './subworkflows/local/SeekDeep/SeekDeep_Nanopore_Clustering.nf'


workflow {

main:
    // Validate inputs
    if (params.input_fastq_dir == null  || params.primers_fnp == null || params.genome_dir == null || params.gff_dir == null || params.outdir == null) {
        error "flags '--input_fastq_dir', '--primers_fnp', '--genome_dir', '--outdir', and '--gff_dir' must be specified!"
    }
    // input_fastq_dir //a directory with all fastqs to be analyzed
    // primers_fnp  // a file with primer info, 3 columsn, target,fwd_primer,rev_primer
    // genome_dir  // a directory with genomes
    // gff_dir // a directory with annotation gff files
    // primers_errors_allowed // errors to allow in primers
    // meta_fnp //a sample meta file

    NANOPORE_AMPLICON_CLUSTERING(file(params.input_fastq_dir), file(params.primers_fnp), file(params.genome_dir), file(params.gff_dir), params.nanopore_primers_errors_allowed, file(params.meta_fnp), file(params.outdir))
}

