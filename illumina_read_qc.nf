#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
    Illumina amplicon read QC only.

    Puts every read pair of each sample into one category (on target, primer dimer, short product,
    mismatched primers etc) and works out where the products that weren't expected came from, without
    any of the extraction or clustering. Doesn't delete anything already in --outdir so it can be
    pointed at the output of illumina_clustering.nf.
*/

include { ILLUMINA_READ_QC } from './subworkflows/local/SeekDeep/SeekDeep_Illumina_Read_QC.nf'


workflow {

main:
    def required_params = [
        input_fastq_dir : '--input_fastq_dir',
        primers_fnp     : '--primers_fnp',
        genome_dir      : '--genome_dir',
        outdir          : '--outdir',
        illumina_read_qc_primary_genome : '--illumina_read_qc_primary_genome',
    ]

    def missing = required_params.findAll { key, _flag ->
        if (params[key] == null) return true
        def s = params[key].toString().trim()
        return s == ""
    }

    if (!missing.isEmpty()) {
        def missing_list = missing.collect { it -> it.value }.join('\n  ')
        error """
Input validation failed:

  Missing required parameters or blank values:
  ${missing_list}

Please provide non-empty values and re-run.
"""
    }

    ILLUMINA_READ_QC(
        params.input_fastq_dir,
        params.primers_fnp,
        params.genome_dir,
        params.gff_dir,
        params.illumina_read_qc_target_info_dir,
        params.nanopore_primers_errors_allowed,
        params.illumina_clustering_paired_end_length,
        params.outdir
    )
}
