#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ILLUMINA_AMPLICON_CLUSTERING } from './subworkflows/local/SeekDeep/SeekDeep_Illumina_Clustering.nf'

workflow {

main:

    //
    // Required params (for presence / non-blank checks only)
    //
    def required_params = [
        input_fastq_dir   : '--input_fastq_dir',
        primers_fnp       : '--primers_fnp',
        genome_dir        : '--genome_dir',
        gff_dir           : '--gff_dir',
        outdir            : '--outdir',
        illumina_clustering_paired_end_length : '--illumina_clustering_paired_end_length',
    ]

    def missing = required_params.findAll { key, _flag ->
        if (params[key] == null) return true
        def s = params[key].toString().trim()
        return s == ""
    }

    if (!missing.isEmpty()) {
        def missing_list = missing.collect {it -> it.value }.join('\n  ')
        error """
Input validation failed:

  Missing required parameters or blank values:
  ${missing_list}

Please provide non-empty values and re-run.
"""
    }

    //
    // Hand off to subworkflow — it can do the file()/existence checks
    //
    ILLUMINA_AMPLICON_CLUSTERING(
        params.input_fastq_dir,
        params.primers_fnp,
        params.genome_dir,
        params.gff_dir,
        params.nanopore_primers_errors_allowed,
        params.illumina_clustering_paired_end_length,
        params.meta_fnp,
        params.outdir
    )
}
