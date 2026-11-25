#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SETUP_PLASMODIUM_GENOMES } from './subworkflows/local/plasmodium_genomes.nf'
include { ILLUMINA_AMPLICON_CLUSTERING } from './subworkflows/local/SeekDeep/SeekDeep_Illumina_Clustering.nf'

workflow {

main:

    //
    // Required params (for presence / non-blank checks only)
    //
    def required_params = [
        input_fastq_dir   : '--input_fastq_dir',
        primers_fnp       : '--primers_fnp',
        outdir            : '--outdir',
        illumina_clustering_paired_end_length : '--illumina_clustering_paired_end_length',
    ]

    def download_genomes = PlasmodiumGenomes.enabled(params)
    if (!download_genomes) {
        required_params.genome_dir = '--genome_dir'
        required_params.gff_dir = '--gff_dir'
    }

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
    // --download_plasmodium_genomes points --genome_dir and --gff_dir at the downloaded genomes
    //
    def genome_dir = params.genome_dir
    def gff_dir = params.gff_dir
    def genomes_ready = channel.value(true)
    if (download_genomes) {
        def genomes_root = PlasmodiumGenomes.rootDir(params)
        if (params.genome_dir || params.gff_dir) {
            log.warn "--genome_dir and --gff_dir are ignored with --download_plasmodium_genomes, using ${genomes_root}"
        }
        SETUP_PLASMODIUM_GENOMES(genomes_root, params.plasmodium_genomes_version)
        genome_dir = PlasmodiumGenomes.genomeDir(genomes_root).toString()
        gff_dir = PlasmodiumGenomes.gffDir(genomes_root).toString()
        genomes_ready = SETUP_PLASMODIUM_GENOMES.out.ready
    }

    //
    // Hand off to subworkflow — it can do the file()/existence checks
    //
    ILLUMINA_AMPLICON_CLUSTERING(
        params.input_fastq_dir,
        params.primers_fnp,
        genome_dir,
        gff_dir,
        genomes_ready,
        params.nanopore_primers_errors_allowed,
        params.illumina_clustering_paired_end_length,
        params.meta_fnp,
        params.outdir
    )
}
