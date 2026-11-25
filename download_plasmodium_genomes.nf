#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SETUP_PLASMODIUM_GENOMES } from './subworkflows/local/plasmodium_genomes.nf'

/*
 * Download, gunzip and index https://github.com/nickjhathaway/plasmodium_genomes into --plasmodium_genomes_dir
 * (or <outdir>/plasmodium_genomes), skipping whatever is already done. The result can be given to the other
 * workflows with --download_plasmodium_genomes --plasmodium_genomes_dir <dir>, or as --genome_dir <dir>/genomes
 * --gff_dir <dir>/info/gff (pathweaver: --genome_fnp <dir>/genomes/<genome>.fasta)
 */
workflow {

main:
    def dir = params.plasmodium_genomes_dir?.toString()?.trim()
    def outdir = params.outdir?.toString()?.trim()
    if (!dir && !outdir) {
        error "flag '--plasmodium_genomes_dir' (or '--outdir', to use <outdir>/plasmodium_genomes) must be specified!"
    }
    def genomes_root = PlasmodiumGenomes.rootDir(params)
    SETUP_PLASMODIUM_GENOMES(genomes_root, params.plasmodium_genomes_version)
    SETUP_PLASMODIUM_GENOMES.out.ready.subscribe { _ready ->
        log.info "Plasmodium genomes ready in ${genomes_root}: ${PlasmodiumGenomes.genomeNames(genomes_root).join(", ")}"
    }
}
