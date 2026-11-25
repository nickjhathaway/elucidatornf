#!/usr/bin/env nextflow

include { INDEX_PLASMODIUM_GENOMES } from '../../modules/local/index_plasmodium_genomes/main.nf'

/*
 * --download_plasmodium_genomes: download https://github.com/nickjhathaway/plasmodium_genomes into
 * genomes_root (skipped when it's already extracted) and index it (skipped when it's already indexed).
 * The download happens while the workflow is being built so the files are there for the input checks;
 * gate the first process that reads the genomes on `ready` so it waits for the indexing.
 */
workflow SETUP_PLASMODIUM_GENOMES {
    take:
    genomes_root // directory to download into, see PlasmodiumGenomes.rootDir()
    version      // git tag/branch of the repo to download

    main:
    def root = genomes_root
    def indexed = PlasmodiumGenomes.isIndexed(root)
    if (indexed) {
        log.info "Using the already set up plasmodium genomes in ${root}"
    } else {
        try {
            PlasmodiumGenomes.download(root, version.toString(), log)
        } catch (Exception e) {
            error "--download_plasmodium_genomes failed: ${e.message}"
        }
        log.info "Indexing the plasmodium genomes in ${root}"
        INDEX_PLASMODIUM_GENOMES(file(root.toString()), params.max_cpus)
    }

    emit:
    // value channel, emits once the genomes are downloaded and indexed
    ready = indexed ? channel.value(true) : INDEX_PLASMODIUM_GENOMES.out.ready
}
