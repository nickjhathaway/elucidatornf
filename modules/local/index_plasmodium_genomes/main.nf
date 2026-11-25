process INDEX_PLASMODIUM_GENOMES {
    label 'process_medium_low_memory'
    cpus   { ncpus }

    // writes into the staged directory (indexes the genomes in place), so it's never cached
    cache false

    input:
    path genomes_root
    val ncpus

    output:
    val true, emit: ready

    script:
    """
    elucidator bioIndexGenomes --genomeDir ${genomes_root}/genomes --numThreads ${ncpus}
    """
}
