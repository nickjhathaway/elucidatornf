process INDEX_GENOMES_DIR {
    label 'process_low'
    cpus   { ncpus }

    input:
    path genome_dir
    val ncpus
    script:

    """
    elucidator bioIndexGenomes --genomeDir ${genome_dir} --numThreads ${ncpus}
    """
}
