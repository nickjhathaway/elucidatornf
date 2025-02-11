process INDEX_GENOMES_DIR {
    label 'process_low'
    cpus   { ncpus }

    input:
    path genome_dir
    val ncpus
    script:

    """
    #unzip any genomes that are currently zipped, keep the zipped in case it was needed
    find ${genome_dir}/ -type f -name "*.gz" -exec gunzip -k {} ";"

    elucidator bioIndexGenomes --genomeDir ${genome_dir} --numThreads ${ncpus}
    """
}
