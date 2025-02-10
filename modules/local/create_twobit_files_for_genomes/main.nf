process CREATE_TWOBIT_FILES_FOR_GENOMES {
    label 'process_single'

    input:
    path genome_dir
    script:

    """
    for x in ${genome_dir}/*.fa*; do elucidator faToTwoBit --in \${x} --out \${x%%.fa*}.2bit --overWrite; done;
    """
}
