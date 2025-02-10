process CREATE_TWOBIT_FILE_FOR_GENOME {
    label 'process_single'

    input:
    path genome_dir
    val genome_fname
    val genome_basename

    output:
    path "${genome_dir}/${genome_basename}.2bit", emit: twobit_fnp
    script:

    def twobit_fnp = file("${genome_dir}/${genome_fname}").baseName + ".2bit"

    """
    elucidator faToTwoBit --in ${genome_dir}/${genome_fname} --out ${genome_dir}/${twobit_fnp} --overWrite
    """
}
