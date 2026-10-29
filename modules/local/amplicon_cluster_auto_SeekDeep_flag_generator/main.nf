process AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR {
    label 'process_low'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true

    cpus   { ncpus }

    input:
    path fastq_dir_fnp
    path primers_fnp
val technology
    val ncpus
    val pub_dir
    val input_file_pat // regex for the files in fastq_dir_fnp to investigate, blank uses SeekDeep's default (fastqs)
    val recursive // search fastq_dir_fnp recursively


    output:
    path "outSeekDeepExtractorFlags.txt", emit: out_seekdeep_extractor_flags

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def input_file_pat_arg = input_file_pat ? "--inputFilePat '${input_file_pat}'" : ""
    def recursive_arg = recursive ? "--recursive" : ""

    """
    SeekDeep gatherInfoOnTargetedAmpliconSeqFile \
            --id ${primers_fnp} \
            --reads_dir ${fastq_dir_fnp} \
            ${input_file_pat_arg} \
            ${recursive_arg} \
            --dout info_dir \
            --technology ${technology} \
            --numThreads ${ncpus}\
            --overWriteDir\
            ${extra_args}
    ln -s info_dir/outSeekDeepExtractorFlags.txt outSeekDeepExtractorFlags.txt
    """
}
