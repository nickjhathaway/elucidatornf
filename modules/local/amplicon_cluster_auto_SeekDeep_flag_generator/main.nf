process AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR {
    label 'process_low'

    publishDir "${pub_dir}", mode: 'copy', overwrite: true

    cpus   { ncpus }

    input:
    path fastq_dir_fnp
    path primers_fnp
    val technology
    val ncpus
    val pub_dir


    output:
    path "outSeekDeepExtractorFlags.txt", emit: out_seekdeep_extractor_flags

    script:
    def extra_args = task.ext.args ? task.ext.args : ''

    """
    SeekDeep gatherInfoOnTargetedAmpliconSeqFile \
            --id ${primers_fnp} \
            --reads_dir ${fastq_dir_fnp} \
            --dout info_dir \
            --technology ${technology} \
            --numThreads ${ncpus}\
            --overWriteDir\
            ${extra_args}
    ln -s info_dir/outSeekDeepExtractorFlags.txt outSeekDeepExtractorFlags.txt
    """
}
