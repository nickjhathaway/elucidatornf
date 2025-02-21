process AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR {
    label 'process_low'

    cpus   { ncpus }

    input:
    path fastq_dir_fnp
    path primers_fnp
    val tehcnology
    val ncpus


    output:
    path "outSeekDeepExtractorFlags.txt", emit: out_seekdeep_extractor_flags
    script:
    
    """
    SeekDeep gatherInfoOnTargetedAmpliconSeqFile \
            --id ${primers_fnp} \
            --reads_dir ${fastq_dir_fnp} \
            --dout info_dir \
            --technology ${tehcnology} \
            --numThreads ${ncpus}
            --overWriteDir
    ln -s info_dir/outSeekDeepExtractorFlags.txt outSeekDeepExtractorFlags.txt
    """
}
            //--testNumber 5000 \
            //--numberOfFilesToInvestigate 5