process AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR {
    label 'process_low'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true

    cpus   { ncpus }

    input:
    path input_dir, stageAs: 'input_dir' // staged as a single link however many files it holds
    path input_files // tsv of a path relative to input_dir and the name to link it as, only the files the run processes (see lib/FlagGeneratorInputs.groovy)
    path primers_fnp
    val technology
    val ncpus
    val pub_dir
    val input_file_pat // regex for the input files to investigate, blank uses SeekDeep's default (fastqs)


    output:
    path "outSeekDeepExtractorFlags.txt", emit: out_seekdeep_extractor_flags

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def input_file_pat_arg = input_file_pat ? "--inputFilePat '${input_file_pat}'" : ""

    // SeekDeep looks at every file in --reads_dir and fails on ones it can't pair or name, so it gets a directory of
    // links to only the listed files; illumina mates have to share a directory for SeekDeep to pair them
    """
    real_input_dir=\$(readlink -f input_dir)
    mkdir reads_dir
    while IFS=\$'\\t' read -r rel_path link_name; do
        [ -n "\$rel_path" ] || continue
        ln -s "\$real_input_dir/\$rel_path" "reads_dir/\$link_name"
    done < ${input_files}

    SeekDeep gatherInfoOnTargetedAmpliconSeqFile \
            --id ${primers_fnp} \
            --reads_dir reads_dir \
            ${input_file_pat_arg} \
            --dout info_dir \
            --technology ${technology} \
            --numThreads ${ncpus}\
            --overWriteDir\
            ${extra_args}
    ln -s info_dir/outSeekDeepExtractorFlags.txt outSeekDeepExtractorFlags.txt
    """
}
