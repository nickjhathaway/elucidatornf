process EXTRACTOR_BY_KMER_MATCHING {
    tag "${sample_name}_extract_by_kmer_matching"
    label 'process_single'


    input:
    // reads_fnp can be a .fastq, .fastq.gz or an unaligned .bam (e.g. from dorado demux)
    tuple path(reads_fnp), val(sample_name), path(primers_fnp), path(kmers_sets), path(length_cut_offs_per_target), path(auto_flags_fnp), val (min_len_cut_off)

    output:
    tuple val(sample_name), path("extraction/*.fastq.gz"), emit: fastqs_per_target, optional : true
    tuple val(sample_name), path("extraction/"), emit: sample_dir

    path "${sample_name}_extractionProfile.tab.txt", emit: extraction_profile_per_target
    path "${sample_name}_extractionStats.tab.txt", emit: extraction_stats

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def reads_arg = reads_fnp.name.endsWith(".bam") ? "--bam ${reads_fnp}"
        : reads_fnp.name.endsWith(".gz") ? "--fastqgz ${reads_fnp}"
        : "--fastq ${reads_fnp}"

    """
    SeekDeep extractorByKmerMatching \
            --dout extraction \
            --uniqueKmersPerTarget ${kmers_sets} \
            --rename \
            ${reads_arg} \
            --id ${primers_fnp} \
            --sampleName ${sample_name} \
            --lenCutOffs ${length_cut_offs_per_target} \
            --minLenCutOff ${min_len_cut_off}\
            ${extra_args} \$(cat ${auto_flags_fnp})
    ln -s extraction/extractionProfile.tab.txt ${sample_name}_extractionProfile.tab.txt
    ln -s extraction/extractionStats.tab.txt ${sample_name}_extractionStats.tab.txt

    """
}
